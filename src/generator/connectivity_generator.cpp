#include "connectivity_generator.h"
#include "curve/hermite_init.h"
#include "curve/curve_utils.h"
#include "constraints/fence_check.h"
#include "utils.h"
#include "optimizer/sdf_field.h"
#include "constraints/infeasibility_detector.h"
#include "utils/quadtree.h"
#include <chrono>
#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <memory>

// ── IntersectionInput helpers ─────────────────────────────────────────────────
const bool IntersectionInput::IsEntryLane(const LaneId& id) const {
    for (auto& lg : lane_groups)
        if (std::find(lg.lanes.begin(), lg.lanes.end(), id) != lg.lanes.end())
            return lg.role == GroupRole::Entry;
    return false;
}

const bool IntersectionInput::IsEntryLaneEdge(const LaneEdgeId& id) const {
    for (auto& lg : lane_groups)
        if (std::find(lg.boundaries.begin(), lg.boundaries.end(), id) != lg.boundaries.end())
            return lg.role == GroupRole::Entry;
    for (auto& lg : lane_groups)
        for (auto& lid : lg.lanes) {
            auto l = findLane(lid);
            if (l && (l->left_edge_id == id || l->right_edge_id == id))
                return lg.role == GroupRole::Entry;
        }
    return false;
}

const Lane* IntersectionInput::findLane(const LaneId& id) const {
    for (auto& l : lanes) if (l.id == id)return &l;
    return nullptr;
}

const LaneGroup* IntersectionInput::findGroup(const LaneGroupId& id) const {
    for (auto& g : lane_groups) if (g.id == id)return &g;
    return nullptr;
}

const LaneEdge* IntersectionInput::findEdge(const LaneEdgeId& id) const {
    for (auto& e : lane_edges) if (e.id == id)return &e;
    return nullptr;
}

bool IntersectionInput::laneGroupExists(const LaneGroupId& id) const {
    return findGroup(id) != nullptr;
}

std::pair<Vec2d, Vec2d> IntersectionInput::entryPtDir(const LaneId& lid) const {
    auto* l = findLane(lid);
    if (!l || l->geometry.points.empty()) {
        std::cout << "[WARN] entrylane:" << lid << " no geometry!\n";
        return {Vec2d(0, 0), Vec2d(1, 0)};
    }
    return {entryLinePoint(l->geometry.points), entryLineTangent(l->geometry.points)};
}

std::pair<Vec2d, Vec2d> IntersectionInput::exitPtDir(const LaneId& lid) const {
    auto* l = findLane(lid);
    if (!l || l->geometry.points.empty()) {
        std::cout << "[WARN] exitlane:" << lid << " no geometry!\n";
        return {Vec2d(10, 0), Vec2d(1, 0)};
    }
    return {exitLinePoint(l->geometry.points), exitLineTangent(l->geometry.points)};
}

// ── GlobalCoordinator ─────────────────────────────────────────────────────────
//  Generation order:
//    Priority 0: Straight (anchor curves; generated first)
//    Priority 1: TurnLeft, TurnRight (adjust around straights)
//    Priority 2: UTurnLeft, UTurnRight (most outer; generated last)
//
//  Within each priority group, sort by CLUSTER SPATIAL ORDER:
//    - For each connectivity, look up its position in the entry group's
//      spatial-sorted list (from ClusterOrderSolver::entryGroupOrder()).
//    - Descending sort puts index 0 = LEFTMOST curve.
//    - "Innermost" = closest to straight = highest index among left-turns,
//      lowest index among right-turns, i.e., CLOSEST TO MIDDLE INDEX.
//    - Generate by ascending |distance from middle| so inner curves come first.
//
//  This ensures that when an outer curve is generated, its inner neighbour
//  is already in `done` and will be passed as a sibling with correct
//  expected_side, giving the optimizer a correct reference from the start.
// ─────────────────────────────────────────────────────────────────────────────
//
//
// "排序时转向信息优先按方位计算为准,而不要强依赖输入数据中的ConnTurnType输入信息数据中转向类型不一定准确".
// Compute turn type geometrically from endpoint tangents, not from conn.turn_type.
//
// Algorithm:
//   t0 = entry lane tangent (pointing INTO the intersection)
//   t1 = exit  lane tangent (pointing OUT of the intersection)
//   • T0·T1 < −0.5  → anti-parallel exit → U-turn  (priority 2)
//   • |cross2d(t0, path_dir)| < 0.25  → near-straight (priority 0)
//   • cross2d(t0, path_dir) > 0 → left turn, < 0 → right turn  (priority 1)
// where path_dir = (exit_pt − entry_pt).normalized()
static int globalTurnPriorityGeometric(
    const Connectivity& conn, const IntersectionInput& inp) {
    // Get endpoint geometry
    auto kv0 = inp.entryPtDir(conn.entry_lane_id);
    auto& p0 = kv0.first;
    auto& t0 = kv0.second;
    auto kv1 = inp.exitPtDir (conn.exit_lane_id);
    auto& p1 = kv1.first;
    auto& t1 = kv1.second;

    if (t0.norm() < 1e-9 || t1.norm() < 1e-9) {
        // Fallback to declared type when geometry is degenerate
        switch (conn.turn_type) {
        case ConnTurnType::Straight:   return 0;
        case ConnTurnType::UTurnLeft:
        case ConnTurnType::UTurnRight: return 2;
        default:                       return 1;
        }
    }
    t0.normalize(); t1.normalize();

    // U-turn: exit tangent anti-parallel to entry tangent
    if (t0.dot(t1) < -0.5) return 2;

    // Distinguish straight / turn by bearing from entry to exit
    Vec2d d = p1 - p0;
    if (d.norm() < 1e-9) return 0;
    d.normalize();
    double c = cross2d(t0, d); // positive → left, negative → right
    if (std::abs(c) < 0.25)   return 0; // straight (< ~15°)
    return 1;                            // left or right turn
}

// Legacy wrapper kept for readability – not used for priority any more.
static int globalTurnPriority(ConnTurnType t) {
    switch (t) {
    case ConnTurnType::Straight: return 0;
    case ConnTurnType::TurnLeft: return 1;
    case ConnTurnType::TurnRight: return 1;
    case ConnTurnType::UTurnLeft: return 2;
    case ConnTurnType::UTurnRight: return 2;
    default: return 0;
    }
}

void GlobalCoordinator::build(
    const std::vector<Connectivity>& conns, const IntersectionInput& inp, const ClusterOrderSolver& cs) {
    // Map connId → spatial index within its entry group (from solver)
    std::unordered_map<ConnId, int> spatial_idx;
    std::unordered_map<ConnId, int> group_size;
    for (auto& kv : cs.entryGroupOrder()) {
        auto& gid=kv.first; auto& sorted=kv.second;
        int sz = (int)sorted.size();
        for (int i = 0; i < sz; ++i) {
            spatial_idx[sorted[i]] = i;
            group_size [sorted[i]] = sz;
        }
    }

    // Group by geometric turn priority (not conn.turn_type which may be wrong)
    std::map<int, std::vector<const Connectivity*>> pm;
    for (auto& c : conns)
        pm[globalTurnPriorityGeometric(c, inp)].push_back(&c);

    groups_.clear();
    for (auto& kv : pm) {
        auto& pri = kv.first; auto& cv = kv.second;
        // Sort within priority: innermost (smallest |spatial_idx - sz/2|) first
        // Ties: sort by enterGroupId (group members together), then by spatial_idx
        std::stable_sort(cv.begin(), cv.end(), [&](const Connectivity* a, const Connectivity* b) {
            int ia = spatial_idx.count(a->id) ? spatial_idx[a->id] : 0;
            int ib = spatial_idx.count(b->id) ? spatial_idx[b->id] : 0;
            int sa = group_size .count(a->id) ? group_size [a->id] : 1;
            int sb = group_size .count(b->id) ? group_size [b->id] : 1;
            // "Distance from middle" = |i - (sz-1)/2.0|
            double da = std::abs(ia - (sa-1)/2.0);
            double db = std::abs(ib - (sb-1)/2.0);
            if (std::abs(da-db) > 0.5) return da < db;   // inner before outer
            if (a->enterGroupId != b->enterGroupId)
                return a->enterGroupId < b->enterGroupId;  // same group together
            return ia < ib;
        });

        OptGroup g; g.priority = pri;
        for (auto* c : cv) g.conn_ids.push_back(c->id);
        groups_.push_back(std::move(g));
    }
}

void GlobalCoordinator::addSoftObstacles(SDFField&,
    const std::vector<ConnectivityCurve>&, double) const {}

// ── ConnectivityGenerator ─────────────────────────────────────────────────────

ConnectivityGenerator::ConnectivityGenerator(const LBFGSConfig& cfg) : solver_(cfg) {}

// ─────────────────────────────────────────────────────────────────────────────
//  buildSiblings
//
//  For each already-generated curve:
//  - Determine cluster membership (same enterGroup OR same exitGroup)
//  - Set exempt_a1: true only if NEITHER enterGroup NOR exitGroup matches
//    (= cross-traffic from a different road arm = structural intersection)
//  - Set expected_side from ClusterOrderSolver::expectedSideOf(cid, id)
//    which now returns CORRECT sign because sort is DESCENDING
//  - Set ref_perp from ClusterOrderSolver::refPerpOf()
// ─────────────────────────────────────────────────────────────────────────────
std::vector<SiblingCurve> ConnectivityGenerator::buildSiblings(
    const ConnId& id, const std::unordered_map<ConnId,BezierCurve>& done,
    const ClusterOrderSolver& cs, const std::vector<Connectivity>& conns) const {
    std::vector<SiblingCurve> sibs;
    for (auto& kv : done) {
        auto& cid = kv.first;
        auto& curve = kv.second;
        if (cid == id)
            continue;

        SiblingCurve s;
        s.curve = curve;

        auto ex = cs.exemptionOf(id, cid);
        bool in_cluster = cs.pairExists(id, cid);
        if (!in_cluster) {
            // No cluster relationship → different arm → structural cross-traffic
            s.exempt_a1 = true;
        } else if (ex == CrossExemption::StructuralCross) {
            // Topological inversion → must cross → exempt
            s.exempt_a1 = true;
        } else {
            // In same cluster, constrained pair → apply penalty
            s.exempt_a1 = false;
        }

        s.exempt_a2_radius = (ex == CrossExemption::ObstacleCross) ? 1.5 : 0.0;

        // expected_side from cid's perspective relative to id:
        // cs.expectedSideOf(cid, id):
        //   +1 → sibling(cid) is LEFT  of current(id)
        //   -1 → sibling(cid) is RIGHT of current(id)
        s.expected_side = cs.expectedSideOf(cid, id);

        // Fixed lateral reference for evalCluster
        s.ref_perp = cs.refPerpOf(id, cid);

        sibs.push_back(std::move(s));
    }
    return sibs;
}

// ─────────────────────────────────────────────────────────────────────────────
//  computeRequiredLateralOffset
//
//  For the initial curve of connectivity `id`, check if the direct Bezier arc
//  from p0→p1 would violate any inner sibling's ordering.  Returns a (side,
//  magnitude) lateral offset to apply if needed:  side = +1 (push LEFT),
//  side = -1 (push RIGHT).
//
//  Only inner siblings are considered (those that should be to the RIGHT of
//  the current curve, i.e., expected_side = -1).  If the direct arc is ALREADY
//  outside (LEFT of) all of them → no offset needed.
//
//  ref_perp: the lateral reference direction (arm left-normal).
// ─────────────────────────────────────────────────────────────────────────────
static bool needsLateralOffset(
    const Vec2d& p0, const Vec2d& p1,
    const std::vector<SiblingCurve>& siblings,
    int& out_side, double& out_offset_m)
{
    out_side = 0; out_offset_m = 0.0;

    // Find the sibling that should be to my RIGHT (expected_side=-1)
    // and check if I'm actually to its LEFT (correct) or not.
    // The most critical one: the inner sibling with the LARGEST lat value
    // (most to the LEFT among inner siblings = the one I must stay outside of).
    double worst_violation = 0.0;
    Vec2d ref_perp(0,0);
    int required_push_side = 0;

    for (auto& s : siblings) {
        if (s.exempt_a1) continue;
        if (s.expected_side != -1) continue;  // only inner siblings (should be to my RIGHT)
        if (s.ref_perp.norm() < 1e-9) continue;

        // Sibling's mean lateral position in ref_perp frame
        auto sib_pts = s.curve.sampleByArcLength(20);
        if (sib_pts.empty()) continue;
        double sib_lat = 0;
        for (auto& pt : sib_pts) sib_lat += pt.dot(s.ref_perp);
        sib_lat /= sib_pts.size();

        // My (current curve) estimated mean lateral position (midpoint of direct arc)
        Vec2d mid = 0.5*(p0+p1);
        double my_lat = mid.dot(s.ref_perp);

        // For expected_side=-1 (sibling should be to my RIGHT = sibling has SMALLER lat):
        // Violation if my_lat < sib_lat + MARGIN (I'm not far enough to the LEFT)
        double margin = 1.0;  // Increased from 0.5 to 1.0m clearance for initial curve to prevent intersections
        double viol = sib_lat + margin - my_lat;  // positive = violation
        if (viol > worst_violation) {
            worst_violation = viol;
            ref_perp = s.ref_perp;
            required_push_side = +1;  // push LEFT (increase lat)
        }
    }

    // Also check siblings that should be to my LEFT (expected_side=+1):
    // These are the outer siblings I must stay inside of.
    for (auto& s : siblings) {
        if (s.exempt_a1) continue;
        if (s.expected_side != +1) continue;  // only outer siblings (should be to my LEFT)
        if (s.ref_perp.norm() < 1e-9) continue;

        auto sib_pts = s.curve.sampleByArcLength(20);
        if (sib_pts.empty()) continue;
        double sib_lat = 0;
        for (auto& pt : sib_pts) sib_lat += pt.dot(s.ref_perp);
        sib_lat /= sib_pts.size();

        Vec2d mid = 0.5*(p0+p1);
        double my_lat = mid.dot(s.ref_perp);

        double margin = 1.0;  // Increased from 0.5 to 1.0m
        double viol = my_lat + margin - sib_lat;  // positive = violation (I'm too far LEFT)
        // Only apply if stronger than current violation from other direction
        if (viol > worst_violation && required_push_side == 0) {
            worst_violation = viol;
            ref_perp = s.ref_perp;
            required_push_side = -1;  // push RIGHT (decrease lat)
        }
    }

    // Additionally, check for potential crossings with same-group siblings
    // by examining both endpoints and considering curve trajectories
    for (auto& s : siblings) {
        if (s.exempt_a1) continue;
        // Check for potential crossings regardless of expected_side for same-group connections
        if (s.ref_perp.norm() < 1e-9) continue;

        // Calculate distances between current path and sibling path at start/end points
        Vec2d sib_start = s.curve.startPt();
        Vec2d sib_end = s.curve.endPt();

        // Check if this is a same-group connection (same entry or exit group)
        // If so, we may need additional separation
        double start_sep = std::abs((p0 - sib_start).dot(s.ref_perp));
        double end_sep = std::abs((p1 - sib_end).dot(s.ref_perp));

        // If already quite close laterally and in same group context, enforce more separation
        if (start_sep < 2.0 && end_sep < 2.0) {
            double required_sep = 2.0; // Enforce 2m separation for potentially conflicting curves
            double current_sep = std::min(start_sep, end_sep);
            if (required_sep > current_sep) {
                double needed_sep = required_sep - current_sep;
                // Determine direction based on relative positions
                Vec2d my_mid = 0.5 * (p0 + p1);
                Vec2d sib_mid = 0.5 * (sib_start + sib_end);

                double my_lat = my_mid.dot(s.ref_perp);
                double sib_lat = sib_mid.dot(s.ref_perp);

                int sep_side = (my_lat > sib_lat) ? -1 : +1;  // Push away from sibling

                if (needed_sep > worst_violation) {
                    worst_violation = needed_sep;
                    ref_perp = s.ref_perp;
                    required_push_side = sep_side;
                }
            }
        }
    }

    if (worst_violation < 1e-3 || required_push_side == 0) return false;

    out_side     = required_push_side;
    out_offset_m = worst_violation + 0.5;  // Increased extra margin from 0.3 to 0.5
    return true;
}

BezierCurve ConnectivityGenerator::postProcess(
    const BezierCurve& c, const SDFField& sdf, const Polygon2d& fence, double kmax,
    const Vec2d& t0_orig, const Vec2d& t1_orig, bool skip_elastic_band,
    const Vec2d* p0_exact, const Vec2d* p1_exact) {
    auto refined = adaptiveRefine(c, sdf, kmax);
    BezierCurve cur = refined.curve;
    if (refined.was_split) {
        PenaltyCost cost2;
        cost2.proto = cur;
        cost2.sdf = &sdf;
        cost2.fence = fence;
        cost2.start_tan_dir = t0_orig.norm() > 1e-8 ? t0_orig.normalized() : c.startTan();
        cost2.end_tan_dir = t1_orig.norm() > 1e-8 ? t1_orig.normalized() : c.endTan();
        cost2.obstacle_clearance = 0.0;
        cost2.full_param_mode = (cur.numSegments() > 1);
        cost2.buildCache();
        cur = optimiseCurve(cost2, solver_, cur, 2);
    }

    Vec2d st = t0_orig.norm() > 1e-8 ? t0_orig.normalized() : cur.startTan();
    Vec2d et = t1_orig.norm() > 1e-8 ? t1_orig.normalized() : cur.endTan();
    Vec2d ep0 = p0_exact ? *p0_exact : cur.startPt();
    Vec2d ep1 = p1_exact ? *p1_exact : cur.endPt();

    if (skip_elastic_band) {
        if (!cur.segs.empty()) {
            cur.segs.front().ctrl[0] = ep0;
            cur.segs.back().ctrl[3] = ep1;
            Vec2d& c1 = cur.segs.front().ctrl[1];
            double sl = (cur.segs.front().ctrl[3] - ep0).norm();
            double lam = std::max((c1 - ep0).dot(st), sl * 0.1);
            c1 = ep0 + lam * st;
            Vec2d& c2 = cur.segs.back().ctrl[2];
            double el = (ep1 - cur.segs.back().ctrl[0]).norm();
            double mu = std::max((ep1 - c2).dot(et), el * 0.1);
            c2 = ep1 - mu * et;
        }
        return cur;
    }

    double arc = cur.arcLength();
    int n_samp = std::max(40, (int)(arc / 0.15));
    auto pts = cur.sampleByArcLength(n_samp);
    if ((int)pts.size() >= 3) {
        pts.front() = ep0;
        pts.back() = ep1;

        double max_k = cur.maxCurvature(20);
        double move_st = std::min(0.15, std::max(0.02, max_k * 0.3));
        auto sm = elasticBandSmooth(pts, sdf, fence, kmax, move_st, 80, 0.1);
        sm.front() = ep0;
        sm.back() = ep1;

        cur = rebuildFromSmoothedPts(sm, st, et);
        if (!cur.segs.empty()) {
            cur.segs.front().ctrl[0] = ep0;
            cur.segs.back().ctrl[3] = ep1;
        }
    }
    return cur;
}

void ConnectivityGenerator::validate(
    ConnectivityCurve& cc, const IntersectionInput& input, const SDFField& sdf) const {
    if (!cc.curve) return;
    auto& c = *cc.curve;
    double ms = minSDFAlongCurve(c, sdf, 20);
    cc.violation.max_obstacle_penetration = std::max(0.0, -ms);
    if (!input.area.geometry.outer.empty()) {
        double ov = 0;
        for (auto& pt : c.sampleByArcLength(30))
            if (!polygonContains(input.area.geometry, pt))
                ov = std::max(ov, pointToPolygonDist(pt, input.area.geometry));
        cc.violation.max_fence_overflow = ov;
    }
    if (cc.violation.max_obstacle_penetration > 0.05)
        cc.status = CurveStatus::Degraded;
    else if (!cc.violation.exempt_crosses.empty())
        cc.status = CurveStatus::WarnA2;
    else
        cc.status = CurveStatus::OK;
}

ConnectivityCurve ConnectivityGenerator::generateOne(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf,
    const SDFField& sdf_coarse, const std::vector<SiblingCurve>& siblings) {
    ConnectivityCurve cc;
    cc.id = conn.id;
    cc.entry_lane_id = conn.entry_lane_id;
    cc.exit_lane_id = conn.exit_lane_id;
    cc.turn_type = conn.turn_type;

    auto _entry = input.entryPtDir(conn.entry_lane_id);
    Vec2d p0 = _entry.first;
    Vec2d t0 = _entry.second;
    auto _exit = input.exitPtDir(conn.exit_lane_id);
    Vec2d p1 = _exit.first;
    Vec2d t1 = _exit.second;

    double lw = 3.5;
    if (auto* l = input.findLane(conn.entry_lane_id))
        lw = l->width;

    auto pre = preCheck(sdf_coarse, input.area.geometry, p0, p1, lw, input.boundaries);
    if (pre.type == ViolationInfo::InfeasibilityType::TopologicalBlock)
        return makeFallbackCurve(pre, conn, p0, p1);

    // ── Build initial curve with sibling topology awareness ──────────────────
    //
    // Step A: Collect sibling polys for obstacle-bypass side selection.
    //   Include ALL non-exempt siblings (same entry OR same exit group).
    //   This prevents choosing a bypass direction that conflicts with existing curves.
    std::vector<std::vector<Vec2d>> sib_polys_for_init;
    for (auto& sib : siblings) {
        if (!sib.exempt_a1)
            sib_polys_for_init.push_back(sib.curve.sampleByArcLength(20));
    }
    BezierCurve initial = buildInitialCurve(
        p0, t0, p1, t1, sdf, input.area.geometry, sib_polys_for_init);

    // Step B: Check if initial curve needs lateral adjustment for sibling ordering.
    //   If the direct arc would place us on the wrong side of an inner sibling,
    //   rebuild with a lateral offset so the optimizer starts on the correct side.
    //   (Topology correction from wrong side requires crossing siblings, which
    //    is hard for gradient-based optimizer to escape.)
    {
        int push_side = 0;
        double push_m = 0.0;
        if (needsLateralOffset(p0, p1, siblings, push_side, push_m) && push_m > 0.1) {
            // Compute offset direction: push perpendicularly from the path p0→p1
            Vec2d along = (p1-p0);
            if (along.norm() > 1e-9) along.normalize();
            Vec2d left_normal{-along[1], along[0]};
            // Use the sibling's ref_perp as the offset direction if available,
            // otherwise use the path left normal.
            Vec2d offset_dir = left_normal;
            for (auto& s : siblings) {
                if (!s.exempt_a1 && s.ref_perp.norm() > 1e-9 && s.expected_side == -1) {
                    offset_dir = s.ref_perp;
                    break;
                }
            }
            double magnitude = std::min(push_m, 4.0);  // cap at 4m
            Vec2d apex_shift = offset_dir * (push_side * magnitude);

            // Build a new initial curve via a 3-point arch passing through
            // an apex offset from the midpoint
            Vec2d mid = 0.5*(p0+p1) + apex_shift;
            // Use hermite_init arch builder: p0 → mid → p1
            std::vector<Vec2d> knots = {p0, mid, p1};
            Vec2d along_norm = (p1-p0).norm()>1e-9 ? (p1-p0).normalized() : Vec2d(1,0);
            Vec2d apex_tan   = along_norm; // approximate apex tangent
            std::vector<Vec2d> tans = {t0.normalized(), apex_tan, t1.normalized()};
            BezierCurve arch = makeCurveFromKnots(knots, tans, 0.35);

            // Use the adjusted arch only if it's obstacle-free
            bool arch_clear = true;
            if (sdf.valid()) {
                for (auto& seg : arch.segs)
                    for (int i=1; i<15; ++i) {
                        auto _q = sdf.queryWithGrad(seg.evaluate((double)i/15));
                        if (_q.first < 0) { arch_clear = false; break; }
                    }
            }
            if (arch_clear) initial = arch;
        }
    }

    // ── Adaptive alpha: if initial curve crosses siblings, try shorter handles ──
    // Root cause: buildInitialCurve uses alpha=0.4 and returns immediately when
    // obstacle-free, ignoring sibling crossings. For same-exit curves (e.g. a right
    // turn whose northward entry tangent arch briefly crosses a straight's path),
    // a shorter handle avoids the crossing from the start, letting the optimizer
    // refine from a correct topology rather than fighting out of a local minimum.
    {
        std::vector<double> try_alphas = {0.2, 0.12};
        for (double alpha : try_alphas) {
            // Check if initial curve crosses any constrained sibling
            bool has_cross = false;
            auto init_pts = initial.sampleByArcLength(32);
            for (auto& sib : siblings) {
                if (sib.exempt_a1) continue;
                auto sib_pts = sib.curve.sampleByArcLength(32);
                if (curvesIntersectBusiness(initial, sib.curve, 1.5)) {
                    has_cross = true; break;
                }
            }
            if (!has_cross) break;   // current initial is already crossing-free
            // Rebuild initial with shorter handles (G1 maintained)
            BezierSegment shorter = makeCubicG1(
                p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                alpha);
            // Verify shorter arc is still obstacle-free
            bool arc_ok = true;
            for (int k=1; k<16; ++k) {
                auto kv = sdf.queryWithGrad(shorter.evaluate((double)k/16));
                if (kv.first < 0) { arc_ok = false; break; }
            }
            if (arc_ok) {
                BezierCurve c_short; c_short.segs.push_back(shorter);
                initial = c_short;
            }
        }
    }

    // ── Optimise ──────────────────────────────────────────────────────────────
    PenaltyCost cost;
    cost.proto = initial;
    cost.sdf = &sdf;
    cost.boundaries = input.boundaries;
    cost.fence = input.area.geometry;
    cost.siblings = siblings;
    cost.crosswalks = input.crosswalks;
    cost.obstacle_clearance = 0.0;
    cost.start_tan_dir = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    cost.end_tan_dir = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    cost.full_param_mode = (initial.numSegments() > 1);

    // raised outer_iters 5 → 7; with initial cluster weight 20.0 the
    // ramp now reaches 96 within 4 iterations, giving 3 extra high-weight steps
    // to resolve any residual adjacent-lane crossings.
    BezierCurve opt = optimiseCurve(cost, solver_, initial, /*outer_iters=*/7);

    bool skip_band = true; // always skip: preserve G1, rely on optimizer
    BezierCurve final_c = postProcess(opt, sdf, input.area.geometry, 0.25, t0, t1, skip_band, &p0, &p1);
    cc.curve = std::make_shared<BezierCurve>(final_c);
    validate(cc, input, sdf);
    if (pre.narrow_passage && cc.status == CurveStatus::OK)
        cc.status = CurveStatus::WarnA2;
    return cc;
}

std::vector<ConnectivityCurve> ConnectivityGenerator::generate(
    const IntersectionInput& input, SDFField& sdf, double* out_ms) {
    auto t0 = std::chrono::steady_clock::now();

    SDFField sdf_coarse;
    auto roi = input.area.geometry.empty() ? BoundingBox2d{} : input.area.geometry.bbox();
    if (roi.width() < 1) {
        for (auto& l : input.lanes) {
            for (auto& p : l.geometry.points)
                roi.expand(p);
        }
        roi.min_pt -= Vec2d(20, 20);
        roi.max_pt += Vec2d(20, 20);
    }
    sdf_coarse.build(roi, input.obstacles, 0.5, 0.4);

    // Build group-based cluster solver
    cluster_solver_.build(input.connectivities, input.lanes, input.lane_groups);

    // Build generation order
    GlobalCoordinator coord;
    coord.build(input.connectivities, input, cluster_solver_);

    std::vector<ConnectivityCurve> results;
    results.reserve(input.connectivities.size());
    std::unordered_map<ConnId, const Connectivity*> cmap;
    for (auto& c : input.connectivities) {
        cmap[c.id] = &c;
    }
    std::unordered_map<ConnId, BezierCurve> done;
    for (auto& group : coord.groups()) {
        for (auto& cid : group.conn_ids) {
            auto* conn = cmap[cid];
            if (!conn) continue;

            auto sibs = buildSiblings(cid, done, cluster_solver_, input.connectivities);

            auto cc = generateOne(*conn, input, sdf, sdf_coarse, sibs);
            if (cc.curve) {
                done[cid] = *cc.curve;
            }

            results.push_back(std::move(cc));
        }
        // After each priority group: mark obstacle-adjacent crossings as soft
        cluster_solver_.checkAndMarkA2(done, sdf, 1.5);
    }

    auto t1 = std::chrono::steady_clock::now();
    if (out_ms) *out_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return results;
}
