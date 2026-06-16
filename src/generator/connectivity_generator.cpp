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
#include <limits>
#include <memory>
#include <unordered_set>

namespace isg {

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

struct LaneDirectionSample {
    LaneId lane_id;
    Vec2d direction{1, 0};
    int lane_order = 0;
    int group_index = 0;
};

static double directionAngleDiff(const Vec2d& a, const Vec2d& b) {
    if (a.norm() < 1e-10 || b.norm() < 1e-10)
        return M_PI;
    double c = a.normalized().dot(b.normalized());
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c);
}

static const Lane* findLaneById(const std::vector<Lane>& lanes, const LaneId& id) {
    for (const auto& lane : lanes)
        if (lane.id == id)
            return &lane;
    return nullptr;
}

static Lane* findLaneById(std::vector<Lane>& lanes, const LaneId& id) {
    for (auto& lane : lanes)
        if (lane.id == id)
            return &lane;
    return nullptr;
}

static Vec2d laneDirectionForRole(const Lane& lane, GroupRole role) {
    Vec2d dir = role == GroupRole::Entry
        ? entryLineTangent(lane.geometry.points)
        : exitLineTangent(lane.geometry.points);
    return dir.norm() > 1e-10 ? dir.normalized() : Vec2d(1, 0);
}

static std::vector<LaneDirectionSample> collectLaneDirections(
    const LaneGroup& group, const std::vector<Lane>& lanes) {
    std::vector<LaneDirectionSample> samples;
    samples.reserve(group.lanes.size());
    std::unordered_set<LaneId> seen_lanes;
    for (int i = 0; i < (int)group.lanes.size(); ++i) {
        if (!seen_lanes.insert(group.lanes[i]).second)
            continue;
        const Lane* lane = findLaneById(lanes, group.lanes[i]);
        if (!lane || lane->geometry.points.empty())
            continue;
        LaneDirectionSample s;
        s.lane_id = lane->id;
        s.direction = laneDirectionForRole(*lane, group.role);
        s.lane_order = lane->laneOrder;
        s.group_index = i;
        samples.push_back(s);
    }
    return samples;
}

static const LaneDirectionSample* innermostLaneDirection(
    const std::vector<LaneDirectionSample>& samples) {
    if (samples.empty())
        return nullptr;
    return &*std::min_element(
        samples.begin(), samples.end(),
        [](const LaneDirectionSample& a, const LaneDirectionSample& b) {
            if (a.lane_order != b.lane_order)
                return a.lane_order < b.lane_order;
            return a.group_index < b.group_index;
        });
}

static Vec2d meanDirectionFromSupport(
    const std::vector<LaneDirectionSample>& samples,
    const std::vector<int>& support,
    const Vec2d& fallback) {
    Vec2d sum(0, 0);
    for (int idx : support) {
        if (idx < 0 || idx >= (int)samples.size())
            continue;
        if (samples[idx].direction.norm() > 1e-10)
            sum += samples[idx].direction.normalized();
    }
    return sum.norm() > 1e-10 ? sum.normalized() : fallback;
}

static Vec2d groupUnifiedDirection(
    const LaneGroup& group, const std::vector<Lane>& lanes, double threshold_deg) {
    auto samples = collectLaneDirections(group, lanes);
    const auto* inner = innermostLaneDirection(samples);
    Vec2d fallback = inner ? inner->direction : Vec2d(1, 0);

    if (samples.size() <= 2)
        return fallback;

    double threshold_rad = std::max(0.0, std::min(180.0, threshold_deg)) * M_PI / 180.0;
    int best_idx = -1;
    int best_count = 0;
    for (int i = 0; i < (int)samples.size(); ++i) {
        int count = 0;
        for (int j = 0; j < (int)samples.size(); ++j) {
            if (i == j)
                continue;
            if (directionAngleDiff(samples[i].direction, samples[j].direction) <= threshold_rad)
                ++count;
        }
        if (count > best_count ||
            (count == best_count && count > 0 &&
             (best_idx < 0 ||
              samples[i].lane_order < samples[best_idx].lane_order ||
              (samples[i].lane_order == samples[best_idx].lane_order &&
               samples[i].group_index < samples[best_idx].group_index)))) {
            best_idx = i;
            best_count = count;
        }
    }

    if (best_idx < 0 || best_count == 0)
        return fallback;

    std::vector<int> support;
    support.push_back(best_idx);
    for (int j = 0; j < (int)samples.size(); ++j) {
        if (j == best_idx)
            continue;
        if (directionAngleDiff(samples[best_idx].direction, samples[j].direction) <= threshold_rad)
            support.push_back(j);
    }
    return meanDirectionFromSupport(samples, support, fallback);
}

static void setEntryLaneDirection(Lane& lane, const Vec2d& direction) {
    if (direction.norm() < 1e-10)
        return;
    Vec2d dir = direction.normalized();
    auto& pts = lane.geometry.points;
    if (pts.empty())
        return;
    Vec2d p = pts.back();
    if (pts.size() == 1)
        pts.insert(pts.begin(), p - dir);
    else
        pts.insert(pts.end() - 1, p - dir);
}

static void setExitLaneDirection(Lane& lane, const Vec2d& direction) {
    if (direction.norm() < 1e-10)
        return;
    Vec2d dir = direction.normalized();
    auto& pts = lane.geometry.points;
    if (pts.empty())
        return;
    Vec2d p = pts.front();
    if (pts.size() == 1)
        pts.push_back(p + dir);
    else
        pts.insert(pts.begin() + 1, p + dir);
}

static void applyGroupUnifiedDirections(
    IntersectionInput& input, const ConnectivityDirectionConfig& cfg) {
    if (cfg.mode != ConnectivityDirectionMode::GroupUnified)
        return;

    std::unordered_map<LaneGroupId, Vec2d> group_dirs;
    for (const auto& group : input.lane_groups)
        group_dirs[group.id] = groupUnifiedDirection(
            group, input.lanes, cfg.group_similarity_angle_deg);

    for (const auto& group : input.lane_groups) {
        auto it = group_dirs.find(group.id);
        if (it == group_dirs.end())
            continue;
        for (const auto& lane_id : group.lanes) {
            Lane* lane = findLaneById(input.lanes, lane_id);
            if (!lane)
                continue;
            if (group.role == GroupRole::Entry)
                setEntryLaneDirection(*lane, it->second);
            else
                setExitLaneDirection(*lane, it->second);
        }
    }
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

ConnectivityGenerator::ConnectivityGenerator(
    const LBFGSConfig& cfg,
    const ConnectivityDirectionConfig& direction_cfg)
    : solver_(cfg), direction_cfg_(direction_cfg) {}

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
    const ClusterOrderSolver& cs, const std::vector<Connectivity>& conns,
    bool constrained_only,
    const std::unordered_set<ConnId>* fixed_shape_ids) const {
    std::vector<SiblingCurve> sibs;
    for (auto& kv : done) {
        auto& cid = kv.first;
        auto& curve = kv.second;
        if (cid == id)
            continue;

        SiblingCurve s;
        s.curve = curve;
        s.fixed_shape = fixed_shape_ids && fixed_shape_ids->count(cid) > 0;

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
        if (constrained_only && s.exempt_a1)
            continue;

        s.exempt_a2_radius = (ex == CrossExemption::ObstacleCross) ? 1.5 : 0.0;

        // expected_side from cid's perspective relative to id:
        // cs.expectedSideOf(cid, id):
        //   +1 → sibling(cid) is LEFT  of current(id)
        //   -1 → sibling(cid) is RIGHT of current(id)
        s.expected_side = cs.expectedSideOf(cid, id);

        // Fixed lateral reference for evalCluster
        s.ref_perp = cs.refPerpOf(id, cid);

        // Propagate shared_endpoint flag so evalCluster uses wider skip zone
        s.shared_endpoint = cs.isSharedEndpoint(id, cid);

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
        double margin = 0.25;
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

        double margin = 0.25;
        double viol = my_lat + margin - sib_lat;  // positive = violation (I'm too far LEFT)
        // Only apply if stronger than current violation from other direction
        if (viol > worst_violation && required_push_side == 0) {
            worst_violation = viol;
            ref_perp = s.ref_perp;
            required_push_side = -1;  // push RIGHT (decrease lat)
        }
    }

    if (worst_violation < 1e-3 || required_push_side == 0) return false;

    out_side     = required_push_side;
    out_offset_m = worst_violation + 0.1;
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

static double minSDFAlongCurveAdaptive(const BezierCurve& curve, const SDFField& sdf);
static bool curveIntersectsObstacles(
    const BezierCurve& curve, const std::vector<Obstacle>& obstacles, Vec2d* out);
static bool curveIntersectsBoundaries(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries);

void ConnectivityGenerator::validate(
    ConnectivityCurve& cc, const IntersectionInput& input, const SDFField& sdf) const {
    if (!cc.curve) return;
    auto& c = *cc.curve;
    double ms = minSDFAlongCurveAdaptive(c, sdf);
    cc.violation.max_obstacle_penetration = std::max(0.0, -ms);
    Vec2d obstacle_hit;
    if (curveIntersectsObstacles(c, input.obstacles, &obstacle_hit)) {
        cc.violation.max_obstacle_penetration =
            std::max(cc.violation.max_obstacle_penetration, 0.01);
        cc.violation.reason = "curve intersects obstacle";
    }
    if (!input.area.geometry.outer.empty()) {
        double ov = 0;
        int n = std::max(40, (int)std::ceil(c.arcLength() / 0.25) + 1);
        n = std::min(n, 240);
        for (auto& pt : c.sampleByArcLength(n))
            if (!polygonContains(input.area.geometry, pt))
                ov = std::max(ov, pointToPolygonDist(pt, input.area.geometry));
        cc.violation.max_fence_overflow = ov;
    }
    bool self_cross = curveSelfIntersectsBusiness(c, 1.0);
    if (self_cross)
        cc.violation.reason = "curve self-intersects away from endpoints";
    bool boundary_cross = curveIntersectsBoundaries(c, input.boundaries);
    if (boundary_cross)
        cc.violation.reason = "curve intersects boundary away from endpoints";
    if (self_cross || boundary_cross || cc.violation.max_obstacle_penetration > 0.05)
        cc.status = CurveStatus::Degraded;
    else if (!cc.violation.exempt_crosses.empty())
        cc.status = CurveStatus::WarnA2;
    else
        cc.status = CurveStatus::OK;
}

static std::unordered_map<ConnId, size_t> resultIndexById(
    const std::vector<ConnectivityCurve>& results) {
    std::unordered_map<ConnId, size_t> idx;
    for (size_t i = 0; i < results.size(); ++i)
        idx[results[i].id] = i;
    return idx;
}

static std::unordered_map<ConnId, BezierCurve> curveMapFromResults(
    const std::vector<ConnectivityCurve>& results) {
    std::unordered_map<ConnId, BezierCurve> curves;
    for (auto& cc : results)
        if (cc.curve)
            curves[cc.id] = *cc.curve;
    return curves;
}

static bool hasFixedGeometry(const Connectivity& conn) {
    return conn.geometry.points.size() >= 2;
}

static BezierCurve fixedGeometryToCurve(
    const Connectivity& conn, const IntersectionInput& input) {
    (void)input;
    BezierCurve curve;
    const auto& pts = conn.geometry.points;
    if (pts.size() < 2)
        return curve;

    for (int i = 0; i + 1 < (int)pts.size(); ++i) {
        Vec2d chord = pts[i + 1] - pts[i];
        double len = chord.norm();
        if (len < 1e-8)
            continue;
        Vec2d seg_dir = chord / len;
        curve.segs.push_back(makeCubicG1(pts[i], seg_dir, pts[i + 1], seg_dir, 1.0 / 3.0));
    }
    return curve;
}

static ConnectivityCurve makeFixedGeometryCurve(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf) {
    ConnectivityCurve cc;
    cc.id = conn.id;
    cc.entry_lane_id = conn.entry_lane_id;
    cc.exit_lane_id = conn.exit_lane_id;
    cc.turn_type = conn.turn_type;
    cc.geometry = conn.geometry;
    BezierCurve curve = fixedGeometryToCurve(conn, input);
    if (!curve.empty())
        cc.curve = std::make_shared<BezierCurve>(curve);

    if (cc.curve) {
        double ms = minSDFAlongCurveAdaptive(*cc.curve, sdf);
        cc.violation.max_obstacle_penetration = std::max(0.0, -ms);
        Vec2d obstacle_hit;
        if (curveIntersectsObstacles(*cc.curve, input.obstacles, &obstacle_hit)) {
            cc.violation.max_obstacle_penetration =
                std::max(cc.violation.max_obstacle_penetration, 0.01);
            cc.violation.reason = "fixed geometry intersects obstacle";
            cc.status = CurveStatus::Degraded;
        }
    }
    return cc;
}

static bool fixedGeometryHitsObstacle(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf) {
    if (!hasFixedGeometry(conn))
        return false;
    BezierCurve curve = fixedGeometryToCurve(conn, input);
    if (curve.empty())
        return false;
    (void)sdf;
    return curveIntersectsObstacles(curve, input.obstacles, nullptr);
}

struct SampledCurve {
    std::vector<Vec2d> pts;
    BoundingBox2d bbox;
    Vec2d start{0, 0};
    Vec2d end{0, 0};
};

struct SampledSiblingCurve {
    SampledCurve sampled;
    bool exempt_a1 = false;
};

static SampledCurve sampleCurveForIntersections(const BezierCurve& curve, int n = 32) {
    SampledCurve s;
    s.start = curve.startPt();
    s.end = curve.endPt();
    int adaptive_n = std::max(n, (int)std::ceil(curve.arcLength() / 0.20) + 1);
    adaptive_n = std::min(adaptive_n, 240);
    s.pts = curve.sampleByArcLength(adaptive_n);
    for (const auto& pt : s.pts)
        s.bbox.expand(pt);
    return s;
}

static double minSDFAlongCurveAdaptive(const BezierCurve& curve, const SDFField& sdf) {
    if (!sdf.valid())
        return 1e18;
    int n = std::max(40, (int)std::ceil(curve.arcLength() / 0.15) + 1);
    n = std::min(n, 320);
    double m = 1e18;
    for (const auto& pt : curve.sampleByArcLength(n))
        m = std::min(m, sdf.queryWithGrad(pt).first);
    return m;
}

static bool pointOnPolygonBoundary(const Vec2d& pt, const Polygon2d& poly, double tol = 1e-6) {
    const auto& ring = poly.outer;
    if (ring.size() < 2)
        return false;
    for (int i = 0; i < (int)ring.size(); ++i) {
        const Vec2d& a = ring[i];
        const Vec2d& b = ring[(i + 1) % ring.size()];
        if (pointToSegment(pt, a, b).first <= tol)
            return true;
    }
    return false;
}

static bool pointInsideOrOnPolygon(const Vec2d& pt, const Polygon2d& poly, double tol = 1e-6) {
    return polygonContains(poly, pt) || pointOnPolygonBoundary(pt, poly, tol);
}

static bool polygonSegmentIntersects(
    const Vec2d& a, const Vec2d& b, const Polygon2d& poly, Vec2d* out = nullptr) {
    const auto& ring = poly.outer;
    if (ring.size() < 2)
        return false;
    for (int i = 0; i < (int)ring.size(); ++i) {
        Vec2d isect;
        if (segmentsIntersect(a, b, ring[i], ring[(i + 1) % ring.size()], &isect)) {
            if (out)
                *out = isect;
            return true;
        }
    }
    return false;
}

static const Polygon2d& obstacleHardGeometry(const Obstacle& obs) {
    return obs.geometry.outer.empty() ? obs.buffered_geometry : obs.geometry;
}

static bool curveIntersectsObstacles(
    const BezierCurve& curve, const std::vector<Obstacle>& obstacles, Vec2d* out = nullptr) {
    if (obstacles.empty())
        return false;

    auto sampled = sampleCurveForIntersections(curve, 64);
    if (sampled.pts.size() < 2)
        return false;

    for (const auto& obs : obstacles) {
        const Polygon2d& poly = obstacleHardGeometry(obs);
        if (poly.outer.size() < 3)
            continue;

        BoundingBox2d obs_box = poly.bbox();
        if (!sampled.bbox.intersects(obs_box))
            continue;

        for (int i = 1; i + 1 < (int)sampled.pts.size(); ++i) {
            if (pointInsideOrOnPolygon(sampled.pts[i], poly)) {
                if (out)
                    *out = sampled.pts[i];
                return true;
            }
        }

        for (int i = 0; i + 1 < (int)sampled.pts.size(); ++i) {
            Vec2d isect;
            if (polygonSegmentIntersects(sampled.pts[i], sampled.pts[i + 1], poly, &isect)) {
                if (out)
                    *out = isect;
                return true;
            }
        }
    }
    return false;
}

static double distToSampledEndpoints(const Vec2d& pt, const SampledCurve& a, const SampledCurve& b) {
    double d = (pt - a.start).norm();
    d = std::min(d, (pt - a.end).norm());
    d = std::min(d, (pt - b.start).norm());
    d = std::min(d, (pt - b.end).norm());
    return d;
}

static bool sampledCurvesIntersectBusiness(
    const SampledCurve& a, const SampledCurve& b, double endpoint_tol, Vec2d* out = nullptr) {
    if (a.pts.size() < 2 || b.pts.size() < 2 || !a.bbox.intersects(b.bbox))
        return false;
    for (int ai = 0; ai + 1 < (int)a.pts.size(); ++ai) {
        Vec2d amid = 0.5 * (a.pts[ai] + a.pts[ai + 1]);
        for (int bi = 0; bi + 1 < (int)b.pts.size(); ++bi) {
            Vec2d bmid = 0.5 * (b.pts[bi] + b.pts[bi + 1]);
            if ((amid - bmid).squaredNorm() > 900.0)
                continue;
            Vec2d isect;
            if (!segmentsIntersect(a.pts[ai], a.pts[ai + 1], b.pts[bi], b.pts[bi + 1], &isect))
                continue;
            if (distToSampledEndpoints(isect, a, b) <= endpoint_tol)
                continue;
            if (out)
                *out = isect;
            return true;
        }
    }
    return false;
}

static std::vector<SampledSiblingCurve> sampleSiblingsForIntersections(
    const std::vector<SiblingCurve>& siblings, int n = 32) {
    std::vector<SampledSiblingCurve> sampled;
    sampled.reserve(siblings.size());
    for (const auto& sib : siblings) {
        SampledSiblingCurve s;
        s.exempt_a1 = sib.exempt_a1;
        s.sampled = sampleCurveForIntersections(sib.curve, n);
        sampled.push_back(std::move(s));
    }
    return sampled;
}

static int sampledSiblingCrossCount(
    const SampledCurve& curve, const std::vector<SampledSiblingCurve>& siblings,
    bool constrained_only, double endpoint_tol = 1.5) {
    int count = 0;
    for (const auto& sib : siblings) {
        if (constrained_only && sib.exempt_a1)
            continue;
        if (sampledCurvesIntersectBusiness(curve, sib.sampled, endpoint_tol))
            ++count;
    }
    return count;
}

static int sampledSiblingCrossCount(
    const BezierCurve& curve, const std::vector<SampledSiblingCurve>& siblings,
    bool constrained_only, double endpoint_tol = 1.5) {
    return sampledSiblingCrossCount(
        sampleCurveForIntersections(curve), siblings, constrained_only, endpoint_tol);
}

static std::unordered_set<ConnId> crossingIdsTouchingSeeds(
    const std::vector<ConnectivityCurve>& results,
    const ClusterOrderSolver& cs, const std::unordered_set<ConnId>& seeds,
    double endpoint_tol = 1.5) {
    std::unordered_set<ConnId> bad;
    if (seeds.empty())
        return bad;
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, SampledCurve> samples;
    samples.reserve(results.size());
    for (const auto& cc : results)
        if (cc.curve)
            samples.emplace(cc.id, sampleCurveForIntersections(*cc.curve));
    for (auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        if (!seeds.count(p.id_a) && !seeds.count(p.id_b))
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        auto sa = samples.find(p.id_a);
        auto sb = samples.find(p.id_b);
        if (sa == samples.end() || sb == samples.end())
            continue;
        if (sampledCurvesIntersectBusiness(sa->second, sb->second, endpoint_tol)) {
            if (seeds.count(p.id_a))
                bad.insert(p.id_a);
            if (seeds.count(p.id_b))
                bad.insert(p.id_b);
        }
    }
    return bad;
}

static void annotateClusterCrossings(
    std::vector<ConnectivityCurve>& results,
    const ClusterOrderSolver& cs, double endpoint_tol = 1.5) {
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, SampledCurve> samples;
    samples.reserve(results.size());
    for (const auto& cc : results)
        if (cc.curve)
            samples.emplace(cc.id, sampleCurveForIntersections(*cc.curve, 48));
    for (auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        auto& ca = results[ia->second];
        auto& cb = results[ib->second];
        if (!ca.curve || !cb.curve)
            continue;
        auto sa = samples.find(p.id_a);
        auto sb = samples.find(p.id_b);
        if (sa == samples.end() || sb == samples.end())
            continue;
        Vec2d pt;
        bool interior = sampledCurvesIntersectBusiness(sa->second, sb->second, endpoint_tol, &pt);
        if (!interior)
            continue;
        ca.violation.exempt_crosses.push_back(pt);
        cb.violation.exempt_crosses.push_back(pt);
        CurveStatus st = (p.exempt == CrossExemption::ObstacleCross)
            ? CurveStatus::WarnA2 : CurveStatus::Degraded;
        if (ca.status == CurveStatus::OK || st == CurveStatus::Degraded)
            ca.status = st;
        if (cb.status == CurveStatus::OK || st == CurveStatus::Degraded)
            cb.status = st;
    }
}

static bool naturalCubicHitsObstacle(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const SDFField& sdf, const std::vector<Obstacle>& obstacles) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    BezierCurve trial;
    trial.segs.push_back(makeCubicG1(p0, T0, p1, T1, 0.4));
    if (curveIntersectsObstacles(trial, obstacles))
        return true;
    return sdf.valid() && minSDFAlongCurveAdaptive(trial, sdf) < 0.0;
}

static bool curveIntersectsBoundaries(
    const BezierCurve& curve, const std::vector<Boundary>& boundaries) {
    if (boundaries.empty())
        return false;
    auto sampled = sampleCurveForIntersections(curve, 64);
    for (const auto& bnd : boundaries) {
        const auto& pts = bnd.geometry.points;
        if (pts.size() < 2)
            continue;
        BoundingBox2d bnd_box;
        for (const auto& pt : pts)
            bnd_box.expand(pt);
        if (!sampled.bbox.intersects(bnd_box))
            continue;
        for (int i = 0; i + 1 < (int)sampled.pts.size(); ++i) {
            for (int j = 0; j + 1 < (int)pts.size(); ++j) {
                Vec2d isect;
                if (!segmentsIntersect(sampled.pts[i], sampled.pts[i + 1], pts[j], pts[j + 1], &isect))
                    continue;
                if ((isect - sampled.start).norm() <= 0.75 ||
                    (isect - sampled.end).norm() <= 0.75)
                    continue;
                if ((isect - pts[j]).norm() <= 0.10 ||
                    (isect - pts[j + 1]).norm() <= 0.10)
                    continue;
                return true;
            }
        }
    }
    return false;
}

static bool curveLeavesFence(const BezierCurve& curve, const Polygon2d& fence) {
    if (fence.outer.empty())
        return false;
    auto pts = curve.sample(24);
    for (int i = 1; i + 1 < (int)pts.size(); ++i) {
        const auto& pt = pts[i];
        if (!polygonContains(fence, pt) && pointToPolygonDist(pt, fence) > 0.10)
            return true;
    }
    return false;
}

struct CurveRisk {
    bool obstacle = false;
    bool boundary = false;
    bool fence = false;
    int sibling_crosses = 0;

    bool physical() const {
        return obstacle || boundary || fence;
    }
};

static void setConnectivityCurveGeometry(
    ConnectivityCurve& cc, const BezierCurve& curve, const LineString2d* fixed_geometry = nullptr) {
    cc.curve = std::make_shared<BezierCurve>(curve);
    cc.geometry.points.clear();
    if (fixed_geometry && fixed_geometry->points.size() >= 2) {
        cc.geometry = *fixed_geometry;
        return;
    }
    if (!curve.empty()) {
        int n = std::max(2, std::min(240, (int)std::ceil(curve.arcLength() / 0.3) + 1));
        cc.geometry.points = curve.sampleByArcLength(n);
    }
}

static CurveRisk assessCurveRisk(
    const BezierCurve& curve, const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence = true) {
    CurveRisk risk;
    double ms = minSDFAlongCurveAdaptive(curve, sdf);
    risk.obstacle = curveIntersectsObstacles(curve, input.obstacles) || (ms < 0.0);
    risk.boundary = curveIntersectsBoundaries(curve, input.boundaries);
    risk.fence = include_fence && (!input.area.is_rough && curveLeavesFence(curve, input.area.geometry));
    risk.sibling_crosses = sampledSiblingCrossCount(curve, sampled_siblings, true, 1.5);
    return risk;
}

static bool tryPhysicalSafeSingleCubic(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, BezierCurve& curve) {
    CurveRisk current = assessCurveRisk(curve, input, sdf, sampled_siblings);
    if (!current.physical())
        return false;

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    bool require_topology_clear = (T0.dot(T1) < -0.5);
    BezierCurve best;
    bool have_best = false;
    int best_cross = std::numeric_limits<int>::max();
    double best_shape = std::numeric_limits<double>::max();

    for (double alpha : {0.50, 0.55, 0.60, 0.46, 0.42, 0.38, 0.34, 0.30, 0.70, 0.80}) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        if (curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings);
        if (risk.physical() || (require_topology_clear && risk.sibling_crosses > 0))
            continue;
        double shape_score = std::abs(alpha - 0.40);
        if (!have_best ||
            risk.sibling_crosses < best_cross ||
            (risk.sibling_crosses == best_cross && shape_score < best_shape)) {
            best = candidate;
            have_best = true;
            best_cross = risk.sibling_crosses;
            best_shape = shape_score;
        }
    }

    if (!have_best)
        return false;
    curve = best;
    return true;
}

static bool tryShapeSafeSingleCubic(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, bool include_fence,
    BezierCurve& curve) {
    double chord_len = (p1 - p0).norm();
    if (chord_len < 1e-6)
        return false;
    bool shape_bad = curve.maxCurvature(20) > 2.0 ||
        curve.arcLength() > std::max(chord_len * 1.8, chord_len + 12.0);
    if (!shape_bad)
        return false;

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : (p1 - p0).normalized();
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : (p1 - p0).normalized();
    BezierCurve best = curve;
    CurveRisk best_risk = assessCurveRisk(best, input, sdf, sampled_siblings, include_fence);
    double best_score = 1000.0 * best_risk.sibling_crosses + best.maxCurvature(20) + 0.02 * best.arcLength();
    bool improved = false;

    for (double alpha : {0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.12, 0.08}) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        if (curveSelfIntersectsBusiness(candidate, 1.0))
            continue;
        CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings, include_fence);
        if (risk.physical())
            continue;
        double score = 1000.0 * risk.sibling_crosses + candidate.maxCurvature(20) + 0.02 * candidate.arcLength();
        if (score + 1e-6 < best_score) {
            best = candidate;
            best_score = score;
            improved = true;
            if (risk.sibling_crosses == 0 && candidate.maxCurvature(20) < 1.0)
                break;
        }
    }

    if (!improved)
        return false;
    curve = best;
    return true;
}

static const Connectivity* findConnectivityById(
    const std::vector<Connectivity>& conns, const ConnId& id) {
    for (const auto& conn : conns)
        if (conn.id == id)
            return &conn;
    return nullptr;
}

static std::unordered_map<ConnId, std::vector<ConnId>> constrainedNeighborMap(
    const std::vector<ConnectivityCurve>& results, const ClusterOrderSolver& cs) {
    auto idx = resultIndexById(results);
    std::unordered_map<ConnId, std::vector<ConnId>> neighbors;
    for (const auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        if (!results[ia->second].curve || !results[ib->second].curve)
            continue;
        neighbors[p.id_a].push_back(p.id_b);
        neighbors[p.id_b].push_back(p.id_a);
    }
    return neighbors;
}

static int constrainedCrossCountForId(
    const ConnId& id, const BezierCurve& curve,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    double endpoint_tol = 1.5) {
    int count = 0;
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return 0;
    for (const auto& other : nit->second) {
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (curvesIntersectBusiness(curve, *other_cc.curve, endpoint_tol))
            ++count;
    }
    return count;
}

static Vec2d dominantConstrainedRefPerpForId(
    const ConnId& id, const BezierCurve& curve,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    const ClusterOrderSolver& cs,
    double endpoint_tol = 1.5) {
    auto nit = neighbors.find(id);
    if (nit == neighbors.end())
        return Vec2d(0, 0);
    for (const auto& other : nit->second) {
        auto it = result_idx.find(other);
        if (it == result_idx.end())
            continue;
        const auto& other_cc = results[it->second];
        if (!other_cc.curve)
            continue;
        if (!curvesIntersectBusiness(curve, *other_cc.curve, endpoint_tol))
            continue;
        Vec2d ref = cs.refPerpOf(id, other);
        if (ref.norm() > 1e-8)
            return ref.normalized();
    }
    return Vec2d(0, 0);
}

static std::unordered_set<ConnId> collectPureGeometryRepairIds(
    const std::vector<ConnectivityCurve>& results, const ClusterOrderSolver& cs,
    const std::unordered_set<ConnId>& preserved_fixed_ids,
    double endpoint_tol = 1.5) {
    std::unordered_set<ConnId> repair_ids;
    auto idx = resultIndexById(results);
    for (const auto& p : cs.pairs()) {
        if (p.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = idx.find(p.id_a), ib = idx.find(p.id_b);
        if (ia == idx.end() || ib == idx.end())
            continue;
        const auto& ca = results[ia->second];
        const auto& cb = results[ib->second];
        if (!ca.curve || !cb.curve)
            continue;
        if (!curvesIntersectBusiness(*ca.curve, *cb.curve, endpoint_tol))
            continue;
        ConnId repair_id = ia->second > ib->second ? p.id_a : p.id_b;
        if (!preserved_fixed_ids.count(repair_id))
            repair_ids.insert(repair_id);
    }
    return repair_ids;
}

static bool tryPureGeometryTopologyRepair(
    const Connectivity& conn, const IntersectionInput& input,
    const std::vector<ConnectivityCurve>& results,
    const std::unordered_map<ConnId, size_t>& result_idx,
    const std::unordered_map<ConnId, std::vector<ConnId>>& neighbors,
    const ClusterOrderSolver& cs,
    BezierCurve& out_curve) {
    auto it = result_idx.find(conn.id);
    if (it == result_idx.end() || !results[it->second].curve)
        return false;

    auto _entry = input.entryPtDir(conn.entry_lane_id);
    Vec2d p0 = _entry.first;
    Vec2d t0 = _entry.second;
    auto _exit = input.exitPtDir(conn.exit_lane_id);
    Vec2d p1 = _exit.first;
    Vec2d t1 = _exit.second;
    Vec2d chord = p1 - p0;
    if (chord.norm() < 1e-6)
        return false;
    const BezierCurve& current = *results[it->second].curve;
    int best_cross = constrainedCrossCountForId(
        conn.id, current, results, result_idx, neighbors, 1.5);
    if (best_cross <= 0)
        return false;

    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : chord.normalized();
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : chord.normalized();
    BezierCurve best = current;
    double best_shape = current.maxCurvature(20) + 0.02 * current.arcLength();
    bool improved = false;

    auto consider = [&](BezierCurve candidate) {
        if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
            return;

        int cross_count = constrainedCrossCountForId(
            conn.id, candidate, results, result_idx, neighbors, 1.5);
        double shape = candidate.maxCurvature(20) + 0.02 * candidate.arcLength();
        if (cross_count < best_cross ||
            (cross_count == best_cross && shape + 1e-6 < best_shape)) {
            best = std::move(candidate);
            best_cross = cross_count;
            best_shape = shape;
            improved = true;
        }
    };

    for (double alpha : {0.10, 0.12, 0.14, 0.16, 0.18, 0.22, 0.26, 0.30, 0.34}) {
        BezierCurve candidate;
        candidate.segs.push_back(makeCubicG1(p0, T0, p1, T1, alpha));
        consider(candidate);
        if (best_cross == 0)
            break;
    }

    Vec2d ref_perp = dominantConstrainedRefPerpForId(
        conn.id, current, results, result_idx, neighbors, cs, 1.5);
    if (best_cross > 0 && ref_perp.norm() > 1e-8) {
        double chord_len = chord.norm();
        for (double alpha : {0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42}) {
            double handle = std::max(0.05, alpha) * chord_len;
            for (double side : {-1.0, 1.0}) {
                for (double offset : {0.75, 1.25, 1.8, 2.5, 3.4, 4.5, 6.0}) {
                    BezierSegment seg;
                    seg.ctrl[0] = p0;
                    seg.ctrl[1] = p0 + T0 * handle + ref_perp * (side * offset);
                    seg.ctrl[2] = p1 - T1 * handle + ref_perp * (side * offset);
                    seg.ctrl[3] = p1;
                    BezierCurve candidate;
                    candidate.segs.push_back(seg);
                    consider(candidate);
                    if (best_cross == 0)
                        break;
                }
                if (best_cross == 0)
                    break;
            }
            if (best_cross == 0)
                break;
        }
    }

    if (!improved)
        return false;
    out_curve = best;
    return true;
}

static bool chordIntersectsObstacle(
    const Vec2d& p0, const Vec2d& p1, const Obstacle& obs) {
    const Polygon2d& poly = obstacleHardGeometry(obs);
    if (poly.outer.size() < 3)
        return false;
    if (pointInsideOrOnPolygon(0.5 * (p0 + p1), poly))
        return true;
    return polygonSegmentIntersects(p0, p1, poly);
}

static BezierCurve buildWaypointCurve(
    const std::vector<Vec2d>& pts, const Vec2d& start_tan, const Vec2d& end_tan) {
    if (pts.size() < 2)
        return {};
    std::vector<Vec2d> tans(pts.size(), Vec2d(1, 0));
    tans.front() = start_tan.norm() > 1e-8 ? start_tan.normalized() : (pts[1] - pts[0]).normalized();
    tans.back() = end_tan.norm() > 1e-8 ? end_tan.normalized() : (pts.back() - pts[pts.size() - 2]).normalized();
    for (int i = 1; i + 1 < (int)pts.size(); ++i) {
        Vec2d d = pts[i + 1] - pts[i - 1];
        tans[i] = d.norm() > 1e-8 ? d.normalized() : tans[i - 1];
    }
    return makeCurveFromKnots(pts, tans, 0.34);
}

static bool tryObstacleBypassCandidate(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings, BezierCurve& curve) {
    if (input.obstacles.empty())
        return false;

    Vec2d along = p1 - p0;
    double len = along.norm();
    if (len < 1e-6)
        return false;
    along /= len;
    Vec2d perp{-along[1], along[0]};
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    bool require_topology_clear = (T0.dot(T1) < -0.5);

    BezierCurve best;
    bool have_best = false;
    int best_cross = std::numeric_limits<int>::max();
    double best_score = std::numeric_limits<double>::max();
    Vec2d mid = 0.5 * (p0 + p1);
    Vec2d center = input.area.geometry.outer.empty() ? Vec2d(0, 0) : Vec2d(0, 0);
    if (!input.area.geometry.outer.empty()) {
        center.setZero();
        int cnt = 0;
        for (const auto& p : input.area.geometry.outer) {
            center += p;
            ++cnt;
        }
        if (cnt > 0)
            center /= cnt;
    }
    double center_lat = (center - mid).dot(perp);

    for (const auto& obs : input.obstacles) {
        const Polygon2d& poly = obstacleHardGeometry(obs);
        if (poly.outer.size() < 3)
            continue;
        if (!chordIntersectsObstacle(p0, p1, obs) &&
            !curveIntersectsObstacles(curve, std::vector<Obstacle>{obs}))
            continue;

        double lon_min = 1e18, lon_max = -1e18;
        double lat_min = 1e18, lat_max = -1e18;
        for (const auto& p : poly.outer) {
            double lon = (p - p0).dot(along);
            double lat = (p - p0).dot(perp);
            lon_min = std::min(lon_min, lon);
            lon_max = std::max(lon_max, lon);
            lat_min = std::min(lat_min, lat);
            lat_max = std::max(lat_max, lat);
        }
        double lon_mid = std::max(len * 0.15, std::min(len * 0.85, 0.5 * (lon_min + lon_max)));
        double obs_lat_mid = 0.5 * (lat_min + lat_max);
        std::vector<int> sides;
        if (std::abs(center_lat) > 0.25 && std::abs(obs_lat_mid) > 0.25 &&
            ((center_lat > 0) != (obs_lat_mid > 0))) {
            sides = {(center_lat > 0) ? +1 : -1, (center_lat > 0) ? -1 : +1};
        } else {
            sides = {+1, -1};
        }

        for (int side : sides) {
            double edge_lat = side > 0 ? lat_max : lat_min;
            for (double clearance : {0.65, 1.0, 1.4, 1.9, 2.5, 3.2}) {
                double apex_lat = edge_lat + side * clearance;
                Vec2d apex = p0 + lon_mid * along + apex_lat * perp;
                double spread = std::max(2.0, std::min(len * 0.28, (lon_max - lon_min) * 0.8 + 1.5));
                double lon_a = std::max(len * 0.08, lon_mid - spread);
                double lon_b = std::min(len * 0.92, lon_mid + spread);
                Vec2d q0 = p0 + lon_a * along + apex_lat * 0.55 * perp;
                Vec2d q1 = p0 + lon_b * along + apex_lat * 0.55 * perp;
                std::vector<Vec2d> pts;
                if ((q0 - p0).norm() > 1.0 && (q1 - p1).norm() > 1.0 && lon_b > lon_a + 0.5)
                    pts = {p0, q0, apex, q1, p1};
                else
                    pts = {p0, apex, p1};

                BezierCurve candidate = buildWaypointCurve(pts, t0, t1);
                if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings);
                if (risk.physical() || (require_topology_clear && risk.sibling_crosses > 0))
                    continue;
                double shape_score = clearance + 0.02 * candidate.arcLength() +
                    (side * center_lat >= 0.0 ? 0.0 : 0.4);
                if (!have_best ||
                    risk.sibling_crosses < best_cross ||
                    (risk.sibling_crosses == best_cross && shape_score < best_score)) {
                    best = candidate;
                    have_best = true;
                    best_cross = risk.sibling_crosses;
                    best_score = shape_score;
                }
            }
        }
    }

    if (!have_best)
        return false;
    curve = best;
    return true;
}

static BezierCurve buildAlignedUTurn(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const Vec2d& offset_dir, double offset_m, double forward_scale = 1.0) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    double bias = 0.0;
    if (offset_m > 0.0 && offset_dir.norm() > 1e-8)
        bias = offset_m * offset_dir.normalized().dot(axis);

    BezierCurve c;
    c.segs.push_back(makeAlignedUTurnCubic(p0, T0, p1, T1, forward_scale, bias));
    return c;
}

static bool curveExceedsUTurnEnvelope(
    const BezierCurve& curve, const BezierCurve& reference,
    const Vec2d& p0, const Vec2d& p1) {
    double ref_len = std::max(reference.arcLength(), (p1 - p0).norm());
    double max_len = std::max(80.0, ref_len * 4.0);
    if (curve.arcLength() > max_len)
        return true;

    Vec2d mid = 0.5 * (p0 + p1);
    double max_radius = std::max(60.0, ref_len * 3.0);
    for (const auto& pt : curve.sampleByArcLength(80)) {
        if ((pt - mid).norm() > max_radius)
            return true;
    }
    return false;
}

static double curveLateralBulge(
    const BezierCurve& curve, const Vec2d& p0, const Vec2d& p1) {
    Vec2d chord = p1 - p0;
    if (chord.norm() < 1e-8)
        return 0.0;
    Vec2d perp{-chord.y(), chord.x()};
    perp.normalize();
    double best = 0.0;
    for (const auto& pt : curve.sampleByArcLength(80)) {
        double v = (pt - p0).dot(perp);
        if (std::abs(v) > std::abs(best))
            best = v;
    }
    return best;
}

static bool curveCollapsesUTurnEnvelope(
    const BezierCurve& curve, const BezierCurve& reference,
    const Vec2d& p0, const Vec2d& p1) {
    double ref_len = reference.arcLength();
    double chord_len = (p1 - p0).norm();
    double ref_bulge = std::abs(curveLateralBulge(reference, p0, p1));
    double cur_bulge = std::abs(curveLateralBulge(curve, p0, p1));

    if (ref_len > 1.0 && curve.arcLength() < std::max(chord_len * 1.08, ref_len * 0.72))
        return true;
    if (ref_bulge > 2.0 && cur_bulge < ref_bulge * 0.50)
        return true;
    if (curve.maxCurvature(20) > std::max(1.0, reference.maxCurvature(20) * 8.0))
        return true;
    return false;
}

static bool tryBoundedUTurnCandidate(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const IntersectionInput& input, const SDFField& sdf,
    const std::vector<SampledSiblingCurve>& sampled_siblings,
    const BezierCurve& reference, BezierCurve& out_curve) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;
    Vec2d axis = T0 - T1;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();

    std::vector<Vec2d> dirs;
    auto add_dir = [&](const Vec2d& raw) {
        if (raw.norm() < 1e-8)
            return;
        Vec2d d = raw.normalized();
        for (const auto& existing : dirs)
            if (existing.dot(d) > 0.98)
                return;
        dirs.push_back(d);
    };
    add_dir(axis);
    add_dir(-axis);
    add_dir(Vec2d(-T0.y(), T0.x()));
    add_dir(Vec2d(T0.y(), -T0.x()));

    bool have_best = false;
    BezierCurve best;
    int best_cross = std::numeric_limits<int>::max();
    double best_score = std::numeric_limits<double>::max();
    double ref_len = reference.arcLength();
    for (double scale : {0.85, 0.70, 0.55, 0.40, 1.0, 0.30}) {
        for (const auto& dir : dirs) {
            for (double mag : {0.0, 1.0, 2.0, 3.0, 4.5, 6.0}) {
                BezierCurve candidate = buildAlignedUTurn(p0, t0, p1, t1, dir, mag, scale);
                if (candidate.empty() || curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                if (curveExceedsUTurnEnvelope(candidate, reference, p0, p1))
                    continue;
                if (curveCollapsesUTurnEnvelope(candidate, reference, p0, p1))
                    continue;
                CurveRisk risk = assessCurveRisk(candidate, input, sdf, sampled_siblings);
                if (risk.physical())
                    continue;
                double score = std::abs(candidate.arcLength() - ref_len)
                    + 3.0 * std::abs(scale - 0.85)
                    + 0.2 * mag;
                if (!have_best ||
                    risk.sibling_crosses < best_cross ||
                    (risk.sibling_crosses == best_cross && score < best_score)) {
                    have_best = true;
                    best = candidate;
                    best_cross = risk.sibling_crosses;
                    best_score = score;
                }
            }
        }
    }

    if (!have_best)
        return false;
    out_curve = best;
    return true;
}

ConnectivityCurve ConnectivityGenerator::generateOne(
    const Connectivity& conn, const IntersectionInput& input, const SDFField& sdf,
    const SDFField& sdf_coarse, const std::vector<SiblingCurve>& siblings,
    bool allow_uturn_search, bool* out_physical_risk) {
    if (out_physical_risk)
        *out_physical_risk = false;
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
    bool is_uturn_geom = (t0.norm() > 1e-8 && t1.norm() > 1e-8 &&
                          angleBetween(t0, t1) > M_PI * 0.85);

    double lw = 3.5;
    if (auto* l = input.findLane(conn.entry_lane_id))
        lw = l->width;
    bool enforce_fence = !(input.obstacles.empty() && input.boundaries.empty());

    // ── Build initial curve with sibling topology awareness ──────────────────
    //
    // Step A: Collect sibling polys for obstacle-bypass side selection.
    //   Include ALL non-exempt siblings (same entry OR same exit group).
    //   This prevents choosing a bypass direction that conflicts with existing curves.
    BezierCurve initial;
    if (is_uturn_geom) {
        Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
        Vec2d offset_dir{-T0[1], T0[0]};
        initial = buildAlignedUTurn(p0, t0, p1, t1, offset_dir, 0.0, 1.0);
    } else {
        std::vector<std::vector<Vec2d>> sib_polys_for_init;
        if (!siblings.empty() &&
            !input.obstacles.empty() &&
            naturalCubicHitsObstacle(p0, t0, p1, t1, sdf, input.obstacles)) {
            sib_polys_for_init.reserve(siblings.size());
            for (auto& sib : siblings) {
                if (!sib.exempt_a1)
                    sib_polys_for_init.push_back(sib.curve.sampleByArcLength(12));
            }
        }
        initial = buildInitialCurve(
            p0, t0, p1, t1, sdf, input.area.geometry, sib_polys_for_init);
    }
    std::vector<SampledSiblingCurve> sampled_siblings;
    auto ensure_sampled_siblings = [&]() -> const std::vector<SampledSiblingCurve>& {
        if (sampled_siblings.empty() && !siblings.empty())
            sampled_siblings = sampleSiblingsForIntersections(siblings);
        return sampled_siblings;
    };

    if (!is_uturn_geom) {
        int current_cross = sampledSiblingCrossCount(initial, ensure_sampled_siblings(), true, 1.5);
        double turn_strength = 0.0;
        bool has_shared_endpoint_sibling = false;
        for (const auto& sib : siblings) {
            if (!sib.exempt_a1 && sib.shared_endpoint) {
                has_shared_endpoint_sibling = true;
                break;
            }
        }
        Vec2d chord = p1 - p0;
        if (chord.norm() > 1e-8 && t0.norm() > 1e-8)
            turn_strength = std::abs(cross2d(t0.normalized(), chord.normalized()));
        if (turn_strength > 0.35 && !has_shared_endpoint_sibling) {
            for (double alpha : {0.50, 0.46, 0.42}) {
                BezierCurve candidate;
                candidate.segs.push_back(makeCubicG1(
                    p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                    p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                    alpha));
                if (curveSelfIntersectsBusiness(candidate, 1.0))
                    continue;
                if (assessCurveRisk(candidate, input, sdf, ensure_sampled_siblings(), enforce_fence).physical())
                    continue;
                int cross_count = sampledSiblingCrossCount(candidate, ensure_sampled_siblings(), true, 1.5);
                if (cross_count < current_cross ||
                    (current_cross == 0 && cross_count == 0)) {
                    initial = candidate;
                    current_cross = cross_count;
                    break;
                }
            }
        }
    }

    // ── Adaptive alpha: if initial curve crosses siblings, try shorter handles ──
    // Root cause: buildInitialCurve uses alpha=0.4 and returns immediately when
    // obstacle-free, ignoring sibling crossings. For same-exit curves (e.g. a right
    // turn whose northward entry tangent arch briefly crosses a straight's path),
    // a shorter handle avoids the crossing from the start, letting the optimizer
    // refine from a correct topology rather than fighting out of a local minimum.
    if (!is_uturn_geom) {
        bool has_shared_endpoint_sibling = false;
        for (const auto& sib : siblings) {
            if (!sib.exempt_a1 && sib.shared_endpoint) {
                has_shared_endpoint_sibling = true;
                break;
            }
        }
        std::vector<double> try_alphas = {0.42, 0.38, 0.34, 0.30, 0.26, 0.22};
        int best_cross = sampledSiblingCrossCount(initial, ensure_sampled_siblings(), true, 1.5);
        BezierCurve best_curve = initial;
        double best_alpha = 0.4;
        for (double alpha : try_alphas) {
            if (best_cross == 0) break;
            // Rebuild initial with shorter handles (G1 maintained)
            BezierSegment shorter = makeCubicG1(
                p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                alpha);
            // Verify shorter arc is still obstacle-free
            bool arc_ok = true;
            if (sdf.valid()) {
                for (int k=1; k<16; ++k) {
                    auto kv = sdf.queryWithGrad(shorter.evaluate((double)k/16));
                    if (kv.first < 0) { arc_ok = false; break; }
                }
            }
            if (arc_ok) {
                BezierCurve c_short; c_short.segs.push_back(shorter);
                if (assessCurveRisk(c_short, input, sdf, ensure_sampled_siblings(), enforce_fence).physical())
                    continue;
                int cross_count = sampledSiblingCrossCount(c_short, ensure_sampled_siblings(), true, 1.5);
                if (cross_count < best_cross ||
                    (cross_count == best_cross && alpha > best_alpha)) {
                    best_cross = cross_count;
                    best_curve = c_short;
                    best_alpha = alpha;
                }
            }
        }
        initial = best_curve;
        if (best_cross > 0 && has_shared_endpoint_sibling) {
            for (double alpha : {0.18, 0.14, 0.10, 0.06}) {
                BezierCurve c_short;
                c_short.segs.push_back(makeCubicG1(
                    p0, t0.norm()>1e-8 ? t0.normalized() : Vec2d(1,0),
                    p1, t1.norm()>1e-8 ? t1.normalized() : Vec2d(1,0),
                    alpha));
                int cross_count = sampledSiblingCrossCount(c_short, ensure_sampled_siblings(), true, 1.5);
                if (cross_count < best_cross) {
                    initial = c_short;
                    best_cross = cross_count;
                    if (best_cross == 0)
                        break;
                }
            }
        }
    } else if (allow_uturn_search &&
               sampledSiblingCrossCount(initial, ensure_sampled_siblings(), true, 1.5) > 0) {
        Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
        Vec2d forward = (t1.norm() > 1e-8) ? (T0 - t1.normalized()) : T0;
        if (forward.norm() < 1e-8)
            forward = T0;
        forward.normalize();
        std::vector<Vec2d> dirs;
        auto add_dir = [&](const Vec2d& raw) {
            if (raw.norm() < 1e-8)
                return;
            Vec2d d = raw.normalized();
            for (const auto& existing : dirs)
                if (existing.dot(d) > 0.98)
                    return;
            dirs.push_back(d);
        };
        add_dir(forward);
        add_dir(-forward);

        BezierCurve best = initial;
        int best_cross = sampledSiblingCrossCount(best, ensure_sampled_siblings(), false, 1.5);
        double best_shape = 0.0;
        for (const auto& dir_raw : dirs) {
            if (dir_raw.norm() < 1e-8)
                continue;
            Vec2d dir = dir_raw.normalized();
            for (double scale : {1.0, 0.85, 0.70, 0.55, 0.40}) {
                for (double mag : {0.0, 1.0, 2.0, 3.0}) {
                    BezierCurve candidate = buildAlignedUTurn(p0, t0, p1, t1, dir, mag, scale);
                    int cross_count = sampledSiblingCrossCount(candidate, ensure_sampled_siblings(), false, 1.5);
                    double shape_score = std::abs(scale - 1.0) + mag * 0.05;
                    if (cross_count < best_cross ||
                        (cross_count == best_cross && shape_score < best_shape)) {
                        best = candidate;
                        best_cross = cross_count;
                        best_shape = shape_score;
                    }
                }
            }
        }
        initial = best;
    }

    const auto& sampled_for_gate = ensure_sampled_siblings();
    CurveRisk risk = assessCurveRisk(initial, input, sdf, sampled_for_gate, enforce_fence);
    bool fixed_shape_cross = false;
    if (risk.sibling_crosses > 0) {
        auto current_sample = sampleCurveForIntersections(initial);
        for (const auto& sib : siblings) {
            if (!sib.fixed_shape || sib.exempt_a1)
                continue;
            auto fixed_sample = sampleCurveForIntersections(sib.curve);
            if (sampledCurvesIntersectBusiness(current_sample, fixed_sample, 1.5)) {
                fixed_shape_cross = true;
                break;
            }
        }
    }
    bool needs_optimization = risk.physical() || fixed_shape_cross ||
        (allow_uturn_search && risk.sibling_crosses > 0);
    if (out_physical_risk)
        *out_physical_risk = false;
    PreCheckResult pre;
    if (needs_optimization && sdf_coarse.valid() && !input.area.is_rough) {
        pre = preCheck(sdf_coarse, input.area.geometry, p0, p1, lw, input.boundaries);
        if (pre.type == ViolationInfo::InfeasibilityType::TopologicalBlock)
            return makeFallbackCurve(pre, conn, p0, p1);
    }

    if (!needs_optimization) {
        if (!is_uturn_geom)
            tryShapeSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, initial);
        setConnectivityCurveGeometry(cc, initial);
        validate(cc, input, sdf);
        return cc;
    }

    if (is_uturn_geom) {
        BezierCurve bounded_uturn = initial;
        if (tryBoundedUTurnCandidate(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, initial, bounded_uturn)) {
            setConnectivityCurveGeometry(cc, bounded_uturn);
            validate(cc, input, sdf);
            return cc;
        }
        setConnectivityCurveGeometry(cc, initial);
        validate(cc, input, sdf);
        if (out_physical_risk)
            *out_physical_risk = true;
        return cc;
    }

    BezierCurve safe_single = initial;
    if (!is_uturn_geom &&
        tryPhysicalSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, safe_single)) {
        setConnectivityCurveGeometry(cc, safe_single);
        validate(cc, input, sdf);
        return cc;
    }

    BezierCurve bypass_candidate = initial;
    if (tryObstacleBypassCandidate(
            p0, t0, p1, t1, input, sdf, sampled_for_gate, bypass_candidate)) {
        setConnectivityCurveGeometry(cc, bypass_candidate);
        validate(cc, input, sdf);
        return cc;
    }

    // ── Optimise ──────────────────────────────────────────────────────────────
    PenaltyCost cost;
    cost.proto = initial;
    cost.sdf = &sdf;
    cost.boundaries = input.boundaries;
    if (enforce_fence && !input.area.is_rough)
        cost.fence = input.area.geometry;
    cost.siblings = siblings;
    cost.obstacle_clearance = 0.0;
    cost.start_tan_dir = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    cost.end_tan_dir = t1.norm() > 1e-8 ? t1.normalized() : Vec2d(1, 0);
    cost.full_param_mode = (initial.numSegments() > 1);

    BezierCurve opt = optimiseCurve(cost, solver_, initial, /*outer_iters=*/3);

    bool skip_band = true; // always skip: preserve G1, rely on optimizer
    BezierCurve final_c = postProcess(opt, sdf, input.area.geometry, 0.25, t0, t1, skip_band, &p0, &p1);
    bool shape_risk = is_uturn_geom &&
        (curveExceedsUTurnEnvelope(final_c, initial, p0, p1) ||
         curveCollapsesUTurnEnvelope(final_c, initial, p0, p1));
    if (allow_uturn_search && shape_risk) {
        BezierCurve fallback = initial;
        if (tryBoundedUTurnCandidate(p0, t0, p1, t1, input, sdf, sampled_for_gate, initial, fallback)) {
            final_c = fallback;
        } else {
            final_c = initial;
        }
        shape_risk = false;
    }
    CurveRisk final_risk = assessCurveRisk(final_c, input, sdf, sampled_for_gate, enforce_fence);
    if (final_risk.physical()) {
        BezierCurve repaired = final_c;
        if (tryObstacleBypassCandidate(
                p0, t0, p1, t1, input, sdf, sampled_for_gate, repaired)) {
            final_c = repaired;
            final_risk = assessCurveRisk(final_c, input, sdf, sampled_for_gate, enforce_fence);
        }
    }
    if (!is_uturn_geom)
        tryShapeSafeSingleCubic(p0, t0, p1, t1, input, sdf, sampled_for_gate, enforce_fence, final_c);
    if (out_physical_risk) {
        *out_physical_risk = final_risk.physical() || shape_risk;
    }
    setConnectivityCurveGeometry(cc, final_c);
    validate(cc, input, sdf);
    if (pre.narrow_passage && cc.status == CurveStatus::OK)
        cc.status = CurveStatus::WarnA2;
    return cc;
}

std::vector<ConnectivityCurve> ConnectivityGenerator::generate(
    const IntersectionInput& raw_input, SDFField& sdf, double* out_ms) {
    auto t0 = std::chrono::steady_clock::now();
    IntersectionInput input = raw_input;
    applyGroupUnifiedDirections(input, direction_cfg_);

    bool pure_geometry = input.obstacles.empty()
        && input.boundaries.empty();

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
    if (!input.obstacles.empty())
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
    std::unordered_set<ConnId> physical_risk_ids;
    std::unordered_set<ConnId> preserved_fixed_ids;
    for (const auto& conn : input.connectivities) {
        if (!hasFixedGeometry(conn) || fixedGeometryHitsObstacle(conn, input, sdf))
            continue;
        auto cc = makeFixedGeometryCurve(conn, input, sdf);
        if (cc.curve) {
            done[conn.id] = *cc.curve;
            preserved_fixed_ids.insert(conn.id);
        }
        results.push_back(std::move(cc));
    }

    for (auto& group : coord.groups()) {
        for (auto& cid : group.conn_ids) {
            if (preserved_fixed_ids.count(cid))
                continue;
            auto* conn = cmap[cid];
            if (!conn) continue;

            auto sibs = buildSiblings(
                cid, done, cluster_solver_, input.connectivities, pure_geometry, &preserved_fixed_ids);

            bool phys_risk = false;
            auto cc = generateOne(*conn, input, sdf, sdf_coarse, sibs, false, &phys_risk);
            if (cc.curve) {
                done[cid] = *cc.curve;
                if (phys_risk)
                    physical_risk_ids.insert(cid);
            }

            results.push_back(std::move(cc));
        }
        // After each priority group: mark obstacle-adjacent crossings as soft
        cluster_solver_.checkAndMarkA2(done, sdf, 1.5);
    }

    if (pure_geometry) {
        auto result_idx = resultIndexById(results);
        auto neighbors = constrainedNeighborMap(results, cluster_solver_);
        for (int pass = 0; pass < 2; ++pass) {
            auto repair_ids = collectPureGeometryRepairIds(
                results, cluster_solver_, preserved_fixed_ids, 1.5);
            if (repair_ids.empty())
                break;
            bool changed = false;
            for (const auto& cid : repair_ids) {
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end())
                    continue;
                const Connectivity* conn = findConnectivityById(input.connectivities, cid);
                if (!conn)
                    continue;
                BezierCurve repaired;
                if (!tryPureGeometryTopologyRepair(
                        *conn, input, results, result_idx, neighbors, cluster_solver_, repaired))
                    continue;
                setConnectivityCurveGeometry(results[ri->second], repaired);
                validate(results[ri->second], input, sdf);
                changed = true;
            }
            if (!changed)
                break;
        }
        annotateClusterCrossings(results, cluster_solver_, 1.5);
        if (out_ms) {
            *out_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0).count();
        }
        return results;
    }

    auto result_idx = resultIndexById(results);
    std::unordered_map<ConnId, BezierCurve> all_done = curveMapFromResults(results);
    for (int repair_pass = 0; repair_pass < 3; ++repair_pass) {
        auto bad = crossingIdsTouchingSeeds(results, cluster_solver_, physical_risk_ids, 1.5);
        if (bad.empty())
            break;
        int repaired_count = 0;
        int max_repairs = std::min(12, (int)input.connectivities.size());
        for (auto git = coord.groups().rbegin(); git != coord.groups().rend(); ++git) {
            for (auto& cid : git->conn_ids) {
                if (!bad.count(cid))
                    continue;
                if (preserved_fixed_ids.count(cid))
                    continue;
                if (repaired_count >= max_repairs)
                    break;
                auto* conn = cmap[cid];
                if (!conn)
                    continue;
                auto sibs = buildSiblings(cid, all_done, cluster_solver_, input.connectivities, false, &preserved_fixed_ids);
                bool phys_risk = false;
                auto cc = generateOne(*conn, input, sdf, sdf_coarse, sibs, true, &phys_risk);
                auto ri = result_idx.find(cid);
                if (ri == result_idx.end())
                    continue;
                results[ri->second] = std::move(cc);
                if (results[ri->second].curve)
                    all_done[cid] = *results[ri->second].curve;
                if (phys_risk)
                    physical_risk_ids.insert(cid);
                else
                    physical_risk_ids.erase(cid);
                ++repaired_count;
            }
            if (repaired_count >= max_repairs)
                break;
        }
        cluster_solver_.checkAndMarkA2(all_done, sdf, 1.5);
    }
    annotateClusterCrossings(results, cluster_solver_, 1.5);

    auto t1 = std::chrono::steady_clock::now();
    if (out_ms) *out_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return results;
}

}
