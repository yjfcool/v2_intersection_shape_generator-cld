#include "hermite_init.h"
#include "optimizer/sdf_field.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <cmath>
#include <algorithm>

namespace isg {

// ─────────────────────────────────────────────────────────────────────────────
//  Shared helpers
// ─────────────────────────────────────────────────────────────────────────────
static bool segmentsIntersect_internal(
    const Vec2d& a, const Vec2d& b, const Vec2d& c, const Vec2d& d, Vec2d* out = nullptr) {
    Vec2d r = b - a, s = d - c;
    double den = cross2d(r, s);
    if (std::abs(den) < 1e-12) return false;
    Vec2d ac = c - a;
    double t = cross2d(ac, s) / den;
    double u = cross2d(ac, r) / den;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        if (out) *out = a + t * r;
        return true;
    }
    return false;
}

static constexpr double MAX_TAN_DEV = 60.0 * M_PI / 180.0;

static Vec2d clampTangent(const Vec2d& tan, const Vec2d& ref, double max_a) {
    Vec2d t = tan.norm() > 1e-10 ? tan.normalized() : ref;
    if (t.dot(ref) >= std::cos(max_a)) return t;
    double sg = cross2d(ref, t) >= 0 ? 1.0 : -1.0;
    double ca = std::cos(max_a), sa = std::sin(max_a);
    return Vec2d(ref[0] * ca - sg * ref[1] * sa,
                 ref[0] * sa * sg + ref[1] * ca);
}

// SDF-sampled straight line clearance check.
// clearance is the minimum acceptable SDF value along the line.
// Use clearance=0.0 to only detect actual physical penetration (SDF<0),
// which avoids false-triggering bypass for lanes that merely pass close
// to an obstacle but do not physically cross it.
static bool straightLineClear(
    const SDFField& sdf, const Vec2d& p0, const Vec2d& p1, double clearance = 0.0, int n = 40) {
    if (!sdf.valid()) return true;
    for (int i = 1; i < n; ++i) {
        double t = (double)i / n;
        std::pair<double, Vec2d> _q = sdf.queryWithGrad((1 - t) * p0 + t * p1);
        if (_q.first < clearance) return false;
    }
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Level-1: Geometric Direct Construction (方案A)
// ─────────────────────────────────────────────────────────────────────────────
// Union-AABB of all obstacle buffered geometries visible to the SDF field.
// We probe the SDF along the direct path to find the obstacle's bounding box
// without needing obstacle polygon data (SDF is already available).
struct ObstacleAABB {
    double x_min, x_max, y_min, y_max;
    bool valid = false;
};

static ObstacleAABB probeObstacleAABB(
    const SDFField& sdf, const Vec2d& p0, const Vec2d& p1, double clearance = 0.0) {
    // Sample dense grid in the neighbourhood of the direct path
    ObstacleAABB box;
    box.x_min = box.y_min = 1e18;
    box.x_max = box.y_max = -1e18;

    int nx = 60, ny = 30;
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) return box;
    along = along * (1.0 / len);
    Vec2d perp{-along[1], along[0]};
    // Sweep width: use a fixed narrow strip (one lane width = 3.5m each side)
    // to avoid capturing obstacles far off the path axis.  A large sweep
    // caused the AABB to span the entire junction when multiple obstacles exist,
    // making obs_lat unreliable for side selection.
    double sweep = 4.0; // ±4m from path centre-line

    for (int ix = 0; ix <= nx; ++ix) {
        for (int iy = -ny / 2; iy <= ny / 2; ++iy) {
            double s = (double)ix / nx * len;
            double t = (double)iy / (ny / 2) * sweep;
            Vec2d pt = p0 + s * along + t * perp;
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(pt);
            if (_q.first < clearance) {
                box.x_min = std::min(box.x_min, pt.x());
                box.x_max = std::max(box.x_max, pt.x());
                box.y_min = std::min(box.y_min, pt.y());
                box.y_max = std::max(box.y_max, pt.y());
                box.valid = true;
            }
        }
    }
    return box;
}

// Compute the bypass apex point for a given side (+1 = left of path, -1 = right).
// Apex is placed at the obstacle's maximum lateral projection + clearance,
// at the longitudinal midpoint of the obstacle along the path.
static Vec2d computeApex(
    const ObstacleAABB& box, const Vec2d& p0, const Vec2d& p1, int side, double clearance) {
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) return 0.5 * (p0 + p1);
    along = along * (1.0 / len);
    Vec2d perp{-along[1], along[0]}; // left normal (unit vector)

    // Longitudinal midpoint of obstacle box (world coords → project onto path axis)
    // The box stores world coordinates of sampled points, so we must project
    // all four corners onto the path tangent and perpendicular axes.
    // Corners of the axis-aligned world box:
    double corners_x[4] = {box.x_min, box.x_max, box.x_max, box.x_min};
    double corners_y[4] = {box.y_min, box.y_min, box.y_max, box.y_max};

    double lon_min = 1e18, lon_max = -1e18;
    double lat_min = 1e18, lat_max = -1e18;
    for (int i = 0; i < 4; ++i) {
        Vec2d c(corners_x[i], corners_y[i]);
        double lon = c.dot(along);
        double lat = c.dot(perp);
        lon_min = std::min(lon_min, lon);
        lon_max = std::max(lon_max, lon);
        lat_min = std::min(lat_min, lat);
        lat_max = std::max(lat_max, lat);
    }

    // Longitudinal midpoint of obstacle, clamped to [20%, 80%] of path
    double path_lon_0 = p0.dot(along);
    double path_lon_1 = p1.dot(along);
    double lon_mid_obs = 0.5 * (lon_min + lon_max);
    double frac = (lon_mid_obs - path_lon_0) / (path_lon_1 - path_lon_0 + 1e-12);
    frac = std::max(0.2, std::min(0.8, frac));

    // Lateral edge of obstacle on chosen side + clearance
    double lat_edge = (side > 0) ? (lat_max + clearance) : (lat_min - clearance);

    // Apex: move from the path along perp to the required lateral position
    Vec2d base = p0 + frac * (p1 - p0);
    double base_lat = base.dot(perp);
    Vec2d apex = base + (lat_edge - base_lat) * perp;
    return apex;
}

// Choose which side (left/right of path) is clear for bypass.
//
// Right-hand traffic rules (RHT) and junction-centre preference:
//
//   Coordinate conventions:
//     along = (p1-p0).normalized()
//     perp  = (-along.y, along.x)   ← LEFT of the path direction
//     side = +1 → bypass to the LEFT  of along (= perp direction)
//     side = -1 → bypass to the RIGHT of along (= -perp direction)
//
//   Turn geometry (cross2d(t0, p1-p0)):
//     > 0  → p1-p0 is LEFT  of t0 → LEFT TURN  → obstacle on left-diagonal of path
//     < 0  → p1-p0 is RIGHT of t0 → RIGHT TURN → obstacle on right-diagonal of path
//
//   RHT avoidance principle: prefer the side closer to the junction CENTRE:
//     - RIGHT TURN: natural arc swings right; obstacle sits on the right-diagonal
//       (inner corner of the right turn, close to junction centre).
//       The correct bypass goes to the LEFT of the path direction (+1),
//       which is the PATH'S left side = geometrically the OUTER side of the
//       right-turn arc. Wait — this is wrong per user requirement.
//
//   CORRECTED reasoning (confirmed by geometry derivation):
//     For N→W right turn: p0=(+1.75,-10), p1=(+10,-1.75), t0=(0,+1)
//       along = (8.25, 8.25).norm ≈ (0.707, 0.707)
//       perp  = (-0.707, +0.707)   ← points toward junction OUTSIDE (NW corner)
//       side=+1 (left of path) = NW = OUTER = WRONG
//       side=-1 (right of path) = SE = junction CENTRE diagonal = CORRECT
//
//     For a right turn, cross2d(t0, p1-p0) < 0.
//     The inner (centre-side) bypass corresponds to side = -1 (right of path).
//
//   Left-turn example N→E: p0=(+1.75,-10), p1=(-10,+1.75), t0=(0,+1)
//       along = (-11.75, 11.75).norm ≈ (-0.707, 0.707)
//       perp  = (-0.707, -0.707)   ← points toward junction CENTRE (SW, inner of left turn)
//       side=+1 (left of path) = SW = junction centre = inner side = CORRECT
//       side=-1 (right of path) = NE = outer side
//     For a left turn, cross2d(t0, p1-p0) > 0.
//     The inner (centre-side) bypass corresponds to side = +1 (left of path).
//
//   RULE: side_inner = +1 when left-turn, -1 when right-turn.
//         Equivalently: side_inner = sign(cross2d(t0, p1-p0)) with fallback to
//         gradient-voting for straight / near-straight paths.
//
// t0_hint: optional entry tangent for turn-type detection. If zero-vector,
//          fall back to gradient voting only.
static int chooseSide(
    const SDFField& sdf, const ObstacleAABB& box, const Vec2d& p0, const Vec2d& p1, double clearance,
    const std::vector<std::vector<Vec2d>>& sibling_polys, const Vec2d& t0_hint = Vec2d(0, 0)) {
    Vec2d along = (p1 - p0).normalized();
    Vec2d perp{-along[1], along[0]}; // LEFT of path direction

    // ── Side selection: RHT centre-priority with obstacle position awareness ──
    //
    // Key principle: prefer the bypass side that keeps the curve close to the
    // junction centre (inner arc / short arc), UNLESS the obstacle is between
    // the path and the centre, in which case we must go outer.
    //
    // Algorithm:
    //   1. Measure lat position of the junction centre relative to the path.
    //   2. Measure lat position of the obstacle AABB centre relative to the path.
    //   3. If centre and obstacle are on OPPOSITE sides of the path:
    //        → obstacle is on the outer side; bypass TOWARD the centre (same side as centre).
    //   4. If centre and obstacle are on the SAME side of the path:
    //        → obstacle blocks the inner (centre) side; must bypass to the outer side.
    //
    // Fallback: if no obstacle detected in AABB or centre is ambiguous, use
    // SDF clearance voting.
    if (t0_hint.norm() > 1e-8 || box.valid) {
        Vec2d junction_centre(0.0, 0.0);
        Vec2d mid_path = 0.5 * (p0 + p1);
        Vec2d to_centre = junction_centre - mid_path;
        double centre_lat = to_centre.dot(perp); // +: centre left, -: centre right

        // Obstacle AABB centre
        double obs_cx = 0.5 * (box.x_min + box.x_max);
        double obs_cy = 0.5 * (box.y_min + box.y_max);
        Vec2d obs_centre_world(obs_cx, obs_cy);
        double obs_lat = (obs_centre_world - mid_path).dot(perp);

        const double LAT_THRESH = 0.5;
        bool centre_clear = std::abs(centre_lat) > LAT_THRESH;
        bool obs_clear = std::abs(obs_lat) > LAT_THRESH;

        if (centre_clear && obs_clear) {
            // Both clearly on defined sides
            bool same_side = (centre_lat > 0) == (obs_lat > 0);
            if (!same_side) {
                // Obstacle on outer side → bypass toward centre
                return (centre_lat > 0) ? +1 : -1;
            } else {
                // Obstacle on inner/centre side → forced outer bypass
                return (centre_lat > 0) ? -1 : +1;
            }
        }
        if (centre_clear && !obs_clear) {
            // Obstacle is nearly on-path (head-on) — prefer inner (centre) bypass
            // since the obstacle doesn't clearly block the inner lane.
            return (centre_lat > 0) ? +1 : -1;
        }
        if (centre_clear) {
            // Only centre direction is clear: default to inner bypass
            return (centre_lat > 0) ? +1 : -1;
        }
    }
    // The SDF gradient at points near the obstacle points AWAY from the obstacle
    // (toward free space).  We accumulate the lateral component of the gradient
    // to determine which side of the path has MORE clearance — that is the side
    // we choose to bypass TOWARD.
    // Note: clearance=0 in normal operation so "d > clearance*1.5" always holds;
    // use a fixed probe threshold instead.
    double weighted_lat = 0.0;
    double total_weight = 0.0;
    constexpr double PROBE = 0.15;
    constexpr double SAMPLE_THRESH = 3.0; // probe within 3m of obstacle
    constexpr int N = 30;
    for (int i = 1; i < N; ++i) {
        double t = (double)i / N;
        Vec2d pt = p0 + t * (p1 - p0);
        std::pair<double, Vec2d> _q = sdf.queryWithGrad(pt);
        double d = _q.first;
        if (d > SAMPLE_THRESH) continue;
        std::pair<double, Vec2d> _ql = sdf.queryWithGrad(pt + PROBE * perp);
        std::pair<double, Vec2d> _qm = sdf.queryWithGrad(pt - PROBE * perp);
        double d_lp = _ql.first;
        double d_lm = _qm.first;
        double lat_grad = (d_lp - d_lm) / (2 * PROBE);
        double w = 1.0 / (d + 0.05);
        weighted_lat += w * lat_grad;
        total_weight += w;
    }

    // side = direction of MORE clearance (toward free space)
    int primary_side = (total_weight > 1e-10 && weighted_lat > 0) ? +1 : -1;

    // Sibling penalty
    Vec2d mid = 0.5 * (p0 + p1);
    int sib_left = 0, sib_right = 0;
    for (auto& poly : sibling_polys) {
        for (auto& pt : poly) {
            double lat = (pt - mid).dot(perp);
            if (lat > 0.2) sib_left++;
            if (lat < -0.2) sib_right++;
        }
    }
    int occ_primary = (primary_side > 0) ? sib_left : sib_right;
    int occ_opposite = (primary_side > 0) ? sib_right : sib_left;
    if (occ_primary > occ_opposite + 5) primary_side = -primary_side;

    return primary_side;
}

// Build a smooth 2-segment arch: p0(t0) → apex(apex_tan) → p1(t1)
// apex_tan is perpendicular to the p0→p1 direction (or aligned with t0→t1 average)
static BezierCurve buildArch(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& apex,
    const Vec2d& p1, const Vec2d& t1, double alpha = 0.4) {
    // Tangent at apex: average of incoming (p0→apex) and outgoing (apex→p1),
    // clamped to be perpendicular to the net p0→p1 direction for a smooth arch.
    Vec2d leg0 = (apex - p0).norm() > 1e-8 ? (apex - p0).normalized() : t0.normalized();
    Vec2d leg1 = (p1 - apex).norm() > 1e-8 ? (p1 - apex).normalized() : t1.normalized();
    Vec2d apex_tan = (leg0 + leg1);
    if (apex_tan.norm() < 1e-8) apex_tan = leg0;
    apex_tan.normalize();

    BezierSegment s0 = makeCubicG1(p0, t0.normalized(), apex, apex_tan, alpha);
    BezierSegment s1 = makeCubicG1(apex, apex_tan, p1, t1.normalized(), alpha);

    BezierCurve c;
    c.segs.push_back(s0);
    c.segs.push_back(s1);
    return c;
}

// Level-1 entry point
static BezierCurve geometricBypass(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const SDFField& sdf, const Polygon2d& fence, double clearance = 0.0,
    const std::vector<std::vector<Vec2d>>& sibling_polys = {}) {
    // Use the same clearance for probing — detect only actual penetrations
    auto box = probeObstacleAABB(sdf, p0, p1, clearance);
    if (!box.valid) {
        // No obstacle detected within clearance → direct single-segment curve
        BezierCurve c;
        c.segs.push_back(makeCubicG1(p0, t0.normalized(), p1, t1.normalized(), 0.4));
        return c;
    }

    // Pass entry tangent so chooseSide can detect turn type (RHT centre-priority)
    int side = chooseSide(sdf, box, p0, p1, clearance, sibling_polys, t0);
    // apex_clearance: place apex beyond obstacle edge + small buffer
    double apex_clearance = clearance + 0.3;
    Vec2d apex = computeApex(box, p0, p1, side, apex_clearance);

    // Verify apex is obstacle-free; if not (apex lands inside obstacle), try other side.
    // Use a meaningful threshold regardless of clearance value.
    std::pair<double, Vec2d> _qa = sdf.queryWithGrad(apex);
    if (_qa.first < -0.1) {
        apex = computeApex(box, p0, p1, -side, apex_clearance);
    }

    // Verify arch clears obstacles (sample each segment)
    BezierCurve arch = buildArch(p0, t0, apex, p1, t1);
    bool arch_clear = true;
    for (auto& seg : arch.segs) {
        for (int i = 1; i < 20; ++i) {
            std::pair<double, Vec2d> _qs = sdf.queryWithGrad(seg.evaluate((double)i / 20));
            if (_qs.first < clearance * 0.5) {
                arch_clear = false;
                break;
            }
        }
        if (!arch_clear) break;
    }

    if (!arch_clear) {
        // Try opposite side
        BezierCurve arch2 = buildArch(p0, t0, computeApex(box, p0, p1, -side, apex_clearance), p1, t1);
        bool arch2_clear = true;
        for (auto& seg : arch2.segs) {
            for (int i = 1; i < 20; ++i) {
                std::pair<double, Vec2d> _qs2 = sdf.queryWithGrad(seg.evaluate((double)i / 20));
                if (_qs2.first < clearance * 0.5) {
                    arch2_clear = false;
                    break;
                }
            }
            if (!arch2_clear) break;
        }
        if (!arch2_clear) return {}; // signal Level-2 fallback
        arch = arch2;
    }

    // Validate no intersection with existing sibling curves.
    // Use a larger endpoint tolerance for same-entry siblings, since curves
    // sharing the same start point have overlapping initial segments by design.
    if (!sibling_polys.empty()) {
        auto arch_pts = arch.sampleByArcLength(20);
        double arch_len = arch.arcLength();
        // Minimum "interior" fraction: skip intersections within the first/last
        // 15% of arc length (where same-entry curves naturally overlap).
        double endpoint_skip_dist = std::max(0.05, arch_len * 0.15);
        for (auto& sp : sibling_polys) {
            for (int ai = 0; ai + 1 < (int)arch_pts.size(); ++ai) {
                for (int si = 0; si + 1 < (int)sp.size(); ++si) {
                    Vec2d isect;
                    if (!segmentsIntersect_internal(
                        arch_pts[ai], arch_pts[ai + 1],
                        sp[si], sp[si + 1], &isect))
                        continue;
                    // Skip intersections close to either endpoint of the arch
                    double de = std::min(
                        (isect - arch_pts.front()).norm(),
                        (isect - arch_pts.back()).norm());
                    if (de < endpoint_skip_dist) continue;
                    // Genuine interior intersection → Level-2 fallback
                    return {};
                }
            }
        }
    }

    return arch;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Level-2: Full control-point geometry initialisation (方案C 初始化)
//
//  Creates a plausible initial arch WITHOUT grid routing.
//  All control points (including the join point) are free for optimisation.
//  Strategy: tangent-extension intersection gives the arch apex analytically.
// ─────────────────────────────────────────────────────────────────────────────
static BezierCurve geometricInitLevel2(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1, const SDFField& sdf) {
    // Step 1: find a safe lateral direction (away from obstacle cloud)
    Vec2d along = (p1 - p0);
    double len = along.norm();
    if (len < 1e-6) {
        BezierCurve c;
        c.segs.push_back(makeCubicG1(p0, t0.normalized(), p1, t1.normalized(), 0.4));
        return c;
    }
    along = along * (1.0 / len);
    Vec2d perp{-along[1], along[0]};

    // RHT centre-priority: prefer bypassing toward the junction centre (0,0),
    // unless the obstacle is ON the centre side (then must go outer).
    // Use probeObstacleAABB to get actual obstacle world-coords, then determine
    // which side of the path the obstacle lies on.
    double side = 0.0;
    {
        Vec2d junction_centre(0.0, 0.0);
        Vec2d mid_path = 0.5 * (p0 + p1);
        double centre_lat = (junction_centre - mid_path).dot(perp);

        // Probe obstacle AABB (reuse Level-1 helper)
        auto box = probeObstacleAABB(sdf, p0, p1, 0.0);

        if (box.valid && std::abs(centre_lat) > 0.5) {
            // Compute obstacle AABB centre's lateral position relative to path
            Vec2d box_world_centre((box.x_min + box.x_max) * 0.5,
                                   (box.y_min + box.y_max) * 0.5);
            double obs_lat = (box_world_centre - mid_path).dot(perp);

            // Determine which side to bypass:
            //   If obs is on the OUTER side (opposite to centre): bypass toward CENTRE
            //   If obs is on the INNER side (same as centre, well separated): forced OUTER
            //   If obs is nearly on-path (|obs_lat| < 1.5m): prefer INNER (centre) side,
            //     since we have a head-on obstacle and the inner arc is shorter.
            const double ON_PATH_THRESH = 1.5;
            if (std::abs(obs_lat) < ON_PATH_THRESH) {
                // Head-on obstacle: prefer inner (centre) side
                side = (centre_lat > 0) ? 1.0 : -1.0;
            } else {
                bool same_side = (centre_lat > 0) == (obs_lat > 0);
                side = (!same_side)
                           ? ((centre_lat > 0) ? 1.0 : -1.0) // outer obstacle → inner bypass
                           : ((centre_lat > 0) ? -1.0 : 1.0); // inner obstacle → outer bypass
            }
        } else if (std::abs(centre_lat) > 0.5) {
            side = (centre_lat > 0) ? 1.0 : -1.0; // default inner bypass
        } else {
            // Near-straight: pick side with more SDF clearance
            double d_left = 0, d_right = 0;
            int n_probe = 10;
            for (int i = 1; i <= n_probe; ++i) {
                double s = (double)i / (n_probe + 1) * len;
                Vec2d mid_pt = p0 + s * along;
                std::pair<double, Vec2d> _ql2 = sdf.queryWithGrad(mid_pt + 1.5 * perp);
                std::pair<double, Vec2d> _qr2 = sdf.queryWithGrad(mid_pt - 1.5 * perp);
                d_left += _ql2.first;
                d_right += _qr2.first;
            }
            side = (d_left >= d_right) ? 1.0 : -1.0;
        }
    }
    // Safety check: if preferred side is inside an obstacle, flip to other side.
    {
        Vec2d mid_pt = p0 + 0.5 * (p1 - p0);
        std::pair<double, Vec2d> _qc = sdf.queryWithGrad(mid_pt + side * 1.5 * perp);
        std::pair<double, Vec2d> _qo = sdf.queryWithGrad(mid_pt - side * 1.5 * perp);
        if (_qc.first < 0.0 && _qo.first > _qc.first) {
            side = -side; // preferred side is inside obstacle, flip
        }
    }

    // Step 2: place 3-waypoint arch
    //   apex at 50% longitudinal, lateral offset = obstacle_radius + clearance
    //   We estimate obstacle_radius from min SDF along direct line
    double min_d = 1e18;
    for (int i = 1; i < 20; ++i) {
        double t = (double)i / 20;
        std::pair<double, Vec2d> _qd = sdf.queryWithGrad((1 - t) * p0 + t * p1);
        min_d = std::min(min_d, _qd.first);
    }
    double lateral_needed = std::max(1.0, -min_d + 1.5); // how far to bypass

    Vec2d apex = p0 + 0.5 * (p1 - p0) + side * lateral_needed * perp;

    // Step 3: build 4-segment arch through: p0 → q1 → apex → q2 → p1
    //   q1 = 25% along, lifted half-way to apex
    //   q2 = 75% along, lifted half-way to apex
    Vec2d q1 = p0 + 0.25 * (p1 - p0) + side * (lateral_needed * 0.6) * perp;
    Vec2d q2 = p0 + 0.75 * (p1 - p0) + side * (lateral_needed * 0.6) * perp;

    // Tangents: Catmull-Rom from waypoints
    std::vector<Vec2d> pts = {p0, q1, apex, q2, p1};
    std::vector<Vec2d> tans = {t0.normalized(), {}, {}, {}, t1.normalized()};
    for (int i = 1; i <= 3; ++i) {
        Vec2d d = 0.5 * (pts[i + 1] - pts[i - 1]);
        tans[i] = d.norm() > 1e-10 ? d.normalized() : (pts[i + 1] - pts[i]).normalized();
        // Clamp against net direction to avoid loops
        tans[i] = clampTangent(tans[i], along, MAX_TAN_DEV);
    }

    return makeCurveFromKnots(pts, tans, 0.35);
}

// ─────────────────────────────────────────────────────────────────────────────
//  Public API
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve buildInitialCurve(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    const SDFField& sdf, const Polygon2d& fence, const std::vector<std::vector<Vec2d>>& sibling_polys) {
    // U-turn: special case
    if (angleBetween(t0, t1) > M_PI * 0.85)
        return buildTwoSegmentUTurn(p0, t0, p1, t1, sdf, fence);

    // Penetration-only clearance: only trigger bypass when the curve actually
    // enters an obstacle (SDF < 0).  A clearance of 1.0 was too aggressive —
    // it caused straight-through lanes (whose centre-line passes 1.25m from an
    // obstacle) to incorrectly trigger the bypass logic, producing wild detours.
    //
    // The optimizer (obstacle_clearance=1.0) still pushes the final curve away
    // from obstacles to maintain safe lateral distance, but the initial curve
    // shape should be the geometrically natural one unless it literally crosses
    // an obstacle.
    constexpr double INIT_CLEARANCE = 0.0; // only bypass on actual penetration

    // ── Fast path: keep the natural single-arc turn when it is already clear.
    //
    // The old gate checked the endpoint chord before checking the actual Bezier.
    // For left/right turns the chord can pass through an obstacle even though the
    // natural turning arc is clear, which incorrectly triggered bypass shaping.
    {
        BezierSegment trial = makeCubicG1(p0, t0.normalized(), p1, t1.normalized(), 0.4);
        bool bezier_clear = true;
        int samples = std::max(40, (int)std::ceil(trial.arcLength(24) / 0.15));
        samples = std::min(samples, 320);
        for (int i = 1; i < samples; ++i) {
            std::pair<double, Vec2d> _qt = sdf.queryWithGrad(trial.evaluate((double)i / samples));
            if (_qt.first < INIT_CLEARANCE) {
                bezier_clear = false;
                break;
            }
        }
        if (bezier_clear) {
            BezierCurve c;
            c.segs.push_back(trial);
            return c;
        }
    }

    // ── Level-1: geometric direct construction ───────────────────────────────
    {
        BezierCurve arch = geometricBypass(p0, t0, p1, t1, sdf, fence, INIT_CLEARANCE, sibling_polys);
        if (!arch.empty()) return arch;
    }

    // ── Level-2: geometry-initialised full-control-point curve ───────────────
    return geometricInitLevel2(p0, t0, p1, t1, sdf);
}

// ─────────────────────────────────────────────────────────────────────────────
//  U-turn: single cubic built from aligned endpoint rays.
//
//  The old implementation inserted an apex knot between p0 and p1. For nearly
//  anti-parallel lanes this can invert the control polygon and create a visible
//  kink/S bend. The cubic below keeps only endpoint controls: first align the
//  longitudinal offset along the entry/reverse-exit rays, then use 2/3 of the
//  aligned lateral gap as the circular-arc handle length.
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve buildTwoSegmentUTurn(
    const Vec2d& p0, const Vec2d& t0,
    const Vec2d& p1, const Vec2d& t1, const SDFField& sdf, const Polygon2d&) {
    (void)sdf;
    BezierCurve c;
    c.segs.push_back(makeAlignedUTurnCubic(p0, t0, p1, t1));
    return c;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Legacy: sdfMaxClearancePath retained for any external callers
// ─────────────────────────────────────────────────────────────────────────────
std::vector<Vec2d> sdfMaxClearancePath(
    const SDFField& sdf, const Polygon2d& fence, const Vec2d& start, const Vec2d& goal, double /*alpha*/) {
    // Simplified: just return start→goal; buildInitialCurve no longer uses this.
    (void)sdf;
    (void)fence;
    return {start, goal};
}

}
