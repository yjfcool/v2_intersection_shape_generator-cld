#include "penalty_cost.h"
#include "constraints/fence_check.h"
#include "curve/curve_utils.h"
#include "utils.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace isg {

// ─── Cache rebuild ────────────────────────────────────────────────────────────
void PenaltyCostCache::rebuild(
        const std::vector<Boundary> &boundaries,
        const std::vector<SiblingCurve> &siblings, const Polygon2d &fence) {
    bnd_segs.clear();
    for (auto& bnd : boundaries) {
        auto& pts = bnd.geometry.points;
        for (int i = 0; i + 1 < (int)pts.size(); ++i)
            bnd_segs.push_back({pts[i], pts[i + 1]});
    }

    sib_polys.clear();
    for (auto& sib : siblings) {
        SiblingPoly sp;
        sp.exempt_a1 = sib.exempt_a1;
        sp.a2_radius = sib.exempt_a2_radius;
        sp.expected_side = sib.expected_side;
        sp.ref_perp = sib.ref_perp;
        sp.shared_endpoint = sib.shared_endpoint; // wider skip zone for convergent pairs
        sp.pts = sib.curve.sampleByArcLength(48); // 48pts for finer sign-change detection
        sib_polys.push_back(std::move(sp));
    }
    fence_empty = fence.outer.empty();
}

void PenaltyCost::buildCache() {
    cache_.rebuild(boundaries, siblings, fence);
}

// ─── Scalar terms ─────────────────────────────────────────────────────────────
double PenaltyCost::evalSmooth(const BezierCurve &c) const {
    double cost = 0;
    constexpr int S = 16;
    for (auto& seg : c.segs) {
        double ds = seg.arcLength(S) / S;
        for (int i = 0; i <= S; ++i) {
            double k = seg.curvature((double)i / S);
            cost += k * k * ds;
        }
    }
    return cost;
}

double PenaltyCost::evalObstacle(const BezierCurve& c) const {
    if (!sdf) return 0;
    double cost = 0;
    constexpr int S = 24;
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i)
            cost += sdf->obstaclePenalty(seg.evaluate((double)i / S), obstacle_clearance);
    return cost;
}

double PenaltyCost::evalBoundary(const BezierCurve& c) const {
    double cost = 0;
    constexpr int S = 20;
    std::vector<Vec2d> cp;
    cp.reserve(c.segs.size() * (S + 1));
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i)
            cp.push_back(seg.evaluate((double)i / S));

    for (auto& bs : cache_.bnd_segs)
        for (int ci = 0; ci + 1 < (int)cp.size(); ++ci)
            if (segmentsIntersect(cp[ci], cp[ci + 1], bs.a, bs.b))
                cost += 1.0;
    return cost;
}

double PenaltyCost::evalFence(const BezierCurve& c) const {
    if (cache_.fence_empty)
        return 0;
    double cost = 0;
    constexpr int S = 20;
    for (auto& seg : c.segs)
        for (int i = 0; i <= S; ++i) {
            Vec2d pt = seg.evaluate((double)i / S);
            if (!polygonContains(fence, pt))
                cost += pointToPolygonDist(pt, fence) * 2.0;
        }
    return cost;
}

// ─────────────────────────────────────────────────────────────────────────────
//  evalCluster — Order-Preserving Topology Constraint
// ─────────────────────────────────────────────────────────────────────────────
double PenaltyCost::evalCluster(const BezierCurve& c) const {
    if (cache_.sib_polys.empty())return 0;

    constexpr double MARGIN = 0.3;   // m lateral separation
    // raised 30.0 → 50.0 — a detected sign-change (actual path crossing)
    // now produces a much stronger gradient signal than a simple margin violation,
    // driving the optimizer to uncross adjacent same-direction turns.
    constexpr double SIGN_CHANGE_COST = 50.0;  // per-crossing base cost
    constexpr double OBS_REDUCE = 0.4;   // weight when obstacle forced (NOT zero)
    constexpr double MAX_DIST = 12.0;  // m nearest-point search radius
    constexpr double SKIP_FRAC = 0.07;  // skip first/last 7%
    constexpr int N_SAMP = 48;    // match sibling cache size

    auto curve_pts = c.sampleByArcLength(N_SAMP);
    int N = (int) curve_pts.size();
    if (N < 4)return 0;

    double cost = 0;

    for (auto& sp : cache_.sib_polys) {
        if (sp.exempt_a1)continue; // StructuralCross or different arm → skip

        double w = (sp.a2_radius > 0) ? OBS_REDUCE : 1.0; // reduced but NOT zero for obstacle
        int M = (int) sp.pts.size();
        if (M < 2 || sp.ref_perp.norm() < 1e-9)continue;

        // For curves sharing the same exit endpoint (shared_endpoint=true) the
        // constraint is meaningless at the convergence point.  Use a 25% skip
        // zone so the penalty is only applied to the 50% central portion of the
        // arc, well away from the shared endpoint.
        double eff_skip = SKIP_FRAC; // double eff_skip = sp.shared_endpoint ? 0.25 : SKIP_FRAC;

        // Pre-compute nearest-sibling-point lateral diffs for all samples
        std::vector<double> diffs(N, std::numeric_limits<double>::quiet_NaN());
        for (int i = 0; i < N; ++i) {
            double frac = (double) i / std::max(N - 1, 1);
            if (frac < eff_skip || frac > 1.0 - eff_skip)continue;
            Vec2d cp = curve_pts[i];
            // Find nearest sibling point within MAX_DIST
            double best_d2 = MAX_DIST * MAX_DIST + 1;
            int best_k = -1;
            for (int k = 0; k < M; ++k) {
                double d2 = (sp.pts[k] - cp).squaredNorm();
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_k = k;
                }
            }
            if (best_k < 0)continue;
            Vec2d sib_pt = sp.pts[best_k];
            double cur_lat = cp.dot(sp.ref_perp);
            double sib_lat = sib_pt.dot(sp.ref_perp);
            // diff = latA - latB (positive = current is LEFT of sibling)
            diffs[i] = cur_lat - sib_lat;

            // [A] Margin penalty
            // expected_side=+1: sibling LEFT of current → diff should be NEGATIVE (sib>cur)
            //   → violation = max(0, diff + MARGIN)
            // expected_side=-1: sibling RIGHT of current → diff should be POSITIVE (cur>sib)
            //   → violation = max(0, -diff + MARGIN)
            double violation = 0;
            if (sp.expected_side == +1)
                violation = std::max(0.0, diffs[i] + MARGIN); // diff should be < -MARGIN
            else
                violation = std::max(0.0, -diffs[i] + MARGIN); // diff should be >  MARGIN
            cost += w * violation * violation;
        }

        // [B] Sign-change (Order Preservation) penalty
        // Walk through consecutive sample pairs; detect sign flip in diffs[]
        {
            double prev_diff = std::numeric_limits<double>::quiet_NaN();
            for (int i = 0; i < N; ++i) {
                if (std::isnan(diffs[i])) {
                    prev_diff = std::numeric_limits<double>::quiet_NaN();
                    continue;
                }
                if (!std::isnan(prev_diff) && std::abs(prev_diff) > 1e-6 && std::abs(diffs[i]) > 1e-6) {
                    if (prev_diff * diffs[i] < 0) {
                        // Sign change → crossing detected
                        double depth = (std::abs(prev_diff) + std::abs(diffs[i])) * 0.5;
                        cost += w * SIGN_CHANGE_COST * depth;
                    }
                }
                prev_diff = diffs[i];
            }
        }

        // [C] Exhaustive segment intersection check (cross-fraction crossing detection)
        // Covers cases where curves cross at very different arc fractions:
        //   e.g. conn32 at t=0.25 crossing conn18 at t=0.84 — nearest-point misses it
        //   because at those positions the nearest sibling point is > MAX_D away.
        // O(N*M) but N=M=48 → 2304 checks, acceptable.
        constexpr double EP_TOL = 1.5; // endpoint exclusion tolerance (m)
        for (int ci = 0; ci + 1 < N; ++ci) {
            double frac_ci = (double) ci / std::max(N - 1, 1);
            if (frac_ci < eff_skip || frac_ci > 1.0 - eff_skip)
                continue;
            Vec2d A0 = curve_pts[ci], A1 = curve_pts[ci + 1];
            Vec2d mid_A = 0.5 * (A0 + A1);
            for (int si = 0; si + 1 < M; ++si) {
                Vec2d B0 = sp.pts[si], B1 = sp.pts[si + 1];
                // Quick bbox-distance reject: skip if midpoints more than MAX_D+2 apart
                if ((mid_A - 0.5 * (B0 + B1)).norm() > MAX_DIST + 2.0)
                    continue;
                Vec2d isect;
                if (!segmentsIntersect(A0, A1, B0, B1, &isect))
                    continue;
                // Skip if crossing is within endpoint zone of CURRENT curve
                double d_start = (isect - curve_pts.front()).norm();
                double d_end = (isect - curve_pts.back()).norm();
                if (d_start < EP_TOL || d_end < EP_TOL)continue;
                // Crossing depth = lateral separation just before crossing
                // Use ref_perp projection difference at A0 vs nearest B point
                Vec2d near_B = B0;
                double bd2 = (B0 - A0).squaredNorm();
                double b1d2 = (B1 - A0).squaredNorm();
                if (b1d2 < bd2) {
                    bd2 = b1d2;
                    near_B = B1;
                }
                double lat_A0 = A0.dot(sp.ref_perp);
                double lat_NB = near_B.dot(sp.ref_perp);
                double depth = std::max(std::abs(lat_A0 - lat_NB), 0.3);
                cost += w * SIGN_CHANGE_COST * depth;
            }
        }
    }
    return cost;
}

double PenaltyCost::evalCrosswalk(const BezierCurve& c) const {
    constexpr double DEG15 = 15.0 * M_PI / 180.0;
    double cost = 0;
    for (auto& cw : crosswalks) {
        for (auto& seg : c.segs)
            for (int i = 1; i < 20; ++i) {
                double t = (double)i / 20;
                Vec2d pt = seg.evaluate(t);
                if (!polygonContains(cw.geometry, pt))
                    continue;
                Vec2d tan = seg.tangent(t).normalized();
                double cosA = std::abs(tan.dot(cw.crossing_direction.normalized()));
                double exc = cosA - std::cos(M_PI / 2 - DEG15);
                if (exc > 0) cost += exc * exc;
            }
    }
    return cost;
}

// ─── G1 endpoint enforcement ──────────────────────────────────────────────────
static void enforceG1(
    BezierCurve& c, const Vec2d& start_dir, const Vec2d& end_dir) {
    if (c.segs.empty())
        return;

    // Start: ctrl[1] must lie on ray ctrl[0] + λ*start_dir, λ > 0
    if (start_dir.norm() > 1e-10) {
        Vec2d& p0 = c.segs.front().ctrl[0];
        Vec2d& c1 = c.segs.front().ctrl[1];
        Vec2d sd = start_dir.normalized();
        // Project current ctrl[1] onto ray to find scale λ
        double lam = (c1 - p0).dot(sd);
        lam = std::max(lam, 0.05); // minimum tension: avoid degenerate segs
        c1 = p0 + lam * sd;
    }

    // End: ctrl[2] must lie on ray ctrl[3] - μ*end_dir, μ > 0
    if (end_dir.norm() > 1e-10) {
        Vec2d& p1 = c.segs.back().ctrl[3];
        Vec2d& c2 = c.segs.back().ctrl[2];
        Vec2d ed = end_dir.normalized();
        double mu = (p1 - c2).dot(ed);
        mu = std::max(mu, 0.05);
        c2 = p1 - mu * ed;
    }
}

// ─── Analytic obstacle gradient ───────────────────────────────────────────────
void PenaltyCost::addObstacleGrad(
    const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const {
    if (!sdf || w < 1e-15)
        return;
    constexpr int S = 24;
    for (int k = 0; k < c.segs.size(); ++k) {
        auto& seg = c.segs[k];
        for (int i = 0; i <= S; ++i) {
            double t = (double)i / S;
            Vec2d pt = seg.evaluate(t);
            auto kv = sdf->queryWithGrad(pt);
            double& d = kv.first;
            Vec2d& gd = kv.second;
            double slack = d - obstacle_clearance;
            if (slack >= 0)
                continue; // no violation → zero gradient

            double coeff = 2.0 * w * slack; // negative (violation)

            // Basis function derivatives at t for ctrl[1] and ctrl[2]
            double B1 = 3 * (1 - t) * (1 - t) * t; // ∂c/∂ctrl[1]
            double B2 = 3 * (1 - t) * t * t; // ∂c/∂ctrl[2]

            int base = 4 * k;
            // p_{4k}   = ctrl[1].x  →  ∂c.x/∂p = B1, ∂c.y/∂p = 0
            grad[base + 0] += coeff * (gd[0] * B1);
            // p_{4k+1} = ctrl[1].y  →  ∂c.x/∂p = 0, ∂c.y/∂p = B1
            grad[base + 1] += coeff * (gd[1] * B1);
            // p_{4k+2} = ctrl[2].x
            grad[base + 2] += coeff * (gd[0] * B2);
            // p_{4k+3} = ctrl[2].y
            grad[base + 3] += coeff * (gd[1] * B2);
        }
    }
}

//  Analytic smooth gradient
void PenaltyCost::addSmoothGrad(
    const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const {
    if (w < 1e-15)
        return;
    constexpr int S = 12;
    const double h = 1e-4;
    for (int k = 0; k < c.segs.size(); ++k) {
        int base = 4 * k;
        double ds = c.segs[k].arcLength(S) / S;
        for (int j = 0; j < 4; ++j) {
            VecXd pp = params, pm = params;
            pp[base + j] += h;
            pm[base + j] -= h;

            auto cp = full_param_mode ? curveFromParamsFull(pp, proto) : curveFromParams(pp, proto);
            auto cm = full_param_mode ? curveFromParamsFull(pm, proto) : curveFromParams(pm, proto);
            enforceG1(cp, start_tan_dir, end_tan_dir);
            enforceG1(cm, start_tan_dir, end_tan_dir);

            double fp = 0, fm = 0;
            for (int i = 0; i <= S; ++i) {
                double t = (double)i / S;
                double kp = cp.segs[k].curvature(t);
                double km = cm.segs[k].curvature(t);
                fp += kp * kp * ds;
                fm += km * km * ds;
            }
            grad[base + j] += w * (fp - fm) / (2 * h);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Numeric gradient for secondary terms (boundary/fence/cluster/crosswalk)
//  ONLY runs when these terms have non-zero cost → skipped in clean iterations.
// ─────────────────────────────────────────────────────────────────────────────
void PenaltyCost::addNumericGrad(
        const VecXd &params, double wb, double wf, double wc, double wx, VecXd &grad) {
    BezierCurve c0 = full_param_mode ? curveFromParamsFull(params, proto) : curveFromParams(params, proto);
    bool nb = (wb > 1e-6 && evalBoundary(c0) > 1e-9);
    bool nf = (wf > 1e-6 && evalFence(c0) > 1e-9);
    bool nc = (wc > 1e-6 && evalCluster(c0) > 1e-9);
    bool nx = (wx > 1e-6 && evalCrosswalk(c0) > 1e-9);
    if (!nb && !nf && !nc && !nx)return;
    const double h = 1e-4;
    const int n = (int)params.size();
    VecXd p = params;
    for (int i = 0; i < n; ++i) {
        double orig = p[i];
        auto eval = [&](VecXd &pv) {
            auto cv = full_param_mode ? curveFromParamsFull(pv, proto) : curveFromParams(pv, proto);
            enforceG1(cv, start_tan_dir, end_tan_dir);
            double f = 0;
            if (nb)f += wb * evalBoundary(cv);
            if (nf)f += wf * evalFence(cv);
            if (nc)f += wc * evalCluster(cv);
            if (nx)f += wx * evalCrosswalk(cv);
            return f;
        };
        p[i] = orig + h;
        double fp = eval(p);
        p[i] = orig - h;
        double fm = eval(p);
        grad[i] += (fp - fm) / (2 * h);
        p[i] = orig;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
//  Main operator(): one forward pass → f + grad
//  Total cost: O(n_segs × S × analytic) instead of O(2n × full_eval)
// ─────────────────────────────────────────────────────────────────────────────
double PenaltyCost::operator()(const VecXd& params, VecXd& grad) {
    const int n = (int)params.size();
    grad.resize(n);
    grad.setZero();
    BezierCurve c = full_param_mode ? curveFromParamsFull(params, proto) : curveFromParams(params, proto);

    enforceG1(c, start_tan_dir, end_tan_dir);

    // ── Cost evaluation ───────────────────────────────────────────────────
    double f_smooth = evalSmooth(c);
    double f_obstacle = evalObstacle(c);
    double f_boundary = evalBoundary(c);
    double f_fence = evalFence(c);
    double f_cluster = evalCluster(c);
    double f_xwalk = evalCrosswalk(c);

    double f = weights.smooth * f_smooth
        + weights.obstacle * f_obstacle
        + weights.boundary * f_boundary
        + weights.fence * f_fence
        + weights.cluster * f_cluster
        + weights.crosswalk * f_xwalk;

    // smooth — sparse analytic (4 params per segment, S samples each)
    addSmoothGrad(c, params, weights.smooth, grad);

    // obstacle — analytic via SDF gradient (dominant + cheap)
    addObstacleGrad(c, params, weights.obstacle, grad);

    // numeric diff, but skipped when cost is zero (common after convergence)
    addNumericGrad(params, weights.boundary, weights.fence, weights.cluster, weights.crosswalk, grad);

    return f;
}

// ─────────────────────────────────────────────────────────────────────────────
//  optimiseCurve — outer adaptive-weight loop
// ─────────────────────────────────────────────────────────────────────────────
BezierCurve optimiseCurve(
    PenaltyCost& cost, LBFGSSolver& solver, const BezierCurve& initial, int outer_iters) {
    cost.proto = initial;
    cost.buildCache();
    VecXd params = cost.full_param_mode ? curveToParamsFull(initial) : curveToParams(initial);
    for (int outer = 0; outer < outer_iters; ++outer) {
        auto res = solver.solve([&](const VecXd& p, VecXd& g) { return cost(p, g); }, params);
        params = res.x;
        BezierCurve c = cost.full_param_mode ?
                        curveFromParamsFull(params, cost.proto) : curveFromParams(params, cost.proto);
        double op = cost.evalObstacle(c);
        double bp = cost.evalBoundary(c);
        double fp = cost.evalFence(c);
        double cp = cost.evalCluster(c);
        cost.weights.update(op, bp, fp, cp);

        if (op + bp + fp + cp < 1e-6)
            break; // all constraints satisfied
    }
    return cost.full_param_mode ? curveFromParamsFull(params, cost.proto) : curveFromParams(params, cost.proto);
}

}