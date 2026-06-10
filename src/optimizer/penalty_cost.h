#pragma once

#include "types.h"
#include "sdf_field.h"
#include "lbfgs_solver.h"
#include <vector>

namespace isg {

struct SiblingCurve {
    BezierCurve curve;
    bool exempt_a1 = false; // true = StructuralCross or cross-arm → skip penalty
    double exempt_a2_radius = 0.0;   // >0 = obstacle forced (soft, 0.4× weight, NOT zero)
    int expected_side = 0;     // +1: sibling LEFT of current; -1: RIGHT; 0: unknown
    Vec2d ref_perp{0, 0};   // pair-specific fixed lateral axis
};

struct PenaltyWeights {
    double smooth = 1.0;
    double obstacle = 30.0;  // highest priority — obstacle avoidance
    double boundary = 10.0;
    double fence = 6.0;
    double cluster = 12.0;  // topology order — lower than obstacle
    double crosswalk = 0.5;

    void update(double op, double bp, double fp, double cp) {
        if (op > 1e-3) obstacle = std::min(obstacle * 1.8, 240.0);
        if (bp > 1e-3) boundary = std::min(boundary * 2.0, 160.0);
        if (fp > 1e-3) fence = std::min(fence * 2.0, 96.0);
        if (cp > 1e-3) cluster = std::min(cluster * 1.5, 96.0);
    }
};

struct BoundarySegment {
    Vec2d a, b;
};

struct PenaltyCostCache {
    struct SiblingPoly {
        std::vector<Vec2d> pts;
        bool exempt_a1;
        double a2_radius;
        int expected_side;
        Vec2d ref_perp;
    };

    std::vector<BoundarySegment> bnd_segs;
    std::vector<SiblingPoly> sib_polys;
    bool fence_empty = true;

    void rebuild(const std::vector<Boundary>& bounds, const std::vector<SiblingCurve>& sibls, const Polygon2d& fence);
};

class PenaltyCost {
public:
    BezierCurve proto;
    PenaltyWeights weights;
    const SDFField* sdf = nullptr;
    std::vector<Boundary> boundaries;
    Polygon2d fence;
    std::vector<SiblingCurve> siblings;
    std::vector<Crosswalk> crosswalks;
    double obstacle_clearance = 0.0;
    bool full_param_mode = false;
    Vec2d start_tan_dir{0, 0}; // set to t0 from Lane::geometry
    Vec2d end_tan_dir{0, 0}; // set to t1 from Lane::geometry

    void buildCache();

    double operator()(const VecXd& params, VecXd& grad);

    double evalSmooth(const BezierCurve& c) const;
    double evalObstacle(const BezierCurve& c) const;
    double evalBoundary(const BezierCurve& c) const;
    double evalFence(const BezierCurve& c) const;
    double evalCluster(const BezierCurve& c) const;
    double evalCrosswalk(const BezierCurve& c) const;

private:
    PenaltyCostCache cache_;

    // Compute analytic gradient contributions from obstacle SDF term
    // d/dp_j [ Σ_i max(0, cl - d(c(t_i)))² ]
    void addObstacleGrad(const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const;

    // Numeric central-diff gradient for boundary+fence+cluster (cheap compared
    // to obstacle since these terms are usually zero after early iterations)
    void addNumericGrad(const VecXd& params, double w_bnd, double w_fence, double w_cluster, double w_xwalk, VecXd& grad);

    // Smooth term analytic gradient via curvature derivative
    void addSmoothGrad(const BezierCurve& c, const VecXd& params, double w, VecXd& grad) const;
};

BezierCurve optimiseCurve(PenaltyCost& cost, LBFGSSolver& solver, const BezierCurve& initial, int outer_iters = 5);

}