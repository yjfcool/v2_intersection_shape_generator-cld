#pragma once
#include "types.h"
#include "optimizer/penalty_cost.h"
#include "optimizer/lbfgs_solver.h"
#include "constraints/cluster_order.h"
#include "utils/quadtree.h"

struct OptGroup {
    std::vector<ConnId> conn_ids;
    int priority = 0;
};

class GlobalCoordinator {
public:
    void build(const std::vector<Connectivity>&, const IntersectionInput&, const ClusterOrderSolver& cs);
    const std::vector<OptGroup>& groups() const { return groups_; }
    void addSoftObstacles(SDFField&, const std::vector<ConnectivityCurve>&, double sc = 0.3) const;

private:
    std::vector<OptGroup> groups_;
};

class ConnectivityGenerator {
public:
    explicit ConnectivityGenerator(const LBFGSConfig& cfg = {});
    std::vector<ConnectivityCurve> generate(const IntersectionInput&, SDFField&, double* out_ms = nullptr);

private:
    LBFGSSolver solver_;
    ClusterOrderSolver cluster_solver_;
    std::unique_ptr<QuadTree> quad_tree_; // Persistent spatial index for faster sibling queries across all curves

private:
    ConnectivityCurve generateOne(
        const Connectivity& conn, const IntersectionInput& input,
        const SDFField& sdf, const SDFField& sdf_coarse, const std::vector<SiblingCurve>& siblings);

    std::vector<SiblingCurve> buildSiblings(
        const ConnId& id, const std::unordered_map<ConnId, BezierCurve>& done,
        const ClusterOrderSolver& cs, const std::vector<Connectivity>& conns) const;

    BezierCurve postProcess(
        const BezierCurve& c, const SDFField& sdf, const Polygon2d& fence, double kmax,
        const Vec2d& t0_orig, const Vec2d& t1_orig, bool skip_elastic_band = false,
        const Vec2d* p0_exact = nullptr, const Vec2d* p1_exact = nullptr);

    void validate(ConnectivityCurve&, const IntersectionInput&, const SDFField&) const;
};
