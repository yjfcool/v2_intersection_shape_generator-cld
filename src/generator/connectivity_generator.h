#pragma once
#include "types.h"
#include "optimizer/penalty_cost.h"
#include "optimizer/lbfgs_solver.h"
#include "constraints/cluster_order.h"

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

    ConnectivityCurve generateOne(
        const Connectivity&, const IntersectionInput&,
        const SDFField&, const SDFField&, const std::vector<SiblingCurve>&);

    std::vector<SiblingCurve> buildSiblings(
        const ConnId&, const std::unordered_map<ConnId, BezierCurve>&,
        const ClusterOrderSolver&, const std::vector<Connectivity>&) const;

    BezierCurve postProcess(
        const BezierCurve&, const SDFField&, const Polygon2d&, double kmax,
        const Vec2d& t0_orig, const Vec2d& t1_orig, bool skip_elastic_band = false,
        const Vec2d* p0_exact = nullptr, const Vec2d* p1_exact = nullptr);

    void validate(ConnectivityCurve&, const IntersectionInput&, const SDFField&) const;
};
