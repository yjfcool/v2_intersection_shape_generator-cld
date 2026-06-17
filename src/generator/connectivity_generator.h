#pragma once
#include "types.h"
#include "optimizer/penalty_cost.h"
#include "optimizer/lbfgs_solver.h"
#include "constraints/cluster_order.h"
#include "utils/quadtree.h"
#include <unordered_set>

namespace isg {

/// 优化分组：同一优先级内的连通关系集合
struct OptGroup {
    std::vector<ConnId> conn_ids;
    int priority = 0;  ///< 0=直行(先), 1=左/右转, 2=U型调头(后)
};

/// 全局协调器：负责按"直行→左右转→U型调头"的优先级及"由内向外"的空间顺序
/// 规划全部连通曲线的生成顺序
class GlobalCoordinator {
public:
    /// 构造生成顺序
    void build(const std::vector<Connectivity>&, const IntersectionInput&, const ClusterOrderSolver& cs);

    /// 返回所有分组（按优先级升序）
    const std::vector<OptGroup>& groups() const { return groups_; }

    /// 添加软障碍物（将已生成曲线作为软障碍加入SDF，引导后续曲线避让）
    void addSoftObstacles(SDFField&, const std::vector<ConnectivityCurve>&, double sc = 0.3) const;

private:
    std::vector<OptGroup> groups_;
};

/// 连通曲线生成器：核心生成流程
class ConnectivityGenerator {
public:
    explicit ConnectivityGenerator(
        const LBFGSConfig& cfg = {},
        const ConnectivityDirectionConfig& direction_cfg = {});

    /// 生成全部连通曲线
    /// @param out_ms 输出优化耗时(毫秒)
    std::vector<ConnectivityCurve> generate(const IntersectionInput&, SDFField&, double* out_ms = nullptr);

private:
    LBFGSSolver solver_;
    ClusterOrderSolver cluster_solver_;
    ConnectivityDirectionConfig direction_cfg_;
    std::unique_ptr<QuadTree> quad_tree_;  ///< 持久化空间索引，加速同簇兄弟曲线查询

private:
    /// 生成单条连通曲线
    ConnectivityCurve generateOne(
        const Connectivity& conn, const IntersectionInput& input,
        const SDFField& sdf, const SDFField& sdf_coarse, const std::vector<SiblingCurve>& siblings,
        bool allow_uturn_search = false, bool* out_physical_risk = nullptr);

    /// 构造同簇兄弟曲线列表
    std::vector<SiblingCurve> buildSiblings(
        const ConnId& id, const std::unordered_map<ConnId, BezierCurve>& done,
        const ClusterOrderSolver& cs, const std::vector<Connectivity>& conns,
        bool constrained_only = false,
        const std::unordered_set<ConnId>* fixed_shape_ids = nullptr) const;

    /// 后处理：弹性带平滑+G1强制+自适应细分
    BezierCurve postProcess(
        const BezierCurve& c, const SDFField& sdf, const Polygon2d& fence, double kmax,
        const Vec2d& t0_orig, const Vec2d& t1_orig, bool skip_elastic_band = false,
        const Vec2d* p0_exact = nullptr, const Vec2d* p1_exact = nullptr);

    /// 校验曲线违规信息
    void validate(ConnectivityCurve&, const IntersectionInput&, const SDFField&) const;
};

}
