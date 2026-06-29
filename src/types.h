#pragma once
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <map>

namespace isg {

using Vec2d = Eigen::Vector2d;
using VecXd = Eigen::VectorXd;

inline std::vector<std::array<double, 2>> toArray(const std::vector<Vec2d>& pts) {
    std::vector<std::array<double, 2>> arrpts;
    for (auto p : pts)
        arrpts.emplace_back(std::array<double, 2>{p[0], p[1]});
    return arrpts;
}

inline std::vector<Vec2d> toArray(const std::vector<std::array<double, 2>>& arrpts) {
    std::vector<Vec2d> pts;
    for (auto p : arrpts)
        pts.emplace_back(Vec2d{p[0], p[1]});
    return pts;
}

// ── 几何基元 ──────────────────────────────────────────────────
/// 二维包围盒
struct BoundingBox2d {
    Vec2d min_pt{1e18, 1e18};
    Vec2d max_pt{-1e18, -1e18};

    void expand(const Vec2d& p) {
        min_pt = min_pt.cwiseMin(p);
        max_pt = max_pt.cwiseMax(p);
    }

    bool intersects(const BoundingBox2d& o) const {
        return max_pt[0] >= o.min_pt[0] && min_pt[0] <= o.max_pt[0] &&
            max_pt[1] >= o.min_pt[1] && min_pt[1] <= o.max_pt[1];
    }

    bool contains(const Vec2d& p) const {
        return p[0] >= min_pt[0] && p[0] <= max_pt[0] &&
            p[1] >= min_pt[1] && p[1] <= max_pt[1];
    }

    double width() const { return max_pt[0] - min_pt[0]; }
    double height() const { return max_pt[1] - min_pt[1]; }
    bool empty() const { return min_pt[0] > max_pt[0]; }
};

/// 二维折线
struct LineString2d {
    std::vector<Vec2d> points;

    BoundingBox2d bbox() const {
        BoundingBox2d b;
        for (auto& p : points)
            b.expand(p);
        return b;
    }
};

/// 二维多边形（含外环与洞）
struct Polygon2d {
    std::vector<Vec2d> outer;
    std::vector<std::vector<Vec2d>> holes;
    bool empty() const { return outer.empty(); }

    BoundingBox2d bbox() const {
        BoundingBox2d b;
        for (auto& p : outer)
            b.expand(p);
        return b;
    }
};

// ── Bezier 类型（在此定义以避免循环依赖）─────────────────────
/// 三次Bezier段（4个控制点）
struct BezierSegment {
    std::array<Vec2d, 4> ctrl;

    Vec2d evaluate(double t) const;  ///< B(t) 求值
    Vec2d tangent(double t) const;
    double curvature(double t) const;  ///< 曲率 κ(t) = |B'×B''| / |B'|³
    std::pair<BezierSegment, BezierSegment> splitAt(double t) const;
    BoundingBox2d bbox() const;
    double arcLength(int samples = 20) const;             ///< 弧长（数值积分 - 辛普森法）
    double arcLengthToParam(double s, int samples = 50) const;
    Vec2d evalDeriv1(double t) const;  ///< B'(t) 一阶导数
    Vec2d evalDeriv2(double t) const;  ///< B''(t) 二阶导数
    double maxCurvature(int samples = 50) const;          ///< 最大曲率（采样估算）
    std::vector<Vec2d> sampleCount(int n) const;          ///< 固定数量采样
    std::vector<Vec2d> sampleBySpacing(double spacing) const;  ///< 固定间距采样（弧长近似）
    std::vector<Vec2d> sampleAdaptive(double maxAngleDeg, double maxSeg, double minSeg) const;  ///< 自适应采样（基于角度误差）
    std::vector<Vec2d> sampleByShape(double straightSpacing, double angularResolution, double minSpacing) const;
    void subdivide(double t0, double t1, double maxAngle, double maxSeg, double minSeg,
                   std::vector<Vec2d>& pts, int depth) const;
    double getAlpha(const Vec2d& T0) const;  ///< 获取 α（P1相对P0的标量，T0方向）
    double getBeta(const Vec2d& T3) const;
};

/// Bezier曲线：由若干段三次Bezier段拼接而成
struct BezierCurve {
    std::vector<BezierSegment> segs;

    bool empty() const { return segs.empty(); }
    int numSegments() const { return (int)segs.size(); }
    Vec2d startPt() const;
    Vec2d endPt() const;
    Vec2d startTan() const;
    Vec2d endTan() const;
    Vec2d evaluate(double u) const;
    Vec2d tangent(double u) const;
    std::vector<Vec2d> sample(int n) const;
    std::vector<Vec2d> sampleByArcLength(int n) const;
    std::vector<Vec2d> sampleByShape(double straightSpacing = 3.0, double angularResolution = 10.0, double minSpacing = 0.05) const;  ///< 形态自适应采样
    double arcLength() const;
    BoundingBox2d bbox() const;
    double maxCurvature(int sps = 20) const;
};

// ── ID 与属性别名 ─────────────────────────────────────────────
using LaneId       = std::string;  ///< 车道ID
using LaneGroupId  = std::string;  ///< 车道组ID
using LaneEdgeId   = std::string;  ///< 车道边线ID
using ConnId       = std::string;  ///< 连通关系ID
using InterId      = std::string;  ///< 路口ID
using AttrMap      = std::map<std::string, std::string>;  ///< 属性映射

/// 车道组角色：进入组或退出组
enum class GroupRole { Entry, Exit };

/// 转向类型（按几何方位计算优先，输入数据的turn_type不一定可信）
enum class ConnTurnType { Unknown = 0, TurnLeft = 1, UTurnLeft = 2, Straight = 3, TurnRight = 4, UTurnRight = 5 };

enum class ConnLaneType { Motorway = 0, NonMotorway = 1 };

/// 车道边线：车道左右侧边界线，可被相邻车道共享
struct LaneEdge {
    LaneEdgeId id;
    LineString2d geometry;
    bool is_shared = false;                              ///< 是否被相邻车道共享
    std::shared_ptr<std::pair<LaneId, LaneId>> shared_by = nullptr;  ///< 共享该边线的两个车道
    AttrMap attrs;
    std::string groupId;
    int lineOrder = 0; ///< 组内横向排序（0=最内侧）
};

/// 车道：包含左右边线及几何中心线
struct Lane {
    LaneId id;
    LaneEdgeId left_edge_id;
    LaneEdgeId right_edge_id;
    double width = 3.5;       ///< 车道宽度（米）
    LineString2d geometry;    ///< 车道中心线，进入车道末点为路口端点，退出车道首点为路口端点
    AttrMap attrs;
    std::string groupId;
    int laneOrder = 0; ///< 组内横向排序（0=最内侧）
};

/// 车道组：同方向同侧的若干车道的集合
struct LaneGroup {
    LaneGroupId id;
    GroupRole role;
    std::vector<LaneId> lanes;
    std::vector<LaneEdgeId> boundaries;  ///< 组内车道间边界线集合
    AttrMap attrs;
};

/// 连通关系：从进入车道到退出车道的一条行驶路径
struct Connectivity {
    ConnId id;
    LaneId entry_lane_id;
    LaneId exit_lane_id;
    ConnTurnType turn_type = ConnTurnType::Straight;

    LaneGroupId enterGroupId;
    LaneGroupId exitGroupId;

    LineString2d geometry;  ///< 固有形态（可选）：通常为已生成好的特定形态曲线非空（坐标数>=2）时除穿障必须避障可重新生成外，
                            ///  其他场景不能修改固有形态；同簇非交约束和避障时都需要额外考虑此固有形态影响
    ConnLaneType lane_type = ConnLaneType::Motorway;
    bool fixed_shape = false;
};

/// 障碍物（带缓冲多边形）
struct Obstacle {
    std::string id;
    Polygon2d geometry, buffered_geometry;
};

/// 边界线（道路边缘、中央隔离带、绿化带等）
struct Boundary {
    enum class Type { RoadEdge, MedianStrip, GreenBelt, Other };
    std::string id;
    Type type = Type::Other;
    LineString2d geometry;
};

/// 停止线：红灯停车位置（与人行横道是不同概念）
struct StopLine {
    std::string id;
    LineString2d geometry;
    LaneGroupId associated_group_id;
    Vec2d normal_direction{0, 1};
};

/// 人行横道（面）：调头需跨越此多边形后才能开始拱弧
struct Crosswalk {
    std::string id;
    Polygon2d geometry;
    Vec2d crossing_direction{0, 1};  ///< 行人通行方向
};

/// 路口面
struct IntersectionArea {
    InterId id;
    Polygon2d geometry;
    bool is_rough = false;  ///< true:粗糙路口面；false:精细路口面
};

/// 曲线状态
enum class CurveStatus { OK, WarnA2, Degraded, Infeasible };

/// 违规信息：记录曲线生成过程中的各类违规
struct ViolationInfo {
    enum class InfeasibilityType { None, NarrowPassage, Sandwich, TopologicalBlock, ForcedCross };

    InfeasibilityType type = InfeasibilityType::None;
    double max_obstacle_penetration = 0, max_fence_overflow = 0, fence_expansion_applied = 0;
    std::vector<Vec2d> exempt_crosses;
    std::string reason;
};

/// 连通曲线：连通关系对应的Bezier曲线及违规信息
struct ConnectivityCurve {
    ConnId id;
    LaneId entry_lane_id;
    LaneId exit_lane_id;
    ConnTurnType turn_type = ConnTurnType::Straight;
    std::shared_ptr<BezierCurve> curve = nullptr;  ///< 生成的Bezier曲线
    LineString2d geometry;
    CurveStatus status = CurveStatus::OK;
    ViolationInfo violation;

    LaneEdgeId left_edge_id = "";
    LaneEdgeId right_edge_id = "";
    ConnLaneType lane_type = ConnLaneType::Motorway;
    bool fixed_shape = false;
};

/// 连通车道边线：可被相邻连通曲线共享
struct ConnectivityLaneEdge {
    LaneEdgeId id;
    LineString2d geometry;
    bool is_shared = false;
    std::shared_ptr<std::pair<LaneId, LaneId>> shared_by = nullptr;  ///< 边线左侧车道，边线右侧车道

    AttrMap attrs;
    std::string groupId;
    int lineOrder = 0;  ///< 组内横向排序（0=最内侧）
};

/// 路口生成输出
struct IntersectionOutput {
    std::vector<ConnectivityCurve> connectivity_curves;  ///< 全部连通曲线
    std::vector<ConnectivityLaneEdge> lane_edges;        ///< 全部车道边线
    IntersectionArea area;                               ///< 精细路口面

    /// 性能统计（毫秒）
    struct PerfStats {
        double sdf_build_ms = 0, precheck_ms = 0, optimize_ms = 0,
               smooth_ms = 0, edge_gen_ms = 0, area_gen_ms = 0;
    } perf;
};


/// 路口生成输入
struct IntersectionInput {
    std::vector<LaneGroup> lane_groups;
    std::vector<Lane> lanes;
    std::vector<LaneEdge> lane_edges;
    std::vector<Connectivity> connectivities;
    std::vector<Obstacle> obstacles;
    std::vector<Boundary> boundaries;
    std::vector<StopLine> stop_lines;
    std::vector<Crosswalk> crosswalks;
    IntersectionArea area;

    std::string id;  ///< 路口ID

    const bool IsEntryLaneEdge(const LaneEdgeId& id) const;  ///< 是否为进入车道边线
    const bool IsEntryLane(const LaneId& id) const;          ///< 是否为进入车道
    const Lane* findLane(const LaneId&) const;
    const LaneGroup* findGroup(const LaneGroupId&) const;
    const LaneEdge* findEdge(const LaneEdgeId&) const;
    bool laneGroupExists(const LaneGroupId&) const;
    std::pair<Vec2d, Vec2d> entryPtDir(const LaneId&) const;  ///< 返回进入车道末点及切向（指向路口内）
    std::pair<Vec2d, Vec2d> exitPtDir(const LaneId&) const;   ///< 返回退出车道首点及切向（指向路口外）
};

// ── 配置相关 ─────────────────────────────────────────────────
struct SDFField;  ///< 前向声明

/// 连通方向统一模式
enum class ConnectivityDirectionMode {
    PerLane = 0,        ///< 每条车道独立计算方向
    GroupUnified = 1    ///< 同组车道统一方向（推荐）
};

/// 连通方向配置
struct ConnectivityDirectionConfig {
    ConnectivityDirectionMode mode = ConnectivityDirectionMode::GroupUnified;
    double group_similarity_angle_deg = 5.0;  ///< 同组车道方向相似性阈值（度）
};

/// LBFGS优化器配置
struct LBFGSConfig {
#ifndef NDEBUG
    // Debug模式: 大幅降低迭代数,避免调试时分钟级耗时
    int max_iter = 10;          ///< Debug最大迭代(性能优化:从30降至10)
    int history_size = 3;       ///< Debug历史大小
    int max_ls_iter = 5;        ///< Debug线搜索迭代
    double grad_tol = 5e-3;     ///< Debug梯度收敛(放宽)
    double func_tol = 1e-4;     ///< Debug函数值收敛(放宽)
#else
    int max_iter = 30;          ///< 最大迭代次数(性能调优:从80降至30)
    int history_size = 5;       ///< 历史记录大小(从10降至5,内存与计算双赢)
    int max_ls_iter = 8;        ///< 线搜索最大迭代次数(从10降至8)
    double grad_tol = 1e-3;     ///< 梯度收敛阈值(从5e-4放宽至1e-3)
    double func_tol = 1e-6;     ///< 函数值收敛阈值(从1e-7放宽至1e-6)
#endif
    double wolfe_c1 = 1e-4;     ///< Wolfe条件c1
    double wolfe_c2 = 0.9;      ///< Wolfe条件c2
};

/// 自适应细分结果
struct AdaptiveRefineResult {
    BezierCurve curve;
    bool was_split = false;
};

}
