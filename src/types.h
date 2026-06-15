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

// ── Geometry primitives ──────────────────────────────────────
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

struct LineString2d {
    std::vector<Vec2d> points;

    BoundingBox2d bbox() const {
        BoundingBox2d b;
        for (auto& p : points)
            b.expand(p);
        return b;
    }
};

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

// ── Bézier types (defined here to avoid circular deps) ───────
struct BezierSegment {
    std::array<Vec2d, 4> ctrl;

    Vec2d evaluate(double t) const; // B(t) 求值
    Vec2d tangent(double t) const;
    double curvature(double t) const; // 曲率 κ(t) = |B'×B''| / |B'|³
    std::pair<BezierSegment, BezierSegment> splitAt(double t) const;
    BoundingBox2d bbox() const;
    double arcLength(int samples = 20) const; // 弧长（数值积分 - 辛普森法）
    double arcLengthToParam(double s, int samples = 50) const;
    Vec2d evalDeriv1(double t) const; // B'(t) 一阶导数
    Vec2d evalDeriv2(double t) const; // B''(t) 二阶导数
    double maxCurvature(int samples = 50) const; // 最大曲率（采样估算）
    std::vector<Vec2d> sampleCount(int n) const; // 固定数量采样
    std::vector<Vec2d> sampleBySpacing(double spacing) const; // 固定间距采样（弧长近似）
    std::vector<Vec2d> sampleAdaptive(double maxAngleDeg, double maxSeg, double minSeg) const; // 自适应采样（基于角度误差）
    std::vector<Vec2d> sampleByShape(double straightSpacing, double angularResolution, double minSpacing) const;
    void subdivide(double t0, double t1, double maxAngle, double maxSeg, double minSeg,
                   std::vector<Vec2d>& pts, int depth) const;
    double getAlpha(const Vec2d& T0) const; // 获取 α（P1相对P0的标量，T0方向）
    double getBeta(const Vec2d& T3) const;
};

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
    std::vector<Vec2d> sampleByShape(double straightSpacing = 3.0, double angularResolution = 10.0, double minSpacing = 0.05) const;  // 形态自适应采样
    double arcLength() const;
    BoundingBox2d bbox() const;
    double maxCurvature(int sps = 20) const;
};

// ── Lane IDs ─────────────────────────────────────────────────
using LaneId = std::string;
using LaneGroupId = std::string;
using LaneEdgeId = std::string;
using ConnId = std::string;
using InterId = std::string;
using AttrMap = std::map<std::string, std::string>;

enum class GroupRole { Entry, Exit };

enum class ConnTurnType { Unknown = 0, TurnLeft = 1, UTurnLeft = 2, Straight = 3, TurnRight = 4, UTurnRight = 5 };

struct LaneEdge {
    LaneEdgeId id;
    LineString2d geometry;
    bool is_shared = false;
    std::shared_ptr<std::pair<LaneId, LaneId>> shared_by = nullptr;
    AttrMap attrs;

    // Vec2d connectionPt;
    // Vec2d tangentDir;
    std::string groupId;
    int lineOrder = 0; // 组内横向排序（0=最内侧）
};

struct Lane {
    LaneId id;
    LaneEdgeId left_edge_id;
    LaneEdgeId right_edge_id;
    double width = 3.5;
    LineString2d geometry;
    AttrMap attrs;

    // Vec2d connectionPt; // 连接点坐标
    // Vec2d tangentDir; // 连接点切线（指向路口内侧）
    std::string groupId;
    int laneOrder = 0; // 组内横向排序（0=最内侧）
};

struct LaneGroup {
    LaneGroupId id;
    GroupRole role;
    std::vector<LaneId> lanes;
    std::vector<LaneEdgeId> boundaries;
    AttrMap attrs;
};

struct Connectivity {
    ConnId id;
    LaneId entry_lane_id;
    LaneId exit_lane_id;
    ConnTurnType turn_type = ConnTurnType::Straight;

    LaneGroupId enterGroupId;
    LaneGroupId exitGroupId;
};

struct Obstacle {
    std::string id;
    Polygon2d geometry, buffered_geometry;
};

struct Boundary {
    enum class Type { RoadEdge, MedianStrip, GreenBelt, Other };

    std::string id;
    Type type = Type::RoadEdge;
    LineString2d geometry;
};

struct StopLine {
    std::string id;
    LineString2d geometry;
    LaneGroupId associated_group_id;
    Vec2d normal_direction{0, 1};
};

struct Crosswalk {
    std::string id;
    Polygon2d geometry;
    Vec2d crossing_direction{0, 1};
};

// 路口面
struct IntersectionArea {
    InterId id;
    Polygon2d geometry;
    bool is_rough = false; //true:粗糙路口面，false:精细路口面
};

enum class CurveStatus { OK, WarnA2, Degraded, Infeasible };

struct ViolationInfo {
    enum class InfeasibilityType { None, NarrowPassage, Sandwich, TopologicalBlock, ForcedCross };

    InfeasibilityType type = InfeasibilityType::None;
    double max_obstacle_penetration = 0, max_fence_overflow = 0, fence_expansion_applied = 0;
    std::vector<Vec2d> exempt_crosses;
    std::string reason;
};

struct ConnectivityCurve {
    ConnId id;
    LaneId entry_lane_id;
    LaneId exit_lane_id;
    ConnTurnType turn_type = ConnTurnType::Straight;
    std::shared_ptr<BezierCurve> curve = nullptr; // BezierCurve is complete above
    CurveStatus status = CurveStatus::OK;
    ViolationInfo violation;

    LaneEdgeId left_edge_id = "";
    LaneEdgeId right_edge_id = "";
};

struct ConnectivityLaneEdge {
    LaneEdgeId id;
    LineString2d geometry;
    bool is_shared = false;
    std::shared_ptr<std::pair<LaneId, LaneId>> shared_by = nullptr; //边线左侧车道，边线右侧车道

    AttrMap attrs;
    // Vec2d connectionPt;
    // Vec2d tangentDir;
    std::string groupId;
    int lineOrder = 0; // 组内横向排序（0=最内侧）
};

struct IntersectionOutput {
    std::vector<ConnectivityCurve> connectivity_curves;
    std::vector<ConnectivityLaneEdge> lane_edges;
    IntersectionArea area;

    struct PerfStats {
        double sdf_build_ms = 0, precheck_ms = 0, optimize_ms = 0,
               smooth_ms = 0, edge_gen_ms = 0, area_gen_ms = 0;
    } perf;
};


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

    std::string id; //路口ID

    const bool IsEntryLaneEdge(const LaneEdgeId& id) const;
    const bool IsEntryLane(const LaneId& id) const;
    const Lane* findLane(const LaneId&) const;
    const LaneGroup* findGroup(const LaneGroupId&) const;
    const LaneEdge* findEdge(const LaneEdgeId&) const;
    bool laneGroupExists(const LaneGroupId&) const;
    std::pair<Vec2d, Vec2d> entryPtDir(const LaneId&) const;
    std::pair<Vec2d, Vec2d> exitPtDir(const LaneId&) const;
};

// ── AdaptiveRefineResult forward-declared here ───────────────
struct SDFField; // forward

enum class ConnectivityDirectionMode {
    PerLane = 0,
    GroupUnified = 1
};

struct ConnectivityDirectionConfig {
    ConnectivityDirectionMode mode = ConnectivityDirectionMode::PerLane;
    double group_similarity_angle_deg = 5.0;
};

struct LBFGSConfig {
    int max_iter = 80;
    int history_size = 10;
    int max_ls_iter = 10;
    double grad_tol = 5e-4;
    double func_tol = 1e-7;
    double wolfe_c1 = 1e-4;
    double wolfe_c2 = 0.9;
};

struct AdaptiveRefineResult {
    BezierCurve curve;
    bool was_split = false;
};

}
