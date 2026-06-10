/*
 * InfeasibilityDetector (无解预检测)
 *  - BFS 可达性（C 类：拓扑死锁）
 *  - 通道宽度检测（A 类：过窄）
 *  - 围栏-障碍夹断检测（B 类）
 */
#pragma once
#include "types.h"

namespace isg {

struct PreCheckResult {
    bool topological_block = false, narrow_passage = false, fence_sandwich = false;
    double min_corridor_width = 0.0;
    ViolationInfo::InfeasibilityType type = ViolationInfo::InfeasibilityType::None;
};

PreCheckResult preCheck(const SDFField& sdf, const Polygon2d& fence,
                        const Vec2d& entry_pt, const Vec2d& exit_pt,
                        double required_width, const std::vector<Boundary>& boundaries);
bool bfsReachable(const SDFField&, const Polygon2d&, const Vec2d&, const Vec2d&, double mc = 0.3);
bool detectSandwich(const SDFField&, const Polygon2d&, double cl = 0.5);

struct FenceRelaxResult {
    bool success = false;
    Polygon2d relaxed_fence;
    double max_expansion = 0;

    FenceRelaxResult(bool success, const Polygon2d& relaxed_fence, double max_expansion)
        : success(success), relaxed_fence(relaxed_fence) {}
};

FenceRelaxResult tryRelaxFence(const Polygon2d&, const std::vector<Boundary>&, const Vec2d&, double mx = 1.5);
ConnectivityCurve makeFallbackCurve(const PreCheckResult&, const Connectivity&, const Vec2d&, const Vec2d&);
bool isSimplePolygon(const Polygon2d&);
bool segmentsIntersect2(const Vec2d&, const Vec2d&, const Vec2d&, const Vec2d&);
// Overload with default empty boundary list
inline PreCheckResult preCheck(
    const SDFField& sdf, const Polygon2d& fence, const Vec2d& ep, const Vec2d& xp, double rw) {
    return preCheck(sdf, fence, ep, xp, rw, {});
}

}