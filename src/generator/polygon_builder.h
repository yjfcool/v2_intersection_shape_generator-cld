#pragma once
#include "types.h"
#include "utils.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>

namespace isg {

/**
 * 精细路口面构建
 * 处理：道路边缘线（独立线段）+ 进入组(或退出组)尾端(或首端)边线回缩线(外扩线) / 停止线回缩线;
 *      如果缩扩延长线与道路边缘相交时需截取有效段作为路口边线;
 * 输出：封闭可凹陷多边形（路口面）
 */
struct EdgeEndpoint {
    Vec2d   pt;
    int     roadEdgeIdx = -1;
    bool    isStart = true;

    EdgeEndpoint(Vec2d _pt, int _idx, bool _start):pt(_pt), roadEdgeIdx(_idx), isStart(_start) {}
};
class IntersectionAreaBuilder {
    // 路口面顶点吸附容差(米),距离小于此值的顶点会被合并。取值范围：0.1~2.0
    double snapTolerance = 0.5;
    // 路口面多边形绕向 (逆时针:"counterclockwise" | 顺时针:"clockwise")
    std::string winding = "clockwise";
    // 进入组端侧边界向路口内回缩、退出组端侧边界向路口外扩的距离(米)
    double groupCutOffset = 0.05;
    // 端侧边界折线两端用于对齐道路边缘线的外延长度(米)
    double cutLineExtension = 0.20;

public:
    explicit IntersectionAreaBuilder(
        double _snapTolerance = 0.5,
        std::string _winding = "clockwise",
        double _groupCutOffset = 0.05,
        double _cutLineExtension = 0.20)
        : snapTolerance(_snapTolerance),
          winding(_winding),
          groupCutOffset(std::max(0.0, _groupCutOffset)),
          cutLineExtension(std::max(0.0, _cutLineExtension)) {}

    IntersectionArea build(
        const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines,
        const std::vector<ConnectivityLaneEdge>& edgelines)
    {
        IntersectionArea junctArea;
        junctArea.id = "int_polygon_" + inp.id;

        // 1. 确定路口中心和范围
        Vec2d center = computeCenter(inp, centerlines);
        double radius = computeRadius(inp, centerlines, center);

        std::vector<Vec2d> finePoly = buildFinePolygon(inp, centerlines, center, radius);
        if (finePoly.size() >= 4) {
            junctArea.geometry.outer = finePoly;
            junctArea.is_rough = false;
            return junctArea;
        }

        // 2. 裁剪道路边缘线，获取路口侧端点集
        std::vector<EdgeEndpoint> endpoints;

        for (int i = 0; i < inp.boundaries.size(); ++i) {
            const auto& re = inp.boundaries[i];
            if (re.geometry.points.size() < 2) continue;

            // 找到路口侧的端点（靠近center的端点）
            double d0 = dist(re.geometry.points.front(), center);
            double d1 = dist(re.geometry.points.back(), center);

            Vec2d candidatePt;
            bool isStart;
            if (d0 < d1) {
                candidatePt = re.geometry.points.front();
                isStart = true;
            } else {
                candidatePt = re.geometry.points.back();
                isStart = false;
            }

            // 检查是否在路口范围内（在 rough_area 内或在 radius 内）
            bool inRange = false;
            if (!inp.area.geometry.empty() && !inp.area.geometry.outer.empty()) {
                inRange = pointInPolygon(candidatePt, inp.area.geometry.outer) ||
                    dist(candidatePt, center) < radius * 1.5;
            } else {
                inRange = dist(candidatePt, center) < radius * 1.5;
            }
            if (inRange) {
                // 吸附到最近的车道边线端点
                Vec2d snapped = snapToEdgePt(candidatePt, edgelines, inp, snapTolerance);
                endpoints.emplace_back(EdgeEndpoint{snapped, i, isStart});
            }
        }

        // 若道路边缘线端点不足，从连通关系连接点的外侧边线端点补充
        if (endpoints.size() < 3) {
            supplementEndpointsFromLaneEdges(endpoints, inp, edgelines, center, radius);
        }

        if (endpoints.size() < 3) {
            // 终极回退：用所有连接点的凸包
            junctArea.geometry.outer = buildConvexHullFallback(inp, centerlines);
            junctArea.is_rough = false;
            return junctArea;
        }

        // 3. 按极角排序
        std::sort(endpoints.begin(), endpoints.end(),
                  [&](const EdgeEndpoint& a, const EdgeEndpoint& b) {
                      double angA = std::atan2(a.pt[1] - center[1], a.pt[0] - center[0]);
                      double angB = std::atan2(b.pt[1] - center[1], b.pt[0] - center[0]);
                      if (winding == "clockwise")
                          return angA > angB;
                      else
                          return angA < angB;
                  });

        // 去重（距离过近的端点合并）
        deduplicateEndpoints(endpoints, 0.1);

        // 4. 连接端点构建多边形
        std::vector<Vec2d> rawPoly;
        for (int i = 0; i < (int)endpoints.size(); ++i) {
            const EdgeEndpoint& cur = endpoints[i];
            const EdgeEndpoint& next = endpoints[(i + 1) % endpoints.size()];

            rawPoly.push_back(cur.pt);

            // 若相邻端点来自同一条道路边缘线，沿边缘线几何连接
            if (cur.roadEdgeIdx >= 0 && cur.roadEdgeIdx == next.roadEdgeIdx) {
                std::vector<Vec2d> seg = extractEdgeSegment(
                    inp.boundaries[cur.roadEdgeIdx],
                    cur.pt, next.pt,
                    snapTolerance);
                for (auto& p : seg) rawPoly.push_back(p);
            }
            // 否则直连（或弧线连，按参数）
        }

        // 5. 闭合多边形
        if (!rawPoly.empty() && dist(rawPoly.front(), rawPoly.back()) > EPS)
            rawPoly.push_back(rawPoly.front());

        // 6. 多边形修复（移除共线点，确保逆/顺时针）
        junctArea.geometry.outer = repairPolygon(rawPoly, center);

        junctArea.is_rough = false;
        return junctArea;
    }

private:
    struct AreaChain {
        std::vector<Vec2d> pts;
    };

    struct BreakLineInfo {
        const StopLine* stop_line = nullptr;
        std::vector<Vec2d> raw_pts;
        std::vector<Vec2d> shifted_pts;
        std::vector<LaneGroupId> cut_entry_group_ids;
    };

    struct LineIntersection {
        Vec2d pt{0, 0};
        double station_a = 0.0;
        double station_b = 0.0;
    };

    struct PolylineProjection {
        Vec2d pt{0, 0};
        double station = 0.0;
        double distance = 1e18;
    };

    struct ClosestPolylineProjection {
        Vec2d pt_on_a{0, 0};
        Vec2d pt_on_b{0, 0};
        double station_a = 0.0;
        double station_b = 0.0;
        double distance = 1e18;
    };

    static std::vector<Vec2d> openRing(std::vector<Vec2d> ring) {
        while (ring.size() > 1 && dist(ring.front(), ring.back()) < 1e-6)
            ring.pop_back();
        return ring;
    }

    void addPointChain(
        std::vector<AreaChain>& chains, const Vec2d& pt, double tol = 0.05) const {
        for (const auto& c : chains) {
            if (c.pts.size() == 1 && dist(c.pts.front(), pt) <= tol)
                return;
        }
        AreaChain c;
        c.pts.push_back(pt);
        chains.push_back(c);
    }

    Vec2d stopLineInwardDir(const StopLine& sl, const Vec2d& center) const {
        const auto& pts = sl.geometry.points;
        if (pts.empty())
            return Vec2d(0, 0);
        Vec2d mid(0, 0);
        for (const auto& pt : pts)
            mid += pt;
        mid /= (double)pts.size();
        Vec2d to_center = center - mid;
        if (to_center.norm() < 1e-8)
            to_center = Vec2d(1, 0);

        Vec2d geom_normal(0, 0);
        if (pts.size() >= 2) {
            Vec2d tangent = pts.back() - pts.front();
            if (tangent.norm() > 1e-8) {
                tangent.normalize();
                geom_normal = rotLeft(tangent);
                if (geom_normal.dot(to_center) < 0.0)
                    geom_normal = -geom_normal;
            }
        }

        if (sl.normal_direction.norm() > 1e-8) {
            Vec2d normal = sl.normal_direction.normalized();
            bool normal_matches_line = true;
            if (pts.size() >= 2) {
                Vec2d tangent = pts.back() - pts.front();
                if (tangent.norm() > 1e-8) {
                    tangent.normalize();
                    normal_matches_line = std::abs(normal.dot(tangent)) < 0.35;
                }
            }
            if (normal_matches_line) {
                if (normal.dot(to_center) < 0.0)
                    normal = -normal;
                return normal;
            }
        }

        if (geom_normal.norm() > 1e-8)
            return geom_normal.normalized();
        return to_center.normalized();
    }

    std::vector<Vec2d> shiftedStopLine(const StopLine& sl, const Vec2d& center) const {
        std::vector<Vec2d> shifted;
        Vec2d dir = stopLineInwardDir(sl, center);
        if (dir.norm() < 1e-8)
            return shifted;
        shifted.reserve(sl.geometry.points.size());
        for (const auto& pt : sl.geometry.points)
            shifted.push_back(pt + groupCutOffset * dir.normalized());
        return shifted;
    }

    bool segmentIntersection(
        const Vec2d& a, const Vec2d& b, const Vec2d& c, const Vec2d& d,
        Vec2d* out = nullptr, double* ta = nullptr, double* tb = nullptr) const {
        Vec2d r = b - a;
        Vec2d s = d - c;
        double den = cross2d(r, s);
        if (std::abs(den) < 1e-12)
            return false;
        Vec2d ac = c - a;
        double t = cross2d(ac, s) / den;
        double u = cross2d(ac, r) / den;
        if (t < -1e-9 || t > 1.0 + 1e-9 || u < -1e-9 || u > 1.0 + 1e-9)
            return false;
        t = std::max(0.0, std::min(1.0, t));
        u = std::max(0.0, std::min(1.0, u));
        if (out)
            *out = a + t * r;
        if (ta)
            *ta = t;
        if (tb)
            *tb = u;
        return true;
    }

    double polylineLength(const std::vector<Vec2d>& pts) const {
        double len = 0.0;
        for (int i = 0; i + 1 < (int)pts.size(); ++i)
            len += dist(pts[i], pts[i + 1]);
        return len;
    }

    double stopLineExtensionLength(const IntersectionInput& inp, const StopLine& sl) const {
        BoundingBox2d box;
        for (const auto& lane : inp.lanes)
            for (const auto& pt : lane.geometry.points)
                box.expand(pt);
        for (const auto& bnd : inp.boundaries)
            for (const auto& pt : bnd.geometry.points)
                box.expand(pt);
        for (const auto& pt : inp.area.geometry.outer)
            box.expand(pt);

        double stop_len = polylineLength(sl.geometry.points);
        double diag = box.empty() ? 0.0 : (box.max_pt - box.min_pt).norm();
        double by_context = diag > 1e-6 ? diag * 0.08 : stop_len * 2.0;
        return std::max(6.0, std::min(30.0, std::max(stop_len, by_context)));
    }

    std::vector<Vec2d> extendedStopLine(const IntersectionInput& inp, const StopLine& sl) const {
        std::vector<Vec2d> line = sl.geometry.points;
        if (line.size() < 2)
            return line;

        double extension = stopLineExtensionLength(inp, sl);
        Vec2d first_dir = line[1] - line[0];
        if (first_dir.norm() > 1e-8) {
            first_dir.normalize();
            line.front() -= first_dir * extension;
        }
        Vec2d last_dir = line.back() - line[line.size() - 2];
        if (last_dir.norm() > 1e-8) {
            last_dir.normalize();
            line.back() += last_dir * extension;
        }
        return line;
    }

    Vec2d pointAtStation(const std::vector<Vec2d>& pts, double station) const {
        if (pts.empty())
            return Vec2d(0, 0);
        if (pts.size() == 1)
            return pts.front();
        double total = polylineLength(pts);
        station = std::max(0.0, std::min(total, station));
        double acc = 0.0;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            double seg_len = dist(pts[i], pts[i + 1]);
            if (seg_len < 1e-9)
                continue;
            if (acc + seg_len >= station) {
                double t = (station - acc) / seg_len;
                return pts[i] * (1.0 - t) + pts[i + 1] * t;
            }
            acc += seg_len;
        }
        return pts.back();
    }

    std::vector<Vec2d> subPolylineByStation(
        const std::vector<Vec2d>& pts, double s0, double s1) const {
        std::vector<Vec2d> out;
        if (pts.size() < 2)
            return pts;
        double total = polylineLength(pts);
        s0 = std::max(0.0, std::min(total, s0));
        s1 = std::max(0.0, std::min(total, s1));
        if (s1 < s0)
            std::swap(s0, s1);
        pushUnique(out, pointAtStation(pts, s0));
        double acc = 0.0;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            double seg_len = dist(pts[i], pts[i + 1]);
            double next = acc + seg_len;
            if (seg_len > 1e-9 && next > s0 + 1e-9 && next < s1 - 1e-9)
                pushUnique(out, pts[i + 1]);
            acc = next;
        }
        pushUnique(out, pointAtStation(pts, s1));
        return out;
    }

    std::vector<Vec2d> subPolylineByStations(
        const std::vector<Vec2d>& pts,
        double s0,
        double s1,
        std::vector<double> stations) const {
        std::vector<Vec2d> out;
        if (pts.size() < 2)
            return pts;

        double total = polylineLength(pts);
        s0 = std::max(0.0, std::min(total, s0));
        s1 = std::max(0.0, std::min(total, s1));
        if (s1 < s0)
            std::swap(s0, s1);

        stations.push_back(s0);
        stations.push_back(s1);
        double acc = 0.0;
        for (int i = 0; i + 1 < (int)pts.size(); ++i) {
            double seg_len = dist(pts[i], pts[i + 1]);
            acc += seg_len;
            if (seg_len > 1e-9 && acc > s0 + 1e-9 && acc < s1 - 1e-9)
                stations.push_back(acc);
        }

        std::sort(stations.begin(), stations.end());
        for (double s : stations) {
            if (s < s0 - 1e-9 || s > s1 + 1e-9)
                continue;
            if (!out.empty() && std::abs(s - stations.front()) < 1e-12)
                continue;
            pushUnique(out, pointAtStation(pts, s));
        }
        return out;
    }

    PolylineProjection projectPointToPolyline(
        const Vec2d& pt, const std::vector<Vec2d>& line) const {
        PolylineProjection best;
        double acc = 0.0;
        for (int i = 0; i + 1 < (int)line.size(); ++i) {
            auto pr = pointToSegment(pt, line[i], line[i + 1]);
            double seg_len = dist(line[i], line[i + 1]);
            Vec2d q = line[i] * (1.0 - pr.second) + line[i + 1] * pr.second;
            if (pr.first < best.distance) {
                best.pt = q;
                best.station = acc + pr.second * seg_len;
                best.distance = pr.first;
            }
            acc += seg_len;
        }
        return best;
    }

    ClosestPolylineProjection closestPolylineProjection(
        const std::vector<Vec2d>& a, const std::vector<Vec2d>& b) const {
        ClosestPolylineProjection best;
        if (a.size() < 2 || b.size() < 2)
            return best;

        for (const auto& hit : polylineIntersections(a, b)) {
            best.pt_on_a = hit.pt;
            best.pt_on_b = hit.pt;
            best.station_a = hit.station_a;
            best.station_b = hit.station_b;
            best.distance = 0.0;
            return best;
        }

        double acc_a = 0.0;
        for (int i = 0; i + 1 < (int)a.size(); ++i) {
            double seg_a = dist(a[i], a[i + 1]);
            if (seg_a < 1e-9) {
                acc_a += seg_a;
                continue;
            }
            double acc_b = 0.0;
            for (int j = 0; j + 1 < (int)b.size(); ++j) {
                double seg_b = dist(b[j], b[j + 1]);
                if (seg_b < 1e-9) {
                    acc_b += seg_b;
                    continue;
                }

                auto check_a_to_b = [&](const Vec2d& pa, double sa) {
                    auto pr = pointToSegment(pa, b[j], b[j + 1]);
                    Vec2d qb = b[j] * (1.0 - pr.second) + b[j + 1] * pr.second;
                    if (pr.first < best.distance) {
                        best.pt_on_a = pa;
                        best.pt_on_b = qb;
                        best.station_a = sa;
                        best.station_b = acc_b + pr.second * seg_b;
                        best.distance = pr.first;
                    }
                };
                auto check_b_to_a = [&](const Vec2d& pb, double sb) {
                    auto pr = pointToSegment(pb, a[i], a[i + 1]);
                    Vec2d qa = a[i] * (1.0 - pr.second) + a[i + 1] * pr.second;
                    if (pr.first < best.distance) {
                        best.pt_on_a = qa;
                        best.pt_on_b = pb;
                        best.station_a = acc_a + pr.second * seg_a;
                        best.station_b = sb;
                        best.distance = pr.first;
                    }
                };

                check_a_to_b(a[i], acc_a);
                check_a_to_b(a[i + 1], acc_a + seg_a);
                check_b_to_a(b[j], acc_b);
                check_b_to_a(b[j + 1], acc_b + seg_b);

                acc_b += seg_b;
            }
            acc_a += seg_a;
        }
        return best;
    }

    std::vector<LineIntersection> polylineIntersections(
        const std::vector<Vec2d>& a, const std::vector<Vec2d>& b) const {
        std::vector<LineIntersection> out;
        double acc_a = 0.0;
        for (int i = 0; i + 1 < (int)a.size(); ++i) {
            double seg_a = dist(a[i], a[i + 1]);
            double acc_b = 0.0;
            for (int j = 0; j + 1 < (int)b.size(); ++j) {
                double seg_b = dist(b[j], b[j + 1]);
                Vec2d pt;
                double ta = 0.0, tb = 0.0;
                if (segmentIntersection(a[i], a[i + 1], b[j], b[j + 1], &pt, &ta, &tb)) {
                    LineIntersection hit;
                    hit.pt = pt;
                    hit.station_a = acc_a + ta * seg_a;
                    hit.station_b = acc_b + tb * seg_b;
                    bool duplicate = false;
                    for (const auto& existing : out) {
                        if (dist(existing.pt, hit.pt) < 0.03) {
                            duplicate = true;
                            break;
                        }
                    }
                    if (!duplicate)
                        out.push_back(hit);
                }
                acc_b += seg_b;
            }
            acc_a += seg_a;
        }
        return out;
    }

    struct SegmentProjection {
        int index = 0;
        double t = 0.0;
        double distance = 1e18;
    };

    SegmentProjection closestRingSegment(const std::vector<Vec2d>& ring, const Vec2d& pt) const {
        SegmentProjection best;
        int n = (int)ring.size();
        if (n < 2)
            return best;
        for (int i = 0; i < n; ++i) {
            const Vec2d& a = ring[i];
            const Vec2d& b = ring[(i + 1) % n];
            auto pr = pointToSegment(pt, a, b);
            if (pr.first < best.distance) {
                best.index = i;
                best.t = pr.second;
                best.distance = pr.first;
            }
        }
        return best;
    }

    void pushUnique(std::vector<Vec2d>& out, const Vec2d& pt, double tol = 0.03) const {
        if (out.empty() || dist(out.back(), pt) > tol)
            out.push_back(pt);
    }

    bool pointOnChain(const Vec2d& pt, const AreaChain& chain, double tol = 0.05) const {
        if (chain.pts.empty())
            return false;
        for (const auto& q : chain.pts)
            if (dist(pt, q) <= tol)
                return true;
        for (int i = 0; i + 1 < (int)chain.pts.size(); ++i)
            if (pointToSegment(pt, chain.pts[i], chain.pts[i + 1]).first <= tol)
                return true;
        return false;
    }

    bool pointOnAnyChain(const Vec2d& pt, const std::vector<AreaChain>& chains, double tol = 0.05) const {
        for (const auto& chain : chains)
            if (pointOnChain(pt, chain, tol))
                return true;
        return false;
    }

    double angleOf(const Vec2d& pt, const Vec2d& center) const {
        return std::atan2(pt[1] - center[1], pt[0] - center[0]);
    }

    double normalizedAngleDelta(double a0, double a1) const {
        double d = a1 - a0;
        while (d > M_PI)
            d -= 2.0 * M_PI;
        while (d < -M_PI)
            d += 2.0 * M_PI;
        return d;
    }

    AreaChain orientChain(AreaChain chain, const Vec2d& center) const {
        if (chain.pts.size() < 2)
            return chain;
        double a0 = angleOf(chain.pts.front(), center);
        double a1 = angleOf(chain.pts.back(), center);
        double delta = normalizedAngleDelta(a0, a1);
        bool should_reverse = (winding == "clockwise") ? (delta > 0.0) : (delta < 0.0);
        if (should_reverse)
            std::reverse(chain.pts.begin(), chain.pts.end());
        return chain;
    }

    Vec2d chainMidpoint(const AreaChain& chain) const {
        Vec2d mid(0, 0);
        if (chain.pts.empty())
            return mid;
        for (const auto& pt : chain.pts)
            mid += pt;
        return mid / (double)chain.pts.size();
    }

    const LaneGroup* findEntryGroupByLane(const IntersectionInput& inp, const LaneId& lane_id) const {
        for (const auto& group : inp.lane_groups) {
            if (group.role != GroupRole::Entry)
                continue;
            if (std::find(group.lanes.begin(), group.lanes.end(), lane_id) != group.lanes.end())
                return &group;
        }
        return nullptr;
    }

    bool laneInGroup(const LaneGroup& group, const LaneId& lane_id) const {
        return std::find(group.lanes.begin(), group.lanes.end(), lane_id) != group.lanes.end();
    }

    bool groupContainsLaneId(const LaneGroup& group, const LaneId& lane_id) const {
        return laneInGroup(group, lane_id);
    }

    std::vector<LaneId> uniqueGroupLaneIds(const LaneGroup& group) const {
        std::vector<LaneId> ids;
        for (const auto& lane_id : group.lanes) {
            if (std::find(ids.begin(), ids.end(), lane_id) == ids.end())
                ids.push_back(lane_id);
        }
        return ids;
    }

    bool containsLaneGroupId(const std::vector<LaneGroupId>& ids, const LaneGroupId& id) const {
        return std::find(ids.begin(), ids.end(), id) != ids.end();
    }

    bool containsLaneId(const std::vector<LaneId>& ids, const LaneId& id) const {
        return std::find(ids.begin(), ids.end(), id) != ids.end();
    }

    void addUniqueLaneId(std::vector<LaneId>& ids, const LaneId& id) const {
        if (!id.empty() && !containsLaneId(ids, id))
            ids.push_back(id);
    }

    bool laneInAnyGroup(const std::vector<const LaneGroup*>& groups, const LaneId& lane_id) const {
        for (const auto* group : groups) {
            if (group && groupContainsLaneId(*group, lane_id))
                return true;
        }
        return false;
    }

    void addLaneGroupId(std::vector<LaneGroupId>& ids, const LaneGroupId& id) const {
        if (!id.empty() && !containsLaneGroupId(ids, id))
            ids.push_back(id);
    }

    std::vector<const LaneGroup*> stopLineEntryGroups(
        const IntersectionInput& inp, const StopLine& sl) const {
        std::vector<const LaneGroup*> groups;
        auto add_group = [&](const LaneGroup* group) {
            if (!group || group->role != GroupRole::Entry)
                return;
            for (const auto* existing : groups)
                if (existing->id == group->id)
                    return;
            groups.push_back(group);
        };

        if (!sl.associated_group_id.empty())
            add_group(inp.findGroup(sl.associated_group_id));

        auto collect_groups = [&](const std::vector<Vec2d>& probe_line, double cover_tol) {
            for (const auto& lane : inp.lanes) {
                if (lane.geometry.points.size() < 2)
                    continue;
                auto hits = polylineIntersections(lane.geometry.points, probe_line);
                bool covered = !hits.empty();
                if (!covered) {
                    covered = closestPolylineProjection(probe_line, lane.geometry.points).distance <= cover_tol;
                }
                if (covered)
                    add_group(findEntryGroupByLane(inp, lane.id));
            }
        };

        collect_groups(sl.geometry.points, 0.45);
        if (groups.empty())
            collect_groups(extendedStopLine(inp, sl), 0.25);
        return groups;
    }

    std::vector<Vec2d> trimBreakLineByBoundaries(
        const std::vector<Vec2d>& line, const std::vector<Boundary>& boundaries) const {
        if (line.size() < 2 || boundaries.empty())
            return line;
        double total = polylineLength(line);
        double min_s = 0.0;
        double max_s = total;
        std::vector<double> stations;
        for (const auto& bnd : boundaries) {
            auto hits = polylineIntersections(line, bnd.geometry.points);
            for (const auto& hit : hits)
                stations.push_back(hit.station_a);
        }
        if (stations.empty())
            return line;
        std::sort(stations.begin(), stations.end());
        if (stations.size() >= 2) {
            min_s = stations.front();
            max_s = stations.back();
        } else {
            double s = stations.front();
            if (s < total * 0.5)
                min_s = s;
            else
                max_s = s;
        }
        if (max_s - min_s < 0.2)
            return line;
        return subPolylineByStation(line, min_s, max_s);
    }

    std::vector<Vec2d> shiftedBreakLine(
        const StopLine& sl, const std::vector<Vec2d>& break_line, const Vec2d& center) const {
        std::vector<Vec2d> shifted;
        Vec2d dir = stopLineInwardDir(sl, center);
        if (dir.norm() < 1e-8)
            return shifted;
        shifted.reserve(break_line.size());
        for (const auto& pt : break_line)
            shifted.push_back(pt + groupCutOffset * dir.normalized());
        return shifted;
    }

    std::vector<Vec2d> buildStopLineBreakLineRaw(
        const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines,
        const StopLine& sl) const {
        if (sl.geometry.points.size() < 2)
            return {};

        constexpr double COVER_TOL = 0.45;
        constexpr double CONN_TOL = 1.2;
        std::vector<Vec2d> cut_line = extendedStopLine(inp, sl);
        if (cut_line.size() < 2)
            return {};

        double total = polylineLength(cut_line);
        double lane_half_width = 0.25;
        std::vector<double> hit_stations;
        std::vector<double> fallback_stations;
        auto add_station = [&](std::vector<double>& out, double s) {
            s = std::max(0.0, std::min(total, s));
            for (double existing : out)
                if (std::abs(existing - s) < 0.05)
                    return;
            out.push_back(s);
        };

        auto groups = stopLineEntryGroups(inp, sl);
        for (const auto* group : groups) {
            for (const auto& lane_id : uniqueGroupLaneIds(*group)) {
                const Lane* lane = inp.findLane(lane_id);
                if (!lane)
                    continue;
                lane_half_width = std::max(lane_half_width, lane->width * 0.5);
                auto lane_hits = polylineIntersections(cut_line, lane->geometry.points);
                for (const auto& hit : lane_hits)
                    add_station(hit_stations, hit.station_a);
                if (lane_hits.empty()) {
                    auto closest = closestPolylineProjection(cut_line, lane->geometry.points);
                    if (closest.distance <= COVER_TOL)
                        add_station(fallback_stations, closest.station_a);
                }
                if (lane->geometry.points.size() >= 2) {
                    Vec2d conn = getConnPoint(lane->geometry.points, true);
                    auto pr = projectPointToPolyline(conn, cut_line);
                    if (pr.distance <= CONN_TOL)
                        add_station(fallback_stations, pr.station);
                }
            }
            for (const auto& cc : centerlines) {
                if (!cc.curve || !groupContainsLaneId(*group, cc.entry_lane_id))
                    continue;
                auto pts = cc.curve->sampleByArcLength(80);
                auto curve_hits = polylineIntersections(cut_line, pts);
                for (const auto& hit : curve_hits)
                    add_station(hit_stations, hit.station_a);
                if (curve_hits.empty()) {
                    auto closest = closestPolylineProjection(cut_line, pts);
                    if (closest.distance <= COVER_TOL)
                        add_station(fallback_stations, closest.station_a);
                }
                auto pr0 = projectPointToPolyline(cc.curve->startPt(), cut_line);
                if (pr0.distance <= CONN_TOL)
                    add_station(fallback_stations, pr0.station);
            }
        }

        std::vector<double> stations = hit_stations;
        if (stations.size() < 2) {
            for (double s : fallback_stations)
                add_station(stations, s);
        }

        if (stations.empty())
            return trimBreakLineByBoundaries(cut_line, inp.boundaries);

        std::sort(stations.begin(), stations.end());
        double pad = hit_stations.size() >= 2 ? 0.0 : std::max(0.25, std::min(2.0, lane_half_width));
        double s0 = std::max(0.0, stations.front() - pad);
        double s1 = std::min(total, stations.back() + pad);
        if (s1 - s0 < 0.5) {
            double mid = 0.5 * (s0 + s1);
            s0 = std::max(0.0, mid - 0.25);
            s1 = std::min(total, mid + 0.25);
        }
        auto break_line = subPolylineByStations(cut_line, s0, s1, stations);
        break_line = trimBreakLineByBoundaries(break_line, inp.boundaries);
        return break_line;
    }

    BreakLineInfo buildStopLineBreakLine(
        const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines,
        const StopLine& sl,
        const Vec2d& center) const {
        BreakLineInfo info;
        info.stop_line = &sl;
        info.raw_pts = buildStopLineBreakLineRaw(inp, centerlines, sl);
        info.shifted_pts = shiftedBreakLine(sl, info.raw_pts, center);
        std::vector<const LaneGroup*> groups = stopLineEntryGroups(inp, sl);
        std::vector<Vec2d> cut_line = extendedStopLine(inp, sl);
        if (cut_line.size() >= 2) {
            for (const auto* group : groups) {
                if (!group)
                    continue;
                bool has_intersection = false;
                for (const auto& lane_id : uniqueGroupLaneIds(*group)) {
                    const Lane* lane = inp.findLane(lane_id);
                    if (lane && !polylineIntersections(cut_line, lane->geometry.points).empty()) {
                        has_intersection = true;
                        break;
                    }
                }
                if (!has_intersection) {
                    for (const auto& cc : centerlines) {
                        if (!cc.curve || !groupContainsLaneId(*group, cc.entry_lane_id))
                            continue;
                        auto pts = cc.curve->sampleByArcLength(100);
                        if (!polylineIntersections(cut_line, pts).empty()) {
                            has_intersection = true;
                            break;
                        }
                    }
                }
                if (has_intersection)
                    addLaneGroupId(info.cut_entry_group_ids, group->id);
            }
        }
        return info;
    }

    std::vector<Vec2d> trimBoundaryByBreakLines(
        const Boundary& bnd, const std::vector<BreakLineInfo>& break_lines, const Vec2d& center) const {
        const auto& pts = bnd.geometry.points;
        if (pts.size() < 2 || break_lines.empty())
            return pts;
        double total = polylineLength(pts);
        std::vector<double> stations;
        for (const auto& br : break_lines) {
            for (const auto& hit : polylineIntersections(pts, br.raw_pts))
                stations.push_back(hit.station_a);
        }
        if (stations.empty())
            return pts;
        std::sort(stations.begin(), stations.end());
        double best_s0 = 0.0;
        double best_s1 = total;
        double best_score = 1e18;
        auto consider_segment = [&](double a, double b) {
            if (b - a < 0.2)
                return;
            Vec2d mid = pointAtStation(pts, 0.5 * (a + b));
            double score = dist(mid, center);
            if (score < best_score) {
                best_score = score;
                best_s0 = a;
                best_s1 = b;
            }
        };
        if (stations.size() >= 2) {
            consider_segment(stations.front(), stations.back());
        } else {
            double s = stations.front();
            consider_segment(0.0, s);
            consider_segment(s, total);
        }
        if (best_score == 1e18)
            return pts;
        return subPolylineByStation(pts, best_s0, best_s1);
    }

    // 当前需求改为按进入/退出组端侧折线生成。旧方案：道路边缘 + 停止线回缩线，注释保留了调用
    std::vector<Vec2d> buildFromBoundaryElements(
        const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines,
        const Vec2d& center) const {
        std::vector<AreaChain> chains;
        std::vector<BreakLineInfo> break_lines;

        for (const auto& sl : inp.stop_lines) {
            BreakLineInfo info = buildStopLineBreakLine(inp, centerlines, sl, center);
            if (info.shifted_pts.size() >= 2)
                break_lines.push_back(info);
        }

        for (const auto& bnd : inp.boundaries) {
            if (bnd.geometry.points.empty())
                continue;
            AreaChain chain;
            chain.pts = trimBoundaryByBreakLines(bnd, break_lines, center);
            chains.push_back(orientChain(chain, center));
        }
        for (const auto& br : break_lines) {
            AreaChain chain;
            chain.pts = br.shifted_pts;
            if (chain.pts.size() >= 2)
                chains.push_back(orientChain(chain, center));
        }
        for (const auto& cc : centerlines) {
            if (!cc.curve)
                continue;
            bool cut_entry_endpoint = false;
            for (const auto& br : break_lines) {
                const LaneGroup* group = findEntryGroupByLane(inp, cc.entry_lane_id);
                if (group && containsLaneGroupId(br.cut_entry_group_ids, group->id)) {
                    cut_entry_endpoint = true;
                    break;
                }
            }
            if (!cut_entry_endpoint && !pointOnAnyChain(cc.curve->startPt(), chains))
                addPointChain(chains, cc.curve->startPt());
            if (!pointOnAnyChain(cc.curve->endPt(), chains))
                addPointChain(chains, cc.curve->endPt());
        }

        if (chains.empty())
            return {};
        std::sort(chains.begin(), chains.end(), [&](const AreaChain& a, const AreaChain& b) {
            double angA = angleOf(a.pts.front(), center);
            double angB = angleOf(b.pts.front(), center);
            return winding == "clockwise" ? angA > angB : angA < angB;
        });

        std::vector<Vec2d> raw;
        for (const auto& chain : chains)
            for (const auto& pt : chain.pts)
                pushUnique(raw, pt);
        if (raw.size() < 3)
            return {};
        return repairPolygon(raw, center);
    }

    struct GroupCutPoint {
        Vec2d pt{0, 0};
        double lateral = 0.0;
        int order = 0;
        bool has_order = false;
    };

    // 统一车道/边线在组端侧的连接方向：进入组指向路口内，退出组指向路口外。
    // 先取连接端切向，再用路口中心校正方向，避免折线整体平移反向。
    Vec2d directedConnDir(
        const std::vector<Vec2d>& pts, GroupRole role, const Vec2d& center) const {
        bool is_entry = (role == GroupRole::Entry);
        Vec2d p = getConnPoint(pts, is_entry);
        Vec2d d = getConnTangent(pts, is_entry);
        Vec2d wanted = is_entry ? (center - p) : (p - center);
        if (wanted.norm() > 1e-8 && d.dot(wanted) < 0.0)
            d = -d;
        if (d.norm() < 1e-8)
            return Vec2d(1, 0);
        return d.normalized();
    }

    // 收集车道边线 id 时去重，避免 group.boundaries 与车道左右边线重复。
    void addUniqueLaneEdgeId(std::vector<LaneEdgeId>& ids, const LaneEdgeId& id) const {
        if (!id.empty() && std::find(ids.begin(), ids.end(), id) == ids.end())
            ids.push_back(id);
    }

    std::vector<LaneGroup> effectiveLaneGroups(
        const IntersectionInput& inp, const std::vector<ConnectivityCurve>& centerlines) const {
        std::vector<LaneGroup> groups = inp.lane_groups;
        std::map<LaneGroupId, size_t> index;
        for (size_t i = 0; i < groups.size(); ++i)
            index[groups[i].id] = i;

        auto ensure_group = [&](const LaneGroupId& id, GroupRole role) -> LaneGroup& {
            auto it = index.find(id);
            if (it != index.end())
                return groups[it->second];
            LaneGroup group;
            group.id = id;
            group.role = role;
            groups.push_back(group);
            index[id] = groups.size() - 1;
            return groups.back();
        };

        auto add_lane_to_group = [&](const LaneGroupId& id, GroupRole role, const LaneId& lane_id) {
            if (id.empty() || lane_id.empty())
                return;
            LaneGroup& group = ensure_group(id, role);
            addUniqueLaneId(group.lanes, lane_id);
            const Lane* lane = inp.findLane(lane_id);
            if (!lane)
                return;
            addUniqueLaneEdgeId(group.boundaries, lane->left_edge_id);
            addUniqueLaneEdgeId(group.boundaries, lane->right_edge_id);
        };

        for (const auto& conn : inp.connectivities) {
            add_lane_to_group(conn.enterGroupId, GroupRole::Entry, conn.entry_lane_id);
            add_lane_to_group(conn.exitGroupId, GroupRole::Exit, conn.exit_lane_id);
        }

        for (const auto& cc : centerlines) {
            const Connectivity* conn = nullptr;
            for (const auto& c : inp.connectivities) {
                if (c.id == cc.id) {
                    conn = &c;
                    break;
                }
            }
            if (conn) {
                add_lane_to_group(conn->enterGroupId, GroupRole::Entry, cc.entry_lane_id);
                add_lane_to_group(conn->exitGroupId, GroupRole::Exit, cc.exit_lane_id);
                continue;
            }
            const Lane* entry_lane = inp.findLane(cc.entry_lane_id);
            const Lane* exit_lane = inp.findLane(cc.exit_lane_id);
            if (entry_lane)
                add_lane_to_group(entry_lane->groupId, GroupRole::Entry, cc.entry_lane_id);
            if (exit_lane)
                add_lane_to_group(exit_lane->groupId, GroupRole::Exit, cc.exit_lane_id);
        }

        return groups;
    }

    // 提取组内参与端侧边界的车道边线：优先组边界，再补充组内车道左右边线。
    std::vector<LaneEdgeId> laneGroupEdgeIds(
        const IntersectionInput& inp, const LaneGroup& group) const {
        std::vector<LaneEdgeId> ids;
        for (const auto& edge_id : group.boundaries)
            addUniqueLaneEdgeId(ids, edge_id);
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane)
                continue;
            addUniqueLaneEdgeId(ids, lane->left_edge_id);
            addUniqueLaneEdgeId(ids, lane->right_edge_id);
        }
        return ids;
    }

    // 将组端侧形点按容差去重后加入候选集，保留排序所需的 line/lane order 信息。
    void addGroupCutPoint(std::vector<GroupCutPoint>& pts, GroupCutPoint p, double tol = 0.03) const {
        for (const auto& existing : pts) {
            if (dist(existing.pt, p.pt) <= tol)
                return;
        }
        pts.push_back(p);
    }

    // 计算组端侧折线的平移方向：汇总边线、中心线和生成线簇端侧切向，
    // 失败时用组端点到中心兜底。返回方向已区分进入组回缩到路口内、退出组外扩到路口外。
    Vec2d laneGroupShiftDir(
        const IntersectionInput& inp,
        const LaneGroup& group,
        const Vec2d& center,
        const std::vector<ConnectivityCurve>* centerlines = nullptr) const {
        Vec2d dir(0, 0);
        int count = 0;
        if (centerlines) {
            for (const auto& cc : *centerlines) {
                if (!cc.curve)
                    continue;
                if (group.role == GroupRole::Entry) {
                    const Lane* lane = inp.findLane(cc.entry_lane_id);
                    if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, cc.entry_lane_id)))
                        continue;
                    Vec2d t = cc.curve->startTan();
                    if (t.norm() > 1e-8) {
                        dir += t.normalized();
                        ++count;
                    }
                } else {
                    const Lane* lane = inp.findLane(cc.exit_lane_id);
                    if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, cc.exit_lane_id)))
                        continue;
                    Vec2d t = cc.curve->endTan();
                    if (t.norm() > 1e-8) {
                        dir += t.normalized();
                        ++count;
                    }
                }
            }
            if (count > 0 && dir.norm() > 1e-8) {
                dir.normalize();
                Vec2d mid(0, 0);
                int mid_count = 0;
                for (const auto& cc : *centerlines) {
                    if (!cc.curve)
                        continue;
                    if (group.role == GroupRole::Entry) {
                        const Lane* lane = inp.findLane(cc.entry_lane_id);
                        if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, cc.entry_lane_id)))
                            continue;
                        mid += cc.curve->startPt();
                    } else {
                        const Lane* lane = inp.findLane(cc.exit_lane_id);
                        if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, cc.exit_lane_id)))
                            continue;
                        mid += cc.curve->endPt();
                    }
                    ++mid_count;
                }
                if (mid_count > 0) {
                    mid /= (double)mid_count;
                    Vec2d wanted = (group.role == GroupRole::Entry) ? (center - mid) : (mid - center);
                    if (wanted.norm() > 1e-8 && dir.dot(wanted) < 0.0)
                        dir = -dir;
                }
                return dir.normalized();
            }
        }

        for (const auto& edge_id : laneGroupEdgeIds(inp, group)) {
            const LaneEdge* edge = inp.findEdge(edge_id);
            if (!edge || edge->geometry.points.size() < 2)
                continue;
            dir += directedConnDir(edge->geometry.points, group.role, center);
            ++count;
        }
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane || lane->geometry.points.size() < 2)
                continue;
            dir += directedConnDir(lane->geometry.points, group.role, center);
            ++count;
        }
        if (count > 0 && dir.norm() > 1e-8)
            return dir.normalized();

        Vec2d mid(0, 0);
        count = 0;
        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane || lane->geometry.points.empty())
                continue;
            mid += getConnPoint(lane->geometry.points, group.role == GroupRole::Entry);
            ++count;
        }
        if (count > 0) {
            mid /= (double)count;
            dir = (group.role == GroupRole::Entry) ? (center - mid) : (mid - center);
            if (dir.norm() > 1e-8)
                return dir.normalized();
        }
        return Vec2d(1, 0);
    }

    // 由已生成的连通曲线端点构建端侧边界折线：
    // 进入组取生成线簇首端向路口内回缩，退出组取生成线簇尾端向路口外扩。
    std::vector<Vec2d> buildLaneGroupCurveCutLine(
        const IntersectionInput& inp,
        const LaneGroup& group,
        const std::vector<ConnectivityCurve>& centerlines,
        const Vec2d& center) const {
        std::vector<GroupCutPoint> cut_pts;
        Vec2d shift_dir = laneGroupShiftDir(inp, group, center, &centerlines);

        for (const auto& cc : centerlines) {
            if (!cc.curve)
                continue;
            const bool is_entry_group = group.role == GroupRole::Entry;
            const Lane* lane = inp.findLane(is_entry_group ? cc.entry_lane_id : cc.exit_lane_id);
            const LaneId& lane_id = is_entry_group ? cc.entry_lane_id : cc.exit_lane_id;
            if (!lane || (lane->groupId != group.id && !groupContainsLaneId(group, lane_id)))
                continue;

            GroupCutPoint cp;
            cp.pt = is_entry_group ? cc.curve->startPt() : cc.curve->endPt();
            cp.order = lane->laneOrder * 2 + 1;
            cp.has_order = true;
            addGroupCutPoint(cut_pts, cp);
        }

        if (cut_pts.size() < 2)
            return {};

        Vec2d lateral_dir = rotLeft(shift_dir);
        for (auto& p : cut_pts)
            p.lateral = p.pt.dot(lateral_dir);

        std::sort(cut_pts.begin(), cut_pts.end(), [](const GroupCutPoint& a, const GroupCutPoint& b) {
            if (a.has_order != b.has_order)
                return a.has_order > b.has_order;
            if (a.has_order && a.order != b.order)
                return a.order < b.order;
            return a.lateral < b.lateral;
        });

        std::vector<Vec2d> out;
        for (const auto& p : cut_pts)
            pushUnique(out, p.pt + groupCutOffset * shift_dir, 0.03);
        return out;
    }

    // 构建进入/退出组的端侧边界折线：收集车道边线端点和中心线端点，按线序/车道序连接。
    // 最终整条折线按接口传入的 groupCutOffset 做进入组回缩或退出组外扩。
    std::vector<Vec2d> buildLaneGroupCutLine(
        const IntersectionInput& inp, const LaneGroup& group, const Vec2d& center) const {
        std::vector<GroupCutPoint> cut_pts;
        Vec2d shift_dir = laneGroupShiftDir(inp, group, center);

        for (const auto& edge_id : laneGroupEdgeIds(inp, group)) {
            const LaneEdge* edge = inp.findEdge(edge_id);
            if (!edge || edge->geometry.points.size() < 2)
                continue;
            GroupCutPoint cp;
            cp.pt = getConnPoint(edge->geometry.points, group.role == GroupRole::Entry);
            cp.order = edge->lineOrder;
            cp.has_order = true;
            addGroupCutPoint(cut_pts, cp);
        }

        for (const auto& lane_id : uniqueGroupLaneIds(group)) {
            const Lane* lane = inp.findLane(lane_id);
            if (!lane || lane->geometry.points.size() < 2)
                continue;
            GroupCutPoint cp;
            cp.pt = getConnPoint(lane->geometry.points, group.role == GroupRole::Entry);
            cp.order = lane->laneOrder * 2 + 1;
            cp.has_order = true;
            addGroupCutPoint(cut_pts, cp);
        }

        if (cut_pts.size() < 2)
            return {};
        Vec2d lateral_dir = rotLeft(shift_dir);
        for (auto& p : cut_pts)
            p.lateral = p.pt.dot(lateral_dir);

        std::sort(cut_pts.begin(), cut_pts.end(), [](const GroupCutPoint& a, const GroupCutPoint& b) {
            if (a.has_order != b.has_order)
                return a.has_order > b.has_order;
            if (a.has_order && a.order != b.order)
                return a.order < b.order;
            return a.lateral < b.lateral;
        });

        std::vector<Vec2d> out;
        for (const auto& p : cut_pts)
            pushUnique(out, p.pt + groupCutOffset * shift_dir, 0.03);
        return out;
    }

    // 将端侧边界折线两端按接口参数外延，用外延段与道路边缘求交来消除未对齐尖角。
    // 仅延长首尾端点，中间形点保持原折线形态。
    std::vector<Vec2d> extendCutLineEnds(const std::vector<Vec2d>& line) const {
        std::vector<Vec2d> out = line;
        if (out.size() < 2)
            return out;
        double len = cutLineExtension;
        Vec2d first_dir = out.front() - out[1];
        if (first_dir.norm() > 1e-8)
            out.front() += len * first_dir.normalized();
        Vec2d last_dir = out.back() - out[out.size() - 2];
        if (last_dir.norm() > 1e-8)
            out.back() += len * last_dir.normalized();
        return out;
    }

    struct BoundaryHit {
        bool found = false;
        int boundary_index = -1;
        Vec2d pt{0, 0};
        double station = 0.0;
        double distance = 1e18;
    };

    // 查找外延探测线与所有道路边缘线的最近交点，优先选择离原端点最近的交点。
    BoundaryHit nearestBoundaryHit(
        const std::vector<Vec2d>& probe,
        const std::vector<std::vector<Vec2d>>& boundaries,
        const Vec2d& prefer) const {
        BoundaryHit best;
        for (int i = 0; i < (int)boundaries.size(); ++i) {
            for (const auto& hit : polylineIntersections(probe, boundaries[i])) {
                double d = dist(hit.pt, prefer);
                if (d < best.distance) {
                    best.found = true;
                    best.boundary_index = i;
                    best.pt = hit.pt;
                    best.station = hit.station_b;
                    best.distance = d;
                }
            }
        }
        return best;
    }

    // 记录道路边缘被端侧边界折线命中的里程位置，后续据此截取有效道路边缘段。
    void addBoundaryHitStation(
        std::vector<std::vector<double>>& hit_stations,
        const BoundaryHit& hit) const {
        if (!hit.found || hit.boundary_index < 0 ||
            hit.boundary_index >= (int)hit_stations.size())
            return;
        for (double s : hit_stations[hit.boundary_index])
            if (std::abs(s - hit.station) < 0.03)
                return;
        hit_stations[hit.boundary_index].push_back(hit.station);
    }

    // 对齐端侧边界折线与道路边缘：端侧折线两端外延求交，并同步截取对应道路边缘。
    // 若某端没有交点，则保留原端点和原道路边缘，避免过度裁剪。
    std::vector<std::vector<Vec2d>> alignCutLinesAndBoundaries(
        std::vector<std::vector<Vec2d>>& boundary_lines,
        std::vector<std::vector<Vec2d>> cut_lines,
        const Vec2d& center) const {
        std::vector<std::vector<double>> hit_stations(boundary_lines.size());

        std::vector<std::vector<Vec2d>> adjusted_cuts;
        adjusted_cuts.reserve(cut_lines.size() * 2);
        for (auto& cut : cut_lines) {
            if (cut.size() < 2)
                continue;
            std::vector<Vec2d> extended = extendCutLineEnds(cut);
            std::vector<Vec2d> first_probe = {extended.front(), cut.front()};
            std::vector<Vec2d> last_probe = {cut.back(), extended.back()};
            BoundaryHit first_hit = nearestBoundaryHit(first_probe, boundary_lines, cut.front());
            BoundaryHit last_hit = nearestBoundaryHit(last_probe, boundary_lines, cut.back());
            if (first_hit.found)
                cut.front() = first_hit.pt;
            if (last_hit.found)
                cut.back() = last_hit.pt;
            addBoundaryHitStation(hit_stations, first_hit);
            addBoundaryHitStation(hit_stations, last_hit);
            adjusted_cuts.push_back(cut);
        }

        for (int i = 0; i < (int)boundary_lines.size(); ++i) {
            if (boundary_lines[i].size() < 2)
                continue;
            if (hit_stations[i].empty())
                continue;
            double total = polylineLength(boundary_lines[i]);
            std::sort(hit_stations[i].begin(), hit_stations[i].end());
            double s0 = 0.0;
            double s1 = total;
            if (hit_stations[i].size() >= 2) {
                s0 = hit_stations[i].front();
                s1 = hit_stations[i].back();
            } else {
                double s = hit_stations[i].front();
                Vec2d mid0 = pointAtStation(boundary_lines[i], 0.5 * s);
                Vec2d mid1 = pointAtStation(boundary_lines[i], 0.5 * (s + total));
                if (dist(mid0, center) < dist(mid1, center)) {
                    s0 = 0.0;
                    s1 = s;
                } else {
                    s0 = s;
                    s1 = total;
                }
            }
            s0 = std::max(0.0, std::min(total, s0));
            s1 = std::max(0.0, std::min(total, s1));
            if (s1 - s0 > 0.05)
                boundary_lines[i] = subPolylineByStation(boundary_lines[i], s0, s1);
        }
        return adjusted_cuts;
    }

    // 组装精细路口面：道路边缘线与所有进入/退出组端侧边界折线共同排序成可凹多边形。
    // 端侧边界先按组回缩/外扩，再与道路边缘对齐，最后统一修复多边形绕向与闭合。
    std::vector<Vec2d> buildFromLaneGroupCuts(
        const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines,
        const Vec2d& center) const {
        std::vector<AreaChain> chains;
        std::vector<std::vector<Vec2d>> boundary_lines;
        std::vector<std::vector<Vec2d>> cut_lines;

        for (const auto& bnd : inp.boundaries) {
            if (bnd.geometry.points.size() < 2)
                continue;
            boundary_lines.push_back(bnd.geometry.points);
        }

        std::vector<LaneGroup> groups = effectiveLaneGroups(inp, centerlines);
        for (const auto& group : groups) {
            std::vector<Vec2d> cut = buildLaneGroupCurveCutLine(inp, group, centerlines, center);
            if (cut.size() < 2)
                cut = buildLaneGroupCutLine(inp, group, center);
            if (cut.size() >= 2)
                cut_lines.push_back(cut);
        }

        std::vector<std::vector<Vec2d>> aligned_cuts =
            alignCutLinesAndBoundaries(boundary_lines, cut_lines, center);

        for (const auto& line : boundary_lines) {
            if (line.size() < 2)
                continue;
            AreaChain chain;
            chain.pts = line;
            chains.push_back(orientChain(chain, center));
        }
        for (const auto& line : aligned_cuts) {
            if (line.size() < 2)
                continue;
            AreaChain chain;
            chain.pts = line;
            chains.push_back(orientChain(chain, center));
        }

        if (chains.size() < 3)
            return {};

        std::sort(chains.begin(), chains.end(), [&](const AreaChain& a, const AreaChain& b) {
            double angA = angleOf(chainMidpoint(a), center);
            double angB = angleOf(chainMidpoint(b), center);
            return winding == "clockwise" ? angA > angB : angA < angB;
        });

        std::vector<Vec2d> raw;
        for (const auto& chain : chains) {
            for (const auto& pt : chain.pts)
                pushUnique(raw, pt);
        }
        if (raw.size() < 3)
            return {};
        return repairPolygon(raw, center);
    }

    // 旧方案：无道路边缘时在 rough_area 上插入停止线 break-line。
    // 当前精细面生成暂不使用该路径，入口调用已注释保留，函数体留作后续回退参考。
    std::vector<Vec2d> buildFromRoughBase(
        const IntersectionInput& inp, const Vec2d& center) const {
        std::vector<Vec2d> skeleton = openRing(inp.area.geometry.outer);
        if (skeleton.size() < 3)
            return {};

        struct StopLineInsert {
            int segment = 0;
            double t = 0.0;
            std::vector<Vec2d> pts;
        };

        std::vector<std::vector<StopLineInsert>> inserts(skeleton.size());
        for (const auto& sl : inp.stop_lines) {
            AreaChain chain;
            chain.pts = buildStopLineBreakLine(inp, {}, sl, center).shifted_pts;
            if (chain.pts.size() < 2)
                continue;
            chain = orientChain(chain, center);
            SegmentProjection pr = closestRingSegment(skeleton, chainMidpoint(chain));
            StopLineInsert ins;
            ins.segment = pr.index;
            ins.t = pr.t;
            ins.pts = chain.pts;
            inserts[pr.index].push_back(ins);
        }

        std::vector<Vec2d> raw;
        for (int i = 0; i < (int)skeleton.size(); ++i) {
            pushUnique(raw, skeleton[i]);
            auto& seg_inserts = inserts[i];
            std::sort(seg_inserts.begin(), seg_inserts.end(), [](const StopLineInsert& a, const StopLineInsert& b) {
                return a.t < b.t;
            });
            for (const auto& ins : seg_inserts)
                for (const auto& pt : ins.pts)
                    pushUnique(raw, pt);
        }
        if (raw.size() < 3)
            return {};
        return repairPolygon(raw, center);
    }

    // 精细路口面入口：优先使用进入/退出组端侧边界折线围成可凹多边形。
    // 道路边缘存在时参与端侧线对齐和边缘截断；不存在时仍保留组端侧线的回缩/外扩结果。
    // 暂未启用的停止线/rough_area 旧路径只保留注释调用，由 build() 的传统端点法兜底。
    std::vector<Vec2d> buildFinePolygon(
        const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& centerlines,
        const Vec2d& center,
        double radius) const {
        (void)radius;
        std::vector<Vec2d> by_group_cuts = buildFromLaneGroupCuts(inp, centerlines, center);
        if (by_group_cuts.size() >= 4)
            return by_group_cuts;

        if (!inp.boundaries.empty()) {
            // 当前精细路口面暂不使用停止线 break-line 旧逻辑，保留代码便于后续对比/回退。
            // return buildFromBoundaryElements(inp, centerlines, center);
        }

        // 当前精细路口面暂不使用 rough_area + 停止线插入旧逻辑，保留代码便于后续对比/回退。
        // return buildFromRoughBase(inp, center);
        return {};
    }

    // 计算路口中心点（连接点的质心）
    Vec2d computeCenter(const IntersectionInput& inp,
                        const std::vector<ConnectivityCurve>& cls) const {
        Vec2d c{0, 0};
        int cnt = 0;
        for (auto& gcl : cls) {
            if(gcl.curve) {
                c += gcl.curve->startPt(); ++cnt;
                c += gcl.curve->endPt();  ++cnt;
            }
        }
        if (cnt > 0) return c * (1.0 / cnt);

        // 回退：用中心线连接点
        for (auto& cl : inp.lanes) {
            c += getConnPoint(cl.geometry.points, inp.IsEntryLane(cl.id)); //cl.connectionPt;
            ++cnt;
        }
        return cnt > 0 ? c * (1.0 / cnt) : Vec2d{0, 0};
    }

    // 估算路口半径（连接点到中心的最大距离）
    double computeRadius(const IntersectionInput& inp,
                         const std::vector<ConnectivityCurve>& cls,
                         const Vec2d& center) const {
        double maxD = 5.0;
        for (auto& gcl : cls) {
            if (gcl.curve) {
                maxD = std::max(maxD, dist(gcl.curve->startPt(), center));
                maxD = std::max(maxD, dist(gcl.curve->endPt(), center));
            }
        }
        for (auto& cl : inp.lanes) {
            auto connPt = getConnPoint(cl.geometry.points, inp.IsEntryLane(cl.id));
            maxD = std::max(maxD, dist(connPt, center));
        }
        return maxD;
    }

    // 将端点吸附到最近的车道边线端点（若在容差范围内）
    Vec2d snapToEdgePt(
        const Vec2d& pt, const std::vector<ConnectivityLaneEdge>& edgelines,
        const IntersectionInput& inp, double tol) const {
        double minD = tol;
        Vec2d best = pt;
        // 检查生成的路口内边线端点
        for (auto& el : edgelines) {
            if (el.geometry.points.empty()) continue;
            double d0 = dist(pt, el.geometry.points.front());
            double d1 = dist(pt, el.geometry.points.back());
            if (d0 < minD) {
                minD = d0;
                best = el.geometry.points.front();
            }
            if (d1 < minD) {
                minD = d1;
                best = el.geometry.points.back();
            }
        }
        // 检查路口外边线连接点
        for (auto& el : inp.lane_edges) {
            if (el.geometry.points.empty()) continue;
            auto connPt = getConnPoint(el.geometry.points, true); //el.connectionPt//TODO:entry|exit
            double d = dist(pt, connPt);
            if (d < minD) {
                minD = d;
                best = connPt;
            }
        }
        return best;
    }

    // 从车道边线端点补充路口面顶点
    void supplementEndpointsFromLaneEdges(
        std::vector<struct EdgeEndpoint>& endpoints,
        const IntersectionInput& inp,
        const std::vector<ConnectivityLaneEdge>& edgelines,
        const Vec2d& center, double radius) const {
        // 收集路口外边线的连接点（这些是路口面的边界顶点）
        for (auto& el : inp.lane_edges) {
            if (el.geometry.points.empty()) continue;
            Vec2d cp = getConnPoint(el.geometry.points, inp.IsEntryLaneEdge(el.id)); //el.connectionPt;
            if (dist(cp, center) > radius * 2.0) continue;

            // 检查是否已有相近的端点
            bool dup = false;
            for (auto& ep : endpoints) {
                if (dist(ep.pt, cp) < 0.5) {
                    dup = true;
                    break;
                }
            }
            if (!dup) endpoints.emplace_back(EdgeEndpoint{cp, -1, true});
        }

        // 补充中心线连接点的法线方向偏移点（近似边界点）
        for (auto& cl : inp.lanes) {
            if (cl.geometry.points.size() < 2) continue;
            Vec2d cp = getConnPoint(cl.geometry.points, inp.IsEntryLane(cl.id)); //cl.connectionPt;
            Vec2d ct = getConnTangent(cl.geometry.points, inp.IsEntryLane(cl.id)); //cl.tangentDir
            Vec2d n = rotLeft(ct);
            // 左侧点
            Vec2d lp = cp + n * snapTolerance * 2;
            Vec2d rp = cp - n * snapTolerance * 2;
            bool dupL = false, dupR = false;
            for (auto& ep : endpoints) {
                if (dist(ep.pt, lp) < 0.5) dupL = true;
                if (dist(ep.pt, rp) < 0.5) dupR = true;
            }
            if (!dupL) endpoints.emplace_back(EdgeEndpoint{lp, -1, true});
            if (!dupR) endpoints.emplace_back(EdgeEndpoint{rp, -1, true});
        }
    }

    // 从道路边缘线提取两点间的几何段
    std::vector<Vec2d> extractEdgeSegment(
        const Boundary& re, const Vec2d& a, const Vec2d& b, double tol) const {
        if (re.geometry.points.size() < 2) return {};

        // 找到 a 和 b 在 re.geometry.points 上最近的索引
        int idxA = -1, idxB = -1;
        double minDA = 1e18, minDB = 1e18;
        for (int i = 0; i < (int)re.geometry.points.size(); ++i) {
            double da = dist(re.geometry.points[i], a);
            double db = dist(re.geometry.points[i], b);
            if (da < minDA) {
                minDA = da;
                idxA = i;
            }
            if (db < minDB) {
                minDB = db;
                idxB = i;
            }
        }

        if (idxA < 0 || idxB < 0 || idxA == idxB) return {};

        std::vector<Vec2d> seg;
        int step = (idxA < idxB) ? 1 : -1;
        for (int i = idxA; i != idxB; i += step) {
            seg.push_back(re.geometry.points[i]);
        }
        // 不包含 b 本身（由下一个端点自己加）
        return seg;
    }

    // 去重端点（合并距离过近的点）
    void deduplicateEndpoints(std::vector<struct EdgeEndpoint>& eps, double tol) const {
        std::vector<bool> remove(eps.size(), false);
        for (size_t i = 0; i < eps.size(); ++i) {
            if (remove[i]) continue;
            for (size_t j = i + 1; j < eps.size(); ++j) {
                if (!remove[j] && dist(eps[i].pt, eps[j].pt) < tol) {
                    remove[j] = true;
                }
            }
        }
        std::vector<struct EdgeEndpoint> out;
        for (size_t i = 0; i < eps.size(); ++i)
            if (!remove[i]) out.push_back(eps[i]);
        eps = out;
    }

    // 多边形修复：移除共线点，确保方向正确
    std::vector<Vec2d> repairPolygon(const std::vector<Vec2d>& raw, const Vec2d& center) const {
        if (raw.size() < 3) return raw;

        // 移除重复点
        std::vector<Vec2d> pts;
        for (auto& p : raw) {
            if (pts.empty() || dist(pts.back(), p) > EPS) pts.push_back(p);
        }
        // 移除首尾重复
        while (pts.size() > 1 && dist(pts.front(), pts.back()) < EPS)
            pts.pop_back();

        if (pts.size() < 3) return raw;

        // 检查方向（有符号面积）
        double area = polygonSignedArea(pts);
        bool isCCW = area > 0;

        if (winding == "clockwise" && isCCW) {
            std::reverse(pts.begin(), pts.end());
        } else if (winding == "counter_clockwise" && !isCCW) {
            std::reverse(pts.begin(), pts.end());
        }

        // 闭合
        if (dist(pts.front(), pts.back()) > EPS)
            pts.push_back(pts.front());

        return pts;
    }

    // 扩展多边形以包含所有中心线点（简单的凸包扩展）
    void expandToContainCenterlines(
        std::vector<Vec2d>& poly, const std::vector<ConnectivityCurve>& cls, const Vec2d& center) const {
        if (poly.size() < 3) return;
        bool expanded = false;
        for (auto& gcl : cls) {
            auto curve = gcl.curve;
            if (!curve || curve->segs.empty()) {
                continue;
            }
            for (auto& p : gcl.curve->sample(30)) {
                if (!pointInPolygon(p, poly)) {
                    // 将该点加入多边形（扩展最近边）
                    expandByPoint(poly, p, center);
                    expanded = true;
                }
            }
        }
        if (expanded) {
            // 重新排序（保证凸包方向）
            poly = repairPolygon(poly, center);
        }
    }

    void expandByPoint(std::vector<Vec2d>& poly, const Vec2d& p, const Vec2d& center) const {
        if (poly.size() < 3) return;

        // 找到最近边，在该边插入点
        double minD = std::numeric_limits<double>::max();
        int insertAfter = 0;
        size_t n = poly.size();
        // 跳过最后一个（它是闭合点，等于第一个）
        size_t end = (dist(poly.front(), poly.back()) < EPS) ? n - 1 : n;
        for (size_t i = 0; i < end; ++i) {
            size_t j = (i + 1) % end;
            std::pair<double, double> _seg = pointToSegment(p, poly[i], poly[j]);
            if (_seg.first < minD) {
                minD = _seg.first;
                insertAfter = (int)i;
            }
        }
        poly.insert(poly.begin() + insertAfter + 1, p);
    }

    // 最终回退：所有连接点的凸包
    std::vector<Vec2d> buildConvexHullFallback(
        const IntersectionInput& inp,
        const std::vector<ConnectivityCurve>& cls) const
    {
        std::vector<Vec2d> pts;
        for (auto& gcl : cls) {
            if (gcl.curve) {
                pts.push_back(gcl.curve->startPt());
                pts.push_back(gcl.curve->endPt());
            }
        }
        for (auto& cl : inp.lanes) {
            Vec2d cp = getConnPoint(cl.geometry.points, inp.IsEntryLane(cl.id)); //cl.connectionPt;
            pts.push_back(cp);
        }
        if (pts.empty()) return {};
        return convexHull(pts);
    }

    // Andrew's monotone chain 凸包算法
    std::vector<Vec2d> convexHull(std::vector<Vec2d> pts) const {
        int n = pts.size();
        if (n < 3) return pts;
        std::sort(pts.begin(), pts.end(), [](const Vec2d& a, const Vec2d& b) {
            return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
        });

        std::vector<Vec2d> hull;
        // 下凸包
        for (int i = 0; i < n; ++i) {
            while (hull.size() >= 2) {
                Vec2d a = hull[hull.size() - 2], b = hull.back();
                if (cross2d((b - a), pts[i] - a) <= 0) hull.pop_back();
                else break;
            }
            hull.push_back(pts[i]);
        }
        // 上凸包
        int lower = hull.size();
        for (int i = n - 2; i >= 0; --i) {
            while ((int)hull.size() > lower) {
                Vec2d a = hull[hull.size() - 2], b = hull.back();
                if (cross2d((b - a), pts[i] - a) <= 0) hull.pop_back();
                else break;
            }
            hull.push_back(pts[i]);
        }
        hull.pop_back();
        if (!hull.empty()) hull.push_back(hull.front()); // 闭合
        return hull;
    }
};

}
