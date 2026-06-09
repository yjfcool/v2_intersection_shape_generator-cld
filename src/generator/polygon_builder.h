#pragma once
#include "types.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

/**
 * 精细路口面构建
 * 输入：道路边缘线（独立线段）+ 生成的车道边线端点
 * 输出：封闭多边形（路口面）
 *
 * 流程：
 *  1. 确定路口范围（rough_area 或 由连接点包围盒扩展）
 *  2. 裁剪/截取每条道路边缘线在路口侧的端点
 *  3. 按极角排序端点，连接成封闭多边形
 *  4. 相邻端点间若属同一边缘线则沿线连接，否则直连
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

public:
    explicit IntersectionAreaBuilder(double _snapTolerance = 0.5, std::string _winding = "clockwise")
        : snapTolerance(_snapTolerance), winding(_winding) {}

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

        // std::cout<<("Polygon: center=(" + std::to_string(center[0]) + "," +
        //              std::to_string(center[1]) + ") radius=" + std::to_string(radius))<<std::endl;

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
                // std::cout<<("Edge endpoint[" + std::to_string(i) + "]: ("
                //     + std::to_string(snapped[0]) + "," + std::to_string(snapped[1]) + ")")<<std::endl;
            }
        }

        // 若道路边缘线端点不足，从连通关系连接点的外侧边线端点补充
        if (endpoints.size() < 3) {
            // std::cout<<("Road edge endpoints insufficient (" +
            //              std::to_string(endpoints.size()) + "), supplementing from lane edges.")<<std::endl;
            supplementEndpointsFromLaneEdges(endpoints, inp, edgelines, center, radius);
        }

        if (endpoints.size() < 3) {
            // 终极回退：用所有连接点的凸包
            // std::cout<<("Still insufficient endpoints, using convex hull of connection points.")<<std::endl;
            junctArea.geometry.outer = buildConvexHullFallback(inp, centerlines);
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

        // 7. 确保多边形包含所有中心线采样点（扩展检查）
        expandToContainCenterlines(junctArea.geometry.outer, centerlines, center);

        double area = polygonArea(junctArea.geometry.outer);
        // std::cout<<("Polygon built: vertices=" + std::to_string(junctArea.geometry.outer.size())
        //              + " area=" + std::to_string(area))<<std::endl;
        //
        // if(area < 1.0){
        //     std::cout<<("Polygon area too small: " + std::to_string(area));
        // }

        return junctArea;
    }

private:
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

    // Andrew's monotone chain 凸包
    std::vector<Vec2d> convexHull(std::vector<Vec2d> pts) const {
        int n = pts.size();
        if (n < 3) return pts;
        std::sort(pts.begin(), pts.end(), [](const Vec2d& a, const Vec2d& b) {
            return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
        });

        std::vector<Vec2d> hull;
        // Lower hull
        for (int i = 0; i < n; ++i) {
            while (hull.size() >= 2) {
                Vec2d a = hull[hull.size() - 2], b = hull.back();
                if (cross2d((b - a), pts[i] - a) <= 0) hull.pop_back();
                else break;
            }
            hull.push_back(pts[i]);
        }
        // Upper hull
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
