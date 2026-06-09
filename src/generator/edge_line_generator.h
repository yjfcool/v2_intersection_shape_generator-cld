#pragma once
#include <iostream>
#include <map>
#include <set>

#include "types.h"
#include "curve/bezier.h"

/**
 * Edge line generator (optimized shared-edge version)
 *
 * Strategy:
 *  1. Group generated centerlines by (enterGroupId, exitGroupId) pairs
 *  2. Sort within each group by enter lane order (inner to outer)
 *  3. Detect "special" connections that cross other lanes in the same group
 *  4. For adjacent non-special lanes that share edges on BOTH enter and exit sides,
 *     generate a single shared midline (cubic Bezier) instead of two independent edges
 *  5. Outermost right, innermost left, and special lane edges use independent generation
 *
 * Shared midline properties:
 *  - G1 smooth at both endpoints
 *  - Endpoints exactly coincide with shared edge connection points
 *  - Uses CubicBezier::fromAlphaBeta with alpha=0.38, beta=0.38
 */
class EdgeLineGenerator {
    double bezier_maxCurvature = 1.0;
    double edgeLine_defaultLaneWidth = 1.0;

public:
    explicit EdgeLineGenerator() {}

    std::vector<ConnectivityLaneEdge> generate(
        const IntersectionInput& inp, std::vector<ConnectivityCurve>& centerlines) {
        std::vector<ConnectivityLaneEdge> result;
        std::map<std::string, std::pair<std::string, std::string>> edgeIdMap;

        // Build connectionId -> ConnectivityCurve mapping
        std::map<std::string, const ConnectivityCurve*> clMap;
        for (auto& gcl : centerlines) clMap[gcl.id] = &gcl;

        // Build connection id -> Connection mapping
        std::map<std::string, const Connectivity*> connMap;
        for (auto& conn : inp.connectivities) connMap[conn.id] = &conn;

        // Group generated centerlines by (enterGroupId, exitGroupId)
        using GroupKey = std::pair<std::string, std::string>;
        std::map<GroupKey, std::vector<const ConnectivityCurve*>> pairGroups;

        for (auto& gcl : centerlines) {
            auto cIt = connMap.find(gcl.id);
            if (cIt == connMap.end()) continue;
            const Connectivity* conn = cIt->second;
            GroupKey key = {conn->enterGroupId, conn->exitGroupId};
            pairGroups[key].push_back(&gcl);
        }

        // Process each (enterGroupId, exitGroupId) group
        for (auto& kv : pairGroups) {
            auto& key = kv.first;
            (void)key;
            auto& groupCls = kv.second;
            // Sort by enter lane order (inner to outer)
            std::vector<const ConnectivityCurve*> sorted = groupCls;
            std::sort(sorted.begin(), sorted.end(),
                      [&](const ConnectivityCurve* a, const ConnectivityCurve* b) {
                          auto aIt = inp.findLane(a->entry_lane_id);
                          auto bIt = inp.findLane(b->entry_lane_id);
                          int aOrder = (aIt) ? aIt->laneOrder : 0;
                          int bOrder = (bIt) ? bIt->laneOrder : 0;
                          return aOrder < bOrder;
                      });

            // Detect special connections (those that cross other lanes in the group)
            std::set<std::string> specialIds;
            if (sorted.size() > 1) {
                for (size_t i = 0; i < sorted.size(); ++i) {
                    for (size_t j = i + 1; j < sorted.size(); ++j) {
                        auto icurve = sorted[i]->curve;
                        auto jcurve = sorted[j]->curve;
                        if (!icurve || icurve->empty()) continue;
                        if (!jcurve || jcurve->empty()) continue;
                        auto ismps = icurve->sample(50);
                        auto jsmps = jcurve->sample(50);
                        if (ismps.size() < 2 || jsmps.size() < 2) continue;
                        if (polylinesIntersectExcludeEndpoints(ismps, jsmps)) {
                            specialIds.insert(sorted[i]->id);
                            specialIds.insert(sorted[j]->id);
                        }
                    }
                }
            }

            // Determine which adjacent pairs share edges on both sides
            // sharedBetween[i] = true means lanes i and i+1 share an edge
            std::vector<bool> sharedBetween(sorted.size() > 0 ? sorted.size() - 1 : 0, false);
            // Store shared edge IDs for each pair
            struct SharedEdgeInfo {
                std::string enterEdgeId;
                std::string exitEdgeId;
            };
            std::vector<SharedEdgeInfo> sharedEdgeInfos(sharedBetween.size());

            for (size_t i = 0; i + 1 < sorted.size(); ++i) {
                if (specialIds.count(sorted[i]->id) || specialIds.count(sorted[i + 1]->id))
                    continue;

                const std::string& enterLineI = sorted[i]->entry_lane_id;
                const std::string& enterLineI1 = sorted[i + 1]->entry_lane_id;
                const std::string& exitLineI = sorted[i]->exit_lane_id;
                const std::string& exitLineI1 = sorted[i + 1]->exit_lane_id;

                // Enter side: lane_i's right edge == lane_{i+1}'s left edge
                // Compare by geometry (connectionPt) since different IDs may represent the same physical edge
                auto enterI = inp.findLane(enterLineI);
                auto enterI1 = inp.findLane(enterLineI1);
                bool enterShared = false;
                std::string sharedEnterEdgeId;
                if (enterI && enterI1) {
                    const std::string& rightOfI = enterI->right_edge_id;
                    const std::string& leftOfI1 = enterI1->left_edge_id;
                    if (!rightOfI.empty() && !leftOfI1.empty()) {
                        // Check if IDs match directly (nCL+1 convention)
                        if (rightOfI == leftOfI1) {
                            enterShared = true;
                            sharedEnterEdgeId = rightOfI;
                        } else {
                            // Check if connection points match (2*nCL convention with identical geometry)
                            auto rightEdgeIt = inp.findEdge(rightOfI);
                            auto leftEdgeIt = inp.findEdge(leftOfI1);
                            if (rightEdgeIt && leftEdgeIt) {
                                auto rEdgePt = getConnPoint(rightEdgeIt->geometry.points, true);
                                auto lEdgePt = getConnPoint(leftEdgeIt->geometry.points, true);
                                double ptDist = dist(rEdgePt, lEdgePt);
                                if (ptDist < 0.01) {
                                    enterShared = true;
                                    sharedEnterEdgeId = rightOfI;
                                }
                            }
                        }
                    }
                }

                // Exit side: lane_i's right edge == lane_{i+1}'s left edge
                // (same pattern as enter side - shared edge is between adjacent lanes)
                // Compare by geometry (connectionPt) since different IDs may represent the same physical edge
                auto exitI = inp.findLane(exitLineI);
                auto exitI1 = inp.findLane(exitLineI1);
                bool exitShared = false;
                std::string sharedExitEdgeId;
                if (exitI && exitI1) {
                    const std::string& rightOfI = exitI->right_edge_id;
                    const std::string& leftOfI1 = exitI1->left_edge_id;
                    if (!rightOfI.empty() && !leftOfI1.empty()) {
                        // Check if IDs match directly (nCL+1 convention)
                        if (rightOfI == leftOfI1) {
                            exitShared = true;
                            sharedExitEdgeId = rightOfI;
                        } else {
                            // Check if connection points match (2*nCL convention with identical geometry)
                            auto rightEdgeIt = inp.findEdge(rightOfI);
                            auto leftEdgeIt = inp.findEdge(leftOfI1);
                            if (rightEdgeIt && leftEdgeIt) {
                                auto rEdgePt = getConnPoint(rightEdgeIt->geometry.points, false);
                                auto lEdgePt = getConnPoint(leftEdgeIt->geometry.points, false);
                                double ptDist = dist(rEdgePt, lEdgePt);
                                if (ptDist < 0.01) {
                                    exitShared = true;
                                    sharedExitEdgeId = rightOfI;
                                }
                            }
                        }
                    }
                }

                if (enterShared && exitShared) {
                    sharedBetween[i] = true;
                    sharedEdgeInfos[i] = {sharedEnterEdgeId, sharedExitEdgeId};
                }
            }

            // Generate edges for each lane in the group
            for (size_t i = 0; i < sorted.size(); ++i) {
                const ConnectivityCurve& gcl = *sorted[i];
                if (!gcl.curve) continue;

                auto connIt = connMap.find(gcl.id);
                if (connIt == connMap.end()) continue;
                const Connectivity* conn = connIt->second;

                bool isSpecial = specialIds.count(gcl.id) > 0;

                // Determine if left edge is shared (from lane i-1's right)
                bool leftShared = (!isSpecial && i > 0 && sharedBetween[i - 1]);
                // Determine if right edge is shared (with lane i+1's left)
                bool rightShared = (!isSpecial && i + 1 < sorted.size() && sharedBetween[i]);

                // Generate left edge
                std::string leftElId;
                if (leftShared) {
                    // Left edge was already generated as the shared midline for pair (i-1, i)
                    leftElId = "gen_el_shared_" + sorted[i - 1]->id + "_" + gcl.id;
                } else {
                    // Generate independent left edge
                    leftElId = "gen_el_left_" + conn->id;
                    std::vector<Vec2d> leftPts = generateIndependentEdge(gcl, conn, true, inp);
                    ConnectivityLaneEdge el;
                    el.id = leftElId;
                    el.geometry.points = leftPts;
                    //el.shared_by->first = gcl.id;
                    result.push_back(el);
                }

                // Generate right edge
                std::string rightElId;
                if (rightShared) {
                    // Generate a shared midline between lane i and lane i+1
                    rightElId = "gen_el_shared_" + gcl.id + "_" + sorted[i + 1]->id;
                    std::vector<Vec2d> sharedPts = generateSharedMidline(
                        gcl, *sorted[i + 1],
                        sharedEdgeInfos[i].enterEdgeId,
                        sharedEdgeInfos[i].exitEdgeId,
                        inp);
                    ConnectivityLaneEdge el;
                    el.id = rightElId;
                    el.geometry.points = sharedPts;
                    // centerlineId records the inner lane (lane I). The outer lane (lane I+1)
                    // references this shared edge via its backfilled left_edge_id string,
                    // so consumers using left_edge_id/right_edge_id will find it correctly.
                    //el.shared_by->second = gcl.id;
                    // el.qualityFlags = checkSharedMidlineCurvature(
                    //     gcl, *sorted[i+1],
                    //     sharedEdgeInfos[i].enterEdgeId,
                    //     sharedEdgeInfos[i].exitEdgeId,
                    //     inp);
                    result.push_back(el);
                } else {
                    // Generate independent right edge
                    rightElId = "gen_el_right_" + conn->id;
                    std::vector<Vec2d> rightPts = generateIndependentEdge(gcl, conn, false, inp);
                    ConnectivityLaneEdge el;
                    el.id = rightElId;
                    el.geometry.points = rightPts;
                    //el.shared_by->second = gcl.id;
                    result.push_back(el);
                }

                edgeIdMap[gcl.id] = {leftElId, rightElId};
            }
        }

        // Backfill ConnectivityCurve left_edge_id / right_edge_id
        for (auto& gcl : centerlines) {
            auto it = edgeIdMap.find(gcl.id);
            if (it != edgeIdMap.end()) {
                gcl.left_edge_id = it->second.first;
                gcl.right_edge_id = it->second.second;
            }
        }

        return result;
    }

private:
    // Generate a shared midline between two adjacent lanes
    std::vector<Vec2d> generateSharedMidline(
        const ConnectivityCurve& gclI,
        const ConnectivityCurve& gclI1,
        const std::string& sharedEnterEdgeId,
        const std::string& sharedExitEdgeId,
        const IntersectionInput& inp) const {
        // Start point = shared enter edge connection point
        Vec2d startPt{0, 0};
        auto enterEdgeIt = inp.findEdge(sharedEnterEdgeId);
        if (enterEdgeIt) {
            startPt = getConnPoint(enterEdgeIt->geometry.points, true); //enterEdgeIt->connectionPt;
        }

        // End point = shared exit edge connection point
        Vec2d endPt{0, 0};
        auto exitEdgeIt = inp.findEdge(sharedExitEdgeId);
        if (exitEdgeIt) {
            endPt = getConnPoint(exitEdgeIt->geometry.points, false); //exitEdgeIt->connectionPt;
        }

        // Start tangent = average of enter tangent directions from both lanes
        auto enterClI = inp.findLane(gclI.entry_lane_id);
        auto enterClI1 = inp.findLane(gclI1.entry_lane_id);
        Vec2d startTangI = (enterClI) ? getConnTangent(enterClI->geometry.points, true) : Vec2d{0, 1};
        Vec2d startTangI1 = (enterClI1) ? getConnTangent(enterClI1->geometry.points, true) : Vec2d{0, 1};
        Vec2d startTangSum = startTangI + startTangI1;
        Vec2d startTang = (startTangSum.norm() > EPS) ? startTangSum.normalized() : startTangI;

        // End tangent = average of exit tangent directions from both lanes (points inward)
        auto exitClI = inp.findLane(gclI.exit_lane_id);
        auto exitClI1 = inp.findLane(gclI1.exit_lane_id);
        Vec2d endTangI = (exitClI) ? getConnTangent(exitClI->geometry.points, false) : Vec2d{0, 1};
        Vec2d endTangI1 = (exitClI1) ? getConnTangent(exitClI1->geometry.points, false) : Vec2d{0, 1};
        Vec2d endTangSum = endTangI + endTangI1;
        Vec2d endTang = (endTangSum.norm() > EPS) ? endTangSum.normalized() : endTangI;

        // Create cubic Bezier with alpha=0.38, beta=0.38
        double d = dist(startPt, endPt);
        if (d < EPS) return {startPt}; // Degenerate case: return a single-point polyline to avoid zero-length segments
        BezierSegment cb = makeCubicG1(startPt, startTang, endPt, endTang, 0.38);

        // Sample adaptively for quality matching existing edges
        std::vector<Vec2d> pts = cb.sampleAdaptive(3.0, 2.0, 0.1);
        if (pts.empty()) {
            pts = cb.sampleCount(30);
        }

        // Ensure exact endpoints
        if (!pts.empty()) {
            pts.front() = startPt;
            pts.back() = endPt;
        }

        return pts;
    }

    // // Check shared midline curvature and return quality flags if too high
    // int checkSharedMidlineCurvature(
    //     const ConnectivityCurve& gclI,
    //     const ConnectivityCurve& gclI1,
    //     const std::string& sharedEnterEdgeId,
    //     const std::string& sharedExitEdgeId,
    //     const IntersectionInput& inp) const
    // {
    //     // Reconstruct the Bezier to check curvature (same logic as generateSharedMidline)
    //     Vec2d startPt{0,0};
    //     auto enterEdgeIt = inp.findEdge(sharedEnterEdgeId);
    //     if(enterEdgeIt) startPt = enterEdgeIt->connectionPt;
    //
    //     Vec2d endPt{0,0};
    //     auto exitEdgeIt = inp.findEdge(sharedExitEdgeId);
    //     if(exitEdgeIt) endPt = exitEdgeIt->connectionPt;
    //
    //     double d = dist(startPt, endPt);
    //     if(d < EPS) return 0;
    //
    //     auto enterClI = inp.findLane(gclI.entry_lane_id);
    //     auto enterClI1 = inp.findLane(gclI1.entry_lane_id);
    //     Vec2d startTangI = (enterClI) ? enterClI->tangentDir : Vec2d{0,1};
    //     Vec2d startTangI1 = (enterClI1) ? enterClI1->tangentDir : Vec2d{0,1};
    //     Vec2d startTangSum = startTangI + startTangI1;
    //     Vec2d startTang = (startTangSum.norm() > EPS) ? startTangSum.normalized() : startTangI;
    //
    //     auto exitClI = inp.findLane(gclI.exit_lane_id);
    //     auto exitClI1 = inp.findLane(gclI1.exit_lane_id);
    //     Vec2d endTangI = (exitClI)? exitClI->tangentDir : Vec2d{0,1};
    //     Vec2d endTangI1 = (exitClI1) ? exitClI1->tangentDir : Vec2d{0,1};
    //     Vec2d endTangSum = endTangI + endTangI1;
    //     Vec2d endTang = (endTangSum.norm() > EPS) ? endTangSum.normalized() : endTangI;
    //
    //     BezierSegment cb = makeCubicG1(startPt, startTang, endPt, endTang, 0.33);
    //     double maxK = cb.maxCurvature(30);
    //     if(maxK > bezier_maxCurvature * 2.0) {
    //         return QF_WARN_CURVATURE_HIGH;
    //     }
    //     return 0;
    // }

    // Generate an independent edge (left or right) using offset+smooth+bezierAlignEnds
    std::vector<Vec2d> generateIndependentEdge(
        const ConnectivityCurve& gcl, const Connectivity* conn, bool isLeft, const IntersectionInput& inp) const {
        // Find the enter group
        auto grpIt = inp.findGroup(conn->enterGroupId);
        if (!grpIt) {
            // Fallback: simple offset
            double hw = edgeLine_defaultLaneWidth * 0.5;
            return offsetPolyline(gcl.curve->sample(50), isLeft ? hw : -hw);
        }
        const LaneGroup& grp = *grpIt;

        // Estimate half-width
        double hw = estimateHalfWidth(conn->entry_lane_id, isLeft, grp, inp, true);

        // Generate offset polyline
        std::vector<Vec2d> pts = offsetPolyline(gcl.curve->sample(50), isLeft ? hw : -hw);

        // Remove duplicates
        removeDuplicates(pts, 0.02);

        if (pts.size() < 2) return pts;

        // Find edge endpoint and tangent at enter/exit
        std::string enterGrpId = conn->enterGroupId;
        std::string exitGrpId = conn->exitGroupId;

        Vec2d startPt, endPt;
        if (isLeft) {
            startPt = findEdgePtInGroup(conn->entry_lane_id, true, enterGrpId, inp, hw, true);
            // Exit: exit tangent points inward, so left from exit perspective = right from curve perspective (flipped)
            endPt = findEdgePtInGroup(conn->exit_lane_id, false, exitGrpId, inp, hw, false);
        } else {
            startPt = findEdgePtInGroup(conn->entry_lane_id, false, enterGrpId, inp, hw, true);
            endPt = findEdgePtInGroup(conn->exit_lane_id, true, exitGrpId, inp, hw, false);
        }
        Vec2d enterTang = getEdgeTangent(conn->entry_lane_id, inp, true);
        Vec2d exitTang = getEdgeTangent(conn->exit_lane_id, inp, false);

        // Smooth middle
        smoothMiddle(pts, 5);

        // Bezier align ends
        int smoothPts = std::max(6, (int)pts.size() / 4);
        bezierAlignEnds(pts, startPt, enterTang, endPt, exitTang, smoothPts);

        return pts;
    }

    // -----------------------------------------------
    // Remove near-duplicate consecutive points
    // -----------------------------------------------
    void removeDuplicates(std::vector<Vec2d>& pts, double minDist) const {
        if (pts.size() < 3) return;
        std::vector<Vec2d> out;
        out.push_back(pts.front());
        for (size_t i = 1; i < pts.size() - 1; ++i) {
            if (dist(pts[i], out.back()) >= minDist)
                out.push_back(pts[i]);
        }
        out.push_back(pts.back());
        pts = out;
    }

    // -----------------------------------------------
    // Offset polyline using bisector normals
    // offset > 0 = left, < 0 = right
    // -----------------------------------------------
    std::vector<Vec2d> offsetPolyline(const std::vector<Vec2d>& pts, double offset) const {
        if (pts.size() < 2) return pts;
        int n = (int)pts.size();
        std::vector<Vec2d> out;
        out.reserve(n);

        for (int i = 0; i < n; ++i) {
            Vec2d normal;
            if (i == 0) {
                Vec2d dir = (pts[1] - pts[0]).normalized();
                normal = rotLeft(dir);
            } else if (i == n - 1) {
                Vec2d dir = (pts[n - 1] - pts[n - 2]).normalized();
                normal = rotLeft(dir);
            } else {
                Vec2d d1 = (pts[i] - pts[i - 1]).normalized();
                Vec2d d2 = (pts[i + 1] - pts[i]).normalized();
                Vec2d n1 = rotLeft(d1);
                Vec2d n2 = rotLeft(d2);
                Vec2d avg = n1 + n2;
                if (avg.norm() < EPS) {
                    normal = n1;
                } else {
                    normal = avg.normalized();
                    double cosHalf = n1.dot(normal);
                    if (cosHalf > 0.3) {
                        double miter = std::min(1.0 / cosHalf, 3.0);
                        out.push_back(pts[i] + normal * (offset * miter));
                        continue;
                    }
                }
            }
            out.push_back(pts[i] + normal * offset);
        }
        return out;
    }

    // -----------------------------------------------
    // Moving average smoothing (preserve head/tail margins)
    // -----------------------------------------------
    void smoothMiddle(std::vector<Vec2d>& pts, int passes) const {
        int n = (int)pts.size();
        if (n < 5) return;
        int margin = std::max(2, n / 8);

        for (int pass = 0; pass < passes; ++pass) {
            std::vector<Vec2d> tmp = pts;
            for (int i = margin; i < n - margin; ++i) {
                Vec2d sum = pts[i] * 4.0;
                int cnt = 4;
                if (i - 1 >= 0) {
                    sum += pts[i - 1] * 2.0;
                    cnt += 2;
                }
                if (i + 1 < n) {
                    sum += pts[i + 1] * 2.0;
                    cnt += 2;
                }
                if (i - 2 >= 0) {
                    sum += pts[i - 2] * 1.0;
                    cnt += 1;
                }
                if (i + 2 < n) {
                    sum += pts[i + 2] * 1.0;
                    cnt += 1;
                }
                tmp[i] = sum / (double)cnt;
            }
            pts = tmp;
        }
    }

    // -----------------------------------------------
    // G1 Bezier refit at both ends + endpoint alignment
    // -----------------------------------------------
    void bezierAlignEnds(
        std::vector<Vec2d>& pts, const Vec2d& startPt, const Vec2d& startTang,
        const Vec2d& endPt, const Vec2d& endTang, int K) const {
        int n = (int)pts.size();
        if (n < 6) return;

        K = std::max(4, std::min(K, n / 3));

        pts.front() = startPt;
        pts.back() = endPt;

        // Head refit [0..K]
        {
            Vec2d p0 = startPt;
            Vec2d p3 = pts[K];
            double d = dist(p0, p3);
            if (d > EPS) {
                Vec2d t3;
                if (K + 1 < n) {
                    t3 = (pts[K + 1] - pts[K - 1]).normalized();
                } else {
                    t3 = (pts[K] - pts[K - 1]).normalized();
                }
                double alpha = 0.38;
                Vec2d p1 = p0 + startTang * (alpha * d);
                Vec2d p2 = p3 - t3 * (alpha * d);
                BezierSegment cb;
                cb.ctrl = {p0, p1, p2, p3};
                for (int i = 1; i < K; ++i) {
                    double t = (double)i / K;
                    pts[i] = cb.evaluate(t);
                }
            }
        }

        // Tail refit [n-1-K..n-1]
        {
            int startIdx = n - 1 - K;
            if (startIdx < K) startIdx = K;
            Vec2d p0 = pts[startIdx];
            Vec2d p3 = endPt;
            double d = dist(p0, p3);
            if (d > EPS) {
                Vec2d t0;
                if (startIdx > 0 && startIdx + 1 < n) {
                    t0 = (pts[startIdx] - pts[startIdx - 1]).normalized();
                } else {
                    t0 = (p3 - p0).normalized();
                }
                double alpha = 0.38;
                Vec2d p1 = p0 + t0 * (alpha * d);
                Vec2d p2 = p3 + endTang * (alpha * d);
                BezierSegment cb;
                cb.ctrl = {p0, p1, p2, p3};
                int count = n - 1 - startIdx;
                for (int i = 1; i < count; ++i) {
                    double t = (double)i / count;
                    pts[startIdx + i] = cb.evaluate(t);
                }
            }
        }

        pts.front() = startPt;
        pts.back() = endPt;
    }

    // -----------------------------------------------
    // Estimate half-width from edge endpoints
    // -----------------------------------------------
    double estimateHalfWidth(const std::string& entry_lane_id, bool isLeft,
                             const LaneGroup& grp, const IntersectionInput& inp, bool is_entryline) const {
        auto clit = inp.findLane(entry_lane_id);
        if (!clit) return edgeLine_defaultLaneWidth * 0.5;

        const Vec2d& clPt = getConnPoint(clit->geometry.points, is_entryline); //clit->connectionPt;
        Vec2d tangent = getConnTangent(clit->geometry.points, is_entryline); //clit->tangentDir;
        Vec2d normal = rotLeft(tangent);

        double bestDist = -1;
        for (auto& eid : grp.boundaries) {
            auto elit = inp.findEdge(eid);
            if (!elit) continue;
            if (!elit->geometry.points.empty()) continue;
            const Vec2d& ep = getConnPoint(elit->geometry.points, is_entryline); //elit->connectionPt;
            double lateral = (ep - clPt).dot(normal);
            if (isLeft && lateral > 0.01 && lateral < 10.0) {
                if (bestDist < 0 || lateral < bestDist) bestDist = lateral;
            } else if (!isLeft && lateral < -0.01 && lateral > -10.0) {
                if (bestDist < 0 || (-lateral) < bestDist) bestDist = -lateral;
            }
        }

        if (bestDist > 0) return bestDist;
        return edgeLine_defaultLaneWidth * 0.5;
    }

    // -----------------------------------------------
    // Find edge endpoint within a specific group
    // -----------------------------------------------
    Vec2d findEdgePtInGroup(
        const std::string& lineId, bool isLeft, const std::string& groupId,
        const IntersectionInput& inp, double hw, bool is_entryline) const {
        auto clit = inp.findLane(lineId);
        if (!clit) return {0, 0};

        const Vec2d& clPt = getConnPoint(clit->geometry.points, is_entryline); //clit->connectionPt;
        const Vec2d& tang = getConnTangent(clit->geometry.points, is_entryline); //clit->tangentDir;
        Vec2d normal = rotLeft(tang);

        auto git = inp.findGroup(groupId);
        if (git) {
            for (auto& eid : git->boundaries) {
                auto elit = inp.findEdge(eid);
                if (!elit) continue;
                if (!elit->geometry.points.empty()) continue;
                const Vec2d& ep = getConnPoint(elit->geometry.points, is_entryline); //elit->connectionPt;
                double lat = (ep - clPt).dot(normal);
                if (isLeft && lat > 0.01 && std::abs(lat - hw) < hw * 0.8) return ep;
                if (!isLeft && lat < -0.01 && std::abs(-lat - hw) < hw * 0.8) return ep;
            }
        }
        return clPt + normal * (isLeft ? hw : -hw);
    }

    // -----------------------------------------------
    // Get tangent direction for a centerline
    // -----------------------------------------------
    Vec2d getEdgeTangent(const std::string& lineId, const IntersectionInput& inp, bool is_entryline) const {
        auto clit = inp.findLane(lineId);
        if (!clit) return {0, 1};
        return getConnTangent(clit->geometry.points, is_entryline); //clit->tangentDir;
    }
};
