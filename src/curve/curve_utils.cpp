#include "curve_utils.h"
#include "optimizer/sdf_field.h"
#include "constraints/fence_check.h"
#include "utils.h"
#include <cmath>
#include <algorithm>

namespace isg {

double localCurvature(const Vec2d& a, const Vec2d& b, const Vec2d& c) {
    double ab = (b - a).norm(), bc = (c - b).norm(), ac = (c - a).norm();
    double area = 0.5 * std::abs(cross2d(b - a, c - a));
    double den = ab * bc * ac;
    return den > 1e-12 ? (2 * area / den) : 0.0;
}

Vec2d circumcenter(const Vec2d& a, const Vec2d& b, const Vec2d& c) {
    double ax = a[0], ay = a[1], bx = b[0], by = b[1], cx = c[0], cy = c[1];
    double D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if (std::abs(D) < 1e-12)return 0.5 * (a + c);
    double ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / D;
    double uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / D;
    return Vec2d(ux, uy);
}

std::vector<Vec2d> elasticBandSmooth(const std::vector<Vec2d>& pts, const SDFField& sdf,
                                     const Polygon2d& fence, double km, double ms, int mi, double msc) {
    auto res = pts;
    int N = (int)res.size();
    if (N < 3) return res;

    for (int iter = 0; iter < mi; ++iter) {
        bool ch = false;
        for (int i = 1; i < N - 1; ++i) {
            double k = localCurvature(res[i - 1], res[i], res[i + 1]);
            if (k <= km) continue;
            Vec2d cen = circumcenter(res[i - 1], res[i], res[i + 1]);
            Vec2d dir = res[i] - cen;
            if (dir.norm() < 1e-10) continue;
            dir.normalize();
            Vec2d cand = res[i] + dir * ms;
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(cand);
            if (_q.first >= msc && polygonContains(fence, cand)) {
                res[i] = cand;
                ch = true;
            }
        }
        if (!ch) break;
    }
    return res;
}

// 自适应版本: 根据曲线特征调整采样数
std::vector<Vec2d> elasticBandSmoothAdaptive(const std::vector<Vec2d>& pts, const SDFField& sdf,
                                             const Polygon2d& fence, double km, double ms, int mi, double msc) {
    if (pts.size() < 3) return pts;

    // 计算曲线长度确定初始采样数
    double curve_length = 0.0;
    for (size_t i = 1; i < pts.size(); ++i) {
        curve_length += (pts[i] - pts[i-1]).norm();
    }

    // 基于长度的基准采样数,带上下界
    int initial_samples = std::max(10, std::min(200, static_cast<int>(curve_length / 0.5)));

    // 估计曲率变化以指导自适应采样
    double total_curvature = 0.0;
    int curvature_samples = std::min(static_cast<int>(pts.size()), 20); // 限制采样数以提高效率
    for (int i = 1; i < curvature_samples - 1; ++i) {
        double k = localCurvature(pts[std::min(i-1, static_cast<int>(pts.size()-1))],
                                  pts[i],
                                  pts[std::min(i+1, static_cast<int>(pts.size()-1))]);
        total_curvature += k;
    }

    // 根据曲率调整采样数
    double avg_curvature = total_curvature / std::max(1, curvature_samples - 2);
    int adaptive_samples = std::max(10, static_cast<int>(initial_samples * (1.0 + avg_curvature * 0.5)));

    // 点过多时抽稀
    std::vector<Vec2d> current_pts = pts;
    if (current_pts.size() > adaptive_samples) {
        current_pts = downsamplePoints(current_pts, adaptive_samples);
    }

    int N = (int)current_pts.size();
    if (N < 3) return current_pts;

    // 基于改进阈值的早停
    double improvement_threshold = 1e-4;
    std::vector<Vec2d> prev_pts = current_pts;

    for (int iter = 0; iter < mi; ++iter) {
        bool ch = false;
        for (int i = 1; i < N - 1; ++i) {
            double k = localCurvature(current_pts[i - 1], current_pts[i], current_pts[i + 1]);
            if (k <= km) continue;
            Vec2d cen = circumcenter(current_pts[i - 1], current_pts[i], current_pts[i + 1]);
            Vec2d dir = current_pts[i] - cen;
            if (dir.norm() < 1e-10) continue;
            dir.normalize();
            Vec2d cand = current_pts[i] + dir * ms;
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(cand);
            if (_q.first >= msc && polygonContains(fence, cand)) {
                current_pts[i] = cand;
                ch = true;
            }
        }

        // Early termination if no significant improvement
        if (!ch) break;

        // Check for convergence: if points didn't move significantly
        double max_displacement = 0.0;
        for (int i = 0; i < N; ++i) {
            double displacement = (current_pts[i] - prev_pts[i]).norm();
            max_displacement = std::max(max_displacement, displacement);
        }

        if (max_displacement < improvement_threshold) {
            break;
        }

        prev_pts = current_pts;
    }

    return current_pts;
}

// Helper function to downsample points evenly
std::vector<Vec2d> downsamplePoints(const std::vector<Vec2d>& points, int target_count) {
    if (points.size() <= target_count) {
        return points;
    }

    std::vector<Vec2d> result;
    result.reserve(target_count);

    double step = static_cast<double>(points.size() - 1) / (target_count - 1);
    for (int i = 0; i < target_count; i++) {
        int idx = std::min(static_cast<int>(i * step), static_cast<int>(points.size() - 1));
        result.push_back(points[idx]);
    }

    return result;
}

BezierCurve rebuildFromSmoothedPts(const std::vector<Vec2d>& s, const Vec2d& st, const Vec2d& et) {
    return fitBezierWithEndTangents(s, st, et);
}

std::vector<Vec2d> midlineSampleByArcLength(const BezierCurve& a, const BezierCurve& b, int n) {
    auto sa = a.sampleByArcLength(n), sb = b.sampleByArcLength(n);
    std::vector<Vec2d> m(n);
    for (int i = 0; i < n; ++i)m[i] = 0.5 * (sa[i] + sb[i]);
    return m;
}

double signedDistToLine(const Vec2d& pt, const Vec2d& p0, const Vec2d& p1) {
    Vec2d d = p1 - p0;
    double len = d.norm();
    return len < 1e-12 ? (pt - p0).norm() : cross2d(d, pt - p0) / len;
}

bool segmentsIntersect(const Vec2d& a, const Vec2d& b, const Vec2d& c, const Vec2d& d, Vec2d* out) {
    Vec2d r = b - a, s = d - c;
    double den = cross2d(r, s);
    if (std::abs(den) < 1e-12)return false;
    Vec2d ac = c - a;
    double t = cross2d(ac, s) / den, u = cross2d(ac, r) / den;
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        if (out)*out = a + t * r;
        return true;
    }
    return false;
}

double minSDFAlongCurve(const BezierCurve& c, const SDFField& sdf, int sps) {
    double m = 1e18;
    for (auto& seg : c.segs)
        for (int i = 0; i <= sps; ++i) {
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(seg.evaluate((double)i / sps));
            m = std::min(m, _q.first);
        }
    return m;
}

double distToAllEndpoints(const Vec2d& pt, const BezierCurve& a, const BezierCurve& b) {
    double d = 1e18;
    if (!a.empty()) {
        d = std::min(d, (pt - a.startPt()).norm());
        d = std::min(d, (pt - a.endPt()).norm());
    }
    if (!b.empty()) {
        d = std::min(d, (pt - b.startPt()).norm());
        d = std::min(d, (pt - b.endPt()).norm());
    }
    return d;
}


// src/constraints/intersection_check.cpp
bool bboxOverlap(const BezierCurve& a, const BezierCurve& b) {
    return a.bbox().intersects(b.bbox());
}

static void segPairCross(const BezierSegment& a, const BezierSegment& b,
                         std::vector<Vec2d>& out, double tol, int depth = 0) {
    if (!a.bbox().intersects(b.bbox()))
        return;
    if (depth > 10) {
        out.push_back(a.evaluate(0.5));
        return;
    }
    if (a.arcLength(4) < tol && b.arcLength(4) < tol) {
        out.push_back(a.evaluate(0.5));
        return;
    }
    std::pair<BezierSegment, BezierSegment> _sa = a.splitAt(0.5);
    std::pair<BezierSegment, BezierSegment> _sb = b.splitAt(0.5);
    BezierSegment al = _sa.first, ar = _sa.second;
    BezierSegment bl = _sb.first, br = _sb.second;
    segPairCross(al, bl, out, tol, depth + 1);
    segPairCross(al, br, out, tol, depth + 1);
    segPairCross(ar, bl, out, tol, depth + 1);
    segPairCross(ar, br, out, tol, depth + 1);
}

std::vector<Vec2d> curveCrossings(const BezierCurve& a, const BezierCurve& b, double tol) {
    std::vector<Vec2d> raw;
    if (!bboxOverlap(a, b))
        return raw;
    for (auto& sa : a.segs)
        for (auto& sb : b.segs)
            segPairCross(sa, sb, raw, tol);
    std::vector<Vec2d> out;
    for (auto& p : raw) {
        bool dup = false;
        for (auto& q : out)
            if ((p - q).norm() < tol) {
                dup = true;
                break;
            }
        if (!dup)
            out.push_back(p);
    }
    return out;
}

bool curvesIntersectBusiness(const BezierCurve& a, const BezierCurve& b, double ep) {
    if (!bboxOverlap(a, b))
        return false;
    auto adaptiveSample = [](const BezierCurve& c) {
        int n = std::max(48, (int)std::ceil(c.arcLength() / 0.20) + 1);
        n = std::min(n, 240);
        return c.sampleByArcLength(n);
    };
    auto pa = adaptiveSample(a);
    auto pb = adaptiveSample(b);
    for (int i = 0; i + 1 < (int)pa.size(); ++i) {
        for (int j = 0; j + 1 < (int)pb.size(); ++j) {
            Vec2d isect;
            if (!segmentsIntersect(pa[i], pa[i + 1], pb[j], pb[j + 1], &isect))
                continue;
            if (distToAllEndpoints(isect, a, b) > ep)
                return true;
        }
    }
    return false;
}

bool curveSelfIntersectsBusiness(const BezierCurve& c, double ep) {
    // 性能优化: 采样从80降至40 (O(N²) → 减少4倍)
    auto pts = c.sample(40);
    if (pts.size() < 4)
        return false;
    for (int i = 0; i + 1 < (int)pts.size(); ++i) {
        for (int j = i + 2; j + 1 < (int)pts.size(); ++j) {
            if (i == 0 && j + 1 == (int)pts.size() - 1)
                continue;
            Vec2d isect;
            if (!segmentsIntersect(pts[i], pts[i + 1], pts[j], pts[j + 1], &isect))
                continue;
            double d0 = std::min((isect - pts.front()).norm(), (isect - pts.back()).norm());
            if (d0 <= ep)
                continue;
            return true;
        }
    }
    return false;
}

}
