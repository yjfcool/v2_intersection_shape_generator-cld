#include "curve_utils.h"
#include "optimizer/sdf_field.h"
#include "constraints/fence_check.h"
#include "utils.h"
#include <cmath>
#include <algorithm>

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
    if (N < 3)return res;
    for (int iter = 0; iter < mi; ++iter) {
        bool ch = false;
        for (int i = 1; i < N - 1; ++i) {
            double k = localCurvature(res[i - 1], res[i], res[i + 1]);
            if (k <= km)continue;
            Vec2d cen = circumcenter(res[i - 1], res[i], res[i + 1]);
            Vec2d dir = res[i] - cen;
            if (dir.norm() < 1e-10)continue;
            dir.normalize();
            Vec2d cand = res[i] + dir * ms;
            std::pair<double, Vec2d> _q = sdf.queryWithGrad(cand);
            if (_q.first >= msc && polygonContains(fence, cand)) {
                res[i] = cand;
                ch = true;
            }
        }
        if (!ch)break;
    }
    return res;
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
    auto cr = curveCrossings(a, b);
    for (auto& pt : cr)
        if (distToAllEndpoints(pt, a, b) > ep)
            return true;
    return false;
}