#include "bezier.h"
#include "optimizer/sdf_field.h"
#include "utils.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>

namespace isg {

Vec2d BezierSegment::evaluate(double t) const {
    double u = 1 - t, u2 = u * u, u3 = u2 * u, t2 = t * t, t3 = t2 * t;
    return u3 * ctrl[0] + 3 * u2 * t * ctrl[1] + 3 * u * t2 * ctrl[2] + t3 * ctrl[3];
}

Vec2d BezierSegment::tangent(double t) const {
    double u = 1 - t;
    return 3 * (u * u * (ctrl[1] - ctrl[0]) + 2 * u * t * (ctrl[2] - ctrl[1]) + t * t * (ctrl[3] - ctrl[2]));
}

double BezierSegment::curvature(double t) const {
    Vec2d d1 = tangent(t);
    double u = 1 - t;
    Vec2d d2 = 6 * (u * (ctrl[2] - 2 * ctrl[1] + ctrl[0]) + t * (ctrl[3] - 2 * ctrl[2] + ctrl[1]));
    double num = std::abs(cross2d(d1, d2)), den = std::pow(d1.norm(), 3.0);
    return den < 1e-12 ? 0.0 : num / den;
}

std::pair<BezierSegment, BezierSegment> BezierSegment::splitAt(double t) const {
    Vec2d p01 = (1 - t) * ctrl[0] + t * ctrl[1], p12 = (1 - t) * ctrl[1] + t * ctrl[2];
    Vec2d p23 = (1 - t) * ctrl[2] + t * ctrl[3];
    Vec2d p012 = (1 - t) * p01 + t * p12, p123 = (1 - t) * p12 + t * p23;
    Vec2d p = (1 - t) * p012 + t * p123;
    BezierSegment L, R;
    L.ctrl = {ctrl[0], p01, p012, p};
    R.ctrl = {p, p123, p23, ctrl[3]};
    return {L, R};
}

BoundingBox2d BezierSegment::bbox() const {
    BoundingBox2d b;
    for (int i = 0; i <= 20; ++i)b.expand(evaluate(i / 20.0));
    return b;
}

double BezierSegment::arcLength(int s) const {
    double len = 0;
    Vec2d prev = ctrl[0];
    for (int i = 1; i <= s; ++i) {
        Vec2d c = evaluate((double)i / s);
        len += (c - prev).norm();
        prev = c;
    }
    return len;
}

double BezierSegment::arcLengthToParam(double s, int samp) const {
    double tot = arcLength(samp), tgt = s * tot, acc = 0;
    Vec2d prev = ctrl[0];
    double dt = 1.0 / samp;
    for (int i = 1; i <= samp; ++i) {
        double t = i * dt;
        Vec2d cur = evaluate(t);
        double sl = (cur - prev).norm();
        if (acc + sl >= tgt) {
            double f = (tgt - acc) / (sl + 1e-12);
            return (i - 1) * dt + f * dt;
        }
        acc += sl;
        prev = cur;
    }
    return 1.0;
}

// B'(t) 一阶导数
Vec2d BezierSegment::evalDeriv1(double t) const {
    double u = 1.0 - t;
    return (ctrl[1] - ctrl[0]) * (3 * u * u) + (ctrl[2] - ctrl[1]) * (6 * u * t) + (ctrl[3] - ctrl[2]) * (3 * t *
        t);
}

// B''(t) 二阶导数
Vec2d BezierSegment::evalDeriv2(double t) const {
    return (ctrl[2] - ctrl[1] * 2.0 + ctrl[0]) * (6 * (1 - t)) + (ctrl[3] - ctrl[2] * 2.0 + ctrl[1]) * (6 * t);
}

// 最大曲率（采样估算）
double BezierSegment::maxCurvature(int samples) const {
    double maxK = 0;
    for (int i = 0; i <= samples; ++i) {
        maxK = std::max(maxK, curvature(i * 1.0 / samples));
    }
    return maxK;
}

// 固定数量采样
std::vector<Vec2d> BezierSegment::sampleCount(int n) const {
    std::vector<Vec2d> pts;
    pts.reserve(n + 1);
    for (int i = 0; i <= n; ++i)
        pts.push_back(evaluate(i * 1.0 / n));
    return pts;
}

// 固定间距采样（弧长近似）
std::vector<Vec2d> BezierSegment::sampleBySpacing(double spacing) const {
    // 先粗采样，再按间距重采样
    std::vector<Vec2d> dense = sampleCount(200);
    return resampleBySpacing(dense, spacing);
}

void BezierSegment::subdivide(double t0, double t1, double maxAngle, double maxSeg, double minSeg,
                              std::vector<Vec2d>& pts, int depth) const {
    if (depth > 20) {
        pts.push_back(evaluate(t1));
        return;
    }
    Vec2d p0 = evaluate(t0);
    Vec2d p1 = evaluate(t1);
    double segLen = dist(p0, p1);
    if (segLen < minSeg) {
        pts.push_back(p1);
        return;
    }
    double tMid = 0.5 * (t0 + t1);
    Vec2d pMid = evaluate(tMid);

    // 检查角度误差
    Vec2d d01 = (pMid - p0).normalized();
    Vec2d d12 = (p1 - pMid).normalized();
    double angle = std::acos(std::max(-1.0, std::min(1.0, d01.dot(d12))));
    bool needSplit = (angle > maxAngle) || (segLen > maxSeg);
    if (!needSplit) {
        pts.push_back(p1);
    } else {
        subdivide(t0, tMid, maxAngle, maxSeg, minSeg, pts, depth + 1);
        subdivide(tMid, t1, maxAngle, maxSeg, minSeg, pts, depth + 1);
    }
};
// 自适应采样（基于角度误差）
std::vector<Vec2d> BezierSegment::sampleAdaptive(double maxAngleDeg, double maxSeg, double minSeg) const {
    double maxAngleRad = maxAngleDeg * DEG2RAD;
    std::vector<Vec2d> pts;
    pts.push_back(ctrl[0]);
    subdivide(0.0, 1.0, maxAngleRad, maxSeg, minSeg, pts, 0);
    // 确保端点精确
    if (pts.empty() || dist(pts.back(), ctrl[3]) > EPS)
        pts.push_back(ctrl[3]);
    return pts;
}

// 形态自适应采样（根据曲率调整采样密度）
std::vector<Vec2d> BezierSegment::sampleByShape(
    double straightSpacing, double angularResolution, double minSpacing) const {
    std::vector<Vec2d> pts;

    // 首先检查是否近似为直线段
    // 计算控制点形成的总弯曲角度
    Vec2d v1 = (ctrl[1] - ctrl[0]).normalized();
    Vec2d v2 = (ctrl[3] - ctrl[2]).normalized();
    double straightAngle = angleBetween(v1, -v2);

    // 如果几乎成一直线（接近180度），使用直线采样策略
    if (std::abs(straightAngle - PI) < 0.1) { // 0.1 radian ≈ 5.7 degrees
        // 对直线段使用固定间隔采样
        double totalLength = arcLength(20); // 估算弧长
        int nPoints = std::max(2, static_cast<int>(std::ceil(totalLength / straightSpacing)));

        for (int i = 0; i < nPoints; i++) {
            double t = static_cast<double>(i) / (nPoints - 1);
            pts.push_back(evaluate(t));
        }
    } else {
        // 对于弯曲段，采用更复杂的方法
        pts.push_back(ctrl[0]); // 添加起点

        double angularResRad = angularResolution * DEG2RAD;
        int samples = 100; // 采样精度
        double dt = 1.0 / samples;

        Vec2d lastPt = evaluate(0.0);
        double accumulatedLength = 0.0;

        for (int i = 1; i <= samples; i++) {
            double t = i * dt;
            Vec2d currPt = evaluate(t);
            double segmentLength = (currPt - lastPt).norm();

            // 计算当前位置的曲率
            double curvature = this->curvature(t);

            double targetSpacing;
            if (curvature < 0.01) { // 几乎是直线
                targetSpacing = straightSpacing;
            } else {
                // 对弯曲部分：使用曲率半径和角度分辨率计算弧长间隔
                double radius = 1.0 / curvature;
                targetSpacing = radius * angularResRad;

                // 不超过直线最大间距限制，且不低于最小间距
                targetSpacing = std::min(targetSpacing, straightSpacing);
                targetSpacing = std::max(targetSpacing, minSpacing);
            }

            // 如果累计长度达到目标间距，则添加点
            accumulatedLength += segmentLength;
            if (accumulatedLength >= targetSpacing) {
                pts.push_back(currPt);
                lastPt = currPt;
                accumulatedLength = 0.0;
            } else {
                lastPt = currPt;
            }
        }
    }

    // 确保终点被包含且满足最小间距
    if (pts.empty() || dist(pts.back(), ctrl[3]) > minSpacing * 0.5) { // 使用一半的最小间距作为容差
        // 检查最后添加的点是否离终点太近，如果太近则替换为终点
        if (!pts.empty() && dist(pts.back(), ctrl[3]) <= minSpacing) {
            pts.back() = ctrl[3]; // 替换为精确终点
        } else {
            pts.push_back(ctrl[3]); // 添加新点
        }
    }

    // 清理重复点（如果有的话）
    std::vector<Vec2d> cleanedPts;
    if (!pts.empty()) {
        cleanedPts.push_back(pts[0]);
        for (size_t i = 1; i < pts.size(); i++) {
            if (dist(cleanedPts.back(), pts[i]) > minSpacing * 0.5) {
                cleanedPts.push_back(pts[i]);
            } else {
                // 更新最后一个点为当前点（可能是终点），确保精度
                if (i == pts.size() - 1) { // 如果是终点
                    cleanedPts.back() = pts[i];
                }
            }
        }
    }

    return cleanedPts;
}

// 获取 α（P1相对P0的标量，T0方向）
double BezierSegment::getAlpha(const Vec2d& T0) const {
    double d = dist(ctrl[0], ctrl[3]);
    if (d < EPS)
        return 0.35;
    return (ctrl[1] - ctrl[0]).dot(T0) / d;
}

double BezierSegment::getBeta(const Vec2d& T3) const {
    double d = dist(ctrl[0], ctrl[3]);
    if (d < EPS)
        return 0.35;
    // P2 = P3 + T3*beta*d → beta = (P2-P3)·T3 / d
    return (ctrl[2] - ctrl[3]).dot(T3.normalized()) / d;
}


Vec2d BezierCurve::startPt() const { return segs.front().ctrl[0]; }
Vec2d BezierCurve::endPt() const { return segs.back().ctrl[3]; }

Vec2d BezierCurve::startTan() const {
    Vec2d d = segs.front().ctrl[1] - segs.front().ctrl[0];
    return d.norm() > 1e-10 ? d.normalized() : Vec2d(1, 0);
}

Vec2d BezierCurve::endTan() const {
    Vec2d d = segs.back().ctrl[3] - segs.back().ctrl[2];
    return d.norm() > 1e-10 ? d.normalized() : Vec2d(1, 0);
}

Vec2d BezierCurve::evaluate(double u) const {
    assert(!segs.empty());
    int n = (int)segs.size();
    double sc = u * n;
    int idx = std::min((int)sc, n - 1);
    return segs[idx].evaluate(sc - idx);
}

Vec2d BezierCurve::tangent(double u) const {
    assert(!segs.empty());
    int n = (int)segs.size();
    double sc = u * n;
    int idx = std::min((int)sc, n - 1);
    Vec2d t = segs[idx].tangent(sc - idx);
    return t.norm() > 1e-10 ? t.normalized() : Vec2d(1, 0);
}

std::vector<Vec2d> BezierCurve::sample(int n) const {
    std::vector<Vec2d> out(n);
    for (int i = 0; i < n; ++i)out[i] = evaluate((double)i / (n - 1));
    return out;
}

double BezierCurve::arcLength() const {
    double t = 0;
    for (auto& s : segs)t += s.arcLength();
    return t;
}

std::vector<Vec2d> BezierCurve::sampleByArcLength(int n) const {
    if (n <= 1) return {startPt()};
    if (n == 2) return {startPt(), endPt()};
    double tot = arcLength(), step = tot / (n - 1);
    std::vector<Vec2d> out;
    out.push_back(startPt());
    double acc = 0, tgt = step;
    for (auto& seg : segs) {
        int sub = 50;
        Vec2d prev = seg.ctrl[0];
        double dt = 1.0 / sub;
        for (int i = 1; i <= sub; ++i) {
            Vec2d cur = seg.evaluate(i * dt);
            double d = (cur - prev).norm();
            while (acc + d >= tgt && (int)out.size() < n - 1) {
                // stop at n-1, reserve last slot
                double f = (tgt - acc) / (d + 1e-12);
                out.push_back(prev + f * (cur - prev));
                tgt += step;
            }
            acc += d;
            prev = cur;
        }
    }
    // Always end exactly at endPt() — fixes floating-point endpoint drift
    while ((int)out.size() >= n) out.pop_back();
    out.push_back(endPt());
    return out;
}

std::vector<Vec2d> BezierCurve::sampleByShape(double straightSpacing, double angularResolution, double minSpacing) const {
    if (segs.empty()) return {};

    std::vector<Vec2d> result;

    // 处理每个曲线段
    for (size_t i = 0; i < segs.size(); ++i) {
        const BezierSegment& seg = segs[i];

        // 用形状自适应方法对当前段进行采样
        std::vector<Vec2d> segSamples = seg.sampleByShape(straightSpacing, angularResolution, minSpacing);

        // 如果是第一段，添加所有点；如果是后续段，跳过起点以避免重复
        if (i == 0) {
            result.insert(result.end(), segSamples.begin(), segSamples.end());
        } else {
            // 跳过重复的连接点
            if (!segSamples.empty()) {
                result.insert(result.end(), segSamples.begin() + 1, segSamples.end());
            }
        }
    }

    // 确保起点和终点都被包含（如果不在结果中的话）
    if (result.empty()) {
        result.push_back(startPt());
        result.push_back(endPt());
    } else {
        if (dist(result.front(), startPt()) > minSpacing) {
            result.insert(result.begin(), startPt());
        }
        if (dist(result.back(), endPt()) > minSpacing) {
            result.push_back(endPt());
        }
    }

    return result;
}

BoundingBox2d BezierCurve::bbox() const {
    BoundingBox2d b;
    for (auto& s : segs) {
        auto sb = s.bbox();
        b.expand(sb.min_pt);
        b.expand(sb.max_pt);
    }
    return b;
}

double BezierCurve::maxCurvature(int sps) const {
    double km = 0;
    for (auto& s : segs)for (int i = 0; i <= sps; ++i)km = std::max(km, s.curvature((double)i / sps));
    return km;
}

BezierSegment makeCubicG1(const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1, double alpha) {
    double d = (p1 - p0).norm();
    BezierSegment s;
    s.ctrl[0] = p0;
    s.ctrl[1] = p0 + alpha * d * t0.normalized();
    s.ctrl[2] = p1 - alpha * d * t1.normalized();
    s.ctrl[3] = p1;
    return s;
}

BezierSegment makeAlignedUTurnCubic(
    const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1,
    double handle_scale, double handle_bias) {
    Vec2d T0 = t0.norm() > 1e-8 ? t0.normalized() : Vec2d(1, 0);
    Vec2d T1 = t1.norm() > 1e-8 ? t1.normalized() : -T0;

    Vec2d exit_back = -T1;
    Vec2d axis = T0 + exit_back;
    if (axis.norm() < 1e-8)
        axis = T0;
    axis.normalize();
    if (axis.dot(T0) < 0.0)
        axis = -axis;

    double s0 = p0.dot(axis);
    double s1 = p1.dot(axis);
    double common_s = std::max(s0, s1);

    double c0 = std::max(0.2, T0.dot(axis));
    double c1 = std::max(0.2, exit_back.dot(axis));
    double lead0 = std::max(0.0, (common_s - s0) / c0);
    double lead1 = std::max(0.0, (common_s - s1) / c1);

    Vec2d q0 = p0 + lead0 * T0;
    Vec2d q1 = p1 + lead1 * exit_back;
    Vec2d sep = q1 - q0;
    Vec2d lateral = sep - sep.dot(axis) * axis;
    double turn_gap = lateral.norm();
    if (turn_gap < 1e-6) {
        Vec2d perp{-axis[1], axis[0]};
        turn_gap = std::abs((p1 - p0).dot(perp));
    }
    turn_gap = std::max(turn_gap, 0.15);

    double gap_scale = 1.0;
    if (turn_gap > 5.0) {
        double u = std::min(1.0, (turn_gap - 5.0) / 9.0);
        gap_scale += 0.50 * (u * u * (3.0 - 2.0 * u));
    }
    double scale = std::max(0.15, handle_scale * gap_scale);
    double arc_handle = (2.0 / 3.0) * turn_gap * scale + handle_bias;
    arc_handle = std::max(0.05, arc_handle);

    BezierSegment s;
    s.ctrl[0] = p0;
    s.ctrl[1] = p0 + (lead0 + arc_handle) * T0;
    s.ctrl[2] = p1 - (lead1 + arc_handle) * T1;
    s.ctrl[3] = p1;
    return s;
}

BezierCurve makeCurveFromKnots(const std::vector<Vec2d>& pts, const std::vector<Vec2d>& tans, double alpha) {
    assert(pts.size()==tans.size()&&pts.size()>=2);
    BezierCurve c;
    for (int i = 0; i + 1 < (int)pts.size(); ++i)
        c.segs.push_back(makeCubicG1(pts[i], tans[i], pts[i + 1], tans[i + 1], alpha));
    return c;
}

static Vec2d catmullTan(const std::vector<Vec2d>& pts, int i) {
    int n = (int)pts.size();
    if (n < 2)return Vec2d(1, 0);
    if (i == 0)return (pts[1] - pts[0]).normalized();
    if (i == n - 1)return (pts[n - 1] - pts[n - 2]).normalized();
    Vec2d t = 0.5 * (pts[i + 1] - pts[i - 1]);
    return t.norm() > 1e-10 ? t.normalized() : (pts[i + 1] - pts[i]).normalized();
}

BezierCurve fitBezierWithEndTangents(const std::vector<Vec2d>& pts, const Vec2d& st, const Vec2d& et) {
    if (pts.size() < 2)return {};
    std::vector<Vec2d> tans(pts.size());
    tans.front() = st.normalized();
    tans.back() = et.normalized();
    for (int i = 1; i + 1 < (int)pts.size(); ++i)tans[i] = catmullTan(pts, i);
    return makeCurveFromKnots(pts, tans);
}

AdaptiveRefineResult adaptiveRefine(const BezierCurve& c, const SDFField& sdf, double kmax, double ptol) {
    AdaptiveRefineResult res;
    res.curve.segs = c.segs;
    for (int iter = 0; iter < 4; ++iter) {
        std::vector<BezierSegment> next;
        bool split = false;
        for (auto& seg : res.curve.segs) {
            double km = 0;
            for (int i = 0; i <= 20; ++i)km = std::max(km, seg.curvature(i / 20.0));
            double ms = 1e18;
            for (int i = 0; i <= 20; ++i) {
                std::pair<double, Vec2d> _q = sdf.queryWithGrad(seg.evaluate(i / 20.0));
                ms = std::min(ms, _q.first);
            }
            if (km > kmax || ms < ptol) {
                double wt = 0.5, wv = -1e18;
                for (int i = 1; i < 20; ++i) {
                    double t = i / 20.0;
                    double v = seg.curvature(t) - 0.5 * sdf.queryWithGrad(seg.evaluate(t)).first;
                    if (v > wv) {
                        wv = v;
                        wt = t;
                    }
                }
                wt = std::max(0.1, std::min(0.9, wt));
                std::pair<BezierSegment, BezierSegment> _sp = seg.splitAt(wt);
                next.push_back(_sp.first);
                next.push_back(_sp.second);
                split = res.was_split = true;
            }
            else next.push_back(seg);
        }
        res.curve.segs = next;
        if (!split)break;
    }
    return res;
}

// ── Compact encoding: only inner ctrl[1],ctrl[2] — join pts frozen ──────────
// params size = 4*N (2 inner ctrl pts × 2 coords per segment)
VecXd curveToParams(const BezierCurve& c) {
    int n = (int)c.segs.size();
    VecXd p(4 * n);
    for (int i = 0; i < n; ++i) {
        p[4 * i + 0] = c.segs[i].ctrl[1][0];
        p[4 * i + 1] = c.segs[i].ctrl[1][1];
        p[4 * i + 2] = c.segs[i].ctrl[2][0];
        p[4 * i + 3] = c.segs[i].ctrl[2][1];
    }
    return p;
}

BezierCurve curveFromParams(const VecXd& params, const BezierCurve& proto) {
    BezierCurve c;
    int n = (int)proto.segs.size();
    c.segs = proto.segs;
    for (int i = 0; i < n; ++i) {
        c.segs[i].ctrl[1] = Vec2d(params[4 * i + 0], params[4 * i + 1]);
        c.segs[i].ctrl[2] = Vec2d(params[4 * i + 2], params[4 * i + 3]);
    }
    for (int i = 1; i < n; ++i)c.segs[i].ctrl[0] = proto.segs[i].ctrl[0];
    return c;
}

// ── Full encoding: join pts + inner ctrl pts all free (Level-2) ───────────
// Layout: seg[0]: ctrl[1],ctrl[2] = 4 vals
//         seg[k>0]: ctrl[0],ctrl[1],ctrl[2] = 6 vals (join pt free)
// params size = 4 + 6*(N-1)
VecXd curveToParamsFull(const BezierCurve& c) {
    int n = (int)c.segs.size();
    int sz = (n == 0) ? 0 : 4 + 6 * std::max(0, n - 1);
    VecXd p(sz);
    if (n == 0)return p;
    // seg 0: only inner pts
    p[0] = c.segs[0].ctrl[1][0];
    p[1] = c.segs[0].ctrl[1][1];
    p[2] = c.segs[0].ctrl[2][0];
    p[3] = c.segs[0].ctrl[2][1];
    for (int i = 1; i < n; ++i) {
        int b = 4 + 6 * (i - 1);
        p[b + 0] = c.segs[i].ctrl[0][0];
        p[b + 1] = c.segs[i].ctrl[0][1]; // join
        p[b + 2] = c.segs[i].ctrl[1][0];
        p[b + 3] = c.segs[i].ctrl[1][1];
        p[b + 4] = c.segs[i].ctrl[2][0];
        p[b + 5] = c.segs[i].ctrl[2][1];
    }
    return p;
}

// Reconstruct from full params.  G1 at joins enforced by linking ctrl[3] of
// seg[k] to ctrl[0] of seg[k+1] AND aligning ctrl[2] of seg[k] with
// ctrl[0..1] of seg[k+1] so tangents are parallel.
BezierCurve curveFromParamsFull(const VecXd& params, const BezierCurve& proto) {
    int n = (int)proto.segs.size();
    BezierCurve c;
    c.segs = proto.segs;
    if (n == 0)return c;
    // seg 0: endpoint and inner ctrl pts
    // ctrl[0] of seg[0] is ALWAYS proto start — never modified
    c.segs[0].ctrl[1] = Vec2d(params[0], params[1]);
    c.segs[0].ctrl[2] = Vec2d(params[2], params[3]);
    c.segs[0].ctrl[0] = proto.segs[0].ctrl[0]; // pin start

    for (int i = 1; i < n; ++i) {
        int b = 4 + 6 * (i - 1);
        Vec2d join(params[b + 0], params[b + 1]);
        Vec2d c1(params[b + 2], params[b + 3]);
        Vec2d c2(params[b + 4], params[b + 5]);
        // Set join point: ctrl[3] of seg[i-1] = ctrl[0] of seg[i] = join
        c.segs[i - 1].ctrl[3] = join;
        c.segs[i].ctrl[0] = join;
        c.segs[i].ctrl[1] = c1;
        c.segs[i].ctrl[2] = c2;
        // Enforce G1 at join: ctrl[2] of seg[i-1] collinear with
        // join→ctrl[1] of seg[i], opposite direction
        Vec2d tan_out = (c1 - join);
        double tn = tan_out.norm();
        if (tn > 1e-10) {
            Vec2d tan_dir = tan_out * (1.0 / tn);
            Vec2d old_c2 = c.segs[i - 1].ctrl[2];
            double scale = (join - old_c2).norm();
            if (scale < 1e-6)scale = tn * 0.33;
            c.segs[i - 1].ctrl[2] = join - tan_dir * scale;
        }
    }
    // Pin final endpoint exactly to proto — prevents floating-point drift
    c.segs.back().ctrl[3] = proto.segs.back().ctrl[3];
    // Pin first endpoint too
    c.segs.front().ctrl[0] = proto.segs.front().ctrl[0];
    return c;
}

}
