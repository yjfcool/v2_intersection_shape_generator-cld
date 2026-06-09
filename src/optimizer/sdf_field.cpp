#include "sdf_field.h"
#include "constraints/fence_check.h"
#include "utils/clipper.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <functional>

Polygon2d SDFField::bufferPolygon(const Polygon2d& poly, double r) {
    if (r <= 0 || poly.outer.empty()) return poly;
    Polygon2d out;
    auto inf = ClipperUtil::InflatePaths({toArray(poly.outer)}, r, ClipperLib::jtMiter, ClipperLib::etClosedPolygon);
    out.outer = toArray(inf[0]);
    return out;
}

Polygon2d SDFField::unionPolygons(const std::vector<Polygon2d>& polys) {
    if (polys.empty()) return {};
    if (polys.size() == 1) return polys[0];
    std::vector<std::vector<std::array<double, 2>>> subs;
    for (auto& p : polys) if (!p.outer.empty()) subs.emplace_back(toArray(p.outer));
    auto sol = ClipperUtil::UnionPaths(subs, ClipperLib::pftNonZero);
    if (sol.empty()) return polys[0];
    Polygon2d out;
    out.outer = toArray(sol[0]);
    return out;
}

static double ptSegDist(const Vec2d& p, const Vec2d& a, const Vec2d& b) {
    Vec2d ab = b - a, ap = p - a;
    double t = ab.dot(ap);
    double l2 = ab.squaredNorm();
    if (l2 < 1e-20) return ap.norm();
    t = std::max(0.0, std::min(1.0, t / l2));
    return (p - (a + t * ab)).norm();
}

double SDFField::distToPolygons(const Vec2d& pt, const std::vector<Polygon2d>& polys) const {
    double m = std::numeric_limits<double>::max();
    for (auto& poly : polys) {
        auto& ring = poly.outer;
        int n = (int)ring.size();
        for (int i = 0; i < n; ++i)
            m = std::min(m, ptSegDist(pt, ring[i], ring[(i + 1) % n]));
    }
    return m;
}

bool SDFField::insideAny(const Vec2d& pt, const std::vector<Polygon2d>& polys) const {
    for (auto& p : polys)
        if (polygonContains(p, pt))
            return true;
    return false;
}

void SDFField::build(const BoundingBox2d& roi, const std::vector<Obstacle>& obs, double cs, double buf) {
    CacheKey key{roi, obs, cs, buf};

    // Try to retrieve from cache
    auto it = cache_map_.find(key);
    if (it != cache_map_.end()) {
        grid_ = it->second;
        cs_ = cs;
        roi_ = roi;
        // Recreate buffered polygons
        std::vector<Polygon2d> buffered;
        for (auto& o : obs) {
            if (o.geometry.outer.empty()) continue;
            buffered.push_back(buf > 0 ? bufferPolygon(o.geometry, buf) : o.geometry);
        }
        if (buffered.size() > 1) {
            std::vector<std::vector<std::array<double,2>>> subs;
            for (auto& p : buffered) if (!p.outer.empty()) subs.push_back(toArray(p.outer));
            auto sol = ClipperUtil::UnionPaths(subs, ClipperLib::pftNonZero);
            if (!sol.empty()) {
                buffered_.clear();
                for (auto& path : sol) {
                    Polygon2d pg;
                    pg.outer = toArray(path);
                    buffered_.push_back(pg);
                }
            } else {
                buffered_ = buffered;
            }
        } else {
            buffered_ = buffered;
        }
        // Calculate rows and cols based on ROI and cell size
        double mg = cs_;
        BoundingBox2d ext;
        ext.min_pt = roi_.min_pt - Vec2d(mg, mg);
        ext.max_pt = roi_.max_pt + Vec2d(mg, mg);
        roi_ = ext;
        cols_ = std::max(2, (int)std::ceil(roi_.width() / cs_));
        rows_ = std::max(2, (int)std::ceil(roi_.height() / cs_));
        return;
    }

    // Not in cache, build normally and store in cache
    buildInternal(roi, obs, cs, buf);

    // Store in cache
    cache_map_[key] = grid_;
}

void SDFField::buildFromPolygons(const BoundingBox2d& roi, const std::vector<Polygon2d>& polys, double cs) {
    cs_ = cs;
    roi_ = roi;
    buffered_ = polys;
    rebuildGrid(polys);
}

void SDFField::rebuildGrid(const std::vector<Polygon2d>& ops) {
    double mg = cs_;
    BoundingBox2d ext;
    ext.min_pt = roi_.min_pt - Vec2d(mg, mg);
    ext.max_pt = roi_.max_pt + Vec2d(mg, mg);
    roi_ = ext;
    cols_ = std::max(2, (int)std::ceil(roi_.width() / cs_));
    rows_ = std::max(2, (int)std::ceil(roi_.height() / cs_));
    grid_.assign(rows_ * cols_, 0.0);
    for (int r = 0; r < rows_; ++r) {
        for (int c = 0; c < cols_; ++c) {
            Vec2d pt = cellToWorld(r, c);
            double d = distToPolygons(pt, ops);
            bool ins = insideAny(pt, ops);
            grid_[idx(r, c)] = ins ? -d : d;
        }
    }
}

void SDFField::updateRegion(const BoundingBox2d& dirty, const std::vector<Obstacle>& obs, double buf) {
    std::vector<Polygon2d> lp;
    for (auto& o : obs)
        lp.push_back(buf > 0 ? bufferPolygon(o.geometry, buf) : o.geometry);
    std::pair<int, int> _mc = worldToCell(dirty.min_pt);
    int rm = _mc.first;
    int cm = _mc.second;
    std::pair<int, int> _xc = worldToCell(dirty.max_pt);
    int rx = _xc.first;
    int cx = _xc.second;
    rm = std::max(0, rm - 1);
    cm = std::max(0, cm - 1);
    rx = std::min(rows_ - 1, rx + 1);
    cx = std::min(cols_ - 1, cx + 1);
    for (int r = rm; r <= rx; ++r)
        for (int c = cm; c <= cx; ++c) {
            Vec2d pt = cellToWorld(r, c);
            double d = distToPolygons(pt, lp);
            bool ins = insideAny(pt, lp);
            grid_[idx(r, c)] = ins ? -d : d;
        }
}

double SDFField::rawAt(int r, int c) const {
    r = std::max(0, std::min(rows_ - 1, r));
    c = std::max(0, std::min(cols_ - 1, c));
    return grid_[idx(r, c)];
}

std::pair<int, int> SDFField::worldToCell(const Vec2d& p) const {
    return {(int)((p[1] - roi_.min_pt.y()) / cs_), (int)((p[0] - roi_.min_pt.x()) / cs_)};
}

Vec2d SDFField::cellToWorld(int r, int c) const {
    return Vec2d(roi_.min_pt[0] + (c + 0.5) * cs_, roi_.min_pt[1] + (r + 0.5) * cs_);
}

std::pair<double, Vec2d> SDFField::queryWithGrad(const Vec2d& pt) const {
    if (grid_.empty()) return {1e18, Vec2d(0, 0)};
    double fx = (pt[0] - roi_.min_pt.x()) / cs_ - 0.5, fy = (pt[1] - roi_.min_pt.y()) / cs_ - 0.5;
    int x0 = (int)std::floor(fx), y0 = (int)std::floor(fy), x1 = x0 + 1, y1 = y0 + 1;
    double tx = fx - x0, ty = fy - y0;
    x0 = std::max(0, std::min(cols_ - 1, x0));
    x1 = std::max(0, std::min(cols_ - 1, x1));
    y0 = std::max(0, std::min(rows_ - 1, y0));
    y1 = std::max(0, std::min(rows_ - 1, y1));
    double v00 = rawAt(y0, x0), v10 = rawAt(y1, x0), v01 = rawAt(y0, x1), v11 = rawAt(y1, x1);
    double val = (1 - tx) * (1 - ty) * v00 + tx * (1 - ty) * v01 + (1 - tx) * ty * v10 + tx * ty * v11;
    double gx = ((1 - ty) * (v01 - v00) + ty * (v11 - v10)) / cs_;
    double gy = ((1 - tx) * (v10 - v00) + tx * (v11 - v01)) / cs_;
    return {val, Vec2d(gx, gy)};
}

bool SDFField::isSafe(const Vec2d& pt, double cl) const {
    std::pair<double,Vec2d> _q = queryWithGrad(pt);
    return _q.first >= cl;
}

double SDFField::obstaclePenalty(const Vec2d& pt, double cl) const {
    std::pair<double,Vec2d> _q2 = queryWithGrad(pt);
    double s = _q2.first - cl;
    return s >= 0 ? 0.0 : s * s;
}

Vec2d SDFField::obstaclePenaltyGrad(const Vec2d& pt, double cl) const {
    std::pair<double, Vec2d> _q3 = queryWithGrad(pt);
    double d = _q3.first;
    Vec2d gd = _q3.second;
    double s = d - cl;
    return s >= 0 ? Vec2d(0, 0) : 2.0 * s * gd;
}

// ─── Static cache implementation ──────────────────────────────────────────────────────

// Define static members
std::unordered_map<SDFField::CacheKey, std::vector<double>, SDFField::CacheKeyHash> SDFField::cache_map_;

bool SDFField::CacheKey::operator==(const CacheKey& other) const {
    // Compare ROI
    if (roi.min_pt.x() != other.roi.min_pt.x() || roi.min_pt.y() != other.roi.min_pt.y() ||
        roi.max_pt.x() != other.roi.max_pt.x() || roi.max_pt.y() != other.roi.max_pt.y()) {
        return false;
    }

    // Compare cell size and buffer
    if (cell_size != other.cell_size || buffer != other.buffer) {
        return false;
    }

    // Compare obstacles - simplified check based on size and positions for efficiency
    if (obstacles.size() != other.obstacles.size()) {
        return false;
    }

    // For each obstacle, check basic properties
    for (size_t i = 0; i < obstacles.size(); ++i) {
        const auto& obs1 = obstacles[i];
        const auto& obs2 = other.obstacles[i];

        if (obs1.geometry.outer.size() != obs2.geometry.outer.size()) {
            return false;
        }

        // Check a few key points to determine similarity
        if (!obs1.geometry.outer.empty() && !obs2.geometry.outer.empty()) {
            if (std::abs(obs1.geometry.outer[0].x() - obs2.geometry.outer[0].x()) > 1e-6 ||
                std::abs(obs1.geometry.outer[0].y() - obs2.geometry.outer[0].y()) > 1e-6) {
                return false;
            }
        }
    }

    return true;
}

std::size_t SDFField::CacheKeyHash::operator()(const SDFField::CacheKey& k) const {
    std::size_t h1 = std::hash<double>{}(k.roi.min_pt.x());
    std::size_t h2 = std::hash<double>{}(k.roi.min_pt.y());
    std::size_t h3 = std::hash<double>{}(k.roi.max_pt.x());
    std::size_t h4 = std::hash<double>{}(k.roi.max_pt.y());
    std::size_t h5 = std::hash<double>{}(k.cell_size);
    std::size_t h6 = std::hash<double>{}(k.buffer);
    std::size_t h7 = std::hash<size_t>{}(k.obstacles.size());

    // Hash some key points from obstacles if they exist
    std::size_t h8 = 0;
    if (!k.obstacles.empty() && !k.obstacles[0].geometry.outer.empty()) {
        h8 = std::hash<double>{}(k.obstacles[0].geometry.outer[0].x() +
                                 k.obstacles[0].geometry.outer[0].y());
    }

    // Combine hashes
    return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4) ^ (h6 << 5) ^ (h7 << 6) ^ (h8 << 7);
}

void SDFField::clearCache() {
    cache_map_.clear();
}

size_t SDFField::getCacheSize() {
    return cache_map_.size();
}

void SDFField::buildInternal(const BoundingBox2d& roi, const std::vector<Obstacle>& obs, double cs, double buf) {
    cs_ = cs;
    roi_ = roi;
    std::vector<Polygon2d> buffered;
    for (auto& o : obs) {
        if (o.geometry.outer.empty()) continue;
        buffered.push_back(buf > 0 ? bufferPolygon(o.geometry, buf) : o.geometry);
    }
    if (buffered.size() > 1) {
        std::vector<std::vector<std::array<double,2>>> subs;
        for (auto& p : buffered) if (!p.outer.empty()) subs.push_back(toArray(p.outer));
        auto sol = ClipperUtil::UnionPaths(subs, ClipperLib::pftNonZero);
        if (!sol.empty()) {
            buffered_.clear();
            for (auto& path : sol) {
                Polygon2d pg;
                pg.outer = toArray(path);
                buffered_.push_back(pg);
            }
        } else {
            buffered_ = buffered;
        }
    } else {
        buffered_ = buffered;
    }
    rebuildGrid(buffered_);
}

// ─── End of file additions ────────────────────────────────────────────────────────────
