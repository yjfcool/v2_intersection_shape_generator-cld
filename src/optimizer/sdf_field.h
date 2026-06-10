#pragma once
#include "types.h"
#include <unordered_map>

namespace isg {

class SDFField {
public:
    void build(const BoundingBox2d&, const std::vector<Obstacle>&, double cs = 0.2, double buf = 0.0);
    void buildFromPolygons(const BoundingBox2d&, const std::vector<Polygon2d>&, double cs = 0.2);
    std::pair<double, Vec2d> queryWithGrad(const Vec2d&) const;
    bool isSafe(const Vec2d&, double cl = 0.1) const;
    double obstaclePenalty(const Vec2d&, double cl = 0.1) const;
    Vec2d obstaclePenaltyGrad(const Vec2d&, double cl = 0.1) const;
    void updateRegion(const BoundingBox2d&, const std::vector<Obstacle>&, double buf = 0.0);
    double cellSize() const { return cs_; }
    int rows() const { return rows_; }
    int cols() const { return cols_; }
    std::pair<int, int> worldToCell(const Vec2d&) const;
    Vec2d cellToWorld(int r, int c) const;
    bool valid() const { return !grid_.empty(); }

    // Static cache for SDF field reuse
    static void clearCache();
    static size_t getCacheSize();

private:
    void rebuildGrid(const std::vector<Polygon2d>&);
    double distToPolygons(const Vec2d&, const std::vector<Polygon2d>&) const;
    bool insideAny(const Vec2d&, const std::vector<Polygon2d>&) const;
    static Polygon2d bufferPolygon(const Polygon2d&, double r);
    static Polygon2d unionPolygons(const std::vector<Polygon2d>&);
    inline int idx(int r, int c) const { return r * cols_ + c; }
    double rawAt(int r, int c) const;

    // Cache-related methods and data
    struct CacheKey {
        BoundingBox2d roi;
        std::vector<Obstacle> obstacles;
        double cell_size;
        double buffer;

        bool operator==(const CacheKey& other) const;
    };

    struct CacheKeyHash {
        std::size_t operator()(const SDFField::CacheKey& k) const;
    };

    void buildInternal(const BoundingBox2d& roi, const std::vector<Obstacle>& obs, double cs, double buf);

private:
    double cs_ = 0.2;
    BoundingBox2d roi_;
    int rows_ = 0, cols_ = 0;
    std::vector<double> grid_;
    std::vector<Polygon2d> buffered_;
    static std::unordered_map<CacheKey, std::vector<double>, CacheKeyHash> cache_map_; // Static cache storage
};

}