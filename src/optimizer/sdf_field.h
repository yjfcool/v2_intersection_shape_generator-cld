#pragma once
#include "types.h"

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

private:
    void rebuildGrid(const std::vector<Polygon2d>&);
    double distToPolygons(const Vec2d&, const std::vector<Polygon2d>&) const;
    bool insideAny(const Vec2d&, const std::vector<Polygon2d>&) const;
    static Polygon2d bufferPolygon(const Polygon2d&, double r);
    static Polygon2d unionPolygons(const std::vector<Polygon2d>&);
    inline int idx(int r, int c) const { return r * cols_ + c; }
    double rawAt(int r, int c) const;

private:
    double cs_ = 0.2;
    BoundingBox2d roi_;
    int rows_ = 0, cols_ = 0;
    std::vector<double> grid_;
    std::vector<Polygon2d> buffered_;
};
