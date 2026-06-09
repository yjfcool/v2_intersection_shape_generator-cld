#pragma once
#include "bezier.h"
struct SDFField;
double localCurvature(const Vec2d& a, const Vec2d& b, const Vec2d& c);
Vec2d circumcenter(const Vec2d& a, const Vec2d& b, const Vec2d& c);
std::vector<Vec2d> elasticBandSmooth(const std::vector<Vec2d>&, const SDFField&,
                                     const Polygon2d&, double kappa_max = 0.25, double move_step = 0.03,
                                     int max_iter = 50, double min_sdf = 0.1);
std::vector<Vec2d> elasticBandSmoothAdaptive(const std::vector<Vec2d>&, const SDFField&,
                                     const Polygon2d&, double kappa_max = 0.25, double move_step = 0.03,
                                     int max_iter = 50, double min_sdf = 0.1);
std::vector<Vec2d> downsamplePoints(const std::vector<Vec2d>&, int target_count);
BezierCurve rebuildFromSmoothedPts(const std::vector<Vec2d>&, const Vec2d&, const Vec2d&);
std::vector<Vec2d> midlineSampleByArcLength(const BezierCurve&, const BezierCurve&, int n = 20);
double signedDistToLine(const Vec2d& pt, const Vec2d& p0, const Vec2d& p1);
bool segmentsIntersect(const Vec2d& a, const Vec2d& b,
                       const Vec2d& c, const Vec2d& d, Vec2d* out = nullptr);
double minSDFAlongCurve(const BezierCurve&, const SDFField&, int sps = 20);
double distToAllEndpoints(const Vec2d&, const BezierCurve&, const BezierCurve&);

// src/constraints/intersection_check.h
std::vector<Vec2d> curveCrossings(const BezierCurve&, const BezierCurve&, double tol = 0.01);
bool bboxOverlap(const BezierCurve&, const BezierCurve&);
bool curvesIntersectBusiness(const BezierCurve&, const BezierCurve&, double ep = 0.01);
