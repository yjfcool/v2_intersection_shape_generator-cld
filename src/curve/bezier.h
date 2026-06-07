#pragma once
#include "types.h"

struct SDFField;

// ── BezierSegment & BezierCurve structs are defined in types.h ──

// ── Free functions ───────────────────────────────────────────
BezierSegment makeCubicG1(const Vec2d& p0, const Vec2d& t0, const Vec2d& p1, const Vec2d& t1, double alpha = 0.33);

BezierCurve makeCurveFromKnots(const std::vector<Vec2d>& pts, const std::vector<Vec2d>& tans, double alpha = 0.33);

BezierCurve fitBezierWithEndTangents(const std::vector<Vec2d>& pts, const Vec2d& start_tan, const Vec2d& end_tan);

AdaptiveRefineResult adaptiveRefine(const BezierCurve& c, const SDFField& sdf, double kappa_max = 0.25, double pen_tol = 0.05);

VecXd curveToParams(const BezierCurve& c);

BezierCurve curveFromParams(const VecXd& params, const BezierCurve& proto);

// all control points (incl. join pts) are optimisation variables
VecXd curveToParamsFull(const BezierCurve& c);

BezierCurve curveFromParamsFull(const VecXd& params, const BezierCurve& proto);
