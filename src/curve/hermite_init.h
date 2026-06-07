#pragma once
#include "bezier.h"
struct SDFField;
std::vector<Vec2d> sdfMaxClearancePath(
    const SDFField&, const Polygon2d&, const Vec2d&, const Vec2d&, double alpha = 2.0);

// sibling_polys: pre-sampled polylines of already-generated curves in this cluster
// Used by Level-1 to avoid choosing the same bypass side as opposing arches.
BezierCurve buildInitialCurve(
    const Vec2d&, const Vec2d&, const Vec2d&, const Vec2d&,
    const SDFField&, const Polygon2d&, const std::vector<std::vector<Vec2d>>& sibling_polys = {});

BezierCurve buildTwoSegmentUTurn(
    const Vec2d&, const Vec2d&, const Vec2d&, const Vec2d&, const SDFField&, const Polygon2d&);
