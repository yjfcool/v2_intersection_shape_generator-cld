#include "intersection_shape_generator.h"
#include "utils.h"
#include "generator/connectivity_generator.h"
#include "constraints/infeasibility_detector.h"
#include "generator/edge_line_generator.h"
#include "generator/polygon_builder.h"
#include <chrono>
#include <cmath>
#include <unordered_map>
namespace isg {

IntersectionShapeGenerator::IntersectionShapeGenerator(): cfg_{} {
}

IntersectionShapeGenerator::IntersectionShapeGenerator(const Config& cfg): cfg_(cfg) {
}

static IntersectionInput normalizeConnectivityGroups(IntersectionInput input) {
    std::unordered_map<LaneId, LaneGroupId> entry_gid;
    std::unordered_map<LaneId, LaneGroupId> exit_gid;
    for (auto& g : input.lane_groups) {
        for (auto& lid : g.lanes) {
            if (g.role == GroupRole::Entry)
                entry_gid[lid] = g.id;
            else
                exit_gid[lid] = g.id;
        }
    }
    for (auto& l : input.lanes) {
        auto ie = entry_gid.find(l.id);
        auto ix = exit_gid.find(l.id);
        if (l.groupId.empty()) {
            if (ie != entry_gid.end())
                l.groupId = ie->second;
            else if (ix != exit_gid.end())
                l.groupId = ix->second;
        }
        if (!l.groupId.empty()) {
            if (!entry_gid.count(l.id))
                entry_gid[l.id] = l.groupId;
            if (!exit_gid.count(l.id))
                exit_gid[l.id] = l.groupId;
        }
    }
    for (auto& c : input.connectivities) {
        if (c.enterGroupId.empty()) {
            auto it = entry_gid.find(c.entry_lane_id);
            if (it != entry_gid.end())
                c.enterGroupId = it->second;
        }
        if (c.exitGroupId.empty()) {
            auto it = exit_gid.find(c.exit_lane_id);
            if (it != exit_gid.end())
                c.exitGroupId = it->second;
        }
    }
    return input;
}

inline const char* turnTypeName(ConnTurnType t) {
    switch (t) {
    case ConnTurnType::UTurnLeft: return "UTurnLeft";
    case ConnTurnType::TurnLeft: return "TurnLeft";
    case ConnTurnType::Straight: return "Straight";
    case ConnTurnType::TurnRight: return "TurnRight";
    case ConnTurnType::UTurnRight: return "UTurnRight";
    default: return "Unknown";
    }
}

static ConnTurnType inferTurn(const Connectivity& conn, const IntersectionInput& inp) {
    std::pair<Vec2d, Vec2d> _entry = inp.entryPtDir(conn.entry_lane_id);
    Vec2d p0 = _entry.first;
    Vec2d t0 = _entry.second;
    std::pair<Vec2d, Vec2d> _exit = inp.exitPtDir(conn.exit_lane_id);
    Vec2d p1 = _exit.first;
    Vec2d t1 = _exit.second;
    Vec2d te = (p1 - p0);
    if (te.norm() < 1e-6)
        return ConnTurnType::Straight;
    te.normalize();
    double c = cross2d(t0, te), dot = t0.dot(te);
    if (dot < -0.5)
        return ConnTurnType::UTurnLeft;
    if (c > 0.5)
        return ConnTurnType::TurnLeft;
    if (c < -0.5)
        return ConnTurnType::TurnRight;
    return ConnTurnType::Straight;
}

ValidationReport validateTopology(const IntersectionInput& input) {
    ValidationReport r;
    for (auto& g : input.lane_groups)
        if (g.boundaries.size() > 0 && g.boundaries.size() != g.lanes.size() + 1)
            r.errors.push_back("Group " + g.id + ": boundary count mismatch (expected "
                + std::to_string(g.lanes.size() + 1) + ", got "
                + std::to_string(g.boundaries.size()) + ")");
    for (auto& conn : input.connectivities) {
        ConnTurnType inf = inferTurn(conn, input);
        if (std::abs((int)conn.turn_type - (int)inf) > 1)
            r.warnings.push_back("Connectivity " + conn.id + ": declared=" + turnTypeName(conn.turn_type)
                + " geometry=" + turnTypeName(inf));
    }
    if (!input.area.geometry.outer.empty() && !isSimplePolygon(input.area.geometry))
        r.errors.push_back("Coarse intersection area is self-intersecting");
    for (auto& sl : input.stop_lines)
        if (!sl.associated_group_id.empty() && !input.laneGroupExists(sl.associated_group_id))
            r.warnings.push_back("StopLine " + sl.id + ": group '" + sl.associated_group_id + "' not found");
    return r;
}

bool IntersectionShapeGenerator::generate(const IntersectionInput& input, IntersectionOutput& output) {
    try {
        IntersectionInput norm_input = normalizeConnectivityGroups(input);
        report_ = validateTopology(norm_input);
        if (!report_.is_valid())
            return false;

        auto t0 = std::chrono::steady_clock::now();
        SDFField sdf;
        BoundingBox2d roi;
        if (!norm_input.area.geometry.outer.empty())
            roi = norm_input.area.geometry.bbox();
        else {
            for (auto& l : norm_input.lanes) {
                for (auto& p : l.geometry.points)
                    roi.expand(p);
            }
            roi.min_pt -= Vec2d(20, 20);
            roi.max_pt += Vec2d(20, 20);
        }
        if (!norm_input.obstacles.empty())
            sdf.build(roi, norm_input.obstacles, cfg_.sdf_cell_size, cfg_.obstacle_buffer);
        output.perf.sdf_build_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

        double opt_ms = 0;
        ConnectivityGenerator cgen(cfg_.lbfgs, cfg_.connectivity_direction);
        output.connectivity_curves = cgen.generate(norm_input, sdf, &opt_ms);
        output.perf.optimize_ms = opt_ms;

        auto te = std::chrono::steady_clock::now();
        //EdgeLineGenerator elgen;
        //output.lane_edges = elgen.generate(norm_input, output.connectivity_curves);
        output.perf.edge_gen_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - te).count();

        auto ta = std::chrono::steady_clock::now();
        IntersectionAreaBuilder areabuilder;
        output.area = areabuilder.build(norm_input, output.connectivity_curves, output.lane_edges);
        output.perf.area_gen_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ta).count();

    } catch (...) {
        output.connectivity_curves.clear();
        output.lane_edges.clear();
        output.area = IntersectionArea();
        return false;
    }
    return true;
}

}
