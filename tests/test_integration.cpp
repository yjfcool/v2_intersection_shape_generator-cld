#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "intersection_shape_generator.h"
#include "constraints/fence_check.h"
#include "optimizer/sdf_field.h"

#include "io/iodata_shapefile.h"

using Catch::Matchers::WithinAbs;

using namespace isg;

// ──────────────────────────────────────────────────────────────
//  Minimal 2-direction intersection builder
//  East-entry → West-exit (straight)
// ──────────────────────────────────────────────────────────────
static IntersectionInput makeTwoDirectionInput(bool add_obstacle = false) {
    IntersectionInput inp;

    // Lanes with geometry:
    //   Entry lane (LE): outside→junction, last point = junction entry
    //   Exit  lane (LR): junction→outside, first point = junction exit
    Lane le; le.id="LE"; le.width=3.5;
    le.geometry.points = {Vec2d(-3,0), Vec2d(0,0)};   // entry at (0,0), tangent=(1,0)
    Lane lr; lr.id="LR"; lr.width=3.5;
    lr.geometry.points = {Vec2d(10,0), Vec2d(13,0)};  // exit at (10,0), tangent=(1,0)
    inp.lanes = {le, lr};

    // Lane edges (stubs — geometry not used by generator directly)
    LaneEdge be0; be0.id="BE0"; be0.geometry.points={{0,-1.75},{0,0}};
    LaneEdge be1; be1.id="BE1"; be1.geometry.points={{0,0},{0,1.75}};
    LaneEdge be2; be2.id="BE2"; be2.geometry.points={{0,1.75},{0,3.5}};
    LaneEdge br0; br0.id="BR0";
    LaneEdge br1; br1.id="BR1";
    LaneEdge br2; br2.id="BR2";
    inp.lane_edges = {be0,be1,be2,br0,br1,br2};

    // Entry group (east side, going west → x direction = (1,0))
    LaneGroup entry;
    entry.id="EG"; entry.role=GroupRole::Entry;
    entry.lanes={"LE"};
    entry.boundaries={"BE0","BE1"};
    //entry.direction=Vec2d(1,0);
    //entry.ref_point=Vec2d(0,0);

    // Exit group (west side, continuing west)
    LaneGroup exitg;
    exitg.id="XG"; exitg.role=GroupRole::Exit;
    exitg.lanes={"LR"};
    exitg.boundaries={"BR0","BR1"};
    //exitg.direction=Vec2d(1,0);
    //exitg.ref_point=Vec2d(10,0);

    inp.lane_groups = {entry, exitg};

    // Connectivity: straight through
    Connectivity c;
    c.id="CONN1"; c.entry_lane_id="LE"; c.exit_lane_id="LR";
    c.turn_type=ConnTurnType::Straight;
    inp.connectivities = {c};

    // Coarse area
    Polygon2d fence;
    fence.outer = {{-1,-5},{11,-5},{11,5},{-1,5}};
    inp.area.geometry = fence;

    // Boundary (road edges)
    Boundary top_edge;
    top_edge.id="TOP"; top_edge.type=Boundary::Type::RoadEdge;
    top_edge.geometry.points={{-1,5},{11,5}};

    Boundary bot_edge;
    bot_edge.id="BOT"; bot_edge.type=Boundary::Type::RoadEdge;
    bot_edge.geometry.points={{-1,-5},{11,-5}};

    inp.boundaries = {top_edge, bot_edge};

    if (add_obstacle) {
        Obstacle obs;
        obs.geometry.outer = {{4,-1},{6,-1},{6,1},{4,1}};
        inp.obstacles = {obs};
    }

    return inp;
}

// ──────────────────────────────────────────────────────────────
TEST_CASE("Integration: simple straight, no obstacles, OK status", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    bool ok = gen.generate(inp, out);

    save(inp, "./straight-noobstacles/input/", "");
    save(out, "./straight-noobstacles/output/", "");

    REQUIRE(ok);
    REQUIRE(out.connectivity_curves.size() == 1);

    auto& cc = out.connectivity_curves[0];
    REQUIRE(cc.curve);
    REQUIRE((cc.status == CurveStatus::OK ||
             cc.status == CurveStatus::WarnA2));

    // Endpoints near entry/exit reference points
    REQUIRE_THAT((cc.curve->startPt()-Vec2d(0,0)).norm(), WithinAbs(0.0, 0.5));
    REQUIRE_THAT((cc.curve->endPt()  -Vec2d(10,0)).norm(),WithinAbs(0.0, 0.5));
}

TEST_CASE("Integration: start/end tangent G1 continuity", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    gen.generate(inp, out);

    save(inp, "./G1/input/", "");
    save(out, "./G1/output/", "");

    REQUIRE(!out.connectivity_curves.empty());

    auto& cc = out.connectivity_curves[0];
    REQUIRE(cc.curve);

    // Tangent at start should point roughly east (+x)
    Vec2d start_tan = cc.curve->startTan();
    REQUIRE(start_tan[0] > 0.7);

    // Tangent at end should also point roughly east
    Vec2d end_tan = cc.curve->endTan();
    REQUIRE(end_tan[0] > 0.7);
}

TEST_CASE("Integration: curve stays inside fence", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    gen.generate(inp, out);

    save(inp, "./fence/input/", "");
    save(out, "./fence/output/", "");

    REQUIRE(!out.connectivity_curves.empty());
    auto& cc = out.connectivity_curves[0];
    REQUIRE(cc.curve);

    auto& fence = inp.area.geometry;
    for (auto& pt : cc.curve->sampleByArcLength(40)) {
        // Allow small numerical overshoot
        double dist = pointToPolygonDist(pt, fence);
        bool inside = polygonContains(fence, pt);
        if (!inside) {
            REQUIRE(dist < 0.3);  // no more than 30cm overshoot
        }
    }
}

TEST_CASE("Integration: obstacle avoidance, no penetration", "[integration]") {
    auto inp = makeTwoDirectionInput(true);   // add central obstacle

    IntersectionShapeGenerator::Config cfg;
    cfg.lbfgs.max_iter = 200;
    IntersectionShapeGenerator gen(cfg);
    IntersectionOutput out;
    bool ok = gen.generate(inp, out);

    save(inp, "./obstacle-nopenetration/input/", "");
    save(out, "./obstacle-nopenetration/output/", "");

    REQUIRE(ok);
    REQUIRE(!out.connectivity_curves.empty());

    auto& cc = out.connectivity_curves[0];
    if (cc.status == CurveStatus::Infeasible) return;  // acceptable

    REQUIRE(cc.curve);

    // Build SDF to check penetration
    SDFField check_sdf;
    check_sdf.build(inp.area.geometry.bbox(), inp.obstacles, 0.2, 0.0);

    for (auto& pt : cc.curve->sampleByArcLength(50)) {
        std::pair<double,Vec2d> _q = check_sdf.queryWithGrad(pt);
        // Allow tiny numerical penetration
        REQUIRE(_q.first > -0.2);
    }
}

TEST_CASE("Integration: topology validation catches mismatched boundaries", "[integration]") {
    auto inp = makeTwoDirectionInput(false);
    // Break boundary count: remove one boundary from entry group
    inp.lane_groups[0].boundaries.pop_back();

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    bool ok = gen.generate(inp, out);

    save(inp, "./topology-validation/input/", "");
    save(out, "./topology-validation/output/", "");

    REQUIRE_FALSE(ok);
    REQUIRE(!gen.lastReport().is_valid());
    REQUIRE(!gen.lastReport().errors.empty());
}

TEST_CASE("Integration: fine area is non-empty polygon", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    gen.generate(inp, out);

    save(inp, "./fine-area/input/", "");
    save(out, "./fine-area/output/", "");

    REQUIRE_FALSE(out.area.geometry.outer.empty());
    REQUIRE(out.area.geometry.outer.size() >= 3);
}

TEST_CASE("Integration: perf stats are populated", "[integration]") {
    auto inp = makeTwoDirectionInput(false);

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    gen.generate(inp, out);

    save(inp, "./perf-stats/input/", "");
    save(out, "./perf-stats/output/", "");

    REQUIRE(out.perf.sdf_build_ms >= 0.0);
    REQUIRE(out.perf.optimize_ms  >= 0.0);
    REQUIRE(out.perf.edge_gen_ms  >= 0.0);
    REQUIRE(out.perf.area_gen_ms  >= 0.0);
}

// ──────────────────────────────────────────────────────────────
//  4-way intersection: 4 entry + 4 exit groups, 8 connectivities
// ──────────────────────────────────────────────────────────────
static IntersectionInput makeFourWayInput() {
    IntersectionInput inp;

    // Coarse area: 20×20m square centred at origin
    Polygon2d fence;
    fence.outer = {{-10,-10},{10,-10},{10,10},{-10,10}};
    inp.area.geometry = fence;

    // 4 directions: North(N), South(S), East(E), West(W)
    // Each direction: 1 entry lane + 1 exit lane
    struct Dir { std::string name; Vec2d ref_in; Vec2d dir_in;
                                   Vec2d ref_out; Vec2d dir_out; };
    std::vector<Dir> dirs = {
        {"N", {0,-10}, {0,1},  {0,10},  {0,1}},
        {"S", {0, 10}, {0,-1}, {0,-10}, {0,-1}},
        {"E", {-10,0}, {1,0},  {10,0},  {1,0}},
        {"W", {10, 0}, {-1,0}, {-10,0}, {-1,0}},
    };

    for (auto& d : dirs) {
        // Entry lane: outside → ref_in (last pt = junction edge)
        // Exit  lane: ref_out → outside (first pt = junction edge)
        Lane le; le.id="LE_"+d.name; le.width=3.5;
        le.geometry.points={d.ref_in - d.dir_in*5.0, d.ref_in};
        Lane lr; lr.id="LR_"+d.name; lr.width=3.5;
        lr.geometry.points={d.ref_out, d.ref_out + d.dir_out*5.0};
        inp.lanes.push_back(le);
        inp.lanes.push_back(lr);

        LaneGroup entry;
        entry.id="EG_"+d.name; entry.role=GroupRole::Entry;
        entry.lanes={"LE_"+d.name};
        entry.boundaries={"B0_"+d.name,"B1_"+d.name};
        //entry.direction=d.dir_in;
        //entry.ref_point=d.ref_in;

        LaneGroup exitg;
        exitg.id="XG_"+d.name; exitg.role=GroupRole::Exit;
        exitg.lanes={"LR_"+d.name};
        exitg.boundaries={"C0_"+d.name,"C1_"+d.name};
        //exitg.direction=d.dir_out;
        //exitg.ref_point=d.ref_out;

        inp.lane_groups.push_back(entry);
        inp.lane_groups.push_back(exitg);
    }

    // Straight-through connectivities for each direction
    std::vector<std::pair<std::string,std::string>> straights = {
        {"N","S"},{"S","N"},{"E","W"},{"W","E"}
    };
    int idx = 0;
    for (auto& kv : straights) {
        auto& from = kv.first; auto& to = kv.second;
        Connectivity c;
        c.id = "C"+std::to_string(idx++);
        c.entry_lane_id = "LE_"+from;
        c.exit_lane_id  = "LR_"+to;
        c.turn_type     = ConnTurnType::Straight;
        inp.connectivities.push_back(c);
    }

    return inp;
}

TEST_CASE("Integration: 4-way intersection generates all curves", "[integration]") {
    auto inp = makeFourWayInput();

    IntersectionShapeGenerator gen;
    IntersectionOutput out;
    bool ok = gen.generate(inp, out);

    save(inp, "./4-way/input/", "");
    save(out, "./4-way/output/", "");

    REQUIRE(ok);

    // All 4 connectivities should have curves
    REQUIRE(out.connectivity_curves.size() == 4);
    for (auto& cc : out.connectivity_curves) {
        bool feasible = (cc.status != CurveStatus::Infeasible);
        if (feasible) REQUIRE(cc.curve);
    }
}
