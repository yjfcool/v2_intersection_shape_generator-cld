#include <catch2/catch_test_macros.hpp>

#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"

#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

static std::unordered_map<ConnId, const ConnectivityCurve*> curveMap(const IntersectionOutput& output) {
    std::unordered_map<ConnId, const ConnectivityCurve*> curves;
    for (const auto& cc : output.connectivity_curves)
        curves[cc.id] = &cc;
    return curves;
}

static std::vector<std::string> avoidableSameClusterCrossings(
    const IntersectionInput& input, const IntersectionOutput& output) {
    auto curves = curveMap(output);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups);

    std::vector<std::string> bad_pairs;
    for (const auto& pair : solver.pairs()) {
        if (pair.exempt == CrossExemption::StructuralCross)
            continue;
        auto ia = curves.find(pair.id_a);
        auto ib = curves.find(pair.id_b);
        if (ia == curves.end() || ib == curves.end())
            continue;
        if (!ia->second->curve || !ib->second->curve)
            continue;
        if (curvesIntersectBusiness(*ia->second->curve, *ib->second->curve, 1.5))
            bad_pairs.push_back(pair.id_a + "-" + pair.id_b);
    }
    return bad_pairs;
}

static std::string joinPairs(const std::vector<std::string>& pairs) {
    std::string s;
    for (const auto& p : pairs) {
        if (!s.empty()) s += ", ";
        s += p;
    }
    return s;
}

TEST_CASE("intersection_input has no avoidable same-cluster interior crossings", "[regression][cluster]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_input.json";
    IntersectionInput input = IntersectionIO::loadFromFile(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    std::vector<std::string> bad_pairs = avoidableSameClusterCrossings(input, output);

    INFO("avoidable same-cluster crossings: " << joinPairs(bad_pairs));
    CHECK(bad_pairs.empty());

    const std::vector<std::pair<ConnId, ConnId>> named_examples = {
        {"18", "32"},
        {"30", "12"},
        {"2", "3"},
        {"19", "21"},
        {"5", "22"},
        {"23", "19"},
        {"23", "21"},
        {"24", "20"},
        {"24", "22"},
    };
    std::vector<std::string> crossed_examples;
    for (const auto& ids : named_examples) {
        auto ia = curves.find(ids.first);
        auto ib = curves.find(ids.second);
        REQUIRE(ia != curves.end());
        REQUIRE(ib != curves.end());
        REQUIRE(ia->second->curve);
        REQUIRE(ib->second->curve);
        if (curvesIntersectBusiness(*ia->second->curve, *ib->second->curve, 1.5))
            crossed_examples.push_back(ids.first + "-" + ids.second);
    }
    INFO("named example crossings: " << joinPairs(crossed_examples));
    CHECK(crossed_examples.empty());
}

TEST_CASE("100000643 obstacle reroute preserves same-cluster U-turn topology", "[regression][cluster][obstacle]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json";
    IntersectionInput input = IntersectionIO::loadFromFile(path);

    ClusterOrderSolver solver;
    solver.build(input.connectivities, input.lanes, input.lane_groups);
    CHECK(solver.exemptionOf("71", "75") == CrossExemption::None);
    CHECK(solver.exemptionOf("71", "77") == CrossExemption::None);
    CHECK(solver.exemptionOf("75", "77") == CrossExemption::None);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    REQUIRE(curves.count("71") == 1);
    REQUIRE(curves.count("75") == 1);
    REQUIRE(curves.count("77") == 1);
    REQUIRE(curves["75"]->curve);
    REQUIRE(curves["77"]->curve);

    INFO("71/75 must not have a same-entry non-endpoint crossing");
    CHECK_FALSE(curvesIntersectBusiness(*curves["71"]->curve, *curves["75"]->curve, 1.5));

    INFO("71/77 must not have a same-entry non-endpoint crossing");
    CHECK_FALSE(curvesIntersectBusiness(*curves["71"]->curve, *curves["77"]->curve, 1.5));

    INFO("75/77 must only meet at the shared entry endpoint");
    CHECK_FALSE(curvesIntersectBusiness(*curves["75"]->curve, *curves["77"]->curve, 1.5));
}

TEST_CASE("100000643 group-unified direction keeps large U-turn bounded", "[regression][cluster][obstacle][direction]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/100000643.json";
    IntersectionInput input = IntersectionIO::loadFromFile(path);

    IntersectionShapeGenerator::Config cfg;
    cfg.connectivity_direction.mode = ConnectivityDirectionMode::GroupUnified;
    cfg.connectivity_direction.group_similarity_angle_deg = 5.0;
    IntersectionShapeGenerator gen(cfg);
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    auto curves = curveMap(output);
    REQUIRE(curves.count("71") == 1);
    REQUIRE(curves.count("75") == 1);
    REQUIRE(curves.count("77") == 1);
    REQUIRE(curves["71"]->curve);
    REQUIRE(curves["75"]->curve);
    REQUIRE(curves["77"]->curve);

    auto pts77 = curves["77"]->curve->sampleByArcLength(80);
    double min_x = 1e18;
    for (const auto& pt : pts77)
        min_x = std::min(min_x, pt.x());

    CHECK(curves["77"]->curve->arcLength() < 80.0);
    CHECK(min_x > 0.0);
    CHECK_FALSE(curvesIntersectBusiness(*curves["71"]->curve, *curves["77"]->curve, 1.5));
    CHECK_FALSE(curvesIntersectBusiness(*curves["75"]->curve, *curves["77"]->curve, 1.5));
}
