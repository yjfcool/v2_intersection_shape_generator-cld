#include <catch2/catch_test_macros.hpp>

#include "constraints/cluster_order.h"
#include "curve/curve_utils.h"
#include "intersection_shape_generator.h"
#include "io/iodata_json.h"

#include <string>
#include <unordered_map>
#include <vector>

using namespace isg;

TEST_CASE("intersection_input has no avoidable same-cluster interior crossings", "[regression][cluster]") {
    const std::string path = std::string(PROJECT_ROOT_DIR) + "/datas/intersection_input.json";
    IntersectionInput input = IntersectionIO::loadFromFile(path);

    IntersectionShapeGenerator gen;
    IntersectionOutput output;
    REQUIRE(gen.generate(input, output));

    std::unordered_map<ConnId, const ConnectivityCurve*> curves;
    for (const auto& cc : output.connectivity_curves)
        curves[cc.id] = &cc;

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

    INFO("avoidable same-cluster crossings: " << [&] {
        std::string s;
        for (const auto& p : bad_pairs) {
            if (!s.empty()) s += ", ";
            s += p;
        }
        return s;
    }());
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
    INFO("named example crossings: " << [&] {
        std::string s;
        for (const auto& p : crossed_examples) {
            if (!s.empty()) s += ", ";
            s += p;
        }
        return s;
    }());
    CHECK(crossed_examples.empty());
}
