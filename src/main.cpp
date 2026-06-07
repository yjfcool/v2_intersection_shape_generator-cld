// ============================================================
// main
// ============================================================
#include "io/iodata_json.h"
#include "io/iodata_shapefile.h"
#include "intersection_shape_generator.h"

int main(int argc, char* argv[]) {
    std::string fpth = std::string(PROJECT_ROOT_DIR) + "/intersection_input.json";
    std::vector<IntersectionInput> inputs = {IntersectionIO::loadFromFile(fpth)};
    for (auto& inp : inputs) {
        IntersectionShapeGenerator gen;
        IntersectionOutput out;
        if (!gen.generate(inp, out))
            continue;

        save(inp, "./actual/input/", "");
        save(out, "./actual/output/", "");
    }
    return 0;
}