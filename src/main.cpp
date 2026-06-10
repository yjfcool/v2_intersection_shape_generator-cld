// ============================================================
// main
// ============================================================
#include "io/iodata_json.h"
#include "io/iodata_shapefile.h"
#include "intersection_shape_generator.h"
#include <chrono>
#include <unistd.h>

int main(int argc, char* argv[]) {
    std::vector<std::string> files;
    if (argc > 1) {
        files.emplace_back(argv[1]);
    } else {
        files = {
            std::string(PROJECT_ROOT_DIR) + "/datas/" + "intersection_input.json",
            std::string(PROJECT_ROOT_DIR) + "/datas/" + "intersection_ds.json",
            std::string(PROJECT_ROOT_DIR) + "/datas/" + "intersection_jd.json",
            std::string(PROJECT_ROOT_DIR) + "/datas/" + "intersection_cross.json",
            std::string(PROJECT_ROOT_DIR) + "/datas/" + "100000598.json",
            std::string(PROJECT_ROOT_DIR) + "/datas/" + "100000643.json",
            std::string(PROJECT_ROOT_DIR) + "/datas/" + "100000699.json",
        };
    }

    for (const auto& fpth : files) {
        ghc::filesystem::path fspath(fpth);
        if (!ghc::filesystem::exists(fspath)) continue;
        //if (access(fpth.c_str(), F_OK) == -1) continue;
        std::string fname = fspath.stem().string();
        std::cout << "\n===== Loading input from: " << fname << " =====" << std::endl;
        std::vector<isg::IntersectionInput> inputs = {isg::IntersectionIO::loadFromFile(fpth)};
        for (auto& inp : inputs) {
            std::cout << "Processing input with " << inp.connectivities.size() << " connectivities and "
                      << inp.lanes.size() << " lanes" << std::endl;

            isg::IntersectionShapeGenerator gen;
            isg::IntersectionOutput out;

            auto start_time = std::chrono::high_resolution_clock::now();
            bool success = gen.generate(inp, out);
            auto end_time = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << "Generation took: " << duration.count() << " ms" << std::endl;
            std::cout << "Performance stats:" << std::endl;
            std::cout << "  SDF build: " << out.perf.sdf_build_ms << " ms" << std::endl;
            std::cout << "  Optimization: " << out.perf.optimize_ms << " ms" << std::endl;
            std::cout << "  Edge generation: " << out.perf.edge_gen_ms << " ms" << std::endl;
            std::cout << "  Area generation: " << out.perf.area_gen_ms << " ms" << std::endl;

            if (!success) {
                std::cout << "Generation failed!" << std::endl;
                std::cout << "Validation errors: " << gen.lastReport().errors.size() << std::endl;
                for (const auto& err : gen.lastReport().errors) {
                    std::cout << "  " << err << std::endl;
                }
                continue;
            }
            std::cout << "Successfully generated " << out.connectivity_curves.size() << " connectivity curves" << std::endl;

            save(inp, "./" + fname + "/input/", "");
            save(out, "./" + fname + "/output/", "");
        }
    }

    return 0;
}