#pragma once
#include "types.h"
namespace isg {

struct ValidationReport {
    std::vector<std::string> errors, warnings;
    bool is_valid() const { return errors.empty(); }
};

ValidationReport validateTopology(const IntersectionInput&);

class IntersectionShapeGenerator {
public:
    struct Config {
        double sdf_cell_size = 0.2, obstacle_buffer = 0.4, kappa_max = 0.25;
        LBFGSConfig lbfgs;
        ConnectivityDirectionConfig connectivity_direction;
    };

    IntersectionShapeGenerator();
    explicit IntersectionShapeGenerator(const Config& cfg);
    bool generate(const IntersectionInput&, IntersectionOutput&);
    const ValidationReport& lastReport() const { return report_; }

private:
    Config cfg_;
    ValidationReport report_;
};

}
