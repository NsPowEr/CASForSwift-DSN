#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include <vector>
#include <string>

namespace cas::numeric {

struct SamplePoint {
    double x;
    double y;
};

class AdaptiveSampler {
public:
    struct Options {
        double tolerance;
        std::uint32_t max_depth;
        std::uint32_t initial_samples;

        static Options default_options() {
            return {1e-3, 12, 20};
        }
    };

    explicit AdaptiveSampler(ExprPtr expr, std::string variable, Options options = Options::default_options());

    [[nodiscard]] Result<std::vector<SamplePoint>> sample(double x_min, double x_max);

private:
    [[nodiscard]] Result<double> f(double x);
    void sample_recursive(double x1, double y1, double x2, double y2, std::vector<SamplePoint>& points, std::uint32_t depth);

    ExprPtr expr_;
    std::string var_;
    Options options_;
};

} // namespace cas::numeric
