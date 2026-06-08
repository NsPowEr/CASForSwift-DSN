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

/// F7.1-T1 — Parametric curve sampler.
///
/// Walks a parameter t ∈ [t_min, t_max] producing (x(t), y(t)) points
/// adaptively refined on geometric curvature: three consecutive sample
/// points whose middle vertex lies above the `angle_tolerance` from
/// collinear with the outer two trigger a midpoint subdivision.
class ParametricSampler {
public:
    struct Options {
        double angle_tolerance;      // radians; cos-based check threshold
        std::uint32_t max_depth;
        std::uint32_t initial_samples;
        static Options default_options() {
            return {0.05, 12, 32};
        }
    };

    ParametricSampler(
        ExprPtr x_expr, ExprPtr y_expr, std::string variable,
        Options options = Options::default_options());

    [[nodiscard]] Result<std::vector<SamplePoint>> sample(double t_min, double t_max);

private:
    [[nodiscard]] Result<SamplePoint> eval_at(double t);
    void refine(
        double t_a, const SamplePoint& pa,
        double t_b, const SamplePoint& pb,
        std::vector<SamplePoint>& out, std::uint32_t depth);

    ExprPtr x_expr_;
    ExprPtr y_expr_;
    std::string var_;
    Options options_;
};

} // namespace cas::numeric
