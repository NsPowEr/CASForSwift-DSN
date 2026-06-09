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

// F7.1-B1 — Implicit curve sampler (marching squares).
//
// Given a bivariate scalar field f(x,y) and a rectangular domain
// [x_min,x_max] × [y_min,y_max] sampled on a regular grid, emits the
// piecewise-linear approximation of the level set { (x,y) : f(x,y) = level }.
// Output: a list of line segments (each segment = pair of SamplePoint).
struct SampleSegment {
    SamplePoint a;
    SamplePoint b;
};

class ImplicitSampler {
public:
    struct Options {
        std::uint32_t nx;
        std::uint32_t ny;
        double level;
        static Options default_options() {
            return {64, 64, 0.0};
        }
    };

    ImplicitSampler(ExprPtr expr, std::string var_x, std::string var_y,
                    Options options = Options::default_options());

    [[nodiscard]] Result<std::vector<SampleSegment>> sample(
        double x_min, double x_max, double y_min, double y_max);

private:
    [[nodiscard]] Result<double> f(double x, double y);

    ExprPtr expr_;
    std::string var_x_;
    std::string var_y_;
    Options options_;
};

// F7.1-B1 — Contour sampler: one ImplicitSampler per requested level.
class ContourSampler {
public:
    ContourSampler(ExprPtr expr, std::string var_x, std::string var_y,
                   std::vector<double> levels,
                   ImplicitSampler::Options options = ImplicitSampler::Options::default_options());

    [[nodiscard]] Result<std::vector<std::vector<SampleSegment>>> sample(
        double x_min, double x_max, double y_min, double y_max);

private:
    ExprPtr expr_;
    std::string var_x_;
    std::string var_y_;
    std::vector<double> levels_;
    ImplicitSampler::Options options_;
};

// F7.1-B1 — Vector field sampler.
// At each (x,y) grid node, evaluates the 2-component field (fx, fy)
// and returns (base, displacement) pairs ready for arrow rendering.
struct VectorArrow {
    SamplePoint base;
    SamplePoint dir;   // unnormalized displacement (fx, fy) at base
};

class VectorFieldSampler {
public:
    struct Options {
        std::uint32_t nx;
        std::uint32_t ny;
        static Options default_options() { return {16, 16}; }
    };

    VectorFieldSampler(ExprPtr fx_expr, ExprPtr fy_expr,
                       std::string var_x, std::string var_y,
                       Options options = Options::default_options());

    [[nodiscard]] Result<std::vector<VectorArrow>> sample(
        double x_min, double x_max, double y_min, double y_max);

private:
    ExprPtr fx_expr_;
    ExprPtr fy_expr_;
    std::string var_x_;
    std::string var_y_;
    Options options_;
};

} // namespace cas::numeric
