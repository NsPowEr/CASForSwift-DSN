#include "cas/numeric/sampler.hpp"
#include "cas/numeric.hpp"
#include <cmath>
#include <algorithm>

namespace cas::numeric {

AdaptiveSampler::AdaptiveSampler(ExprPtr expr, std::string variable, Options options)
    : expr_(expr), var_(std::move(variable)), options_(options) {}

Result<double> AdaptiveSampler::f(double x) {
    return eval(expr_, {{var_, x}});
}

void AdaptiveSampler::sample_recursive(double x1, double y1, double x2, double y2, std::vector<SamplePoint>& points, std::uint32_t depth) {
    if (depth == 0) return;

    double xm = (x1 + x2) / 2.0;
    auto ym_res = f(xm);
    if (ym_res.is_error()) return;
    double ym = ym_res.value();

    if (std::isinf(ym) || std::isnan(ym)) return;

    // Controllo linearità: se il punto centrale è vicino al segmento (x1,y1)-(x2,y2)
    // Usiamo la distanza punto-retta o più semplicemente la differenza di pendenza
    double y_linear = (y1 + y2) / 2.0;
    double dy = std::abs(ym - y_linear);
    
    // Tolleranza pesata sulla scala locale
    double local_tol = options_.tolerance * (1.0 + std::abs(y1) + std::abs(y2));

    if (dy > local_tol) {
        sample_recursive(x1, y1, xm, ym, points, depth - 1);
        points.push_back({xm, ym});
        sample_recursive(xm, ym, x2, y2, points, depth - 1);
    }
}

Result<std::vector<SamplePoint>> AdaptiveSampler::sample(double x_min, double x_max) {
    if (x_min >= x_max) return ok(std::vector<SamplePoint>{});

    std::vector<SamplePoint> points;
    double step = (x_max - x_min) / options_.initial_samples;

    for (std::uint32_t i = 0; i <= options_.initial_samples; ++i) {
        double x1 = x_min + i * step;
        auto y1_res = f(x1);
        if (y1_res.is_error()) continue;
        double y1 = y1_res.value();
        
        if (i > 0) {
            double x0 = points.back().x;
            double y0 = points.back().y;
            sample_recursive(x0, y0, x1, y1, points, options_.max_depth);
        }
        
        points.push_back({x1, y1});
    }

    // Ordiniamo i punti (la ricorsione potrebbe averli aggiunti fuori ordine)
    std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) {
        return a.x < b.x;
    });

    return ok(std::move(points));
}

// ── F7.1-T1: Parametric sampler ──────────────────────────────────────────────

ParametricSampler::ParametricSampler(
    ExprPtr x_expr, ExprPtr y_expr, std::string variable, Options options)
    : x_expr_(x_expr), y_expr_(y_expr), var_(std::move(variable)),
      options_(options) {}

Result<SamplePoint> ParametricSampler::eval_at(double t) {
    auto xv = eval(x_expr_, {{var_, t}});
    if (xv.is_error()) return fail<SamplePoint>(xv.error());
    auto yv = eval(y_expr_, {{var_, t}});
    if (yv.is_error()) return fail<SamplePoint>(yv.error());
    return ok(SamplePoint{xv.value(), yv.value()});
}

void ParametricSampler::refine(
    double t_a, const SamplePoint& pa,
    double t_b, const SamplePoint& pb,
    std::vector<SamplePoint>& out, std::uint32_t depth)
{
    if (depth == 0) return;
    const double t_mid = 0.5 * (t_a + t_b);
    auto pm_res = eval_at(t_mid);
    if (pm_res.is_error()) return;
    const SamplePoint pm = pm_res.value();
    if (!std::isfinite(pm.x) || !std::isfinite(pm.y)) return;

    // Cosine of the turning angle at pm using the chord vectors pa→pm and
    // pm→pb.  Two collinear chords yield cos = +1; large turns drive cos
    // toward −1.  We refine when cos < (1 − angle_tolerance) so the
    // tolerance maps to "max permitted angle deviation from straight" in
    // small-angle approximation.
    const double v1x = pm.x - pa.x, v1y = pm.y - pa.y;
    const double v2x = pb.x - pm.x, v2y = pb.y - pm.y;
    const double n1 = std::hypot(v1x, v1y);
    const double n2 = std::hypot(v2x, v2y);
    const double cos_angle = (n1 == 0.0 || n2 == 0.0)
        ? 1.0
        : (v1x * v2x + v1y * v2y) / (n1 * n2);
    const double threshold = 1.0 - options_.angle_tolerance;

    if (cos_angle < threshold) {
        refine(t_a, pa, t_mid, pm, out, depth - 1);
        out.push_back(pm);
        refine(t_mid, pm, t_b, pb, out, depth - 1);
    }
}

Result<std::vector<SamplePoint>> ParametricSampler::sample(
    double t_min, double t_max)
{
    if (t_min >= t_max) return ok(std::vector<SamplePoint>{});
    std::vector<SamplePoint> out;
    const double step = (t_max - t_min) /
        static_cast<double>(options_.initial_samples);
    auto p0 = eval_at(t_min);
    if (p0.is_error()) return fail<std::vector<SamplePoint>>(p0.error());
    out.push_back(p0.value());
    for (std::uint32_t i = 1; i <= options_.initial_samples; ++i) {
        const double t = t_min + step * static_cast<double>(i);
        auto p_res = eval_at(t);
        if (p_res.is_error()) continue;
        const SamplePoint pi = p_res.value();
        refine(t_min + step * static_cast<double>(i - 1),
               out.back(), t, pi, out, options_.max_depth);
        out.push_back(pi);
    }
    return ok(std::move(out));
}

} // namespace cas::numeric
