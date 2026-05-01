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

} // namespace cas::numeric
