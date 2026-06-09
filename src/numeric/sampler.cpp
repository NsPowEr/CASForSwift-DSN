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

// ─── F7.1-B1 Implicit / Contour / VectorField samplers ─────────────────

ImplicitSampler::ImplicitSampler(ExprPtr expr, std::string var_x,
                                 std::string var_y, Options options)
    : expr_(expr), var_x_(std::move(var_x)),
      var_y_(std::move(var_y)), options_(options) {}

Result<double> ImplicitSampler::f(double x, double y) {
    return eval(expr_, {{var_x_, x}, {var_y_, y}});
}

namespace {

// Linear interpolation along a cell edge: returns the parameter t ∈ [0,1]
// where (1-t)*va + t*vb == 0 (level=0 baseline applied to caller).
double edge_zero(double va, double vb) {
    if (vb == va) return 0.5;
    return va / (va - vb);
}

}  // namespace

Result<std::vector<SampleSegment>> ImplicitSampler::sample(
    double x_min, double x_max, double y_min, double y_max) {
    if (x_min >= x_max || y_min >= y_max || options_.nx < 2 || options_.ny < 2) {
        return ok(std::vector<SampleSegment>{});
    }

    const std::uint32_t nx = options_.nx;
    const std::uint32_t ny = options_.ny;
    const double dx = (x_max - x_min) / static_cast<double>(nx - 1);
    const double dy = (y_max - y_min) / static_cast<double>(ny - 1);
    const double level = options_.level;

    // Pre-evaluate grid values; record NaN/inf as level (skip cells).
    std::vector<double> grid(static_cast<std::size_t>(nx) * ny);
    for (std::uint32_t j = 0; j < ny; ++j) {
        const double y = y_min + dy * static_cast<double>(j);
        for (std::uint32_t i = 0; i < nx; ++i) {
            const double x = x_min + dx * static_cast<double>(i);
            auto v = f(x, y);
            const double s = v.is_ok() ? v.value() : std::nan("");
            grid[static_cast<std::size_t>(j) * nx + i] =
                std::isfinite(s) ? (s - level) : std::nan("");
        }
    }

    std::vector<SampleSegment> out;
    auto sign_bit = [](double v) -> int { return (v >= 0.0) ? 1 : 0; };

    auto edge_point = [&](double x_a, double y_a, double x_b, double y_b,
                          double va, double vb) {
        const double t = edge_zero(va, vb);
        return SamplePoint{x_a + (x_b - x_a) * t, y_a + (y_b - y_a) * t};
    };

    for (std::uint32_t j = 0; j + 1 < ny; ++j) {
        for (std::uint32_t i = 0; i + 1 < nx; ++i) {
            const double v00 = grid[static_cast<std::size_t>(j) * nx + i];
            const double v10 = grid[static_cast<std::size_t>(j) * nx + (i + 1)];
            const double v11 = grid[static_cast<std::size_t>(j + 1) * nx + (i + 1)];
            const double v01 = grid[static_cast<std::size_t>(j + 1) * nx + i];
            if (std::isnan(v00) || std::isnan(v10) || std::isnan(v11) || std::isnan(v01)) {
                continue;
            }

            const int idx =
                (sign_bit(v00) << 0) |
                (sign_bit(v10) << 1) |
                (sign_bit(v11) << 2) |
                (sign_bit(v01) << 3);

            const double x0 = x_min + dx * static_cast<double>(i);
            const double x1 = x0 + dx;
            const double y0 = y_min + dy * static_cast<double>(j);
            const double y1 = y0 + dy;

            // Edge points labelled by adjacent vertices.
            auto eb = [&]() { return edge_point(x0, y0, x1, y0, v00, v10); };
            auto er = [&]() { return edge_point(x1, y0, x1, y1, v10, v11); };
            auto et = [&]() { return edge_point(x0, y1, x1, y1, v01, v11); };
            auto el = [&]() { return edge_point(x0, y0, x0, y1, v00, v01); };

            // Marching-squares case table (16 entries). Cases 5/10 = saddle
            // disambiguation; here we resolve by cell-average sign.
            switch (idx) {
                case 0:  case 15: break;                              // empty
                case 1:  case 14: out.push_back({el(), eb()}); break; // bottom-left corner
                case 2:  case 13: out.push_back({eb(), er()}); break; // bottom-right
                case 4:  case 11: out.push_back({er(), et()}); break; // top-right
                case 8:  case 7:  out.push_back({et(), el()}); break; // top-left
                case 3:  case 12: out.push_back({el(), er()}); break; // horizontal cut
                case 6:  case 9:  out.push_back({eb(), et()}); break; // vertical cut
                case 5: {
                    const double center = 0.25 * (v00 + v10 + v11 + v01);
                    if (center >= 0.0) {
                        out.push_back({el(), et()});
                        out.push_back({eb(), er()});
                    } else {
                        out.push_back({el(), eb()});
                        out.push_back({et(), er()});
                    }
                    break;
                }
                case 10: {
                    const double center = 0.25 * (v00 + v10 + v11 + v01);
                    if (center >= 0.0) {
                        out.push_back({el(), eb()});
                        out.push_back({et(), er()});
                    } else {
                        out.push_back({el(), et()});
                        out.push_back({eb(), er()});
                    }
                    break;
                }
            }
        }
    }
    return ok(std::move(out));
}

ContourSampler::ContourSampler(ExprPtr expr, std::string var_x,
                               std::string var_y,
                               std::vector<double> levels,
                               ImplicitSampler::Options options)
    : expr_(expr), var_x_(std::move(var_x)), var_y_(std::move(var_y)),
      levels_(std::move(levels)), options_(options) {}

Result<std::vector<std::vector<SampleSegment>>> ContourSampler::sample(
    double x_min, double x_max, double y_min, double y_max) {
    std::vector<std::vector<SampleSegment>> out;
    out.reserve(levels_.size());
    for (double lv : levels_) {
        ImplicitSampler::Options o = options_;
        o.level = lv;
        ImplicitSampler is(expr_, var_x_, var_y_, o);
        auto segs = is.sample(x_min, x_max, y_min, y_max);
        if (segs.is_error()) return fail<std::vector<std::vector<SampleSegment>>>(segs.error());
        out.push_back(std::move(segs.value()));
    }
    return ok(std::move(out));
}

VectorFieldSampler::VectorFieldSampler(ExprPtr fx_expr, ExprPtr fy_expr,
                                       std::string var_x, std::string var_y,
                                       Options options)
    : fx_expr_(fx_expr), fy_expr_(fy_expr),
      var_x_(std::move(var_x)), var_y_(std::move(var_y)),
      options_(options) {}

Result<std::vector<VectorArrow>> VectorFieldSampler::sample(
    double x_min, double x_max, double y_min, double y_max) {
    if (x_min >= x_max || y_min >= y_max || options_.nx < 1 || options_.ny < 1) {
        return ok(std::vector<VectorArrow>{});
    }
    const std::uint32_t nx = options_.nx;
    const std::uint32_t ny = options_.ny;
    const double dx = (x_max - x_min) / static_cast<double>(nx);
    const double dy = (y_max - y_min) / static_cast<double>(ny);
    std::vector<VectorArrow> out;
    out.reserve(static_cast<std::size_t>(nx) * ny);
    for (std::uint32_t j = 0; j < ny; ++j) {
        const double y = y_min + dy * (static_cast<double>(j) + 0.5);
        for (std::uint32_t i = 0; i < nx; ++i) {
            const double x = x_min + dx * (static_cast<double>(i) + 0.5);
            auto vx = eval(fx_expr_, {{var_x_, x}, {var_y_, y}});
            auto vy = eval(fy_expr_, {{var_x_, x}, {var_y_, y}});
            if (vx.is_error() || vy.is_error()) continue;
            const double fx = vx.value();
            const double fy = vy.value();
            if (!std::isfinite(fx) || !std::isfinite(fy)) continue;
            out.push_back(VectorArrow{
                .base = {x, y},
                .dir  = {fx, fy},
            });
        }
    }
    return ok(std::move(out));
}

} // namespace cas::numeric
