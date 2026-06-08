#include "cas/numeric.hpp"
#include <cmath>

namespace cas::numeric {

namespace {

Result<double> simpson_step(
    double a,
    double b,
    double fa,
    double fb,
    double fm) {
    
    return ok((b - a) / 6.0 * (fa + 4.0 * fm + fb));
}

Result<double> adaptive_simpson(
    const std::function<Result<double>(double)>& f,
    double a,
    double b,
    double tolerance,
    double whole,
    double fa,
    double fb,
    double fm,
    std::uint32_t depth) {

    double mid = (a + b) / 2.0;
    double left_mid = (a + mid) / 2.0;
    double right_mid = (mid + b) / 2.0;

    auto f_lm = f(left_mid);
    if (f_lm.is_error()) return f_lm;
    auto f_rm = f(right_mid);
    if (f_rm.is_error()) return f_rm;

    auto left_res = simpson_step(a, mid, fa, fm, f_lm.value());
    if (left_res.is_error()) return left_res;
    auto right_res = simpson_step(mid, b, fm, fb, f_rm.value());
    if (right_res.is_error()) return right_res;

    double left = left_res.value();
    double right = right_res.value();

    if (depth == 0 || std::abs(left + right - whole) <= 15.0 * tolerance) {
        return ok(left + right + (left + right - whole) / 15.0);
    }

    auto left_total = adaptive_simpson(f, a, mid, tolerance / 2.0, left, fa, fm, f_lm.value(), depth - 1);
    if (left_total.is_error()) return left_total;
    auto right_total = adaptive_simpson(f, mid, b, tolerance / 2.0, right, fm, fb, f_rm.value(), depth - 1);
    if (right_total.is_error()) return right_total;

    return ok(left_total.value() + right_total.value());
}

// ── F7.3-T1: Gauss-Kronrod 15-point quadrature on [-1, 1] ────────────────────
//
// Kronrod extends a 7-point Gauss-Legendre rule with 8 additional nodes so a
// 15-point estimate (G_15) and a nested 7-point estimate (K_7) share function
// values.  The two estimates allow a cheap a-posteriori error estimate
// |G_15 - K_7| at the cost of 15 function evaluations per panel.
//
// Tables: standard QUADPACK constants (Piessens et al. 1983, QUADPACK §4),
// rounded to 18 significant digits.

constexpr double kGK15_Nodes[15] = {
    -0.991455371120812639206854697526329,
    -0.949107912342758524526189684047851,
    -0.864864423359769072789712788640926,
    -0.741531185599394439863864773280788,
    -0.586087235467691130294144838258730,
    -0.405845151377397166906606412076961,
    -0.207784955007898467600689403773245,
     0.0,
     0.207784955007898467600689403773245,
     0.405845151377397166906606412076961,
     0.586087235467691130294144838258730,
     0.741531185599394439863864773280788,
     0.864864423359769072789712788640926,
     0.949107912342758524526189684047851,
     0.991455371120812639206854697526329};

constexpr double kGK15_Weights[15] = {
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714,
    0.204432940075298892414161999234649,
    0.190350578064785409913256402421014,
    0.169004726639267902826583426598550,
    0.140653259715525918745189590510238,
    0.104790010322250183839876322541518,
    0.063092092629978553290700663189204,
    0.022935322010529224963732008058970};

// 7-point Gauss-Legendre weights mapping onto the 7 odd-index Kronrod nodes
// (indices 1, 3, 5, 7, 9, 11, 13).  Even indices are the Kronrod-only nodes.
constexpr double kG7_Weights[7] = {
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327,
    0.381830050505118944950369775488975,
    0.279705391489276667901467771423780,
    0.129484966168869693270611432679082};

constexpr int kG7_NodeMap[7] = {1, 3, 5, 7, 9, 11, 13};

[[nodiscard]] Result<std::pair<double, double>> gauss_kronrod_15(
    const std::function<Result<double>(double)>& f,
    double a,
    double b)
{
    const double half = 0.5 * (b - a);
    const double mid = 0.5 * (a + b);
    double kronrod = 0.0;
    double gauss = 0.0;
    double f_vals[15];
    for (int i = 0; i < 15; ++i) {
        const double x = mid + half * kGK15_Nodes[i];
        auto fv = f(x);
        if (fv.is_error()) return fail<std::pair<double, double>>(fv.error());
        f_vals[i] = fv.value();
        kronrod += kGK15_Weights[i] * f_vals[i];
    }
    for (int i = 0; i < 7; ++i) {
        gauss += kG7_Weights[i] * f_vals[kG7_NodeMap[i]];
    }
    return ok(std::make_pair(half * kronrod, half * gauss));
}

[[nodiscard]] Result<double> adaptive_gauss_kronrod(
    const std::function<Result<double>(double)>& f,
    double a,
    double b,
    double tolerance,
    std::uint32_t depth)
{
    auto est = gauss_kronrod_15(f, a, b);
    if (est.is_error()) return fail<double>(est.error());
    const double kronrod = est.value().first;
    const double gauss   = est.value().second;
    const double err = std::abs(kronrod - gauss);
    // QUADPACK error scaling: err ≈ (200·|G-K|)^{1.5}.  We use the direct
    // |G-K| as the conservative bound; it is always ≥ the QUADPACK form.
    if (depth == 0 || err <= tolerance) {
        return ok(kronrod);
    }
    const double mid = 0.5 * (a + b);
    auto left  = adaptive_gauss_kronrod(f, a, mid, tolerance * 0.5, depth - 1);
    if (left.is_error()) return left;
    auto right = adaptive_gauss_kronrod(f, mid, b, tolerance * 0.5, depth - 1);
    if (right.is_error()) return right;
    return ok(left.value() + right.value());
}

} // namespace

Result<double> integrate_numeric(
    ExprPtr expr,
    const std::string& variable,
    double a,
    double b,
    const IntegrationOptions& options) {

    auto f = [&](double x) -> Result<double> {
        return eval(expr, {{variable, x}});
    };

    if (options.scheme == IntegrationScheme::GaussKronrod) {
        return adaptive_gauss_kronrod(f, a, b, options.tolerance, options.max_depth);
    }

    // Default: adaptive Simpson (legacy).
    auto fa = f(a); if (fa.is_error()) return fa;
    auto fb = f(b); if (fb.is_error()) return fb;
    auto fm = f((a + b) / 2.0); if (fm.is_error()) return fm;

    auto whole = simpson_step(a, b, fa.value(), fb.value(), fm.value());
    if (whole.is_error()) return whole;

    return adaptive_simpson(f, a, b, options.tolerance, whole.value(), fa.value(), fb.value(), fm.value(), options.max_depth);
}

} // namespace cas::numeric
