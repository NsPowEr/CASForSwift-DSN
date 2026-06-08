// F7.2-T4 — Continuous distributions: χ², Student-t, Fisher-Snedecor F.
//
// Each distribution exposes a numerical PDF and CDF.  The CDFs reuse the
// incomplete-gamma and incomplete-beta core functions provided here as
// private helpers:
//   - `regularised_lower_gamma(s, x)` via series for x < s+1 and
//     continued fraction otherwise (Numerical Recipes §6.2).
//   - `regularised_beta_inc(x, a, b)` via Lentz continued fraction
//     (Numerical Recipes §6.4) — the symmetry transformation
//     I_x(a, b) = 1 − I_{1−x}(b, a) is applied to keep convergence in the
//     fast-converging half-plane.
//
// All routines use `lgamma` for stability and short-circuit boundary
// cases (x ≤ 0 on χ²/F, x = 0 on F, etc.) before the iterative core.

#include "cas/statistics.hpp"

#include <cmath>
#include <numbers>

namespace cas::statistics {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    CASError err{};
    err.kind = kind;
    err.message = std::move(message);
    return err;
}

// Iteration caps chosen to give ~14 decimal digits on the entire numeric
// range supported by IEEE-754 doubles.
constexpr int kIterMax = 200;
constexpr double kEps = 1.0e-15;

// Series form of the regularised lower incomplete gamma function for
// x < s + 1 (Numerical Recipes §6.2).
[[nodiscard]] double gamma_inc_series(double s, double x) {
    double ap = s;
    double sum = 1.0 / s;
    double term = sum;
    for (int n = 0; n < kIterMax; ++n) {
        ap += 1.0;
        term *= x / ap;
        sum += term;
        if (std::abs(term) < std::abs(sum) * kEps) break;
    }
    return sum * std::exp(-x + s * std::log(x) - std::lgamma(s));
}

// Continued fraction form of the regularised upper incomplete gamma
// function for x ≥ s + 1 (Numerical Recipes §6.2, Lentz).
[[nodiscard]] double gamma_inc_cf(double s, double x) {
    double b = x + 1.0 - s;
    double c = 1.0 / std::numeric_limits<double>::min();
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i <= kIterMax; ++i) {
        const double an = -i * (i - s);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
        c = b + an / c;
        if (std::abs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::abs(delta - 1.0) < kEps) break;
    }
    return std::exp(-x + s * std::log(x) - std::lgamma(s)) * h;
}

[[nodiscard]] double regularised_lower_gamma(double s, double x) {
    if (x < 0.0 || s <= 0.0) return 0.0;
    if (x == 0.0) return 0.0;
    if (x < s + 1.0) return gamma_inc_series(s, x);
    return 1.0 - gamma_inc_cf(s, x);
}

// Continued fraction for the regularised incomplete Beta function.
[[nodiscard]] double betacf(double a, double b, double x) {
    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= kIterMax; ++m) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
        c = 1.0 + aa / c;
        if (std::abs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
        c = 1.0 + aa / c;
        if (std::abs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < kEps) break;
    }
    return h;
}

[[nodiscard]] double regularised_beta_inc(double x, double a, double b) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double front_log = std::lgamma(a + b)
                           - std::lgamma(a)
                           - std::lgamma(b)
                           + a * std::log(x)
                           + b * std::log1p(-x);
    const double front = std::exp(front_log);
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return front * betacf(a, b, x) / a;
    }
    return 1.0 - front * betacf(b, a, 1.0 - x) / b;
}

}  // namespace

// ── Chi-squared ──────────────────────────────────────────────────────────────

Result<double> chi_squared_pdf(double x, double k) {
    if (!(k > 0.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "chi_squared_pdf: degrees of freedom k must be positive"));
    }
    if (x < 0.0) return ok(0.0);
    if (x == 0.0) return ok(k < 2.0 ? std::numeric_limits<double>::infinity()
                                     : (k == 2.0 ? 0.5 : 0.0));
    const double half_k = 0.5 * k;
    const double log_pdf = (half_k - 1.0) * std::log(x)
                         - 0.5 * x
                         - half_k * std::log(2.0)
                         - std::lgamma(half_k);
    return ok(std::exp(log_pdf));
}

Result<double> chi_squared_cdf(double x, double k) {
    if (!(k > 0.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "chi_squared_cdf: degrees of freedom k must be positive"));
    }
    if (x <= 0.0) return ok(0.0);
    return ok(regularised_lower_gamma(0.5 * k, 0.5 * x));
}

// ── Student-t ────────────────────────────────────────────────────────────────

Result<double> student_t_pdf(double x, double nu) {
    if (!(nu > 0.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "student_t_pdf: degrees of freedom nu must be positive"));
    }
    const double log_norm = std::lgamma(0.5 * (nu + 1.0))
                          - std::lgamma(0.5 * nu)
                          - 0.5 * std::log(nu * std::numbers::pi);
    const double log_kernel = -0.5 * (nu + 1.0) * std::log1p(x * x / nu);
    return ok(std::exp(log_norm + log_kernel));
}

Result<double> student_t_cdf(double x, double nu) {
    if (!(nu > 0.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "student_t_cdf: degrees of freedom nu must be positive"));
    }
    // F_t(x) = 1 − ½ I_{nu/(nu+x²)}(nu/2, 1/2)   for x ≥ 0
    //        =       ½ I_{nu/(nu+x²)}(nu/2, 1/2)   for x < 0
    const double t2 = x * x;
    const double w = nu / (nu + t2);
    const double half_inc = 0.5 * regularised_beta_inc(w, 0.5 * nu, 0.5);
    return ok(x >= 0.0 ? 1.0 - half_inc : half_inc);
}

// ── Fisher-Snedecor F ────────────────────────────────────────────────────────

Result<double> f_pdf(double x, double d1, double d2) {
    if (!(d1 > 0.0) || !(d2 > 0.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "f_pdf: degrees of freedom must be positive"));
    }
    if (x <= 0.0) return ok(0.0);
    const double a = 0.5 * d1;
    const double b = 0.5 * d2;
    const double log_pdf = a * std::log(d1)
                         + b * std::log(d2)
                         + (a - 1.0) * std::log(x)
                         - (a + b) * std::log(d2 + d1 * x)
                         - std::lgamma(a)
                         - std::lgamma(b)
                         + std::lgamma(a + b);
    return ok(std::exp(log_pdf));
}

Result<double> f_cdf(double x, double d1, double d2) {
    if (!(d1 > 0.0) || !(d2 > 0.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "f_cdf: degrees of freedom must be positive"));
    }
    if (x <= 0.0) return ok(0.0);
    const double w = (d1 * x) / (d1 * x + d2);
    return ok(regularised_beta_inc(w, 0.5 * d1, 0.5 * d2));
}

}  // namespace cas::statistics
