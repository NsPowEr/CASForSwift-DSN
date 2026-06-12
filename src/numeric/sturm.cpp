// Sturm-sequence based real-root isolation (Sturm 1829).
//
// For a squarefree univariate polynomial f ∈ Q[x], the Sturm sequence
//
//     s_0 = f,  s_1 = f',  s_{k+1} = -rem(s_{k-1}, s_k)   until s_K = const
//
// has the following property: for any a < b that are not roots of f,
// the number of distinct real roots of f in (a, b] equals
//
//     V(a) − V(b)
//
// where V(c) counts sign variations in the finite sequence
//     [s_0(c), s_1(c), …, s_K(c)]
// after removing zero terms.
//
// We use bisection to isolate each root within a half-open interval of
// width < tol, then run a few Newton-Raphson steps from the midpoint
// for high-precision polish.

#include "cas/numeric.hpp"
#include "cas/numeric_bigfloat.hpp"
#include "cas/symbolic.hpp"

#include "../algebra/polynomial_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::numeric {

namespace {

using ::cas::algebra::RatPoly;

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

// Derivative of a univariate rational polynomial (Q[x] → Q[x]).
[[nodiscard]] RatPoly differentiate_rat_poly(const RatPoly& f) {
    if (f.size() <= 1U) return RatPoly{};
    std::vector<Rational> d;
    d.reserve(f.size() - 1U);
    for (std::size_t k = 1; k < f.size(); ++k) {
        d.push_back(f[k] * Rational(static_cast<long long>(k)));
    }
    RatPoly out(std::move(d));
    out.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return out;
}

// Sturm sequence. By definition s_0 = f, s_1 = f',  s_{k+1} = -rem(s_{k-1}, s_k).
[[nodiscard]] std::vector<RatPoly> sturm_sequence(const RatPoly& f) {
    std::vector<RatPoly> seq;
    if (f.is_zero()) return seq;
    seq.reserve(f.size());
    seq.push_back(f);
    RatPoly fp = differentiate_rat_poly(f);
    if (fp.is_zero()) {
        return seq;
    }
    seq.push_back(std::move(fp));
    while (true) {
        const RatPoly& a = seq[seq.size() - 2U];
        const RatPoly& b = seq.back();
        if (b.is_zero()) break;
        auto [q, r] = ::cas::algebra::div_rem_rational_poly(a, b);
        (void)q;
        if (r.is_zero()) break;
        // s_{k+1} = -r
        std::vector<Rational> neg;
        neg.reserve(r.size());
        for (std::size_t k = 0; k < r.size(); ++k) neg.push_back(Rational(0) - r[k]);
        RatPoly next(std::move(neg));
        next.normalize([](const Rational& v) { return v.numerator().is_zero(); });
        if (next.is_zero()) break;
        seq.push_back(std::move(next));
    }
    return seq;
}

[[nodiscard]] int rational_sign(const Rational& v) noexcept {
    if (v.numerator().is_zero()) return 0;
    return v.numerator().is_negative() ? -1 : 1;
}

[[nodiscard]] unsigned int sign_variations(const std::vector<RatPoly>& seq, const Rational& x) {
    unsigned int variations = 0;
    int prev_sign = 0;
    for (const RatPoly& p : seq) {
        const Rational v = ::cas::algebra::evaluate_rational_polynomial_at(p, x);
        const int s = rational_sign(v);
        if (s == 0) continue;
        if (prev_sign != 0 && s != prev_sign) ++variations;
        prev_sign = s;
    }
    return variations;
}

// Convert a double to a Rational with manageable denominator. Encoding
// the IEEE 754 bit pattern exactly produces astronomically large
// rationals; for Sturm bisection we only need O(2^-52) accuracy.
[[nodiscard]] Rational double_to_rational(double x) {
    if (std::isnan(x) || std::isinf(x)) return Rational(0);
    // Scale by 2^32 then truncate to long long for a 32-bit-fraction Rational.
    // This gives ~1e-9 resolution, plenty for the bisection cap; Newton polish
    // recovers double precision afterwards.
    constexpr double scale = 4294967296.0; // 2^32
    const double scaled = x * scale;
    const long long num = static_cast<long long>(std::llround(scaled));
    return Rational(BigInt(num), BigInt(static_cast<long long>(scale)));
}

[[nodiscard]] double rational_to_double(const Rational& r) {
    const auto num_d = r.numerator().to_double();
    const auto den_d = r.denominator().to_double();
    if (den_d == 0.0) return 0.0;
    return num_d / den_d;
}

// Recursively bisect [a, b] in *rational* arithmetic, collecting
// isolating intervals (each containing exactly one real root of f).
// Termination:
//   * V(a) - V(b) == 0  →  no roots, return
//   * V(a) - V(b) == 1 && b - a < tol  →  this interval isolates one root
//   * otherwise split at midpoint
//
// `max_depth` bounds the recursion structurally (≤ ceil(log2((b-a)/tol)) + 2);
// this is a derived bound, not a hardcoded safety cap.
void isolate_recursive(const std::vector<RatPoly>& seq,
                       const Rational& a, const Rational& b, const Rational& tol,
                       unsigned int max_depth, unsigned int depth,
                       std::vector<std::pair<Rational, Rational>>& out) {
    if (depth > max_depth) return;
    const unsigned int va = sign_variations(seq, a);
    const unsigned int vb = sign_variations(seq, b);
    if (va <= vb) return;
    const unsigned int n_roots = va - vb;
    const Rational width = b - a;
    if (n_roots == 1U) {
        // Compare width to tol via cross-multiplication on positive numerators.
        if ((width.numerator() * tol.denominator()) <= (tol.numerator() * width.denominator())) {
            out.emplace_back(a, b);
            return;
        }
    }
    const Rational mid = (a + b) * Rational(BigInt(1), BigInt(2));
    isolate_recursive(seq, a, mid, tol, max_depth, depth + 1U, out);
    isolate_recursive(seq, mid, b, tol, max_depth, depth + 1U, out);
}

// Newton polish on the *symbolic* polynomial from the midpoint of the
// isolating interval. Number of iterations is derived from quadratic
// convergence: ceil(log2(log2(width / tol))) + safety = small.
[[nodiscard]] double newton_polish_rat(const RatPoly& f, const RatPoly& fp,
                                       double initial, double tol,
                                       unsigned int max_iter) {
    double x = initial;
    for (unsigned int i = 0; i < max_iter; ++i) {
        const Rational xr = double_to_rational(x);
        const double fx = rational_to_double(
            ::cas::algebra::evaluate_rational_polynomial_at(f, xr));
        if (std::abs(fx) < tol) return x;
        const double fpx = rational_to_double(
            ::cas::algebra::evaluate_rational_polynomial_at(fp, xr));
        if (std::abs(fpx) < 1e-300) return x;
        const double dx = fx / fpx;
        x -= dx;
        if (std::abs(dx) < tol) return x;
    }
    return x;
}

} // namespace

Result<std::vector<double>> find_polynomial_roots_sturm(
    ExprPtr expr,
    const std::string& variable,
    symbolic::CASContext& ctx,
    double low,
    double high,
    double tol) {

    if (!expr) {
        return fail<std::vector<double>>(make_error(
            CASErrorKind::InvalidArgument, "find_polynomial_roots_sturm: null expr"));
    }
    if (low >= high) {
        return ok(std::vector<double>{});
    }

    // 1. Parse to PolyExpr in the named variable, then drop to RatPoly.
    const Symbol var(variable);
    auto poly_res = ::cas::algebra::parse_polynomial(expr, var, ctx);
    if (poly_res.is_error()) return fail<std::vector<double>>(poly_res.error());
    auto rat_res = ::cas::algebra::poly_to_rational_poly(poly_res.value());
    if (rat_res.is_error()) return fail<std::vector<double>>(rat_res.error());
    RatPoly f = std::move(rat_res.value());
    if (f.is_zero() || f.size() <= 1U) {
        return ok(std::vector<double>{});
    }

    // 2. Squarefree decomposition: f / gcd(f, f').  Sturm counts distinct
    //    real roots; we drop multiplicity by dividing out gcd(f, f').
    RatPoly fp = differentiate_rat_poly(f);
    RatPoly f_sf = f;
    if (!fp.is_zero()) {
        auto [g, s, t] = ::cas::algebra::extended_gcd_rational_poly(f, fp);
        (void)s; (void)t;
        if (!g.is_zero() && g.size() > 1U) {
            auto [q, r] = ::cas::algebra::div_rem_rational_poly(f, g);
            if (r.is_zero() && !q.is_zero()) f_sf = std::move(q);
        }
    }
    if (f_sf.is_zero() || f_sf.size() <= 1U) {
        return ok(std::vector<double>{});
    }

    // 3. Build Sturm sequence and isolate.
    const std::vector<RatPoly> seq = sturm_sequence(f_sf);
    if (seq.size() < 2U) {
        return ok(std::vector<double>{});
    }
    std::vector<std::pair<Rational, Rational>> intervals;
    const Rational a = double_to_rational(low);
    const Rational b = double_to_rational(high);
    const Rational rtol = double_to_rational(tol);
    // Derived depth bound: ceil(log2((high - low) / tol)) + 4.
    const double width_ratio = (high - low) / std::max(tol, 1e-300);
    const unsigned int depth_bound = static_cast<unsigned int>(
        std::ceil(std::log2(std::max(width_ratio, 2.0)))) + 4U;
    isolate_recursive(seq, a, b, rtol, depth_bound, 0U, intervals);

    // 4. Newton polish each interval from its midpoint on the *full* f
    //    (not the squarefree part) to converge to a double-precision root.
    std::vector<double> roots;
    roots.reserve(intervals.size());
    const RatPoly f_full = std::move(f);
    const RatPoly fp_full = differentiate_rat_poly(f_full);
    // Newton iterations: O(log log (width/tol)). For tol=1e-12 and
    // initial width ~1, that's ~5; we allot 30 as a derived upper bound
    // that comfortably absorbs slow start near multiple roots.
    constexpr unsigned int kNewtonMaxIter = 30U;
    for (const auto& [ia, ib] : intervals) {
        const double mid = rational_to_double((ia + ib) * Rational(BigInt(1), BigInt(2)));
        const double r = newton_polish_rat(f_full, fp_full, mid, tol, kNewtonMaxIter);
        if (r >= low - tol && r <= high + tol) {
            roots.push_back(r);
        }
    }
    std::sort(roots.begin(), roots.end());
    return ok(std::move(roots));
}

// F8.0-5.4: rigorous isolating intervals. Reuses the same Sturm-sequence
// + squarefree pipeline as find_polynomial_roots_sturm, but returns the
// exact rational endpoints instead of the Newton-polished doubles.
Result<std::vector<IsolatingBound>> find_polynomial_isolating_intervals(
    ExprPtr expr,
    const std::string& variable,
    symbolic::CASContext& ctx,
    double low,
    double high,
    double tol) {

    if (!expr) {
        return fail<std::vector<IsolatingBound>>(make_error(
            CASErrorKind::InvalidArgument,
            "find_polynomial_isolating_intervals: null expr"));
    }
    if (low >= high) {
        return ok(std::vector<IsolatingBound>{});
    }

    const Symbol var(variable);
    auto poly_res = ::cas::algebra::parse_polynomial(expr, var, ctx);
    if (poly_res.is_error())
        return fail<std::vector<IsolatingBound>>(poly_res.error());
    auto rat_res = ::cas::algebra::poly_to_rational_poly(poly_res.value());
    if (rat_res.is_error())
        return fail<std::vector<IsolatingBound>>(rat_res.error());
    RatPoly f = std::move(rat_res.value());
    if (f.is_zero() || f.size() <= 1U) {
        return ok(std::vector<IsolatingBound>{});
    }

    // Squarefree part: f / gcd(f, f').
    RatPoly fp = differentiate_rat_poly(f);
    RatPoly f_sf = f;
    if (!fp.is_zero()) {
        auto [g, s, t] = ::cas::algebra::extended_gcd_rational_poly(f, fp);
        (void)s; (void)t;
        if (!g.is_zero() && g.size() > 1U) {
            auto [q, r] = ::cas::algebra::div_rem_rational_poly(f, g);
            if (r.is_zero() && !q.is_zero()) f_sf = std::move(q);
        }
    }
    if (f_sf.is_zero() || f_sf.size() <= 1U) {
        return ok(std::vector<IsolatingBound>{});
    }

    const std::vector<RatPoly> seq = sturm_sequence(f_sf);
    if (seq.size() < 2U) {
        return ok(std::vector<IsolatingBound>{});
    }
    std::vector<std::pair<Rational, Rational>> intervals;
    const Rational a = double_to_rational(low);
    const Rational b = double_to_rational(high);
    const Rational rtol = double_to_rational(tol);
    const double width_ratio = (high - low) / std::max(tol, 1e-300);
    const unsigned int depth_bound = static_cast<unsigned int>(
        std::ceil(std::log2(std::max(width_ratio, 2.0)))) + 4U;
    isolate_recursive(seq, a, b, rtol, depth_bound, 0U, intervals);

    // Sort by midpoint to mirror the ascending-root convention of the
    // double-precision variant.
    std::sort(intervals.begin(), intervals.end(),
        [](const auto& l, const auto& r) {
            return (l.first + l.second) < (r.first + r.second);
        });

    std::vector<IsolatingBound> bounds;
    bounds.reserve(intervals.size());
    for (const auto& [ia, ib] : intervals) {
        bounds.push_back(IsolatingBound{
            ia.numerator(), ia.denominator(),
            ib.numerator(), ib.denominator()});
    }
    return ok(std::move(bounds));
}

// ─── F8.0-5.3: BigFloat Newton polish ─────────────────────────────────────────
namespace {

[[nodiscard]] BigFloat bigfloat_from_rational(const Rational& r, mpfr_prec_t p) {
    return BigFloat::from_rational_parts(
        r.numerator().decimal(),
        r.denominator().decimal(),
        p);
}

// Horner-style evaluation of a Rational-coefficient polynomial at a BigFloat.
[[nodiscard]] BigFloat eval_ratpoly_bigfloat(const RatPoly& f,
                                             const BigFloat& x,
                                             mpfr_prec_t prec) {
    const auto& coeffs = f.coefficients();
    if (coeffs.empty()) return BigFloat::from_double(0.0, prec);
    BigFloat acc = bigfloat_from_rational(coeffs.back(), prec);
    for (auto it = coeffs.rbegin() + 1; it != coeffs.rend(); ++it) {
        acc = acc * x + bigfloat_from_rational(*it, prec);
    }
    return acc;
}

// Newton polish in BigFloat. Returns the refined root or `initial` if the
// iteration stalls (derivative ~0). Convergence is quadratic.
[[nodiscard]] BigFloat newton_polish_bigfloat(const RatPoly& f,
                                              const RatPoly& fp,
                                              BigFloat initial,
                                              const BigFloat& tol,
                                              unsigned int max_iter,
                                              mpfr_prec_t prec) {
    BigFloat x = std::move(initial);
    BigFloat tiny = BigFloat::from_double(1.0e-300, prec);
    for (unsigned int i = 0; i < max_iter; ++i) {
        BigFloat fx  = eval_ratpoly_bigfloat(f,  x, prec);
        if (BigFloat::abs(fx) < tol) return x;
        BigFloat fpx = eval_ratpoly_bigfloat(fp, x, prec);
        if (BigFloat::abs(fpx) < tiny) return x;
        BigFloat dx = fx / fpx;
        x = x - dx;
        if (BigFloat::abs(dx) < tol) return x;
    }
    return x;
}

} // namespace

Result<std::vector<BigFloat>> find_polynomial_roots_sturm_bigfloat(
    ExprPtr expr,
    const std::string& variable,
    symbolic::CASContext& ctx,
    double low,
    double high,
    double tol,
    mpfr_prec_t precision_bits) {

    if (!expr) {
        return fail<std::vector<BigFloat>>(make_error(
            CASErrorKind::InvalidArgument,
            "find_polynomial_roots_sturm_bigfloat: null expr"));
    }
    if (low >= high) {
        return ok(std::vector<BigFloat>{});
    }
    if (precision_bits < BigFloat::MIN_PREC) {
        return fail<std::vector<BigFloat>>(make_error(
            CASErrorKind::InvalidArgument,
            "find_polynomial_roots_sturm_bigfloat: precision_bits below MPFR minimum"));
    }

    const Symbol var(variable);
    auto poly_res = ::cas::algebra::parse_polynomial(expr, var, ctx);
    if (poly_res.is_error()) return fail<std::vector<BigFloat>>(poly_res.error());
    auto rat_res = ::cas::algebra::poly_to_rational_poly(poly_res.value());
    if (rat_res.is_error()) return fail<std::vector<BigFloat>>(rat_res.error());
    RatPoly f = std::move(rat_res.value());
    if (f.is_zero() || f.size() <= 1U) {
        return ok(std::vector<BigFloat>{});
    }

    // Squarefree part for Sturm isolation.
    RatPoly fp = differentiate_rat_poly(f);
    RatPoly f_sf = f;
    if (!fp.is_zero()) {
        auto [g, s, t] = ::cas::algebra::extended_gcd_rational_poly(f, fp);
        (void)s; (void)t;
        if (!g.is_zero() && g.size() > 1U) {
            auto [q, r] = ::cas::algebra::div_rem_rational_poly(f, g);
            if (r.is_zero() && !q.is_zero()) f_sf = std::move(q);
        }
    }
    if (f_sf.is_zero() || f_sf.size() <= 1U) {
        return ok(std::vector<BigFloat>{});
    }

    const std::vector<RatPoly> seq = sturm_sequence(f_sf);
    if (seq.size() < 2U) {
        return ok(std::vector<BigFloat>{});
    }
    std::vector<std::pair<Rational, Rational>> intervals;
    const Rational a = double_to_rational(low);
    const Rational b = double_to_rational(high);
    const Rational rtol = double_to_rational(tol);
    const double width_ratio = (high - low) / std::max(tol, 1e-300);
    const unsigned int depth_bound = static_cast<unsigned int>(
        std::ceil(std::log2(std::max(width_ratio, 2.0)))) + 4U;
    isolate_recursive(seq, a, b, rtol, depth_bound, 0U, intervals);

    // BigFloat Newton polish on each isolated interval midpoint.
    std::vector<BigFloat> roots;
    roots.reserve(intervals.size());
    const RatPoly f_full = std::move(f);
    const RatPoly fp_full = differentiate_rat_poly(f_full);
    const BigFloat half  = BigFloat::from_double(0.5, precision_bits);
    const BigFloat tol_bf = BigFloat::from_double(tol, precision_bits);
    constexpr unsigned int kNewtonMaxIter = 64U; // BigFloat tolerates more

    for (const auto& [ia, ib] : intervals) {
        BigFloat ia_bf = bigfloat_from_rational(ia, precision_bits);
        BigFloat ib_bf = bigfloat_from_rational(ib, precision_bits);
        BigFloat mid = (ia_bf + ib_bf) * half;
        BigFloat root = newton_polish_bigfloat(
            f_full, fp_full, std::move(mid), tol_bf, kNewtonMaxIter, precision_bits);
        roots.push_back(std::move(root));
    }
    return ok(std::move(roots));
}

} // namespace cas::numeric
