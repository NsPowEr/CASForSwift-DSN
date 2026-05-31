// algebraic_tower_primitive_internal.hpp — internal helpers for F3.4.
// Shared between algebraic_tower_primitive.cpp (public API) only.
// All functions are in the anonymous namespace via the including translation unit.
// Include only from algebraic_tower_primitive.cpp.
#pragma once

#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/algebraic_number.hpp"
#include "cas/ast.hpp"
#include "cas/algebra.hpp"
#include "cas/bigint.hpp"
#include "cas/symbolic.hpp"

#include "polynomial_internal.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace cas {
namespace algebra {
namespace primitive_internal {

// ── RatPoly / CoeffVec conversion ────────────────────────────────────────────

// Convert std::vector<Rational> (CoeffVec) to RatPoly.
[[nodiscard]] inline RatPoly vec_to_ratpoly(const std::vector<Rational>& v) {
    RatPoly p(v);
    p.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return p;
}

// Make a std::vector<Rational> monic (divide all by leading coeff).
[[nodiscard]] inline std::vector<Rational> vec_make_monic(std::vector<Rational> v) {
    while (!v.empty() && v.back().numerator().is_zero()) v.pop_back();
    if (v.empty()) return v;
    const Rational lc = v.back();
    if (lc == Rational(BigInt(1))) return v;
    const Rational inv{lc.denominator(), lc.numerator()};
    for (Rational& c : v) c = c * inv;
    while (!v.empty() && v.back().numerator().is_zero()) v.pop_back();
    return v;
}

// ── Low-level RatPoly helpers ─────────────────────────────────────────────────

// Differentiate a RatPoly (ascending coefficients).
[[nodiscard]] inline RatPoly ratpoly_deriv(const RatPoly& f) {
    if (f.size() <= 1U) return RatPoly{};
    RatPoly d;
    d.reserve(f.size() - 1U);
    for (std::size_t i = 1U; i < f.size(); ++i) {
        d.push_back(f[i] * Rational(BigInt(static_cast<std::int64_t>(i))));
    }
    d.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return d;
}

// Scale all coefficients of f by scalar k.
[[nodiscard]] inline RatPoly ratpoly_scale(const RatPoly& f, const Rational& k) {
    if (k.numerator().is_zero()) return RatPoly{};
    RatPoly r;
    r.reserve(f.size());
    for (const Rational& c : f.coefficients()) r.push_back(c * k);
    r.normalize([](const Rational& r2) { return r2.numerator().is_zero(); });
    return r;
}

// Make a RatPoly monic.
[[nodiscard]] inline RatPoly ratpoly_make_monic(RatPoly f) {
    if (f.is_zero()) return f;
    const Rational lc = f.leading_coeff();
    if (lc.numerator().is_zero()) return f;
    if (lc == Rational(BigInt(1))) return f;
    const Rational inv{lc.denominator(), lc.numerator()};
    return ratpoly_scale(f, inv);
}

// Compute gcd of two RatPoly.
[[nodiscard]] inline RatPoly ratpoly_gcd(const RatPoly& a, const RatPoly& b) {
    auto [g, s, t] = extended_gcd_rational_poly(a, b);
    return ratpoly_make_monic(g);
}

// Modular fast-path: test if f is NON-squarefree by reducing mod p and
// checking deg(gcd(f mod p, f' mod p)) > 0 in Fp[y].  Returns:
//   - true  → f is non-squarefree (mod p witnesses a repeated factor over Q);
//             repeated factor in Fp[y] lifts uniquely to Q only when p is
//             "good" (does not divide leading coeffs or discriminants), but a
//             positive answer is reliable up to extremely rare unlucky primes.
//   - false → mod p says squarefree.  Could be a lucky reduction; caller may
//             still verify via full gcd, but for the redundant-generator
//             detection path the modular answer is sufficient (we'll reject
//             non-squarefree, accept squarefree-after-modular and verify
//             posteriori via the algebraic identity).
[[nodiscard]] inline bool ratpoly_definitely_non_squarefree_mod_p(
    const RatPoly& f, std::uint64_t p) {
    if (f.degree() <= 1U) return false;
    // Build dense f mod p and f' mod p (as int64 coefficients in [0, p)).
    auto reduce = [&](const RatPoly& g) -> std::vector<std::int64_t> {
        std::vector<std::int64_t> v(g.size(), 0);
        for (std::size_t i = 0; i < g.size(); ++i) {
            const Rational& r = g[i];
            const BigInt& num = r.numerator();
            const BigInt& den = r.denominator();
            // Reduce numerator mod p: take low limb if fits, else fall back to
            // string-based modular reduction.
            BigInt num_mod = num % BigInt(static_cast<std::int64_t>(p));
            BigInt den_mod = den % BigInt(static_cast<std::int64_t>(p));
            if (den_mod.is_zero()) return std::vector<std::int64_t>{};  // unlucky prime
            // num / den mod p = num * den^{-1} mod p.
            std::int64_t n = static_cast<std::int64_t>(num_mod.to_u64());
            if (num_mod.is_negative()) n = (n + static_cast<std::int64_t>(p)) % static_cast<std::int64_t>(p);
            std::int64_t d = static_cast<std::int64_t>(den_mod.to_u64());
            if (den_mod.is_negative()) d = (d + static_cast<std::int64_t>(p)) % static_cast<std::int64_t>(p);
            // Compute d^{-1} mod p via Fermat (p prime, gcd(d,p)=1 unless d≡0).
            if (d == 0) return std::vector<std::int64_t>{};
            std::int64_t d_inv = 1;
            std::int64_t base = d % static_cast<std::int64_t>(p);
            std::uint64_t e = static_cast<std::uint64_t>(p - 2U);
            while (e > 0U) {
                if (e & 1U) d_inv = static_cast<std::int64_t>((static_cast<__int128>(d_inv) * base) % p);
                base = static_cast<std::int64_t>((static_cast<__int128>(base) * base) % p);
                e >>= 1U;
            }
            v[i] = static_cast<std::int64_t>((static_cast<__int128>(n) * d_inv) % p);
        }
        while (!v.empty() && v.back() == 0) v.pop_back();
        return v;
    };
    auto fp = reduce(f);
    if (fp.empty()) return false;  // unlucky prime
    if (fp.size() <= 1U) return false;
    // Derivative mod p.
    std::vector<std::int64_t> dfp(fp.size() - 1U, 0);
    for (std::size_t i = 1U; i < fp.size(); ++i) {
        dfp[i - 1U] = static_cast<std::int64_t>(
            (static_cast<__int128>(fp[i]) * static_cast<std::int64_t>(i)) % p);
    }
    while (!dfp.empty() && dfp.back() == 0) dfp.pop_back();
    if (dfp.empty()) return true;  // f' ≡ 0 mod p → non-squarefree (or char issue)
    // Euclidean gcd in Fp.
    auto a = fp;
    auto b = dfp;
    auto fp_modinv = [p](std::int64_t x) -> std::int64_t {
        std::int64_t r = 1;
        std::int64_t base = ((x % static_cast<std::int64_t>(p)) + static_cast<std::int64_t>(p)) % static_cast<std::int64_t>(p);
        std::uint64_t e = p - 2U;
        while (e > 0U) {
            if (e & 1U) r = static_cast<std::int64_t>((static_cast<__int128>(r) * base) % p);
            base = static_cast<std::int64_t>((static_cast<__int128>(base) * base) % p);
            e >>= 1U;
        }
        return r;
    };
    while (!b.empty()) {
        std::int64_t inv_lc = fp_modinv(b.back());
        while (a.size() >= b.size()) {
            std::int64_t factor = static_cast<std::int64_t>(
                (static_cast<__int128>(a.back()) * inv_lc) % p);
            std::size_t shift = a.size() - b.size();
            for (std::size_t i = 0; i < b.size(); ++i) {
                std::int64_t prod = static_cast<std::int64_t>(
                    (static_cast<__int128>(factor) * b[i]) % p);
                a[i + shift] = ((a[i + shift] - prod) % static_cast<std::int64_t>(p)
                    + static_cast<std::int64_t>(p)) % static_cast<std::int64_t>(p);
            }
            while (!a.empty() && a.back() == 0) a.pop_back();
        }
        std::swap(a, b);
    }
    return a.size() > 1U;  // deg(gcd) > 0 → non-squarefree mod p.
}

// Test if f ∈ Q[y] is squarefree: gcd(f, f') has degree 0.
[[nodiscard]] inline bool ratpoly_is_squarefree(const RatPoly& f) {
    if (f.degree() == 0U) return true;
    RatPoly df = ratpoly_deriv(f);
    if (df.is_zero()) return false;
    RatPoly g = ratpoly_gcd(f, df);
    return (g.degree() == 0U);
}

// Reduce p mod m in Q[y].
[[nodiscard]] inline RatPoly ratpoly_mod(const RatPoly& p, const RatPoly& m) {
    auto [q, r] = div_rem_rational_poly(p, m);
    return r;
}

// Multiply two RatPoly mod m.
[[nodiscard]] inline RatPoly ratpoly_mulmod(const RatPoly& a, const RatPoly& b, const RatPoly& m) {
    return ratpoly_mod(mul_rational_poly(a, b), m);
}

// ── Shift-resultant: R_s(y) = Res_x(m_k(x), q_prev(y − s·x)) ───────────────
//
// Uses evaluation-interpolation + Newton interpolation.
using Deadline = std::optional<std::chrono::steady_clock::time_point>;

[[nodiscard]] inline bool deadline_exceeded(const Deadline& dl) {
    return dl.has_value() && std::chrono::steady_clock::now() >= *dl;
}

[[nodiscard]] inline Result<RatPoly> compute_shift_resultant(
    const RatPoly& q_prev,
    const RatPoly& m_k,
    const BigInt& s,
    const Deadline& deadline = std::nullopt) {

    if (q_prev.is_zero() || m_k.is_zero()) {
        return fail<RatPoly>(CASError{
            CASErrorKind::InvalidArgument,
            "compute_shift_resultant: zero input polynomial",
            std::nullopt});
    }

    const std::size_t deg_R = q_prev.degree() * m_k.degree();
    const std::size_t num_pts = deg_R + 1U;

    std::vector<Rational> eval_pts;
    eval_pts.reserve(num_pts);
    for (std::size_t j = 1U; j <= num_pts; ++j) {
        eval_pts.push_back(Rational(BigInt(static_cast<std::int64_t>(j))));
    }

    std::vector<Rational> values;
    values.reserve(num_pts);

    const Rational s_rat(s);
    const Rational neg_s_rat = Rational(BigInt(0)) - s_rat;
    const std::size_t dq = q_prev.degree();

    // Precompute binomial coefficients C(dq, k).
    std::vector<std::vector<BigInt>> binom_table(dq + 1U, std::vector<BigInt>(dq + 1U, BigInt(0)));
    for (std::size_t row = 0U; row <= dq; ++row) {
        binom_table[row][0] = BigInt(1);
        for (std::size_t col = 1U; col <= row; ++col) {
            binom_table[row][col] = binom_table[row - 1U][col - 1U] + binom_table[row - 1U][col];
        }
    }

    for (const Rational& y_j : eval_pts) {
        if (deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                CASErrorKind::Unimplemented,
                "compute_shift_resultant: ctx.timeout() exceeded during "
                "evaluation-interpolation",
                std::nullopt});
        }
        // Build Q_j(x) = q_prev(y_j − s·x).
        std::vector<Rational> Q_j_coeffs(q_prev.size(), Rational(BigInt(0)));

        std::vector<Rational> y_pows(q_prev.size());
        std::vector<Rational> neg_s_pows(q_prev.size());
        y_pows[0] = Rational(BigInt(1));
        neg_s_pows[0] = Rational(BigInt(1));
        for (std::size_t p = 1U; p < q_prev.size(); ++p) {
            y_pows[p] = y_pows[p - 1U] * y_j;
            neg_s_pows[p] = neg_s_pows[p - 1U] * neg_s_rat;
        }

        for (std::size_t i = 0U; i < q_prev.size(); ++i) {
            const Rational& ci = q_prev[i];
            if (ci.numerator().is_zero()) continue;
            for (std::size_t k = 0U; k <= i; ++k) {
                const Rational term = ci
                    * Rational(binom_table[i][k])
                    * y_pows[i - k]
                    * neg_s_pows[k];
                Q_j_coeffs[k] = Q_j_coeffs[k] + term;
            }
        }
        RatPoly Q_j(Q_j_coeffs);
        Q_j.normalize([](const Rational& r) { return r.numerator().is_zero(); });

        auto res_scalar = resultant_generic<Rational>(
            m_k.coefficients(),
            Q_j.coefficients(),
            nullptr,
            deadline);
        if (res_scalar.is_error()) return fail<RatPoly>(res_scalar.error());
        values.push_back(res_scalar.value());
    }

    // Newton interpolation: divided differences.
    std::vector<std::vector<Rational>> table(num_pts, std::vector<Rational>(num_pts, Rational(BigInt(0))));
    for (std::size_t i = 0U; i < num_pts; ++i) table[i][0] = values[i];
    for (std::size_t j = 1U; j < num_pts; ++j) {
        if (deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                CASErrorKind::Unimplemented,
                "compute_shift_resultant: ctx.timeout() exceeded during "
                "Newton interpolation",
                std::nullopt});
        }
        for (std::size_t i = j; i < num_pts; ++i) {
            const Rational denom = eval_pts[i] - eval_pts[i - j];
            if (denom.numerator().is_zero()) {
                return fail<RatPoly>(CASError{
                    CASErrorKind::InternalError,
                    "compute_shift_resultant: duplicate evaluation points",
                    std::nullopt});
            }
            const Rational denom_inv{denom.denominator(), denom.numerator()};
            table[i][j] = (table[i][j - 1U] - table[i - 1U][j - 1U]) * denom_inv;
        }
    }

    // Expand Newton form into monomial basis.
    std::vector<Rational> result_coeffs(num_pts, Rational(BigInt(0)));
    std::vector<Rational> newton_poly = {Rational(BigInt(1))};

    for (std::size_t j = 0U; j < num_pts; ++j) {
        const Rational coeff_j = table[j][j];
        if (!coeff_j.numerator().is_zero()) {
            for (std::size_t d = 0U; d < newton_poly.size(); ++d) {
                result_coeffs[d] = result_coeffs[d] + coeff_j * newton_poly[d];
            }
        }
        if (j + 1U < num_pts) {
            std::vector<Rational> new_newton(newton_poly.size() + 1U, Rational(BigInt(0)));
            for (std::size_t d = 0U; d < newton_poly.size(); ++d) {
                new_newton[d + 1U] = new_newton[d + 1U] + newton_poly[d];
                new_newton[d] = new_newton[d] - eval_pts[j] * newton_poly[d];
            }
            while (!new_newton.empty() && new_newton.back().numerator().is_zero())
                new_newton.pop_back();
            newton_poly = std::move(new_newton);
        }
    }

    RatPoly result(result_coeffs);
    result.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return ok(std::move(result));
}

// ── Back-substitute previous α_i via polynomial composition ──────────────────
//
// Compose P_i(T(θ_k)) mod q_k where T = prev_theta_in_theta_k.
// Uses Horner's method.
[[nodiscard]] inline std::vector<Rational> reexpress_in_new_theta(
    const std::vector<Rational>& alpha_i_in_prev_theta,
    const std::vector<Rational>& prev_theta_in_theta_k,
    const RatPoly& q_k) {
    const std::size_t deg_qk = q_k.degree();

    const std::size_t d = alpha_i_in_prev_theta.size();
    if (d == 0U) return std::vector<Rational>(deg_qk, Rational(BigInt(0)));

    RatPoly T_poly = vec_to_ratpoly(prev_theta_in_theta_k);

    RatPoly acc;
    acc.push_back(alpha_i_in_prev_theta[d - 1U]);

    for (std::size_t idx = d - 1U; idx > 0U; --idx) {
        acc = ratpoly_mulmod(acc, T_poly, q_k);
        const Rational& ak = alpha_i_in_prev_theta[idx - 1U];
        if (!ak.numerator().is_zero()) {
            if (acc.is_zero()) {
                acc.push_back(ak);
            } else {
                acc[0] = acc[0] + ak;
                acc.normalize([](const Rational& r) { return r.numerator().is_zero(); });
            }
        }
    }
    acc = ratpoly_mod(acc, q_k);

    std::vector<Rational> result(deg_qk, Rational(BigInt(0)));
    for (std::size_t i = 0U; i < acc.size() && i < deg_qk; ++i)
        result[i] = acc[i];
    return result;
}

// Compute the Trager discriminant collision bound = 2·∏deg(m_i)+1.
[[nodiscard]] inline std::size_t trager_primitive_bound(
    const std::vector<AlgebraicNumber::CoeffVec>& min_polys) {
    std::size_t product = 1U;
    for (const auto& mp : min_polys) {
        const std::size_t d = (mp.size() >= 2U) ? (mp.size() - 1U) : 1U;
        product *= d;
    }
    return 2U * product + 1U;
}

// ── F3.4-DEBT-01 fix: absolute resultant for nested RootOf ───────────────────
//
// Given:
//   g(y)        ∈ Q[y]                — rational min-poly of inner β
//   bivariate   = vector<RatPoly>      — f(x, y) as Σ_i (poly_y_i(y)) · x^i,
//                                        where each poly_y_i ∈ Q[y].
// Compute R(x) = Res_y( g(y), f(x, y) ) ∈ Q[x].
//
// Reference: Cohen "A Course in Computational Algebraic Number Theory"
// §3.6.1 (absolute minimal polynomial via resultant). Trager 1976.
//
// Algorithm (evaluation-interpolation in x):
//   D = deg_y(g) · deg_x(f)
//   For j = 1..D+1:
//     f_y(y) = Σ_i (x_j)^i · poly_y_i(y) ∈ Q[y]
//     v_j    = Res_y(g(y), f_y(y))      (scalar via resultant_generic<Rational>)
//   Newton-interpolate {v_j} → R(x) ∈ Q[x].
//
// Honors deadline (F3.5-DEBT-01).
[[nodiscard]] inline Result<RatPoly> compute_absolute_resultant_xy(
    const RatPoly& g_y,
    const std::vector<RatPoly>& f_xy,
    const Deadline& deadline = std::nullopt) {

    if (g_y.is_zero()) {
        return fail<RatPoly>(CASError{
            CASErrorKind::InvalidArgument,
            "compute_absolute_resultant_xy: zero inner minimal polynomial",
            std::nullopt});
    }
    if (f_xy.empty()) {
        return fail<RatPoly>(CASError{
            CASErrorKind::InvalidArgument,
            "compute_absolute_resultant_xy: empty bivariate polynomial",
            std::nullopt});
    }
    const std::size_t deg_x = f_xy.size() - 1U;
    const std::size_t deg_y_g = g_y.degree();
    const std::size_t deg_R = deg_x * deg_y_g;
    const std::size_t num_pts = deg_R + 1U;

    std::vector<Rational> eval_pts;
    eval_pts.reserve(num_pts);
    for (std::size_t j = 1U; j <= num_pts; ++j) {
        eval_pts.push_back(Rational(BigInt(static_cast<std::int64_t>(j))));
    }

    std::vector<Rational> values;
    values.reserve(num_pts);

    for (const Rational& x_j : eval_pts) {
        if (deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                CASErrorKind::Unimplemented,
                "compute_absolute_resultant_xy: ctx.timeout() exceeded during "
                "evaluation-interpolation",
                std::nullopt});
        }
        // Build f_y(y) = Σ_i (x_j)^i · f_xy[i](y) ∈ Q[y].
        // Find max degree across f_xy[i].
        std::size_t max_deg_y = 0U;
        for (const auto& c_y : f_xy) {
            if (!c_y.is_zero() && c_y.degree() > max_deg_y) max_deg_y = c_y.degree();
        }
        std::vector<Rational> f_y_coeffs(max_deg_y + 1U, Rational(BigInt(0)));

        Rational x_pow = Rational(BigInt(1));
        for (std::size_t i = 0U; i <= deg_x; ++i) {
            const RatPoly& c_y = f_xy[i];
            for (std::size_t d2 = 0U; d2 < c_y.size(); ++d2) {
                f_y_coeffs[d2] = f_y_coeffs[d2] + x_pow * c_y[d2];
            }
            x_pow = x_pow * x_j;
        }
        RatPoly f_y(f_y_coeffs);
        f_y.normalize([](const Rational& r) { return r.numerator().is_zero(); });
        if (f_y.is_zero()) {
            // f(x_j, y) ≡ 0 ⇒ Res = 0 (degenerate evaluation point).
            values.push_back(Rational(BigInt(0)));
            continue;
        }

        auto res_scalar = resultant_generic<Rational>(
            g_y.coefficients(),
            f_y.coefficients(),
            nullptr,
            deadline);
        if (res_scalar.is_error()) return fail<RatPoly>(res_scalar.error());
        values.push_back(res_scalar.value());
    }

    // Newton interpolation (same divided-differences scheme as compute_shift_resultant).
    std::vector<std::vector<Rational>> table(num_pts, std::vector<Rational>(num_pts, Rational(BigInt(0))));
    for (std::size_t i = 0U; i < num_pts; ++i) table[i][0] = values[i];
    for (std::size_t j = 1U; j < num_pts; ++j) {
        if (deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                CASErrorKind::Unimplemented,
                "compute_absolute_resultant_xy: ctx.timeout() exceeded during "
                "Newton interpolation",
                std::nullopt});
        }
        for (std::size_t i = j; i < num_pts; ++i) {
            const Rational denom = eval_pts[i] - eval_pts[i - j];
            if (denom.numerator().is_zero()) {
                return fail<RatPoly>(CASError{
                    CASErrorKind::InternalError,
                    "compute_absolute_resultant_xy: duplicate evaluation points",
                    std::nullopt});
            }
            const Rational denom_inv{denom.denominator(), denom.numerator()};
            table[i][j] = (table[i][j - 1U] - table[i - 1U][j - 1U]) * denom_inv;
        }
    }

    std::vector<Rational> result_coeffs(num_pts, Rational(BigInt(0)));
    std::vector<Rational> newton_poly = {Rational(BigInt(1))};
    for (std::size_t j = 0U; j < num_pts; ++j) {
        const Rational coeff_j = table[j][j];
        if (!coeff_j.numerator().is_zero()) {
            for (std::size_t d = 0U; d < newton_poly.size(); ++d) {
                result_coeffs[d] = result_coeffs[d] + coeff_j * newton_poly[d];
            }
        }
        if (j + 1U < num_pts) {
            std::vector<Rational> new_newton(newton_poly.size() + 1U, Rational(BigInt(0)));
            for (std::size_t d = 0U; d < newton_poly.size(); ++d) {
                new_newton[d + 1U] = new_newton[d + 1U] + newton_poly[d];
                new_newton[d] = new_newton[d] - eval_pts[j] * newton_poly[d];
            }
            while (!new_newton.empty() && new_newton.back().numerator().is_zero())
                new_newton.pop_back();
            newton_poly = std::move(new_newton);
        }
    }

    RatPoly result(result_coeffs);
    result.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return ok(std::move(result));
}

// ── F3.5-DEBT-01: irreducible Q-factorization of a resultant poly ────────────
//
// R_s from compute_shift_resultant is squarefree but may be REDUCIBLE over Q
// when generators are algebraically dependent (e.g. Q(√2,√3,√6) where √6=√2·√3).
// Q[y]/(reducible) is NOT a field → ring_inv produces wrong results.
// Solution: factor R_s over Q, return its irreducible factors as candidates for
// q_current.  The calling loop picks the factor for which the ring-GCD succeeds
// and m_k(α_k) ≡ 0 (guaranteed by construction in a genuine field Q[y]/(irred)).
//
// On any error from factor_over_integers, falls back to {f} (treat as irreducible)
// so the existing behaviour is preserved for the common case where R_s is already
// irreducible (no factorization call overhead beyond the single cheaply-verified
// squarefree check that already passed).
[[nodiscard]] inline std::vector<RatPoly>
collect_irred_factors_over_q(const RatPoly& f, symbolic::CASContext& ctx) {
    // Clear denominators: LCM of all coefficient denominators → Z[y].
    BigInt lcm_den(1);
    for (std::size_t i = 0U; i < f.size(); ++i) {
        const BigInt& d = f[i].denominator();
        BigInt g = gcd(lcm_den, d);
        lcm_den = (lcm_den / g) * d;
    }
    PolyExpr int_poly;
    int_poly.reserve(f.size());
    for (std::size_t i = 0U; i < f.size(); ++i) {
        Rational sc = f[i] * Rational(lcm_den);
        int_poly.push_back(sc.numerator().is_zero()
            ? ExprPtr{}
            : ctx.arena().make<IntegerLit>(sc.numerator()));
    }
    normalize_poly(int_poly);
    // Use a fresh non-user-facing symbol as the polynomial variable.
    Symbol pv = ctx.make_fresh_symbol("pv");
    auto e_res = polynomial_to_expr(int_poly, pv, ctx);
    if (e_res.is_error()) return {f};
    auto fr = factor_over_integers(e_res.value(), pv, ctx);
    if (fr.is_error() || fr.value().factors.size() <= 1U) return {f};
    std::vector<RatPoly> out;
    out.reserve(fr.value().factors.size());
    for (const auto& pf : fr.value().factors) {
        auto pp = parse_polynomial(pf.factor, pv, ctx);
        if (pp.is_error()) return {f};
        auto rp = poly_to_rational_poly(pp.value());
        if (rp.is_error()) return {f};
        out.push_back(ratpoly_make_monic(std::move(rp.value())));
    }
    return out;
}

// ── F3.5-DEBT-01: symbolic-consistency filter on candidate min-polys ─────────
//
// When R_s is squarefree but REDUCIBLE on Q, collect_irred_factors_over_q
// returns multiple irreducible factors.  Several may yield a valid algebraic
// ring-GCD inside compute_primitive_element, but only the unique factor
// vanishing at the SYMBOLIC value of new_theta_expr = α_k + s·θ_{k-1}
// corresponds to the user's intended embedding.  Picking the wrong candidate
// causes Q[y]/(cand_q) ↦ Q(θ_expr) projection mismatch: the algebraic α_k
// extracted in the wrong ring, when rendered via y → theta_expr, yields a
// Galois conjugate of α_k instead of α_k itself (e.g. Q(√2,√8) with
// cand_q=y²-2 instead of y²-18: rendered α_2 = 6√2 ≠ 2√2 = √8).
//
// Filter: evaluate cand_q(new_theta_expr) symbolically via Horner; ctx.simplify
// must reduce to literal 0 (IntegerLit/RationalLit with numerator zero).
// When the simplifier cannot conclusively reduce (e.g. nested radicals it
// does not normalise), accept ALL candidates as fallback to preserve
// termination — the algebraic outcome still satisfies C1/C2 invariants; only
// the symbolic rendering may pick a conjugate embedding.
[[nodiscard]] inline bool cand_vanishes_at_theta_expr(
    const RatPoly& cand,
    ExprPtr theta_expr,
    symbolic::CASContext& ctx) {
    if (cand.size() == 0U) return false;
    ExprPtr val = ctx.arena().make<RationalLit>(
        cand[cand.size() - 1U].numerator(),
        cand[cand.size() - 1U].denominator());
    for (std::size_t i = cand.size() - 1U; i > 0U; --i) {
        ExprPtr mul = ctx.arena().make<Binary>(BinaryOp::Mul, val, theta_expr);
        ExprPtr ci = ctx.arena().make<RationalLit>(
            cand[i - 1U].numerator(),
            cand[i - 1U].denominator());
        val = ctx.arena().make<Binary>(BinaryOp::Add, mul, ci);
    }
    auto simp = ctx.simplify(val);
    if (simp.is_error()) return false;
    ExprPtr sv = simp.value();
    if (auto* il = expr_cast<IntegerLit>(sv); il != nullptr) {
        return il->value.is_zero();
    }
    if (auto* rl = expr_cast<RationalLit>(sv); rl != nullptr) {
        return rl->numerator.is_zero();
    }
    return false;
}

[[nodiscard]] inline std::vector<RatPoly> select_candidates_by_theta_expr(
    const std::vector<RatPoly>& all_candidates,
    ExprPtr theta_expr,
    symbolic::CASContext& ctx) {
    std::vector<RatPoly> filtered;
    filtered.reserve(all_candidates.size());
    for (const auto& cq : all_candidates) {
        if (cand_vanishes_at_theta_expr(cq, theta_expr, ctx)) filtered.push_back(cq);
    }
    // No candidate confirmed by symbolic check → fall back to original (any-valid)
    // order. Preserves termination on cases where the simplifier cannot prove
    // vanishing (e.g. deeply nested radicals).
    if (filtered.empty()) return all_candidates;
    return filtered;
}

}  // namespace primitive_internal
}  // namespace algebra
}  // namespace cas
