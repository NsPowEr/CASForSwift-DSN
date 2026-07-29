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
#include "algebraic_tower_resultant.hpp"

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
