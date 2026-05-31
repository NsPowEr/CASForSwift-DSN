// algebraic_tower_primitive.cpp — F3.4 Primitive Element Theorem (Trager 1976).
//
// Public API: compute_primitive_element, detect_tower_n_level.
// Internal helpers (shift-resultant, ring-GCD, back-substitution) are in
// algebraic_tower_primitive_internal.hpp.
//
// Algorithm (Cohen §3.6.2; GCL §8.7):
//
//   Base: θ₁ = α₁,  K₁ = Q(θ₁),  q₁(y) = min-poly of α₁.
//
//   Step k ∈ {2, …, n}:
//     For s = 1, 2, 3, … (bound from ctx.max_trager_tower_shift_attempts):
//       1. SHIFT-RESULTANT: R_s(y) = Res_x(m_k(x), q_{k-1}(y − s·x)).
//          Via evaluation-interpolation + Newton interpolation.
//       2. SQUAREFREE TEST: gcd(R_s, R_s') == constant  ⟺  R_s squarefree.
//       3. If squarefree: θ_k = θ_{k-1} + s·α_k,  q_k = R_s (made monic).
//       4. EXPRESS GENERATORS via ring-GCD in Q[y]/(q_k)[t]:
//            h(t) = q_{k-1}(y − s·t), gcd(m_k(t), h(t)) = (t − α_k).
//            → α_k extracted; θ_{k-1} = y − s·α_k back-expressed.
//       5. BACK-SUBSTITUTE all previous α_i via composition modulo q_k.
//
// CONSTRAINTS (CLAUDE.md):
//   - Shift budget: ctx.max_trager_tower_shift_attempts() > 0 → use that;
//     0 → Trager bound = 2·∏deg(m_i)+1. Cap → Unimplemented with diagnostics.
//   - All arithmetic: RatPoly/Rational/BigInt. No double, no int64.

#include "cas/algebraic_tower_bridge.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include "polynomial_internal.hpp"
#include "algebra_internal.hpp"
#include "algebraic_tower_primitive_internal.hpp"
#include "algebraic_tower_primitive_nested.hpp"

#include <chrono>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cas {
namespace algebra {

using namespace primitive_internal;

// ── compute_primitive_element ────────────────────────────────────────────────

Result<PrimitiveElementResult> compute_primitive_element(
    const std::vector<ExprPtr>& alphas,
    const std::vector<AlgebraicNumber::CoeffVec>& min_polys,
    symbolic::CASContext& ctx) {

    if (alphas.empty() || min_polys.empty()) {
        return fail<PrimitiveElementResult>(CASError{
            CASErrorKind::InvalidArgument,
            "compute_primitive_element: empty generator list",
            std::nullopt});
    }
    if (alphas.size() != min_polys.size()) {
        return fail<PrimitiveElementResult>(CASError{
            CASErrorKind::InvalidArgument,
            "compute_primitive_element: alphas and min_polys must have the same size",
            std::nullopt});
    }
    for (std::size_t i = 0U; i < min_polys.size(); ++i) {
        if (min_polys[i].size() < 2U) {
            return fail<PrimitiveElementResult>(CASError{
                CASErrorKind::InvalidArgument,
                "compute_primitive_element: min_poly[" + std::to_string(i) +
                "] must have degree ≥ 1 (size ≥ 2)",
                std::nullopt});
        }
    }

    const std::size_t n = alphas.size();

    // Trivial case: n == 1 → Q(α₁) is already a simple extension.
    if (n == 1U) {
        std::vector<Rational> monic_mp = vec_make_monic(min_polys[0]);
        const std::size_t deg = monic_mp.size() - 1U;
        // α₁ = θ itself → [0, 1, 0, ..., 0].
        std::vector<Rational> alpha_in_theta(deg, Rational(BigInt(0)));
        if (deg >= 1U) alpha_in_theta[1] = Rational(BigInt(1));
        PrimitiveElementResult res{
            .theta_expr = alphas[0],
            .shifts = {},
            .min_poly_theta = std::move(monic_mp),
            .alphas_in_theta = {std::move(alpha_in_theta)},
        };
        return ok(std::move(res));
    }

    // Determine shift budget.
    const std::size_t user_bound = ctx.max_trager_tower_shift_attempts();
    const std::size_t bound = (user_bound > 0U)
        ? user_bound
        : trager_primitive_bound(min_polys);

    // F3.5-DEBT-01 fix: wall-clock deadline for redundant-generator detection.
    const primitive_internal::Deadline deadline =
        std::chrono::steady_clock::now() + ctx.timeout();

    // Initialize from α₁.
    RatPoly q_current = vec_to_ratpoly(vec_make_monic(min_polys[0]));
    ExprPtr theta_expr = alphas[0];
    std::vector<BigInt> shifts_so_far;

    // alphas_in_current[i] = poly for α_i in Q[y]/(q_current).
    const std::size_t init_deg = q_current.degree();
    std::vector<std::vector<Rational>> alphas_in_current;
    // α₁ = θ itself = [0, 1, 0, ..., 0] (the monomial y).
    std::vector<Rational> alpha1_poly(init_deg, Rational(BigInt(0)));
    if (init_deg >= 2U) alpha1_poly[1] = Rational(BigInt(1));
    alphas_in_current.push_back(std::move(alpha1_poly));

    // Merge generators k = 1..n-1.
    for (std::size_t k = 1U; k < n; ++k) {
        const RatPoly m_k = vec_to_ratpoly(vec_make_monic(min_polys[k]));
        const std::size_t expected_deg = q_current.degree() * m_k.degree();

        bool found = false;
        std::size_t shift_attempts = 0U;

        for (std::size_t s_val = 1U; s_val <= bound && !found; ++s_val) {
            if (primitive_internal::deadline_exceeded(deadline)) {
                return fail<PrimitiveElementResult>(CASError{
                    CASErrorKind::Unimplemented,
                    "compute_primitive_element: ctx.timeout() exceeded during "
                    "shift search at merge step k=" + std::to_string(k + 1U) +
                    " (tried " + std::to_string(shift_attempts) + " shifts). "
                    "Generators may be algebraically dependent (e.g. Q(√2,√3,√6) "
                    "with √6=√2·√3). Increase ctx.set_timeout() or remove "
                    "redundant generators.",
                    std::nullopt});
            }
            ++shift_attempts;
            const BigInt s(static_cast<std::int64_t>(s_val));

            // 1. SHIFT-RESULTANT.
            auto R_res = compute_shift_resultant(q_current, m_k, s, deadline);
            if (R_res.is_error()) continue;
            RatPoly R_s = std::move(R_res.value());
            normalize_rational_coefficients(R_s);
            if (R_s.degree() != expected_deg) continue;

            // 2. SQUAREFREE TEST.
            // Fast modular pre-filter: when R_s is NON-squarefree (typical of
            // algebraically dependent generators such as Q(√2,√3,√6) with
            // √6=√2·√3), the full Q-gcd inside ratpoly_is_squarefree explodes
            // in coefficient size.  A reduction mod p (any prime that doesn't
            // divide the leading coefficient denominator) reveals the repeated
            // factor in milliseconds.  We try a few primes; if ALL of them
            // report squarefree, we proceed to the full check (which now is
            // expected to succeed cheaply because the poly is genuinely
            // squarefree).
            {
                const std::uint64_t test_primes[] = {1000003U, 1000033U, 1000037U};
                bool any_non_sqfree = false;
                for (std::uint64_t pp : test_primes) {
                    if (primitive_internal::ratpoly_definitely_non_squarefree_mod_p(R_s, pp)) {
                        any_non_sqfree = true;
                        break;
                    }
                }
                if (any_non_sqfree) continue;
            }
            if (!ratpoly_is_squarefree(R_s)) continue;
            R_s = ratpoly_make_monic(std::move(R_s));
            normalize_rational_coefficients(R_s);

            // F3.5-DEBT-01: factor R_s over Q; use irreducible factors as q_current
            // candidates so Q[y]/(cand_q) is a field and ring-GCD is valid.
            const auto all_q_candidates = collect_irred_factors_over_q(R_s, ctx);

            // 3. Build ExprPtr for θ_k = θ_{k-1} + s · α_k.
            ExprPtr s_expr = ctx.arena().make<IntegerLit>(s);
            ExprPtr s_alpha_k = ctx.arena().make<Binary>(BinaryOp::Mul, s_expr, alphas[k]);
            ExprPtr new_theta_expr = ctx.arena().make<Binary>(BinaryOp::Add, theta_expr, s_alpha_k);
            {
                auto simplified = ctx.simplify(new_theta_expr);
                if (simplified.is_ok()) new_theta_expr = simplified.value();
            }

            // F3.5-DEBT-01 semantic-consistency filter: see select_candidates_by_theta_expr.
            std::vector<RatPoly> q_candidates =
                select_candidates_by_theta_expr(all_q_candidates, new_theta_expr, ctx);

            // 4. For each irreducible Q-factor of R_s, try to extract α_k via
            //    ring-GCD of m_k(t) and h(t) = q_current(y−s·t) in Q[y]/(cand_q)[t].
            for (const RatPoly& cand_q : q_candidates) {
                const std::size_t new_deg = cand_q.degree();
                if (new_deg == 0U) continue;
            const std::size_t dq = q_current.degree();
            const Rational s_rat(s);
            const Rational neg_s_rat = Rational(BigInt(0)) - s_rat;

            // Precompute binomials and powers for expansion of q_current(y−s·t).
            std::vector<std::vector<BigInt>> binom_tbl(dq + 1U,
                std::vector<BigInt>(dq + 1U, BigInt(0)));
            for (std::size_t row = 0U; row <= dq; ++row) {
                binom_tbl[row][0] = BigInt(1);
                for (std::size_t col = 1U; col <= row; ++col) {
                    binom_tbl[row][col] = binom_tbl[row-1U][col-1U] + binom_tbl[row-1U][col];
                }
            }

            // Precompute y^i mod R_s for i = 0..dq.
            std::vector<std::vector<Rational>> y_pow_mod(dq + 1U,
                std::vector<Rational>(new_deg, Rational(BigInt(0))));
            y_pow_mod[0][0] = Rational(BigInt(1));
            if (dq >= 1U && new_deg >= 2U) {
                y_pow_mod[1][1] = Rational(BigInt(1));
                for (std::size_t pp = 2U; pp <= dq; ++pp) {
                    RatPoly prev_p = vec_to_ratpoly(y_pow_mod[pp - 1U]);
                    std::vector<Rational> shifted_coeffs(prev_p.size() + 1U, Rational(BigInt(0)));
                    for (std::size_t d2 = 0U; d2 < prev_p.size(); ++d2) {
                        shifted_coeffs[d2 + 1U] = prev_p[d2];
                    }
                    RatPoly shifted = vec_to_ratpoly(shifted_coeffs);
                    RatPoly reduced = ratpoly_mod(shifted, cand_q);
                    for (std::size_t d2 = 0U; d2 < new_deg; ++d2) {
                        y_pow_mod[pp][d2] = (d2 < reduced.size())
                            ? reduced[d2] : Rational(BigInt(0));
                    }
                }
            }

            // Precompute (-s)^p for p = 0..dq.
            std::vector<Rational> neg_s_pows_r(dq + 1U, Rational(BigInt(0)));
            neg_s_pows_r[0] = Rational(BigInt(1));
            for (std::size_t pp = 1U; pp <= dq; ++pp)
                neg_s_pows_r[pp] = neg_s_pows_r[pp-1U] * neg_s_rat;

            // Build coeffs_in_t[j] = coeff of t^j in h(t), as vec[new_deg].
            std::vector<std::vector<Rational>> coeffs_in_t(dq + 1U,
                std::vector<Rational>(new_deg, Rational(BigInt(0))));
            for (std::size_t i = 0U; i <= dq; ++i) {
                const Rational& qi = q_current[i];
                if (qi.numerator().is_zero()) continue;
                for (std::size_t j = 0U; j <= i; ++j) {
                    const Rational contrib = qi
                        * Rational(binom_tbl[i][j])
                        * neg_s_pows_r[j];
                    if (contrib.numerator().is_zero()) continue;
                    for (std::size_t d2 = 0U; d2 < new_deg; ++d2) {
                        coeffs_in_t[j][d2] = coeffs_in_t[j][d2]
                            + contrib * y_pow_mod[i - j][d2];
                    }
                }
            }

            // Lift m_k to Q[y]/(cand_q) coefficients (constants).
            std::vector<std::vector<Rational>> mk_in_ring(m_k.size(),
                std::vector<Rational>(new_deg, Rational(BigInt(0))));
            for (std::size_t i = 0U; i < m_k.size(); ++i)
                mk_in_ring[i][0] = m_k[i];

            // ── Ring-GCD in Q[y]/(cand_q)[t] ─────────────────────────────────
            using RingElem = std::vector<Rational>;
            using RingPoly = std::vector<RingElem>;

            auto ring_zero_elem = [&]() -> RingElem {
                return RingElem(new_deg, Rational(BigInt(0)));
            };
            auto ring_is_zero_elem = [&](const RingElem& v) -> bool {
                for (const Rational& r : v)
                    if (!r.numerator().is_zero()) return false;
                return true;
            };
            auto ring_add_elem = [&](const RingElem& a, const RingElem& b) -> RingElem {
                RingElem r(new_deg, Rational(BigInt(0)));
                for (std::size_t d2 = 0U; d2 < new_deg; ++d2)
                    r[d2] = a[d2] + b[d2];
                return r;
            };
            auto ring_sub_elem = [&](const RingElem& a, const RingElem& b) -> RingElem {
                RingElem r(new_deg, Rational(BigInt(0)));
                for (std::size_t d2 = 0U; d2 < new_deg; ++d2)
                    r[d2] = a[d2] - b[d2];
                return r;
            };
            auto ring_mul_elem = [&](const RingElem& a_e, const RingElem& b_e) -> RingElem {
                RatPoly pa = vec_to_ratpoly(a_e);
                RatPoly pb = vec_to_ratpoly(b_e);
                RatPoly prod = ratpoly_mulmod(pa, pb, cand_q);
                RingElem r(new_deg, Rational(BigInt(0)));
                for (std::size_t d2 = 0U; d2 < prod.size() && d2 < new_deg; ++d2)
                    r[d2] = prod[d2];
                return r;
            };
            // Extended Euclidean in Q[y]/(cand_q). cand_q is irreducible here,
            // so Q[y]/(cand_q) is a field and the inverse always exists for non-zero elements.
            auto ring_inv_elem = [&](const RingElem& a_e) -> std::optional<RingElem> {
                RatPoly pa = vec_to_ratpoly(a_e);
                if (pa.is_zero()) return std::nullopt;
                RatPoly r0 = pa;
                RatPoly r1 = cand_q;
                RatPoly s0 = vec_to_ratpoly({Rational(BigInt(1))});
                RatPoly s1 = RatPoly{};
                const std::size_t max_iter = (R_s.degree() + pa.degree()) * 2U + 16U;
                std::size_t iters = 0U;
                // Bit-length cap: when R_s is reducible, Q[y]/(R_s) is not a
                // field and Euclidean step coefficients grow without bound
                // even when iteration count is small.  Cap at 8192 bits per
                // coefficient (= ~2466 decimal digits) — far beyond any
                // legitimate inverse computation in the irreducible (field)
                // case, but tight enough to catch the degenerate path fast.
                const std::size_t max_bits = 8192U;
                auto max_coef_bits = [](const RatPoly& g) -> std::size_t {
                    std::size_t mb = 0U;
                    for (std::size_t i = 0; i < g.size(); ++i) {
                        const Rational& r = g[i];
                        const std::size_t bn = r.numerator().bit_length();
                        const std::size_t bd = r.denominator().bit_length();
                        if (bn > mb) mb = bn;
                        if (bd > mb) mb = bd;
                    }
                    return mb;
                };
                while (!r1.is_zero()) {
                    if (iters++ >= max_iter) return std::nullopt;
                    if (primitive_internal::deadline_exceeded(deadline)) return std::nullopt;
                    if (max_coef_bits(r0) > max_bits || max_coef_bits(r1) > max_bits) {
                        return std::nullopt;
                    }
                    auto [q, r2] = div_rem_rational_poly(r0, r1);
                    RatPoly s2 = sub_rational_poly(s0, mul_rational_poly(q, s1));
                    r0 = std::move(r1);
                    s0 = std::move(s1);
                    r1 = std::move(r2);
                    s1 = std::move(s2);
                }
                if (r0.is_zero() || r0.degree() > 0U) return std::nullopt;
                const Rational g0 = r0[0];
                if (g0.numerator().is_zero()) return std::nullopt;
                const Rational g0_inv{g0.denominator(), g0.numerator()};
                RatPoly result_p = ratpoly_scale(s0, g0_inv);
                result_p = ratpoly_mod(result_p, cand_q);
                RingElem r(new_deg, Rational(BigInt(0)));
                for (std::size_t d2 = 0U; d2 < result_p.size() && d2 < new_deg; ++d2)
                    r[d2] = result_p[d2];
                return r;
            };

            auto rpoly_degree_t = [&](const RingPoly& p) -> std::size_t {
                for (std::size_t i = p.size(); i > 0U; --i)
                    if (!ring_is_zero_elem(p[i - 1U])) return i - 1U;
                return 0U;
            };
            auto rpoly_is_zero_t = [&](const RingPoly& p) -> bool {
                for (const auto& c : p)
                    if (!ring_is_zero_elem(c)) return false;
                return true;
            };

            // Euclidean division in Q[y]/(cand_q)[t]; bit-cap as numeric safety belt.
            auto ring_elem_bits = [](const RingElem& e) -> std::size_t {
                std::size_t mb = 0U;
                for (const Rational& r : e) {
                    const std::size_t bn = r.numerator().bit_length();
                    const std::size_t bd = r.denominator().bit_length();
                    if (bn > mb) mb = bn;
                    if (bd > mb) mb = bd;
                }
                return mb;
            };
            const std::size_t kMaxRingBits = 8192U;

            auto rpoly_divmod_t = [&](RingPoly A, const RingPoly& B)
                -> std::optional<std::pair<RingPoly, RingPoly>> {
                const std::size_t degB = rpoly_degree_t(B);
                auto inv_lc = ring_inv_elem(B[degB]);
                if (!inv_lc) return std::nullopt;
                if (ring_elem_bits(*inv_lc) > kMaxRingBits) return std::nullopt;
                RingPoly Q_poly;
                const std::size_t max_div_steps =
                    (rpoly_degree_t(A) + degB + 4U) * 4U + 16U;
                std::size_t div_steps = 0U;
                while (!rpoly_is_zero_t(A) && rpoly_degree_t(A) >= degB) {
                    if (div_steps++ >= max_div_steps) return std::nullopt;
                    if (primitive_internal::deadline_exceeded(deadline)) return std::nullopt;
                    const std::size_t degA = rpoly_degree_t(A);
                    if (ring_elem_bits(A[degA]) > kMaxRingBits) return std::nullopt;
                    const std::size_t exp_diff = degA - degB;
                    auto q_coeff = ring_mul_elem(A[degA], *inv_lc);
                    if (Q_poly.size() < exp_diff + 1U)
                        Q_poly.resize(exp_diff + 1U, ring_zero_elem());
                    Q_poly[exp_diff] = ring_add_elem(Q_poly[exp_diff], q_coeff);
                    for (std::size_t jj = 0U; jj <= degB; ++jj)
                        A[exp_diff + jj] = ring_sub_elem(A[exp_diff + jj],
                                                         ring_mul_elem(q_coeff, B[jj]));
                    while (!A.empty() && ring_is_zero_elem(A.back())) A.pop_back();
                }
                return std::make_pair(std::move(Q_poly), std::move(A));
            };

            // Euclidean GCD over field[t]. Honors deadline (F3.5-DEBT-01).
            auto rpoly_gcd_t = [&](RingPoly A, RingPoly B) -> std::optional<RingPoly> {
                const std::size_t max_iter = (dq + m_k.degree()) * 4U + 16U;
                std::size_t iters = 0U;
                while (!rpoly_is_zero_t(B) && iters < max_iter) {
                    if (primitive_internal::deadline_exceeded(deadline)) return std::nullopt;
                    ++iters;
                    auto dm = rpoly_divmod_t(A, B);
                    if (!dm) return std::nullopt;
                    A = std::move(B);
                    B = std::move(dm->second);
                }
                return A;
            };

            auto gcd_rpoly = rpoly_gcd_t(coeffs_in_t, mk_in_ring);
            if (!gcd_rpoly) continue;
            if (rpoly_degree_t(*gcd_rpoly) != 1U) continue;
            if (gcd_rpoly->size() < 2U) continue;

            auto inv_c1 = ring_inv_elem((*gcd_rpoly)[1]);
            if (!inv_c1) continue;

            // α_k = -gcd[0] / gcd[1] in Q[y]/(cand_q).
            RingElem neg_c0 = ring_sub_elem(ring_zero_elem(), (*gcd_rpoly)[0]);
            RingElem alpha_k_poly = ring_mul_elem(neg_c0, *inv_c1);
            {
                RatPoly tmp_poly = vec_to_ratpoly(alpha_k_poly);
                tmp_poly = ratpoly_mod(tmp_poly, cand_q);
                alpha_k_poly.assign(new_deg, Rational(BigInt(0)));
                for (std::size_t d2 = 0U; d2 < tmp_poly.size() && d2 < new_deg; ++d2)
                    alpha_k_poly[d2] = tmp_poly[d2];
            }

            // θ_{k-1} = y − s · α_k in Q[y]/(cand_q).
            std::vector<Rational> y_poly(new_deg, Rational(BigInt(0)));
            if (new_deg >= 2U) y_poly[1] = Rational(BigInt(1));
            std::vector<Rational> theta_prev_poly(new_deg, Rational(BigInt(0)));
            for (std::size_t d2 = 0U; d2 < new_deg; ++d2)
                theta_prev_poly[d2] = y_poly[d2] - s_rat * alpha_k_poly[d2];

            // 5. BACK-SUBSTITUTE all previous α_i.
            std::vector<std::vector<Rational>> new_alphas_in_current;
            new_alphas_in_current.reserve(k + 1U);
            for (std::size_t i = 0U; i < k; ++i) {
                new_alphas_in_current.push_back(
                    reexpress_in_new_theta(alphas_in_current[i], theta_prev_poly, cand_q));
            }
            new_alphas_in_current.push_back(std::move(alpha_k_poly));

            // Commit.
            q_current = cand_q;
            theta_expr = new_theta_expr;
            shifts_so_far.push_back(s);
            alphas_in_current = std::move(new_alphas_in_current);
            found = true;
            break;
            }  // end candidate loop
        }  // shift loop

        if (!found) {
            return fail<PrimitiveElementResult>(CASError{
                CASErrorKind::Unimplemented,
                "compute_primitive_element: primitive element search exhausted "
                "ctx.max_trager_tower_shift_attempts() = " +
                std::to_string(bound) + " without finding a squarefree resultant "
                "at merge step k=" + std::to_string(k + 1U) +
                " (tried shifts 1.." + std::to_string(shift_attempts) + "). "
                "Increase ctx.set_max_trager_tower_shift_attempts() or verify "
                "that the extension is separable.",
                std::nullopt});
        }
    }  // generator loop

    PrimitiveElementResult res;
    res.theta_expr = theta_expr;
    res.shifts = std::move(shifts_so_far);
    res.min_poly_theta = vec_make_monic(q_current.coefficients());
    res.alphas_in_theta = std::move(alphas_in_current);
    return ok(std::move(res));
}

// ── detect_tower_n_level moved to algebraic_tower_primitive_nested.cpp ──────
// (was inline here; relocated to satisfy 500-LOC anti-monolith limit and
// because it now depends on the F3.4-DEBT-01 nested-lift helper.)

}  // namespace algebra
}  // namespace cas
