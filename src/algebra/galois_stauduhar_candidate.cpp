// A6 Brick 3.5 — the per-candidate retry protocol of the Stauduhar
// descent, shared by the first layer (galois_stauduhar_descent.cpp) and
// the below-first-layer walk (galois_stauduhar_below.cpp): identity model
// first, precision raises on demand (always applied to the ORIGINAL
// splitting — the only one Newton can safely refine, since a Tschirnhaus
// model g may have p | disc(g)), and a Tschirnhaus GRID sweep when a
// certified multiple integer root makes the resolvent criterion
// inconclusive.
//
// Why a grid and not a one-parameter curve. A single-parameter family
// P_t(x) = Σ t^{m−1}x^m has c_a·c_b = t^{a+b−2}: products of coefficients
// collide whenever a+b agrees, and that algebraic degeneracy can zero a
// separating polynomial IDENTICALLY in t (measured on x⁵−2 with the
// quadratic pentagon invariant below F₂₀: the pair difference is
// proportional to (c₁c₄ − c₂c₃)·√5 ≡ 0 on the curve). The sound family is
// the full coefficient GRID c = (c₀, …, c_{n−1}) ∈ {0..Δ}ⁿ, β = c₀ +
// Σ_{m≥1} c_m·α^m:
//   • for every pair of distinct coset polynomials F^σ ≠ F^τ the
//     difference D_στ(c) = (F^σ − F^τ)(P_c(r)) is NOT the zero polynomial
//     in the free coefficients — Lagrange interpolation reaches any
//     target point β (n free coefficients, n roots), and two distinct
//     polynomials cannot agree on the Zariski-dense set of
//     distinct-coordinate points;
//   • the same holds for every root difference β_a − β_b (linear in c);
//   • the product of all these nonzero polynomials has per-variable
//     degree ≤ Δ = deg(F)·C(m,2) + C(n,2), so by the Combinatorial
//     Nullstellensatz (Alon) it has a nonzero point on the grid
//     {0, …, Δ}ⁿ: the sweep, enumerated by growing max-norm shells,
//     PROVABLY terminates — exhausting the grid is an InternalError, not
//     a retry-budget guess.

#include "galois_stauduhar_internal.hpp"

#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "galois_setpoly_internal.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cas::algebra::galois_stauduhar {

namespace {

using galois_padic::PadicSplitting;
using galois_padic::RingElem;
using permgrp::BsgsGroup;
using primitive_internal::Deadline;

[[nodiscard]] bool rat_is_zero(const Rational& r) {
    return r.numerator().is_zero();
}

[[nodiscard]] RatPoly to_ratpoly(const IntPoly& p) {
    std::vector<Rational> c;
    c.reserve(p.size());
    for (const auto& v : p.coefficients()) c.emplace_back(v);
    RatPoly r(std::move(c));
    r.normalize(rat_is_zero);
    return r;
}

// g(y − s) by exact binomial expansion (Taylor shift, O(n²) — n ≤ 10).
[[nodiscard]] IntPoly taylor_shift(const IntPoly& g, const BigInt& s) {
    const std::size_t d = g.degree();
    std::vector<BigInt> out(d + 1U, BigInt(0));
    for (std::size_t j = 0U; j <= d; ++j) {
        // g_j·(y − s)^j: binom accumulates C(j,k)·(−s)^{j−k}.
        BigInt binom(1);
        for (std::size_t k = j + 1U; k-- > 0U;) {
            out[k] = out[k] + g[j] * binom;
            if (k > 0U) {
                binom = binom * BigInt(static_cast<std::int64_t>(k)) /
                        BigInt(static_cast<std::int64_t>(j - k + 1U));
                binom = binom * (BigInt(0) - s);
            }
        }
    }
    IntPoly r(std::move(out));
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// The Tschirnhaus model for a coefficient vector c = (c_0, …, c_{n−1}):
// g = minpoly of β = c_0 + Σ_{m≥1} c_m·α^m, with the β_i evaluated on the
// SAME lifted roots (order preserved: σ(β_i) = P(σ(r_i)) = β_{frob(i)} is
// an algebraic identity). nullopt = degenerate c (collision among the β
// or non-squarefree g): the caller advances the sweep.
struct DerivedModel {
    IntPoly g;
    PadicSplitting split;
};
[[nodiscard]] Result<std::optional<DerivedModel>> derive_model(
    const IntPoly& f, const PadicSplitting& base,
    const std::vector<BigInt>& coeffs, const Deadline& dl) {
    const std::size_t n = f.degree();
    // Degrees 1..n−1 for tschirnhaus_general; the constant c_0 is applied
    // afterwards as the exact Taylor shift g(y − c_0).
    std::vector<BigInt> pc(coeffs.begin() + 1, coeffs.end());
    auto g_rat = galois_setpoly::tschirnhaus_general(to_ratpoly(f), pc, dl);
    if (g_rat.is_error()) {
        return fail<std::optional<DerivedModel>>(g_rat.error());
    }
    auto sf = galois_setpoly::is_squarefree_q(g_rat.value());
    if (sf.is_error()) return fail<std::optional<DerivedModel>>(sf.error());
    if (!sf.value() || g_rat.value().degree() != n) {
        return ok(std::optional<DerivedModel>{});
    }
    // Monic with integer coefficients is a theorem here (f monic ∈ Z[x],
    // P ∈ Z[x] ⇒ Res_x(f(x), y − P(x)) ∈ Z[y] monic).
    std::vector<BigInt> gc(g_rat.value().size());
    for (std::size_t i = 0U; i < g_rat.value().size(); ++i) {
        const Rational& c = g_rat.value()[i];
        if (!(c.denominator() == BigInt(1))) {
            return fail<std::optional<DerivedModel>>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "stauduhar: Tschirnhaus model not integral"});
        }
        gc[i] = c.numerator();
    }
    IntPoly g(std::move(gc));
    if (!coeffs[0].is_zero()) g = taylor_shift(g, coeffs[0]);
    const IntPoly P{std::vector<BigInt>(coeffs)};
    const auto& R = base.ring;
    std::vector<RingElem> beta;
    beta.reserve(base.roots.size());
    for (const auto& r : base.roots) {
        beta.push_back(R.eval_int_poly(P, r));
    }
    for (std::size_t i = 0U; i < beta.size(); ++i) {
        if (!R.is_zero(R.eval_int_poly(g, beta[i]))) {
            return fail<std::optional<DerivedModel>>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "stauduhar: g(P(root)) != 0 — Tschirnhaus "
                           "identity violated"});
        }
        for (std::size_t j = 0U; j < i; ++j) {
            if (R.equal(beta[i], beta[j])) {
                return ok(std::optional<DerivedModel>{});  // degenerate c
            }
        }
    }
    PadicSplitting split{base.ring, g, std::move(beta), base.frobenius};
    return ok(std::optional<DerivedModel>{
        DerivedModel{std::move(g), std::move(split)}});
}

// Odometer over the max-norm shell {c ∈ {0..radius}ⁿ : max c_i = radius},
// skipping vectors whose degree-≥1 part is all zero (constant P) and the
// identity transform (tested before the sweep starts). Returns false when
// the shell is exhausted.
[[nodiscard]] bool next_in_shell(std::vector<std::size_t>& c,
                                 std::size_t radius) {
    const std::size_t n = c.size();
    while (true) {
        std::size_t i = n;
        while (i-- > 0U) {
            if (c[i] < radius) {
                ++c[i];
                for (std::size_t j = i + 1U; j < n; ++j) c[j] = 0U;
                break;
            }
            if (i == 0U) return false;
        }
        std::size_t maxc = 0U;
        std::size_t tail = 0U;
        for (std::size_t j = 0U; j < n; ++j) {
            if (c[j] > maxc) maxc = c[j];
            if (j >= 1U) tail += c[j];
        }
        if (maxc != radius || tail == 0U) continue;
        const bool is_identity =
            c[0] == 0U && c[1] == 1U && tail == 1U;
        if (!is_identity) return true;
    }
}

}  // namespace

Result<std::optional<DescentHit>> test_candidate_with_retries(
    const BsgsGroup& current, const BsgsGroup& H,
    const galois_invariant::RelativeInvariant& inv, const IntPoly& f_monic,
    PadicSplitting& base, symbolic::CASContext& ctx, const Deadline& dl) {
    const std::size_t n = f_monic.degree();
    const std::size_t m = inv.coset_reps.size();
    // Grid side: per-variable degree of the product of all separating
    // polynomials (see file header) — deg(F)·C(m,2) for the coset-value
    // pairs plus C(n,2) for the root differences.
    const std::size_t delta =
        inv.total_degree * (m * (m - 1U) / 2U) + n * (n - 1U) / 2U;
    // Sweep state: radius-0 sentinel = the identity model; afterwards the
    // grid shells of growing max-norm.
    std::size_t radius = 0U;
    std::vector<std::size_t> cvec(n, 0U);
    bool shell_fresh = false;  // cvec valid and not yet consumed
    while (true) {
        if (auto chk = ctx.check_interrupt(); chk.is_error()) {
            return fail<std::optional<DescentHit>>(chk.error());
        }
        std::optional<PadicSplitting> model_split;
        if (radius == 0U) {
            model_split = base;  // identity transform
        } else {
            if (!shell_fresh) {
                // Advance within the shell; on exhaustion open the next.
                while (!next_in_shell(cvec, radius)) {
                    ++radius;
                    if (radius > delta) {
                        return fail<std::optional<DescentHit>>(CASError{
                            .kind = CASErrorKind::InternalError,
                            .message =
                                "stauduhar: Tschirnhaus grid exhausted — "
                                "impossible past the Combinatorial "
                                "Nullstellensatz bound"});
                    }
                    for (auto& v : cvec) v = 0U;
                }
                shell_fresh = true;
            }
            std::vector<BigInt> coeffs;
            coeffs.reserve(n);
            for (const std::size_t v : cvec) {
                coeffs.emplace_back(static_cast<std::int64_t>(v));
            }
            auto dm = derive_model(f_monic, base, coeffs, dl);
            if (dm.is_error()) {
                return fail<std::optional<DescentHit>>(dm.error());
            }
            if (!dm.value().has_value()) {
                shell_fresh = false;  // degenerate vector
                continue;
            }
            model_split = std::move(dm.value()->split);
        }
        auto out = stauduhar_step(current, H, inv, *model_split, &ctx, dl);
        if (out.is_error()) {
            return fail<std::optional<DescentHit>>(out.error());
        }
        switch (out.value().status) {
            case StepStatus::NeedPrecision: {
                if (out.value().required_precision <=
                    base.ring.precision()) {
                    return fail<std::optional<DescentHit>>(CASError{
                        .kind = CASErrorKind::InternalError,
                        .message = "stauduhar: precision request without "
                                   "progress"});
                }
                auto raised = galois_padic::raise_splitting_precision(
                    base, out.value().required_precision, &ctx, dl);
                if (raised.is_error()) {
                    return fail<std::optional<DescentHit>>(raised.error());
                }
                base = std::move(raised.value());
                break;  // retry the same model at the new precision
            }
            case StepStatus::Collision: {
                // Certified multiple integer root: only a Tschirnhaus
                // change of model can separate it.
                if (radius == 0U) {
                    radius = 1U;
                    for (auto& v : cvec) v = 0U;
                }
                shell_fresh = false;
                break;
            }
            case StepStatus::Descended: {
                return ok(std::optional<DescentHit>{DescentHit{
                    std::move(*out.value().conjugator),
                    std::move(*out.value().conjugated_subgroup)}});
            }
            case StepStatus::NotContained: {
                return ok(std::optional<DescentHit>{});
            }
        }
    }
}

}  // namespace cas::algebra::galois_stauduhar
