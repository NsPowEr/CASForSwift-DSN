// A6 Brick 3c — one certified Stauduhar step (the resolvent test proper).
// Mathematics and certificates: galois_stauduhar_internal.hpp. The
// first-layer driver (ambient, Tschirnhaus retry, parity) lives in
// galois_stauduhar_descent.cpp.

#include "galois_stauduhar_internal.hpp"

#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "perm_group_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cas::algebra::galois_stauduhar {

namespace {

using galois_invariant::RelativeInvariant;
using galois_padic::PadicRing;
using galois_padic::PadicSplitting;
using galois_padic::RingElem;
using permgrp::BsgsGroup;
using permgrp::compose;
using permgrp::inverse;
using permgrp::Perm;
using primitive_internal::Deadline;
using primitive_internal::deadline_exceeded;

[[nodiscard]] Result<void> poll(symbolic::CASContext* ctx,
                                const Deadline& dl) {
    if (ctx) {
        if (auto chk = ctx->check_interrupt(); chk.is_error()) return chk;
    }
    if (deadline_exceeded(dl)) {
        return fail<void>(
            CASError{.kind = CASErrorKind::Unimplemented,
                     .message = "stauduhar_step: deadline exceeded"});
    }
    return ok();
}

// Archimedean bound B ≥ |embedding of any coset value|: every root of the
// monic f satisfies |α| ≤ ρ = 1 + max|f_j| (Cauchy), every monomial of F
// has total degree = inv.total_degree, so |F^σ(α)| ≤ #terms·ρ^deg.
[[nodiscard]] BigInt coset_value_bound(const IntPoly& f,
                                       const RelativeInvariant& inv) {
    BigInt maxc(0);
    for (const auto& c : f.coefficients()) {
        BigInt a = c;
        if (a.is_negative()) a = -a;
        if (a > maxc) maxc = a;
    }
    const BigInt rho = maxc + BigInt(1);
    return BigInt(static_cast<std::int64_t>(inv.monomials.size())) *
           bigint_pow_nonnegative(rho, inv.total_degree);
}

// Smallest safe precision k with p^k > 2·max(|R(c)|, |R′(c)|) for any
// integer candidate c with |c| ≤ B: |R(c)| ≤ ∏(|c|+|v_i|) ≤ (2B)^m and
// |R′(c)| ≤ m·(2B)^{m−1}. Bit-length over-approximation — a larger k is
// always sound.
[[nodiscard]] std::size_t required_precision_k(const BigInt& p,
                                               const BigInt& B,
                                               std::size_t m) {
    const std::size_t bits = m * (B + B + BigInt(1)).bit_length() +
                             (BigInt(static_cast<std::int64_t>(m)) +
                              BigInt(1)).bit_length() +
                             2U;
    const std::size_t denom =
        p.bit_length() > 1U ? p.bit_length() - 1U : 1U;
    return bits / denom + 1U;
}

// Proven separation cap: distinct coset values differ by a nonzero
// algebraic integer whose conjugates are bounded by 2B and whose degree
// is at most m², so |Norm| ≤ (2B)^{m²}; residues at precision beyond
// this cap cannot keep colliding unless the values are equal.
[[nodiscard]] std::size_t separation_cap_k(const BigInt& p, const BigInt& B,
                                           std::size_t m) {
    const std::size_t bits =
        m * m * (B + B + BigInt(1)).bit_length() + 2U;
    const std::size_t denom =
        p.bit_length() > 1U ? p.bit_length() - 1U : 1U;
    return bits / denom + 1U;
}

// v = (σ·F)(r_1..r_n) using a table of root powers.
[[nodiscard]] RingElem eval_coset_value(
    const PadicRing& R, const std::vector<std::vector<RingElem>>& powtab,
    const Perm& sigma, const RelativeInvariant& inv) {
    RingElem acc = R.zero();
    for (const auto& mon : inv.monomials) {
        RingElem term = R.one();
        for (std::size_t i = 0U; i < mon.size(); ++i) {
            if (mon[i] == 0U) continue;
            term = R.mul(term, powtab[sigma[i]][mon[i]]);
        }
        acc = R.add(acc, term);
    }
    return acc;
}

}  // namespace

Result<StepOutcome> stauduhar_step(const BsgsGroup& G, const BsgsGroup& H,
                                   const RelativeInvariant& inv,
                                   const PadicSplitting& split,
                                   symbolic::CASContext* ctx,
                                   const Deadline& deadline) {
    const std::size_t n = split.f.degree();
    if (G.degree() != n || H.degree() != n ||
        inv.coset_reps.size() < 2U) {
        return fail<StepOutcome>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "stauduhar_step: degree mismatch or "
                                "trivial transversal"});
    }
    // Descent invariant tripwire: Frobenius ∈ G_f ≤ G.
    if (!G.contains(split.frobenius)) {
        return fail<StepOutcome>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "stauduhar_step: Frobenius outside the current "
                       "group — descent invariant G_f <= G violated"});
    }
    const std::size_t m = inv.coset_reps.size();
    const BigInt B = coset_value_bound(split.f, inv);
    // Two-phase precision: screening only needs p^k > 2B (so that a TRUE
    // integer value lifts exactly) — impurity or an out-of-bound lift at
    // ANY precision already certifies "not an integer", and most
    // candidates die here. The full bound k_cert (covering |R(c)|,
    // |R′(c)|) is requested only when a pure candidate survives.
    const std::size_t k_screen =
        required_precision_k(split.ring.prime(), B, 1U);
    const std::size_t k_cert =
        required_precision_k(split.ring.prime(), B, m);
    if (split.ring.precision() < k_screen) {
        StepOutcome o;
        o.status = StepStatus::NeedPrecision;
        o.required_precision = k_screen;
        return ok(std::move(o));
    }
    const PadicRing& R = split.ring;

    // Root-power table up to the maximal exponent appearing in F.
    std::size_t max_e = 0U;
    for (const auto& mon : inv.monomials) {
        for (const auto e : mon) max_e = std::max<std::size_t>(max_e, e);
    }
    std::vector<std::vector<RingElem>> powtab(n);
    for (std::size_t i = 0U; i < n; ++i) {
        powtab[i].resize(max_e + 1U, R.one());
        for (std::size_t e = 1U; e <= max_e; ++e) {
            powtab[i][e] = R.mul(powtab[i][e - 1U], split.roots[i]);
        }
    }

    // Coset values.
    std::vector<RingElem> vals;
    vals.reserve(m);
    for (const auto& sigma : inv.coset_reps) {
        if (auto pr = poll(ctx, deadline); pr.is_error()) {
            return fail<StepOutcome>(pr.error());
        }
        vals.push_back(eval_coset_value(R, powtab, sigma, inv));
    }

    // Any integer root of R(y) = ∏(y − v_i) equals some v_i, whose residue
    // is then pure — so only pure residues need testing, and R is never
    // expanded: R(c) and R′(c) are evaluated directly from the values in
    // O(m) ring operations (prefix/suffix products for the derivative).
    // Both are integers whenever G_f ≤ G (G-invariance), so purity and the
    // archimedean bound act as tripwires; |c| ≤ B for a true value, so
    // larger candidates are certified non-roots without any product.
    bool any_pure = false;
    for (std::size_t i = 0U; i < m && !any_pure; ++i) {
        const auto c = R.integer_residue(vals[i]);
        if (!c.has_value()) continue;
        BigInt cmag = *c;
        if (cmag.is_negative()) cmag = -cmag;
        any_pure = !(cmag > B);
    }
    if (!any_pure) {
        // No value can be an integer root at ANY precision: certified.
        StepOutcome o;
        o.status = StepStatus::NotContained;
        return ok(std::move(o));
    }
    if (split.ring.precision() < k_cert) {
        StepOutcome o;
        o.status = StepStatus::NeedPrecision;
        o.required_precision = k_cert;
        return ok(std::move(o));
    }
    bool multiple_integer_root = false;
    std::optional<std::size_t> descend_at;
    for (std::size_t i = 0U; i < m && !descend_at; ++i) {
        const auto c = R.integer_residue(vals[i]);
        if (!c.has_value()) continue;
        BigInt cmag = *c;
        if (cmag.is_negative()) cmag = -cmag;
        if (cmag > B) continue;  // exceeds the bound on every |v_j|
        if (auto pr = poll(ctx, deadline); pr.is_error()) {
            return fail<StepOutcome>(pr.error());
        }
        // R(c) and R′(c) via prefix/suffix products of (c − v_j).
        const RingElem ce = R.from_int(*c);
        std::vector<RingElem> pre(m + 1U, R.one());
        for (std::size_t j = 0U; j < m; ++j) {
            pre[j + 1U] = R.mul(pre[j], R.sub(ce, vals[j]));
        }
        RingElem suf = R.one();
        RingElem deriv = R.zero();
        for (std::size_t j = m; j-- > 0U;) {
            deriv = R.add(deriv, R.mul(pre[j], suf));
            suf = R.mul(suf, R.sub(ce, vals[j]));
        }
        const RingElem rc = pre[m];
        const auto rc_int = R.integer_residue(rc);
        const auto rd_int = R.integer_residue(deriv);
        if (!rc_int.has_value() || !rd_int.has_value()) {
            return fail<StepOutcome>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "stauduhar_step: impure resolvent evaluation — "
                           "G_f <= G violated upstream"});
        }
        if (!rc_int->is_zero()) continue;  // not a root of R
        if (rd_int->is_zero()) {
            // Multiple integer root (equal coset VALUES are perfectly
            // possible for distinct cosets): the criterion is
            // inconclusive for this model — a Tschirnhaus transform must
            // separate it.
            multiple_integer_root = true;
            continue;
        }
        // Simple integer root c: exactly ONE value equals c. Another
        // value sharing the residue is therefore a p-adic coincidence
        // that more precision provably separates (past the norm cap a
        // persistent collision would force equal values — contradicting
        // simplicity).
        bool ambiguous = false;
        for (std::size_t j = 0U; j < m; ++j) {
            if (j != i && R.equal(vals[j], vals[i])) {
                ambiguous = true;
                break;
            }
        }
        if (ambiguous) {
            const std::size_t cap = separation_cap_k(R.prime(), B, m);
            if (split.ring.precision() >= cap) {
                return fail<StepOutcome>(CASError{
                    .kind = CASErrorKind::InternalError,
                    .message = "stauduhar_step: residue collision beyond "
                               "the separation cap at a SIMPLE root — "
                               "impossible"});
            }
            StepOutcome o;
            o.status = StepStatus::NeedPrecision;
            o.required_precision =
                std::min(2U * split.ring.precision(), cap);
            return ok(std::move(o));
        }
        descend_at = i;
    }
    if (!descend_at) {
        StepOutcome o;
        o.status = multiple_integer_root ? StepStatus::Collision
                                         : StepStatus::NotContained;
        return ok(std::move(o));
    }
    {
        const std::size_t i = *descend_at;
        // Certified: v_i ∈ Z is a SIMPLE root of R ⇒ G_f ≤ σ_i·H·σ_i⁻¹.
        const Perm& sigma = inv.coset_reps[i];
        const Perm sigma_inv = inverse(sigma);
        std::vector<Perm> conj_gens;
        conj_gens.reserve(H.generators().size());
        for (const auto& h : H.generators()) {
            conj_gens.push_back(compose(compose(sigma, h), sigma_inv));
        }
        auto conj = BsgsGroup::build(n, std::move(conj_gens));
        if (conj.is_error()) return fail<StepOutcome>(conj.error());
        if (conj.value().order() != H.order() ||
            !conj.value().contains(split.frobenius)) {
            return fail<StepOutcome>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "stauduhar_step: conjugated subgroup failed "
                           "the order/Frobenius certificate"});
        }
        StepOutcome o;
        o.status = StepStatus::Descended;
        o.conjugator = sigma;
        o.conjugated_subgroup = std::move(conj.value());
        return ok(std::move(o));
    }
}

}  // namespace cas::algebra::galois_stauduhar
