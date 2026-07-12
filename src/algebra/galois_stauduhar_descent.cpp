// A6 Brick 3c — first-layer Stauduhar descent driver: discriminant parity
// picks the ambient (S_n / A_n), Brick-2 supplies the maximal transitive
// candidates, Brick-3b the exact invariants/transversals, Brick-3a the
// certified p-adic splitting, galois_stauduhar_step.cpp the resolvent
// test. Collisions retry through a precision doubling and then a
// Tschirnhaus sweep on the moment curve whose length is DERIVED from the
// degree bound of the separating polynomials (never a magic constant);
// the transformed roots β = P_t(α) are computed on the SAME lifted roots,
// preserving the root order the descent state lives on.

#include "galois_stauduhar_internal.hpp"

#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "galois_internal.hpp"
#include "galois_setpoly_internal.hpp"
#include "perm_group_internal.hpp"
#include "perm_maximal_internal.hpp"
#include "polynomial_internal.hpp"
#include "polynomial_resultant_generic.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra::galois_stauduhar {

namespace {

using galois_padic::PadicSplitting;
using galois_padic::RingElem;
using permgrp::AmbientGroup;
using permgrp::BsgsGroup;
using permgrp::factorial_u64;
using permgrp::identity;
using permgrp::Perm;
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

// disc(f) = (−1)^{n(n−1)/2}·Res(f, f′) for monic f.
[[nodiscard]] Result<Rational> discriminant_q(const IntPoly& f,
                                              const Deadline& dl) {
    const RatPoly fr = to_ratpoly(f);
    std::vector<Rational> dr;
    for (std::size_t i = 1U; i < fr.size(); ++i) {
        dr.push_back(fr[i] * Rational(BigInt(static_cast<std::int64_t>(i))));
    }
    auto res = resultant_generic<Rational>(fr.coefficients(), dr, nullptr,
                                           dl);
    if (res.is_error()) return fail<Rational>(res.error());
    Rational d = res.value();
    const std::size_t n = f.degree();
    if ((n * (n - 1U) / 2U) % 2U == 1U) d = Rational(BigInt(0)) - d;
    return ok(std::move(d));
}

// The Tschirnhaus model at parameter t ≥ 1: g = minpoly of β = P_t(α)
// with P_t(x) = Σ_{m=1..n−1} t^{m−1}·x^m, and the β_i evaluated on the
// SAME lifted roots (order preserved). nullopt = degenerate t (collision
// among the β or non-squarefree g): the caller advances the sweep.
struct DerivedModel {
    IntPoly g;
    PadicSplitting split;
};
[[nodiscard]] Result<std::optional<DerivedModel>> derive_model(
    const IntPoly& f, const PadicSplitting& base, std::size_t t,
    const Deadline& dl) {
    const std::size_t n = f.degree();
    std::vector<BigInt> pc;
    pc.reserve(n - 1U);
    BigInt tp(1);
    const BigInt tb(static_cast<std::int64_t>(t));
    for (std::size_t mth = 1U; mth < n; ++mth) {
        pc.push_back(tp);
        tp = tp * tb;
    }
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
    // β_i = P_t(r_i) in the ring; order (and hence the Frobenius) is
    // preserved: σ(β_i) = P(σ(r_i)) = β_{frob(i)} is an algebraic identity.
    std::vector<BigInt> p_coeffs(pc.size() + 1U, BigInt(0));
    for (std::size_t i = 0U; i < pc.size(); ++i) p_coeffs[i + 1U] = pc[i];
    const IntPoly P(std::move(p_coeffs));
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
                return ok(std::optional<DerivedModel>{});  // degenerate t
            }
        }
    }
    PadicSplitting split{base.ring, g, std::move(beta), base.frobenius};
    return ok(std::optional<DerivedModel>{
        DerivedModel{std::move(g), std::move(split)}});
}

}  // namespace

Result<FirstLayerDescent> stauduhar_first_layer(const IntPoly& f_monic,
                                                symbolic::CASContext& ctx,
                                                const Deadline& deadline) {
    const std::size_t n = f_monic.degree();
    if (f_monic.is_zero() || !(f_monic.leading_coeff() == BigInt(1))) {
        return fail<FirstLayerDescent>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "stauduhar_first_layer: monic f required"});
    }
    if (n < 5U || n > 10U) {
        return fail<FirstLayerDescent>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "stauduhar_first_layer: maximal-candidate coverage "
                       "is certified for degrees 5..10 only (Brick 2)"});
    }
    auto disc = discriminant_q(f_monic, deadline);
    if (disc.is_error()) return fail<FirstLayerDescent>(disc.error());
    if (rat_is_zero(disc.value())) {
        return fail<FirstLayerDescent>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "stauduhar_first_layer: zero discriminant "
                                "(not squarefree)"});
    }
    const bool disc_square = is_rational_square_q(disc.value());

    // Certified splitting (precision grows on demand).
    auto sp = galois_padic::choose_splitting_prime(
        f_monic, ctx.max_galois_frobenius_primes(), &ctx, deadline);
    if (sp.is_error()) return fail<FirstLayerDescent>(sp.error());
    auto split0 =
        galois_padic::build_padic_splitting(f_monic, sp.value(), 1U, &ctx,
                                            deadline);
    if (split0.is_error()) return fail<FirstLayerDescent>(split0.error());
    PadicSplitting base = std::move(split0.value());
    // Parity cross-certificate: an odd Frobenius forbids G_f ⊆ A_n.
    if (disc_square && permgrp::is_odd(base.frobenius)) {
        return fail<FirstLayerDescent>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "stauduhar_first_layer: square discriminant with "
                       "odd Frobenius — parity certificate violated"});
    }

    // Ambient group and its maximal transitive candidates (Brick 2).
    const AmbientGroup amb =
        disc_square ? AmbientGroup::Alternating : AmbientGroup::Symmetric;
    auto cands = permgrp::maximal_transitive_candidates(amb, n);
    if (cands.is_error()) return fail<FirstLayerDescent>(cands.error());
    Perm ncyc = identity(n);
    for (std::size_t i = 0U; i < n; ++i) {
        ncyc[i] = static_cast<std::uint8_t>((i + 1U) % n);
    }
    Perm swap01 = identity(n);
    swap01[0] = 1U;
    swap01[1] = 0U;
    Result<BsgsGroup> ambient_res = [&]() -> Result<BsgsGroup> {
        if (!disc_square) return BsgsGroup::build(n, {swap01, ncyc});
        // A_n is the FIRST entry of the Symmetric candidate list
        // (perm_maximal contract), not of the Alternating list.
        auto sym = permgrp::maximal_transitive_candidates(
            AmbientGroup::Symmetric, n);
        if (sym.is_error()) return fail<BsgsGroup>(sym.error());
        return BsgsGroup::build(
            n, sym.value().front().group.generators());
    }();
    if (ambient_res.is_error()) {
        return fail<FirstLayerDescent>(ambient_res.error());
    }
    const BsgsGroup& ambient = ambient_res.value();
    const std::uint64_t expected_order =
        disc_square ? factorial_u64(n) / 2U : factorial_u64(n);
    if (ambient.order() != expected_order) {
        return fail<FirstLayerDescent>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "stauduhar_first_layer: ambient order mismatch"});
    }
    const std::string ambient_name =
        (disc_square ? "A" : "S") + std::to_string(n);

    for (const auto& cand : cands.value()) {
        // A_n below S_n is decided by the discriminant, not a resolvent.
        if (!disc_square && cand.group.order() == factorial_u64(n) / 2U &&
            !cand.group.has_odd_element()) {
            continue;
        }
        auto inv = galois_invariant::relative_invariant(
            ambient, cand.group, ctx.galois_lattice_max_ops(), &ctx);
        if (inv.is_error()) return fail<FirstLayerDescent>(inv.error());
        const std::size_t m = inv.value().coset_reps.size();
        // Sweep bound: each of the C(m,2)+C(n,2) separating polynomials
        // in t has degree ≤ deg(F)·(n−1); one more value must be good.
        const std::size_t t_bound =
            inv.value().total_degree * (n - 1U) *
                (m * (m - 1U) / 2U + n * (n - 1U) / 2U) +
            1U;
        bool decided = false;
        std::size_t t = 0U;
        while (!decided && t <= t_bound) {
            std::optional<PadicSplitting> model_split;
            if (t == 0U) {
                model_split = base;  // identity transform
            } else {
                auto dm = derive_model(f_monic, base, t, deadline);
                if (dm.is_error()) {
                    return fail<FirstLayerDescent>(dm.error());
                }
                if (!dm.value().has_value()) {
                    ++t;  // degenerate parameter
                    continue;
                }
                model_split = std::move(dm.value()->split);
            }
            auto out = stauduhar_step(ambient, cand.group, inv.value(),
                                      *model_split, &ctx, deadline);
            if (out.is_error()) return fail<FirstLayerDescent>(out.error());
            switch (out.value().status) {
                case StepStatus::NeedPrecision: {
                    if (out.value().required_precision <=
                        base.ring.precision()) {
                        return fail<FirstLayerDescent>(CASError{
                            .kind = CASErrorKind::InternalError,
                            .message = "stauduhar_first_layer: precision "
                                       "request without progress"});
                    }
                    auto raised = galois_padic::raise_splitting_precision(
                        base, out.value().required_precision, &ctx,
                        deadline);
                    if (raised.is_error()) {
                        return fail<FirstLayerDescent>(raised.error());
                    }
                    base = std::move(raised.value());
                    break;  // retry same t at the new precision
                }
                case StepStatus::Collision: {
                    // Certified multiple integer root: only a Tschirnhaus
                    // change of model can separate it.
                    ++t;
                    break;
                }
                case StepStatus::Descended: {
                    FirstLayerDescent r;
                    r.disc_square = disc_square;
                    r.is_ambient = false;
                    r.ambient_name = ambient_name;
                    r.subgroup = std::move(out.value().conjugated_subgroup);
                    r.provenance = cand.provenance;
                    r.splitting = std::move(base);
                    return ok(std::move(r));
                }
                case StepStatus::NotContained: {
                    decided = true;
                    break;
                }
            }
        }
        if (!decided) {
            return fail<FirstLayerDescent>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "stauduhar_first_layer: Tschirnhaus sweep "
                           "exhausted without separating the resolvent "
                           "for candidate " +
                           cand.provenance});
        }
    }
    // Coverage contract (Brick 2): no maximal transitive candidate
    // contains G_f ⇒ G_f is the ambient group, exactly.
    FirstLayerDescent r;
    r.disc_square = disc_square;
    r.is_ambient = true;
    r.ambient_name = ambient_name;
    r.splitting = std::move(base);
    return ok(std::move(r));
}

}  // namespace cas::algebra::galois_stauduhar
