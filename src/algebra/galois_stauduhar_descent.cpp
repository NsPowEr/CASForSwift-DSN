// A6 Brick 3c — first-layer Stauduhar descent driver: discriminant parity
// picks the ambient (S_n / A_n), Brick-2 supplies the maximal transitive
// candidates, Brick-3b the exact invariants/transversals, Brick-3a the
// certified p-adic splitting; the per-candidate retry protocol
// (Tschirnhaus sweep, precision raises) lives in
// galois_stauduhar_candidate.cpp and is shared with the below-first-layer
// walk of Brick 3.5.

#include "galois_stauduhar_internal.hpp"

#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "galois_internal.hpp"
#include "perm_group_internal.hpp"
#include "perm_maximal_internal.hpp"
#include "polynomial_internal.hpp"
#include "polynomial_resultant_generic.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra::galois_stauduhar {

namespace {

using galois_padic::PadicSplitting;
using permgrp::AmbientGroup;
using permgrp::BsgsGroup;
using permgrp::factorial_u64;
using permgrp::identity;
using permgrp::Perm;
using primitive_internal::Deadline;

[[nodiscard]] bool rat_is_zero(const Rational& r) {
    return r.numerator().is_zero();
}

// disc(f) = (−1)^{n(n−1)/2}·Res(f, f′) for monic f.
[[nodiscard]] Result<Rational> discriminant_q(const IntPoly& f,
                                              const Deadline& dl) {
    std::vector<Rational> fr;
    fr.reserve(f.size());
    for (const auto& v : f.coefficients()) fr.emplace_back(v);
    std::vector<Rational> dr;
    for (std::size_t i = 1U; i < fr.size(); ++i) {
        dr.push_back(fr[i] * Rational(BigInt(static_cast<std::int64_t>(i))));
    }
    auto res = resultant_generic<Rational>(fr, dr, nullptr, dl);
    if (res.is_error()) return fail<Rational>(res.error());
    Rational d = res.value();
    const std::size_t n = f.degree();
    if ((n * (n - 1U) / 2U) % 2U == 1U) d = Rational(BigInt(0)) - d;
    return ok(std::move(d));
}

}  // namespace

Result<BsgsGroup> ambient_bsgs(std::size_t n, bool alternating) {
    Perm ncyc = identity(n);
    for (std::size_t i = 0U; i < n; ++i) {
        ncyc[i] = static_cast<std::uint8_t>((i + 1U) % n);
    }
    Perm swap01 = identity(n);
    swap01[0] = 1U;
    swap01[1] = 0U;
    Result<BsgsGroup> r = [&]() -> Result<BsgsGroup> {
        if (!alternating) return BsgsGroup::build(n, {swap01, ncyc});
        // A_n is the FIRST entry of the SYMMETRIC candidate list
        // (perm_maximal contract), not of the Alternating list.
        auto sym = permgrp::maximal_transitive_candidates(
            AmbientGroup::Symmetric, n);
        if (sym.is_error()) return fail<BsgsGroup>(sym.error());
        return BsgsGroup::build(n,
                                sym.value().front().group.generators());
    }();
    if (r.is_error()) return r;
    const std::uint64_t expected =
        alternating ? factorial_u64(n) / 2U : factorial_u64(n);
    if (r.value().order() != expected) {
        return fail<BsgsGroup>(
            CASError{.kind = CASErrorKind::InternalError,
                     .message = "ambient_bsgs: order mismatch"});
    }
    return r;
}

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
    auto split0 = galois_padic::build_padic_splitting(f_monic, sp.value(),
                                                      1U, &ctx, deadline);
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
    // Scan by DECREASING candidate order (= increasing index [ambient:H]):
    // the resolvent degree is the index and the certified p-adic precision
    // bound grows with it, so cheap resolvents go first and a hit on a big
    // candidate spares every deeper (more expensive) miss. Descent order
    // never changes the identified group — every step certifies
    // G_f ≤ node and the terminal test is order-independent — this is
    // resolvent economy only. Stable: preserves the Brick-2 relative
    // order (A_n stays first: its order dominates every proper candidate).
    std::stable_sort(cands.value().begin(), cands.value().end(),
                     [](const permgrp::MaximalCandidate& a,
                        const permgrp::MaximalCandidate& b) {
                         return a.group.order() > b.group.order();
                     });
    auto ambient_res = ambient_bsgs(n, disc_square);
    if (ambient_res.is_error()) {
        return fail<FirstLayerDescent>(ambient_res.error());
    }
    const BsgsGroup& ambient = ambient_res.value();
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
        auto hit = test_candidate_with_retries(ambient, cand.group,
                                               inv.value(), f_monic, base,
                                               ctx, deadline);
        if (hit.is_error()) return fail<FirstLayerDescent>(hit.error());
        if (hit.value().has_value()) {
            FirstLayerDescent r;
            r.disc_square = disc_square;
            r.is_ambient = false;
            r.ambient_name = ambient_name;
            r.subgroup = std::move(hit.value()->subgroup);
            r.provenance = cand.provenance;
            r.splitting = std::move(base);
            return ok(std::move(r));
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
