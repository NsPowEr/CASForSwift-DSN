// A6 Brick 3.5 — the below-first-layer Stauduhar walk: full certified
// identification of G_f (float-free Fieker-Klüners route; spec
// Galois_Groups.md read per REGOLA 0.1, with the approved exact p-adic
// course correction over the spec's MPFR sketch).
//
// Mathematics of the walk. After the first layer (Brick 3c) either
// G_f = ambient EXACTLY (Brick-2 coverage) or G_f ≤ current is certified
// for a proper first-layer node. Then, repeatedly:
//
//   • enumerate the maximal transitive subgroup classes of `current` up
//     to current-conjugacy (galois_sublattice — exhaustive
//     single-generator-extension lattice, completeness by construction);
//   • run the certified Stauduhar test on each class (Brick 3b invariant
//     + Brick 3c step with precision raises and the Tschirnhaus sweep);
//   • a hit descends: G_f ≤ σKσ⁻¹ strictly smaller — the p-adic
//     splitting (and its root indexing) carries over unchanged, so every
//     later step keeps certifying against the SAME roots;
//   • no hit terminates: a proper transitive G_f < current would lie in
//     some maximal subgroup M of current, and M ⊇ G_f is transitive, so
//     M is inside one of the scanned classes — hence G_f = current,
//     EXACTLY. No table is consulted anywhere.
//
// Termination is a theorem: each descent strictly decreases |current|.
// The Frobenius element belongs to G_f, so `frobenius ∈ current` is
// re-checked after every descent (tripwire, never silenced). The naming
// of the final group and the driver wiring are Brick 4.

#include "galois_stauduhar_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "galois_sublattice_internal.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <string>
#include <utility>

namespace cas::algebra::galois_stauduhar {

namespace {

using galois_padic::PadicSplitting;
using permgrp::BsgsGroup;
using primitive_internal::Deadline;

}  // namespace

Result<GaloisIdentification> stauduhar_identify(const IntPoly& f_monic,
                                                symbolic::CASContext& ctx,
                                                const Deadline& deadline) {
    auto fl = stauduhar_first_layer(f_monic, ctx, deadline);
    if (fl.is_error()) return fail<GaloisIdentification>(fl.error());
    const bool disc_square = fl.value().disc_square;
    std::string ambient_name = std::move(fl.value().ambient_name);

    if (fl.value().is_ambient) {
        auto amb = ambient_bsgs(f_monic.degree(), disc_square);
        if (amb.is_error()) return fail<GaloisIdentification>(amb.error());
        return ok(GaloisIdentification{disc_square, std::move(ambient_name),
                                       std::move(amb.value()), 0U, {}});
    }

    BsgsGroup current = std::move(*fl.value().subgroup);
    PadicSplitting base = std::move(*fl.value().splitting);
    std::string provenance = std::move(fl.value().provenance);
    std::size_t steps = 1U;

    while (true) {
        if (auto chk = ctx.check_interrupt(); chk.is_error()) {
            return fail<GaloisIdentification>(chk.error());
        }
        auto classes = permgrp::transitive_subgroup_classes_in(
            current, ctx.galois_lattice_max_ops(),
            ctx.galois_sublattice_max_bytes(), &ctx);
        if (classes.is_error()) {
            return fail<GaloisIdentification>(classes.error());
        }
        bool descended = false;
        for (const BsgsGroup& cand : classes.value()) {
            auto inv = galois_invariant::relative_invariant(
                current, cand, ctx.galois_lattice_max_ops(), &ctx);
            if (inv.is_error()) {
                return fail<GaloisIdentification>(inv.error());
            }
            auto hit = test_candidate_with_retries(
                current, cand, inv.value(), f_monic, base, ctx, deadline);
            if (hit.is_error()) {
                return fail<GaloisIdentification>(hit.error());
            }
            if (hit.value().has_value()) {
                if (!hit.value()->subgroup.contains(base.frobenius)) {
                    return fail<GaloisIdentification>(CASError{
                        .kind = CASErrorKind::InternalError,
                        .message =
                            "stauduhar_identify: Frobenius left the "
                            "descended group — containment certificate "
                            "violated"});
                }
                current = std::move(hit.value()->subgroup);
                ++steps;
                descended = true;
                break;
            }
        }
        if (!descended) break;  // G_f = current, exactly (see header)
    }

    return ok(GaloisIdentification{disc_square, std::move(ambient_name),
                                   std::move(current), steps,
                                   std::move(provenance)});
}

}  // namespace cas::algebra::galois_stauduhar
