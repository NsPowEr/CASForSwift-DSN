// A6 Brick 4 — public driver for the Galois group of an irreducible
// f ∈ Q[x] of degree 8, 9 or 10, plus the structural naming of the group.
//
// The identification itself is the fully exact Stauduhar descent built across
// Bricks 1-3.75 (galois_stauduhar::stauduhar_identify): a scalable BSGS engine,
// on-demand maximal transitive candidates, a certified p-adic containment step,
// and the below-first-layer walk (dense sublattice route + structural wreath-
// preimage route). No floating point, no transcribed group tables.
//
// Brick 4 adds two things on top of that machinery:
//   1. expr → monic integer model (same splitting field / Galois group);
//   2. `structural_transitive_group_name`: a canonical label for the certified
//      group derived ONLY from its own invariants — order, parity, primitivity,
//      block system, abelian exponent. Never a lookup table (CLAUDE.md §8): the
//      famous cases (S_n, A_n, C_n, S_s ≀ S_b) are recognised by *exact order
//      certificates*, everything else falls back to an invariant descriptor.

#include "galois_internal.hpp"
#include "galois_stauduhar_internal.hpp"
#include "perm_blocks_internal.hpp"
#include "perm_bsgs_internal.hpp"
#include "perm_group_internal.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cas::algebra {

namespace {

using permgrp::BsgsGroup;
using permgrp::Perm;

[[nodiscard]] std::uint64_t lcm_u64(std::uint64_t a, std::uint64_t b) {
    if (a == 0U || b == 0U) return 0U;
    std::uint64_t x = a, y = b;
    while (y != 0U) {
        std::uint64_t t = x % y;
        x = y;
        y = t;
    }
    return (a / x) * b;  // a/gcd · b
}

// Order of a permutation = lcm of its cycle lengths.
[[nodiscard]] std::uint64_t perm_order(const Perm& p) {
    std::uint64_t o = 1U;
    for (std::size_t len : permgrp::cycle_type(p)) {
        o = lcm_u64(o, static_cast<std::uint64_t>(len));
    }
    return o;
}

// ⟨gens⟩ is abelian ⟺ the generators commute pairwise.
[[nodiscard]] bool generators_commute(const std::vector<Perm>& gens) {
    for (std::size_t i = 0U; i < gens.size(); ++i) {
        for (std::size_t j = i + 1U; j < gens.size(); ++j) {
            if (permgrp::compose(gens[i], gens[j]) !=
                permgrp::compose(gens[j], gens[i])) {
                return false;
            }
        }
    }
    return true;
}

[[nodiscard]] std::uint64_t ipow_u64(std::uint64_t base, std::size_t e) {
    std::uint64_t r = 1U;
    for (std::size_t i = 0U; i < e; ++i) r *= base;
    return r;
}

// Monic integer model F(x) = lc^{n-1}·f(x/lc): monic, integer coefficients,
// roots lc·αᵢ ⇒ identical splitting field and Galois group. (Same construction
// as galois_deg6.cpp's private monic_integer_model, minus the root_scale that
// the Dedekind sieve there needs — the Stauduhar descent recomputes its own
// splitting primes.)
[[nodiscard]] Result<IntPoly> monic_int_model_from_rationals(
    const std::vector<Rational>& coeffs) {
    const std::size_t n = coeffs.empty() ? 0U : coeffs.size() - 1U;
    auto gcd_ = [](BigInt a, BigInt b) {
        if (a.is_negative()) a = -a;
        if (b.is_negative()) b = -b;
        while (!b.is_zero()) {
            BigInt r = a % b;
            a = b;
            b = r;
        }
        return a;
    };
    // Clear denominators: z = lcm_den · f ∈ Z[x].
    BigInt lcm_den(1);
    for (const auto& c : coeffs) {
        BigInt d = c.denominator();
        if (d.is_negative()) d = -d;
        BigInt g = gcd_(lcm_den, d);
        if (g.is_zero()) g = BigInt(1);
        lcm_den = (lcm_den / g) * d;
    }
    std::vector<BigInt> z(coeffs.size());
    for (std::size_t i = 0U; i < coeffs.size(); ++i) {
        z[i] = coeffs[i].numerator() * (lcm_den / coeffs[i].denominator());
    }
    // Root scaling: F_i = z_i · lc^{n-1-i} for i < n, F_n = 1.
    const BigInt lc = z[n];
    if (lc.is_zero()) {
        return fail<IntPoly>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_deg8: zero leading coefficient"});
    }
    std::vector<BigInt> out(n + 1U);
    out[n] = BigInt(1);
    BigInt pw(1);
    for (std::size_t i = n; i-- > 0U;) {
        out[i] = z[i] * pw;
        pw = pw * lc;
    }
    IntPoly F(std::move(out));
    F.normalize([](const BigInt& v) { return v.is_zero(); });
    if (F.degree() != n) {
        return fail<IntPoly>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "galois_deg8: monic integer model degree mismatch"});
    }
    return ok(std::move(F));
}

}  // namespace

std::string structural_transitive_group_name(const BsgsGroup& g) {
    const std::size_t n = g.degree();
    const std::uint64_t ord = g.order();
    const bool even = !g.has_odd_element();
    const std::uint64_t fact_n = permgrp::factorial_u64(n);
    const std::string tag = even ? "+" : "";

    // Full symmetric / alternating — exact certificate (order + parity).
    if (ord == fact_n) return "S" + std::to_string(n);
    if (even && ord == fact_n / 2U) return "A" + std::to_string(n);

    const auto& gens = g.generators();
    const bool abelian = generators_commute(gens);

    // A transitive group of order n acts regularly on its n points.
    if (ord == static_cast<std::uint64_t>(n)) {
        std::uint64_t exponent = 1U;
        for (const auto& p : gens) exponent = lcm_u64(exponent, perm_order(p));
        if (abelian) {
            if (exponent == ord) return "C" + std::to_string(n);
            // (|G|, exponent) is a complete invariant of abelian groups of
            // order ≤ 15 ⇒ unambiguous for the regular degree-≤10 case.
            return "Ab" + std::to_string(n) + "_exp" + std::to_string(exponent);
        }
        return "Reg" + std::to_string(n) + "_" + std::to_string(ord) + tag;
    }

    // Primitive but not the full S_n / A_n.
    auto prim = permgrp::is_primitive(n, gens);
    if (prim.is_ok() && prim.value()) {
        return "P" + std::to_string(n) + "_" + std::to_string(ord) + tag;
    }

    // Imprimitive: name via a block system. Full wreath is an exact certificate:
    // G always stabilises its own block system B, and the full stabiliser of B
    // with b blocks of size s is S_s ≀ S_b of order (s!)^b·b!; equality of order
    // ⇒ G = that stabiliser. Otherwise fall back to an invariant descriptor.
    auto bss = permgrp::minimal_block_systems(n, gens);
    if (bss.is_ok() && !bss.value().empty()) {
        const permgrp::BlockSystem* coarsest = nullptr;
        std::size_t coarsest_blocksize = 0U;
        for (const auto& sys : bss.value()) {
            const std::size_t b = sys.num_blocks;
            const std::size_t s = n / b;
            const std::uint64_t wreath =
                ipow_u64(permgrp::factorial_u64(s), b) * permgrp::factorial_u64(b);
            if (ord == wreath) {
                return "S" + std::to_string(s) + "wrS" + std::to_string(b);
            }
            if (s > coarsest_blocksize) {
                coarsest_blocksize = s;
                coarsest = &sys;
            }
        }
        if (coarsest != nullptr) {
            const std::size_t b = coarsest->num_blocks;
            const std::size_t s = n / b;
            return "Im" + std::to_string(n) + "[" + std::to_string(s) + "^" +
                   std::to_string(b) + "]_" + std::to_string(ord) + tag;
        }
    }
    // Certified order + parity descriptor (never a guess).
    return "T" + std::to_string(n) + "_" + std::to_string(ord) + tag;
}

Result<std::string> galois_group_irreducible_deg8_to_10(
    ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx) {
    auto parsed = parse_polynomial(poly, var, ctx);
    if (parsed.is_error()) return fail<std::string>(parsed.error());
    auto rc = poly_to_rational_coefficients(parsed.value());
    if (rc.is_error()) return fail<std::string>(rc.error());
    const std::vector<Rational>& coeffs = rc.value();
    const std::size_t n = coeffs.empty() ? 0U : coeffs.size() - 1U;
    if (n < 8U || n > 10U) {
        return fail<std::string>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "galois_deg8: exact Stauduhar driver wired for degree "
                       "8, 9, 10 only"});
    }
    auto F = monic_int_model_from_rationals(coeffs);
    if (F.is_error()) return fail<std::string>(F.error());
    auto id = galois_stauduhar::stauduhar_identify(F.value(), ctx);
    if (id.is_error()) return fail<std::string>(id.error());
    return ok(structural_transitive_group_name(id.value().group));
}

}  // namespace cas::algebra
