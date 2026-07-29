// A6 Brick 3a — public entry points of the p-adic splitting engine:
// unramified prime selection, certified root lifting into GR(p^k, L),
// the Frobenius permutation and the precision raise. The mod-p layer
// (Φ selection, residue-field roots) lives in galois_padic_roots.cpp.

#include "galois_padic_detail.hpp"

#include "algebraic_tower_resultant.hpp"
#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/numtheory.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <optional>
#include <random>
#include <utility>
#include <vector>

namespace cas::algebra::galois_padic {

namespace {

using detail::derivative;
using detail::newton_lift;
using detail::poll;
using primitive_internal::Deadline;

// Cycle-length multiset of a permutation.
[[nodiscard]] std::vector<std::size_t> cycle_type(const permgrp::Perm& p) {
    std::vector<std::size_t> cyc;
    std::vector<bool> seen(p.size(), false);
    for (std::size_t i = 0U; i < p.size(); ++i) {
        if (seen[i]) continue;
        std::size_t len = 0U;
        std::size_t j = i;
        while (!seen[j]) {
            seen[j] = true;
            j = p[j];
            ++len;
        }
        cyc.push_back(len);
    }
    std::sort(cyc.begin(), cyc.end());
    return cyc;
}

// Certification tripwires (see galois_padic_internal.hpp). Violation is
// always an InternalError — never a silent wrong result. `factor_degrees`
// is the (sorted) degree multiset of f mod p for the Dedekind check.
[[nodiscard]] Result<void> certify(
    const PadicSplitting& s, const std::vector<std::size_t>& factor_degrees) {
    const PadicRing& R = s.ring;
    const std::size_t n = s.f.degree();
    if (s.roots.size() != n) {
        return fail<void>(
            CASError{.kind = CASErrorKind::InternalError,
                     .message = "galois_padic: root count != degree"});
    }
    RingElem sum = R.zero();
    RingElem prod = R.one();
    for (const auto& r : s.roots) {
        if (!R.is_zero(R.eval_int_poly(s.f, r))) {
            return fail<void>(
                CASError{.kind = CASErrorKind::InternalError,
                         .message = "galois_padic: f(root) != 0 mod p^k"});
        }
        sum = R.add(sum, r);
        prod = R.mul(prod, r);
    }
    for (std::size_t i = 0U; i < n; ++i) {
        for (std::size_t j = i + 1U; j < n; ++j) {
            if (R.equal_mod_p(s.roots[i], s.roots[j])) {
                return fail<void>(CASError{
                    .kind = CASErrorKind::InternalError,
                    .message = "galois_padic: roots collide mod p"});
            }
        }
    }
    // Newton identities at the exact ends: e_1 = −f_{n−1}, e_n = (−1)^n f_0.
    const auto& fc = s.f.coefficients();
    if (!R.equal(sum, R.from_int(-fc[n - 1U]))) {
        return fail<void>(
            CASError{.kind = CASErrorKind::InternalError,
                     .message = "galois_padic: e_1 certificate failed"});
    }
    BigInt e_n = fc[0];
    if (n % 2U == 1U) e_n = -e_n;
    if (!R.equal(prod, R.from_int(e_n))) {
        return fail<void>(
            CASError{.kind = CASErrorKind::InternalError,
                     .message = "galois_padic: e_n certificate failed"});
    }
    // Frobenius cycle type = factor-degree multiset (Dedekind).
    if (cycle_type(s.frobenius) != factor_degrees) {
        return fail<void>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "galois_padic: Frobenius cycle type does not match "
                       "the factor degrees mod p"});
    }
    return ok();
}

// The Frobenius permutation on the lifted roots: σ = the canonical ring
// automorphism with σ(x) ≡ x^p mod p, computed by Newton-lifting t^p as a
// root of Φ and evaluating coordinates. The match against the root list is
// unique because the roots are pairwise distinct mod p.
[[nodiscard]] Result<permgrp::Perm> frobenius_perm(
    const PadicRing& R, const std::vector<RingElem>& roots) {
    const std::size_t n = roots.size();
    permgrp::Perm perm(n);
    if (R.ext_degree() == 1U) {
        for (std::size_t i = 0U; i < n; ++i) {
            perm[i] = static_cast<std::uint8_t>(i);
        }
        return ok(std::move(perm));
    }
    auto ring1 = R.with_precision(1U);
    if (ring1.is_error()) return fail<permgrp::Perm>(ring1.error());
    RingElem s0 =
        ring1.value().pow(ring1.value().basis_power(1U), R.prime());
    const IntPoly& phi = R.phi();
    auto s_res = newton_lift(phi, derivative(phi), R, std::move(s0), 1U,
                             R.precision());
    if (s_res.is_error()) return fail<permgrp::Perm>(s_res.error());
    // Powers of σ(t); σ(Σ c_j t^j) = Σ c_j σ(t)^j.
    std::vector<RingElem> spow(R.ext_degree(), R.one());
    for (std::size_t j = 1U; j < R.ext_degree(); ++j) {
        spow[j] = R.mul(spow[j - 1U], s_res.value());
    }
    for (std::size_t i = 0U; i < n; ++i) {
        RingElem img = R.zero();
        for (std::size_t j = 0U; j < R.ext_degree(); ++j) {
            if (roots[i][j].is_zero()) continue;
            img = R.add(img, R.mul(R.from_int(roots[i][j]), spow[j]));
        }
        std::size_t match = n;
        for (std::size_t j = 0U; j < n; ++j) {
            if (R.equal(img, roots[j])) {
                match = j;
                break;
            }
        }
        if (match == n) {
            return fail<permgrp::Perm>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "galois_padic: Frobenius image is not a root"});
        }
        perm[i] = static_cast<std::uint8_t>(match);
    }
    return ok(std::move(perm));
}

}  // namespace

Result<SplittingPrime> choose_splitting_prime(const IntPoly& f_monic,
                                              std::size_t prime_budget,
                                              symbolic::CASContext* ctx,
                                              const Deadline& deadline) {
    const std::size_t n = f_monic.degree();
    if (n < 1U || !(f_monic.leading_coeff() == BigInt(1))) {
        return fail<SplittingPrime>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "choose_splitting_prime: monic f of degree "
                                ">= 1 required"});
    }
    std::optional<SplittingPrime> best;
    BigInt p(2);
    const IntPoly fprime = derivative(f_monic);
    for (std::size_t tried = 0U; tried < prime_budget; ++tried) {
        if (auto pr = poll(ctx, deadline); pr.is_error()) {
            return fail<SplittingPrime>(pr.error());
        }
        // p unramified ⇔ f squarefree mod p ⇔ gcd(f, f′) constant over F_p.
        // Tested EXPLICITLY: the mod-p factorizer's output shape on a
        // non-squarefree input is not a reliable ramification signal.
        auto fp_ring = PadicRing::make(
            p, 1U, IntPoly(std::vector<BigInt>{BigInt(0), BigInt(1)}));
        if (fp_ring.is_error()) return fail<SplittingPrime>(fp_ring.error());
        const PadicRing& F = fp_ring.value();
        auto fbar = detail::to_field_poly(F, f_monic);
        if (fbar.is_error()) return fail<SplittingPrime>(fbar.error());
        auto dbar = detail::to_field_poly(F, fprime);
        bool squarefree = false;
        if (dbar.is_ok()) {  // f′ ≢ 0 mod p (error = vanishes ⇒ ramified)
            auto g = rp_gcd_monic(F, fbar.value(), dbar.value());
            if (g.is_error()) return fail<SplittingPrime>(g.error());
            squarefree = g.value().size() == 1U;
        }
        if (squarefree) {
            auto fact = factor_polynomial_mod_p(f_monic, p, ctx);
            if (fact.is_error()) return fail<SplittingPrime>(fact.error());
            std::size_t total = 0U;
            std::size_t L = 1U;
            for (const auto& g : fact.value()) {
                total += g.degree();
                L = std::lcm(L, g.degree());
            }
            if (total != n) {
                return fail<SplittingPrime>(CASError{
                    .kind = CASErrorKind::InternalError,
                    .message = "choose_splitting_prime: mod-p factor "
                               "degrees do not sum to deg f on a "
                               "squarefree reduction"});
            }
            if (!best || L < best->ext_degree) {
                best = SplittingPrime{p, fact.value(), L};
                if (L == 1U) break;  // cannot improve
            }
        }
        auto np = numtheory::next_prime(p);
        if (np.is_error()) return fail<SplittingPrime>(np.error());
        p = np.value();
    }
    if (!best) {
        return fail<SplittingPrime>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "choose_splitting_prime: no unramified prime within "
                       "the configured budget (raise "
                       "max_galois_frobenius_primes)"});
    }
    return ok(std::move(*best));
}

Result<PadicSplitting> build_padic_splitting(const IntPoly& f_monic,
                                             const SplittingPrime& sp,
                                             std::size_t precision_k,
                                             symbolic::CASContext* ctx,
                                             const Deadline& deadline) {
    if (precision_k == 0U) {
        return fail<PadicSplitting>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "build_padic_splitting: precision >= 1"});
    }
    auto phi =
        detail::find_phi(sp.p, sp.factors, sp.ext_degree, ctx, deadline);
    if (phi.is_error()) return fail<PadicSplitting>(phi.error());
    auto field = PadicRing::make(sp.p, 1U, phi.value());
    if (field.is_error()) return fail<PadicSplitting>(field.error());
    const PadicRing& K = field.value();

    // Input-derived seed (never a literal): ctx RNG when available
    // (REGOLA hardcode cat. 6; same idiom as equal_degree_factorization).
    std::size_t seed = f_monic.size();
    for (const auto& c : f_monic.coefficients()) {
        seed ^= c.to_u64() + 0x9e3779b9ULL + (seed << 6U) + (seed >> 2U);
    }
    seed ^= sp.p.to_u64();
    std::mt19937 local_rng(static_cast<std::uint32_t>(seed));
    std::mt19937& rng = ctx ? ctx->rng() : local_rng;

    std::vector<RingElem> roots0;
    std::vector<std::size_t> factor_degrees;
    factor_degrees.reserve(sp.factors.size());
    for (const auto& g : sp.factors) {
        if (auto pr = poll(ctx, deadline); pr.is_error()) {
            return fail<PadicSplitting>(pr.error());
        }
        factor_degrees.push_back(g.degree());
        auto rr = detail::roots_in_field(K, g, rng, ctx, deadline);
        if (rr.is_error()) return fail<PadicSplitting>(rr.error());
        for (auto& r : rr.value()) roots0.push_back(std::move(r));
    }
    std::sort(factor_degrees.begin(), factor_degrees.end());

    auto ring_k = K.with_precision(precision_k);
    if (ring_k.is_error()) return fail<PadicSplitting>(ring_k.error());
    const IntPoly fprime = derivative(f_monic);
    std::vector<RingElem> roots;
    roots.reserve(roots0.size());
    for (auto& r0 : roots0) {
        if (auto pr = poll(ctx, deadline); pr.is_error()) {
            return fail<PadicSplitting>(pr.error());
        }
        auto lifted = newton_lift(f_monic, fprime, K, std::move(r0), 1U,
                                  precision_k);
        if (lifted.is_error()) return fail<PadicSplitting>(lifted.error());
        roots.push_back(std::move(lifted.value()));
    }

    PadicSplitting out{ring_k.value(), f_monic, std::move(roots), {}};
    auto frob = frobenius_perm(out.ring, out.roots);
    if (frob.is_error()) return fail<PadicSplitting>(frob.error());
    out.frobenius = std::move(frob.value());
    if (auto cert = certify(out, factor_degrees); cert.is_error()) {
        return fail<PadicSplitting>(cert.error());
    }
    return ok(std::move(out));
}

Result<PadicSplitting> raise_splitting_precision(const PadicSplitting& s,
                                                 std::size_t k2,
                                                 symbolic::CASContext* ctx,
                                                 const Deadline& deadline) {
    if (k2 <= s.ring.precision()) {
        return ok(PadicSplitting{s.ring, s.f, s.roots, s.frobenius});
    }
    const IntPoly fprime = derivative(s.f);
    auto ring2 = s.ring.with_precision(k2);
    if (ring2.is_error()) return fail<PadicSplitting>(ring2.error());
    std::vector<RingElem> roots;
    roots.reserve(s.roots.size());
    for (const auto& r : s.roots) {
        if (auto pr = poll(ctx, deadline); pr.is_error()) {
            return fail<PadicSplitting>(pr.error());
        }
        auto lifted =
            newton_lift(s.f, fprime, s.ring, r, s.ring.precision(), k2);
        if (lifted.is_error()) return fail<PadicSplitting>(lifted.error());
        roots.push_back(std::move(lifted.value()));
    }
    PadicSplitting out{ring2.value(), s.f, std::move(roots), {}};
    auto frob = frobenius_perm(out.ring, out.roots);
    if (frob.is_error()) return fail<PadicSplitting>(frob.error());
    out.frobenius = std::move(frob.value());
    // Each root refines in place, so the Frobenius must be unchanged; the
    // input's Frobenius cycle type was already certified against Dedekind.
    if (out.frobenius != s.frobenius) {
        return fail<PadicSplitting>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "galois_padic: Frobenius changed under precision "
                       "raise"});
    }
    if (auto cert = certify(out, cycle_type(out.frobenius));
        cert.is_error()) {
        return fail<PadicSplitting>(cert.error());
    }
    return ok(std::move(out));
}

}  // namespace cas::algebra::galois_padic
