// A6 / CAS-L3-18 — Galois group of an irreducible f ∈ Q[x] of degree 6 or 7.
//
// Fully exact pipeline (no floats, no transcribed group tables), n-generic
// throughout — only the candidate lattice and the naming depend on n:
//
//   1. Normalize to a monic integer model (root scaling x → lc·x preserves
//      the splitting field, hence the group and its natural action).
//   2. Discriminant parity: disc square ⇔ G ⊆ Aₙ (eliminates half the
//      transitive classes deterministically).
//   3. Dedekind/Frobenius sieve: for budgeted primes p ∤ lc·disc, the
//      factor-degree multiset of f mod p is the cycle type of Frob_p in G —
//      a candidate class lacking that cycle type is eliminated (positive
//      witnesses only: sound eliminations, never probabilistic guesses).
//   4. k-set resolvents (k = 2, then 3 while still ambiguous): the factor-
//      degree multiset of the squarefree R_k(y) = ∏_{|T|=k}(y − Σ_{i∈T} α_i)
//      over Q equals the orbit-length multiset of G on k-subsets (Soicher-
//      McKay). Root-sum collisions are resolved by a Tschirnhaus power-basis
//      transform (bounded sweep, bound derived below).
//   5. Candidate set = exhaustively generated transitive lattice of Sₙ
//      (galois_transitive_lattice.cpp — completeness by construction).
//      Unique survivor → structural name. Several survivors → structured
//      Unimplemented listing them (REGOLA ZERO: never guess); the
//      Stauduhar relative-resolvent descent is the follow-up increment.
//
// Spec: MISSING_FEATURES_SPECS/Galois_Groups.md (REGOLA 0.1 letta).

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/numtheory.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "galois_internal.hpp"
#include "galois_setpoly_internal.hpp"
#include "perm_group_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra {

namespace {

using permgrp::LatticeBudget;
using permgrp::PermGroup;

[[nodiscard]] bool rat_is_zero(const Rational& r) {
    return r.numerator().is_zero();
}

// Memoized transitive lattice of S_n, keyed by degree. Pure mathematics
// (independent of any CASContext), so a process-wide memo is sound; the ops
// budget of the *first* caller for a given n is used to build it (anti-runaway
// belt only — termination is proven, the subgroup lattice is finite).
[[nodiscard]] Result<std::vector<PermGroup>> transitive_lattice(
    unsigned int n, std::uint64_t max_ops) {
    static std::mutex mu;
    static std::map<unsigned int, std::vector<PermGroup>> memo;
    std::lock_guard<std::mutex> lock(mu);
    auto it = memo.find(n);
    if (it == memo.end()) {
        auto r = permgrp::transitive_subgroup_classes(
            n, LatticeBudget{.max_degree = n, .max_ops = max_ops});
        if (r.is_error()) return r;  // not cached: caller may raise budget
        it = memo.emplace(n, std::move(r.value())).first;
    }
    return ok(it->second);
}

// Exact binomial C(n, k).
[[nodiscard]] std::size_t binomial(std::size_t n, std::size_t k) {
    if (k > n) return 0U;
    std::size_t r = 1U;
    for (std::size_t i = 0U; i < k; ++i) {
        r = r * (n - i) / (i + 1U);
    }
    return r;
}

// Monic integer model with the same splitting field: for f = Σ aᵢxⁱ of
// degree n, F(x) = lcⁿ⁻¹·f(x/lc) is monic with integer coefficients and
// roots lc·αᵢ. Also returns the root scale lc: disc(F) = lc^{n(n-1)-(2n-2)}
// ·disc(f), so Dedekind primes must additionally avoid p | lc.
struct MonicModel {
    IntPoly F;
    BigInt root_scale;
};

[[nodiscard]] Result<MonicModel> monic_integer_model(const RatPoly& f) {
    const std::size_t n = f.degree();
    // Clear denominators first: f_z = d·f ∈ Z[x].
    BigInt lcm_den(1);
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
    for (const auto& c : f.coefficients()) {
        BigInt d = c.denominator();
        if (d.is_negative()) d = -d;
        BigInt g = gcd_(lcm_den, d);
        if (g.is_zero()) g = BigInt(1);
        lcm_den = (lcm_den / g) * d;
    }
    std::vector<BigInt> z(f.size());
    for (std::size_t i = 0U; i < f.size(); ++i) {
        z[i] = f[i].numerator() * (lcm_den / f[i].denominator());
    }
    // Root scaling: F_i = a_i · lc^{n-1-i} for i < n, F_n = 1.
    const BigInt lc = z[n];
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
        return fail<MonicModel>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "galois_resolvent: monic integer model degree mismatch"});
    }
    return ok(MonicModel{std::move(F), lc});
}

[[nodiscard]] RatPoly to_ratpoly(const IntPoly& p) {
    std::vector<Rational> c;
    c.reserve(p.size());
    for (const auto& v : p.coefficients()) c.emplace_back(v);
    RatPoly r(std::move(c));
    r.normalize(rat_is_zero);
    return r;
}

// Does the class contain an element with the given cycle type?
[[nodiscard]] bool contains_cycle_type(
    const PermGroup& g, const std::vector<std::size_t>& type) {
    for (const auto& [ct, cnt] : g.cycle_type_distribution()) {
        (void)cnt;
        if (ct == type) return true;
    }
    return false;
}

// Structural, deterministic label for a transitive class of degree n. Famous
// groups are named via order + parity + n-cycle content — all derived from
// the generated group itself (their orders are mathematical facts, e.g.
// |PSL(3,2)| = 168, not transcribed classifications); everything else gets an
// explicit descriptive label carrying the invariants used to reach it.
[[nodiscard]] std::string structural_name(const PermGroup& g, std::size_t n) {
    const std::uint64_t o = g.order();
    const bool odd = g.has_odd_element();
    const bool has_ncyc = contains_cycle_type(g, {n});
    if (n == 6U) {
        if (o == 720U) return "S6";
        if (o == 360U) return "A6";
        if (o == 120U) return "PGL(2,5)";
        if (o == 60U) return "PSL(2,5)";
        if (o == 6U) return has_ncyc ? "C6" : "S3(6)";
        if (o == 12U) return has_ncyc ? "D6" : "A4(6)";
        if (o == 18U) return "F18";
    } else if (n == 7U) {
        // The seven transitive groups of degree 7 all have distinct orders,
        // so the order alone determines the class.
        if (o == 5040U) return "S7";
        if (o == 2520U) return "A7";
        if (o == 168U) return "PSL(3,2)";
        if (o == 42U) return "F42";
        if (o == 21U) return "F21";
        if (o == 14U) return "D7";
        if (o == 7U) return "C7";
    }
    std::string label = std::to_string(n) + "T(o=" + std::to_string(o);
    label += odd ? ",odd" : ",even";
    if (has_ncyc) label += "," + std::to_string(n) + "cyc";
    label += ")";
    return label;
}

}  // namespace

Result<std::string> galois_group_irreducible_resolvent(
    ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx) {
    // ── Parse to RatPoly ────────────────────────────────────────────────────
    auto pe = parse_polynomial(poly, var, ctx);
    if (pe.is_error()) return fail<std::string>(pe.error());
    auto rc = poly_to_rational_coefficients(pe.value());
    if (rc.is_error()) return fail<std::string>(rc.error());
    RatPoly f_in(rc.value());
    f_in.normalize(rat_is_zero);
    const std::size_t n = f_in.degree();
    // The exact resolvent pipeline is wired for n = 6 and n = 7 (the
    // transitive lattices S₆, S₇ are generable in time; n ≥ 8 is the Stauduhar
    // maximal-descent increment). Everything below is n-generic.
    if (n != 6U && n != 7U) {
        return fail<std::string>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "galois_resolvent: exact pipeline wired for degree 6 "
                       "and 7 only (n >= 8 is the Stauduhar increment)"});
    }

    // ── Discriminant parity ─────────────────────────────────────────────────
    auto disc_res = polynomial_discriminant(poly, var, ctx);
    if (disc_res.is_error()) return fail<std::string>(disc_res.error());
    auto disc_simp = ctx.simplify(disc_res.value());
    if (disc_simp.is_error()) return fail<std::string>(disc_simp.error());
    auto disc_rat = as_rational_q(disc_simp.value());
    if (!disc_rat) {
        return fail<std::string>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "galois_resolvent: discriminant did not reduce to a "
                       "rational"});
    }
    if (disc_rat->numerator().is_zero()) {
        return fail<std::string>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_resolvent: zero discriminant (not separable)"});
    }
    const bool disc_square = is_rational_square_q(*disc_rat);

    // ── Monic integer model ─────────────────────────────────────────────────
    auto F_res = monic_integer_model(f_in);
    if (F_res.is_error()) return fail<std::string>(F_res.error());
    const IntPoly& F = F_res.value().F;
    const BigInt& root_scale = F_res.value().root_scale;

    // ── Candidate classes: exhaustive lattice + parity filter ───────────────
    auto lat = transitive_lattice(static_cast<unsigned int>(n),
                                  ctx.galois_lattice_max_ops());
    if (lat.is_error()) return fail<std::string>(lat.error());
    std::vector<const PermGroup*> cands;
    for (const auto& g : lat.value()) {
        if (g.has_odd_element() == !disc_square) cands.push_back(&g);
    }

    // ── Dedekind/Frobenius sieve (positive witnesses eliminate) ─────────────
    const std::size_t prime_budget = ctx.max_galois_frobenius_primes();
    const std::size_t bad_prime_bound = F.leading_coeff().bit_length() +
                                        disc_rat->numerator().bit_length() +
                                        disc_rat->denominator().bit_length();
    const std::size_t max_candidates_p = prime_budget + bad_prime_bound;
    BigInt p_walker(2);
    std::size_t primes_used = 0U, primes_tried = 0U;
    while (primes_used < prime_budget && primes_tried < max_candidates_p &&
           cands.size() > 1U) {
        ++primes_tried;
        if (p_walker > BigInt(std::int64_t(1) << 30)) break;
        const BigInt p = p_walker;
        auto np = numtheory::next_prime(p_walker);
        if (np.is_error()) break;
        p_walker = np.value();
        if ((disc_rat->numerator() % p).is_zero() ||
            (disc_rat->denominator() % p).is_zero() ||
            (root_scale % p).is_zero()) {
            continue;
        }
        auto fact = factor_polynomial_mod_p(F, p, &ctx);
        if (fact.is_error()) continue;
        std::vector<std::size_t> type;
        std::size_t total = 0U;
        for (const auto& gpoly : fact.value()) {
            const std::size_t d = gpoly.degree();
            if (d == 0U) continue;
            type.push_back(d);
            total += d;
        }
        if (total != n) continue;  // repeated factors mod p — skip prime
        std::sort(type.begin(), type.end(), std::greater<std::size_t>());
        ++primes_used;
        std::vector<const PermGroup*> next;
        for (const auto* g : cands) {
            if (contains_cycle_type(*g, type)) next.push_back(g);
        }
        // The true group always survives (Dedekind); an empty result means
        // an upstream bug, caught below.
        cands = std::move(next);
    }

    // ── k-set resolvent orbit patterns (k = 2, then 3 if still ambiguous) ──
    // Factor-degree multiset of the squarefree k-set resolvent over Q equals
    // the orbit-length multiset of G on k-subsets (Soicher-McKay).
    const RatPoly F_rat = to_ratpoly(F);
    auto resolvent_pattern =
        [&](std::size_t k) -> std::optional<std::vector<std::size_t>> {
        const std::size_t n_deg = F_rat.degree();
        const std::size_t rdeg = binomial(n_deg, k);  // C(n,k)
        // Tschirnhaus sweep on the moment curve (derived bound, not magic):
        // P_t(x) = x + t·x² + t²·x³ + t³·x⁴ + t⁴·x⁵. For each pair of
        // distinct k-subsets (T,T′) the separating condition
        //     Σ_m t^{m−1}·(Σ_T α^m − Σ_{T′} α^m) ≠ 0
        // is a NONZERO polynomial of degree ≤ n−2 in t — nonzero because
        // the power-sum vector (m = 1..n−1) determines the subset
        // (Vandermonde is nonsingular on the distinct roots). Same for root
        // injectivity (pairs of 1-subsets). Hence at most
        // (n−2)·[C(rdeg,2) + C(n,2)] values of t are degenerate and a sweep
        // one longer provably contains a good t. t = 0 is the identity
        // transform (P_0 = x), so the original polynomial is attempt 0.
        const std::size_t max_attempts =
            (n_deg - 2U) *
                (rdeg * (rdeg - 1U) / 2U + n_deg * (n_deg - 1U) / 2U) +
            1U;
        for (std::size_t attempt = 0U; attempt <= max_attempts; ++attempt) {
            RatPoly g_model = F_rat;
            if (attempt > 0U) {
                const BigInt t(static_cast<std::int64_t>(attempt));
                std::vector<BigInt> pc;
                pc.reserve(n_deg - 1U);
                BigInt tp(1);
                for (std::size_t m = 1U; m < n_deg; ++m) {
                    pc.push_back(tp);
                    tp = tp * t;
                }
                auto tr = galois_setpoly::tschirnhaus_general(F_rat, pc);
                if (tr.is_error()) continue;
                g_model = std::move(tr.value());
                auto gsf = galois_setpoly::is_squarefree_q(g_model);
                if (gsf.is_error() || !gsf.value()) continue;
            }
            auto res = (k == 2U)
                           ? galois_setpoly::two_set_resolvent(g_model)
                           : galois_setpoly::three_set_resolvent(g_model);
            if (res.is_error()) continue;
            auto sf = galois_setpoly::is_squarefree_q(res.value());
            if (sf.is_error() || !sf.value()) continue;
            // Factor over Z: degrees = orbit lengths of G on k-subsets.
            AstArena& arena = ctx.arena();
            Symbol y_sym = ctx.make_fresh_symbol("ygal");
            PolyExpr r_pe;
            r_pe.reserve(res.value().size());
            for (const auto& c : res.value().coefficients()) {
                if (c.denominator() == BigInt(1)) {
                    r_pe.push_back(arena.make<IntegerLit>(c.numerator()));
                } else {
                    r_pe.push_back(arena.make<RationalLit>(
                        c.numerator(), c.denominator()));
                }
            }
            auto r_expr = polynomial_to_expr(r_pe, y_sym, ctx);
            if (r_expr.is_error()) continue;
            auto fac = factor_over_integers(r_expr.value(), y_sym, ctx);
            if (fac.is_error()) continue;
            std::vector<std::size_t> degs;
            std::size_t total = 0U;
            for (const auto& pf : fac.value().factors) {
                auto pp = parse_polynomial(pf.factor, y_sym, ctx);
                if (pp.is_error()) break;
                const std::size_t d = poly_degree(pp.value());
                if (d == 0U) continue;
                for (unsigned int m = 0U; m < pf.multiplicity; ++m) {
                    degs.push_back(d);
                    total += d;
                }
            }
            if (total != rdeg) continue;
            std::sort(degs.begin(), degs.end(), std::greater<std::size_t>());
            return degs;
        }
        return std::nullopt;
    };
    for (const std::size_t k : {std::size_t{2U}, std::size_t{3U}}) {
        if (cands.size() <= 1U) break;
        const auto pattern = resolvent_pattern(k);
        if (!pattern.has_value()) continue;
        std::vector<const PermGroup*> next;
        for (const auto* g : cands) {
            if (g->orbit_lengths_on_ksubsets(k) == *pattern) {
                next.push_back(g);
            }
        }
        cands = std::move(next);
    }

    if (cands.size() == 1U) return ok(structural_name(*cands[0], n));
    if (cands.empty()) {
        return fail<std::string>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "galois_resolvent: all transitive candidates eliminated "
                       "— invariant violated (Dedekind/resolvent "
                       "inconsistency)"});
    }
    std::string diag =
        "galois_resolvent: ambiguous after discriminant, Dedekind sieve and "
        "2-/3-set resolvents; remaining candidates:";
    for (const auto* g : cands) diag += " " + structural_name(*g, n);
    diag += " — Stauduhar relative-resolvent descent is the follow-up "
            "increment (A6)";
    return fail<std::string>(
        CASError{.kind = CASErrorKind::Unimplemented, .message = diag});
}

}  // namespace cas::algebra
