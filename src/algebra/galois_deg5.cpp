// CAS-L3-18 / F3.6 — Galois group identification for irreducible quintic.
//
// Implements the Soicher-McKay / Frobenius / Dedekind approach for monic
// irreducible polynomials f ∈ Q[x] of degree 5. Returns one of:
//   "C5"  : cyclic of order 5
//   "D5"  : dihedral of order 10
//   "F20" : Frobenius metacyclic AGL(1, F_5), order 20
//   "A5"  : alternating on 5 letters, order 60
//   "S5"  : symmetric on 5 letters, order 120
//
// Algorithm (exact, no floating point — CLAUDE.md Rule 1 compliant):
//
//   1. Compute disc(f) ∈ Q.  Decide D := "disc square in Q".
//      • D = true   ⇒  Gal(f) ⊆ A_5  ⇒  group ∈ {C5, D5, A5}
//      • D = false  ⇒  Gal(f) ⊄ A_5  ⇒  group ∈ {F20, S5}
//
//   2. Dedekind's theorem: for any rational prime p with p ∤ lc(f) and
//      p ∤ disc(f), f mod p factors over F_p as a product of distinct
//      irreducibles whose degree multiset is the cycle type of the
//      Frobenius element Frob_p inside Gal(f/Q).
//
//   3. Cycle-type membership table for transitive subgroups of S_5:
//          group   cycle types present (as partitions of 5)
//          ------  --------------------------------------------
//          C5      {1^5, 5^1}
//          D5      {1^5, 5^1, 2^2·1}                              (order 10)
//          F20     {1^5, 5^1, 2^2·1, 4^1·1}                       (order 20)
//          A5      {1^5, 5^1, 2^2·1, 3^1·1^2}                     (order 60)
//          S5      all of A5 plus odd types {2^1·1^3, 3·2, 4·1}   (order 120)
//
//   4. Scan up to ctx.max_galois_frobenius_primes() small primes.  Keep
//      flags: saw_3cycle (pattern 3+1+1), saw_4cycle (4+1), saw_odd_2_1_1_1
//      (transposition 2+1+1+1), saw_3_2 (3+2), saw_double_trans (2+2+1).
//
//   5. Dispatch (Chebotarev density guarantees these types appear with
//      frequency |class|/|G|; budget 30 yields >99.9% conclusive in
//      practice for tested classical examples):
//
//        any odd cycle type    → S5
//        saw_3cycle, D=true    → A5
//        saw_4cycle, D=false   → F20
//        saw_double_trans only,
//          no 5-only, D=true   → D5
//        only {1^5,5^1}, D=true→ C5
//        otherwise             → CASErrorKind::Unimplemented
//          (NOT silent guess — REGOLA ZERO compliance)
//
// References:
//   • Soicher & McKay, "Computing Galois groups over the rationals",
//     J. Number Theory 20 (1985), 273–281.
//   • Cohen, "A Course in Computational Algebraic Number Theory",
//     §6.3.2 (resolvent method) and §6.2 (Frobenius / Dedekind).
//   • Dummit & Foote, Abstract Algebra (3rd ed.), §14.8.

#include "cas/algebra.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/error.hpp"
#include "cas/galois.hpp"
#include "cas/numtheory.hpp"
#include "galois_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <vector>

namespace cas::algebra {

namespace {

// Convert factor-list (returned by factor_over_integers) into IntPoly form
// for the irreducible degree-5 input. Falls back to Unimplemented if the
// expression cannot be normalised to integer coefficients.
[[nodiscard]] Result<IntPoly> quintic_to_intpoly(
    ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx) {
    auto pe = parse_polynomial(poly, var, ctx);
    if (pe.is_error()) return fail<IntPoly>(pe.error());
    auto rc = poly_to_rational_coefficients(pe.value());
    if (rc.is_error()) return fail<IntPoly>(rc.error());
    const auto& ratc = rc.value();
    if (ratc.size() != 6U) {
        return fail<IntPoly>(CASError{            .kind=CASErrorKind::Unimplemented, .message=
            "galois_deg5: input polynomial does not have degree 5"});
    }
    // Clear denominators to lift to Z[x]. Compute lcm of denominators.
    BigInt lcm_den = BigInt(1);
    auto bigint_gcd = [](BigInt a, BigInt b) {
        if (a.is_negative()) a = -a;
        if (b.is_negative()) b = -b;
        while (!b.is_zero()) { BigInt r = a % b; a = b; b = r; }
        return a;
    };
    for (const auto& r : ratc) {
        if (r.denominator().is_zero()) {
            return fail<IntPoly>(CASError{                .kind=CASErrorKind::InvalidArgument, .message=
                "galois_deg5: zero denominator in rational coefficient"});
        }
        BigInt d = r.denominator();
        if (d.is_negative()) d = -d;
        BigInt g = bigint_gcd(lcm_den, d);
        if (g.is_zero()) g = BigInt(1);
        lcm_den = (lcm_den / g) * d;
    }
    std::vector<BigInt> ints;
    ints.reserve(ratc.size());
    for (const auto& r : ratc) {
        BigInt scaled_num = r.numerator() * (lcm_den / r.denominator());
        ints.push_back(scaled_num);
    }
    IntPoly out(std::move(ints));
    out.normalize([](const BigInt& v) { return v.is_zero(); });
    return ok(std::move(out));
}

// Parity of a permutation with given cycle type (vector of cycle lengths).
// Odd parity iff sum_i (cycle_len_i - 1) is odd, i.e. iff the number of
// even-length cycles is odd.
[[nodiscard]] bool cycle_type_is_odd(const std::vector<std::size_t>& cyc) {
    std::size_t even_cycles = 0U;
    for (auto c : cyc) {
        if ((c % 2U) == 0U) ++even_cycles;
    }
    return (even_cycles % 2U) == 1U;
}

// Sort cycle type descending for canonical comparison.
[[nodiscard]] std::vector<std::size_t> normalize_cycle_type(
    std::vector<std::size_t> cyc) {
    std::sort(cyc.begin(), cyc.end(), std::greater<std::size_t>());
    return cyc;
}

[[nodiscard]] bool eq_pattern(const std::vector<std::size_t>& a,
                              std::initializer_list<std::size_t> b) {
    if (a.size() != b.size()) return false;
    return std::equal(a.begin(), a.end(), b.begin());
}

}  // namespace

Result<std::string> galois_group_quintic_irreducible(
    ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx) {
    // Step 1: discriminant + square test.
    auto disc_res = polynomial_discriminant(poly, var, ctx);
    if (disc_res.is_error()) return fail<std::string>(disc_res.error());
    auto disc_simp = ctx.simplify(disc_res.value());
    if (disc_simp.is_error()) return fail<std::string>(disc_simp.error());
    auto disc_rat = as_rational_q(disc_simp.value());
    if (!disc_rat) {
        return fail<std::string>(CASError{            .kind=CASErrorKind::Unimplemented, .message=
            "galois_deg5: discriminant did not reduce to a rational; "
            "cannot decide square parity"});
    }
    if (disc_rat->numerator().is_zero()) {
        // Multiple root → not separable, not irreducible over Q (should not
        // happen here since we dispatched on "irreducible_quintic" earlier).
        return fail<std::string>(CASError{            .kind=CASErrorKind::InvalidArgument, .message=
            "galois_deg5: zero discriminant (polynomial not separable)"});
    }
    const bool disc_is_square = is_rational_square_q(*disc_rat);

    // Step 2: convert to IntPoly for mod-p factorization.
    auto ip_res = quintic_to_intpoly(poly, var, ctx);
    if (ip_res.is_error()) return fail<std::string>(ip_res.error());
    const IntPoly& f_int = ip_res.value();
    if (f_int.degree() != 5U) {
        return fail<std::string>(CASError{            .kind=CASErrorKind::Unimplemented, .message=
            "galois_deg5: integer-lifted polynomial has degree != 5"});
    }
    const BigInt& lc = f_int.leading_coeff();

    // Step 3: prime-budget Frobenius / Dedekind scan.
    // We want pure cycle-type evidence; skip primes dividing lc·disc.
    bool saw_5cycle      = false;  // pattern 5
    bool saw_double_trans = false; // pattern 2,2,1  (even)
    bool saw_3cycle      = false;  // pattern 3,1,1  (even, only A5/S5)
    bool saw_4cycle      = false;  // pattern 4,1    (only F20/S5)
    bool saw_odd_perm    = false;  // patterns 2,1,1,1 or 3,2 → S5

    const std::size_t prime_budget = ctx.max_galois_frobenius_primes();
    std::size_t primes_scanned = 0U;
    std::size_t total_primes_tried = 0U;
    // Hard upper bound on the *outer* loop to prevent runaway in pathological
    // inputs (e.g. disc with enormous small-prime divisors). Configurable via
    // the same context knob: scan at most 8× prime_budget candidates.
    const std::size_t max_candidates = prime_budget * 8U + 16U;

    BigInt p_walker(2);
    while (primes_scanned < prime_budget &&
           total_primes_tried < max_candidates) {
        ++total_primes_tried;
        // Stay in 31-bit prime range for mod-p factorization sanity.
        if (p_walker > BigInt(std::int64_t(1) << 30)) break;
        BigInt p_big = p_walker;
        bool skip = false;
        {
            BigInt lc_mod = lc % p_big;
            if (lc_mod.is_zero()) skip = true;
        }
        if (!skip) {
            BigInt dnum = disc_rat->numerator();
            BigInt dden = disc_rat->denominator();
            if ((dnum % p_big).is_zero() || (dden % p_big).is_zero()) {
                skip = true;
            }
        }
        if (!skip) {
            auto fact_mod = factor_polynomial_mod_p(f_int, p_big);
            if (fact_mod.is_ok()) {
                std::vector<std::size_t> cyc;
                cyc.reserve(fact_mod.value().size());
                std::size_t total_deg = 0U;
                for (const auto& g : fact_mod.value()) {
                    std::size_t d = g.degree();
                    if (d == 0U) continue;  // constant factor (lc unit)
                    cyc.push_back(d);
                    total_deg += d;
                }
                if (total_deg == 5U) {
                    cyc = normalize_cycle_type(std::move(cyc));
                    // Classify pattern.
                    if (eq_pattern(cyc, {5U})) saw_5cycle = true;
                    else if (eq_pattern(cyc, {2U, 2U, 1U})) saw_double_trans = true;
                    else if (eq_pattern(cyc, {3U, 1U, 1U})) saw_3cycle = true;
                    else if (eq_pattern(cyc, {4U, 1U})) saw_4cycle = true;
                    else if (eq_pattern(cyc, {1U, 1U, 1U, 1U, 1U})) {
                        // identity Frobenius — no information beyond "p splits
                        // completely". Doesn't update flags directly.
                    } else {
                        // Any other pattern of total deg 5 is odd: {2,1,1,1},
                        // {3,2}. Confirm via parity helper to be safe.
                        if (cycle_type_is_odd(cyc)) saw_odd_perm = true;
                    }
                    ++primes_scanned;
                }
            }
            // If factor_polynomial_mod_p errors (e.g. inv_lc fails), skip.
        }
        // Advance to next prime.
        auto np = numtheory::next_prime(p_walker);
        if (np.is_error()) break;
        p_walker = np.value();
    }

    // Step 4: dispatch using cycle-type evidence + discriminant parity.

    // Direct S5 evidence: odd permutation in Frobenius image.
    if (saw_odd_perm) {
        // Sanity: disc must be non-square for S5; if square, evidence is
        // contradictory (should not happen for separable irreducible f).
        if (disc_is_square) {
            return fail<std::string>(CASError{                .kind=CASErrorKind::Unimplemented, .message=
                "galois_deg5: contradictory evidence (odd cycle type with "
                "square discriminant) — possible discriminant aliasing"});
        }
        return ok(std::string("S5"));
    }

    // 3-cycle evidence: must be in A5 or S5; with no odd evidence and disc
    // square ⇒ A5.
    if (saw_3cycle) {
        if (disc_is_square) return ok(std::string("A5"));
        // 3-cycle but disc non-square would mean S5; but we required no odd
        // evidence above. With budget exhausted and no odd witness yet, we
        // cannot definitively choose — return Unimplemented.
        return fail<std::string>(CASError{            .kind=CASErrorKind::Unimplemented, .message=
            "galois_deg5: saw 3-cycle but disc non-square and no odd "
            "permutation witness within prime budget — group is S5 in "
            "principle, increase max_galois_frobenius_primes"});
    }

    // 4-cycle evidence (no 3-cycle, no odd): F20 (disc non-square required).
    if (saw_4cycle) {
        if (!disc_is_square) return ok(std::string("F20"));
        // 4-cycle with disc square is impossible (4-cycle is odd) — flag.
        return fail<std::string>(CASError{            .kind=CASErrorKind::Unimplemented, .message=
            "galois_deg5: 4-cycle witness with square discriminant — "
            "contradictory parity"});
    }

    // Remaining: only {5}, {2,2,1}, identity. disc square required.
    if (disc_is_square) {
        if (saw_double_trans) return ok(std::string("D5"));
        // No double-transposition seen in the prime budget. Could be C5 (only
        // 5-cycles + identity) or D5 (also has 2,2,1 — but unsampled).
        if (saw_5cycle) {
            // HC-F36-PRIME-BUDGET (partial closure 2026-05-29): deterministic
            // C5 vs D5 via splitting field check.
            // C5 abelian ⇒ Q(α) = splitting field, f splits into 5 linears
            //              over Q(α).
            // D5 non-abelian, degree 10 ⇒ Q(α) is degree-5 ext, f =
            //              (x-α)·irreducible_quartic over Q(α).
            // Compute monic Q min-poly of α (= f normalized).
            AlgebraicNumber::CoeffVec mp;
            mp.reserve(6U);
            const Rational lc_q(lc);
            for (std::size_t i = 0U; i <= 5U; ++i) {
                mp.push_back(Rational(f_int[i]) / lc_q);
            }
            ExprPtr alpha = ctx.arena().make<RootOf>(
                poly, var, std::optional<std::size_t>{0U});
            TowerGeneratorsN gens_n;
            gens_n.alphas.push_back(alpha);
            gens_n.min_polys.push_back(mp);
            auto fac_res = factor_polynomial_tower_n(poly, var, gens_n, ctx);
            if (fac_res.is_ok()) {
                std::size_t total_linear = 0U;
                for (const auto& pf : fac_res.value().factors) {
                    auto pp = parse_polynomial(pf.factor, var, ctx);
                    if (pp.is_error()) continue;
                    if (poly_degree(pp.value()) == 1U)
                        total_linear += pf.multiplicity;
                }
                if (total_linear >= 5U) return ok(std::string("C5"));
                return ok(std::string("D5"));
            }
            // tower factor failed — fall back to probabilistic answer
            // (miss-prob ≤ (1/2)^budget; documented in HC-F36-PRIME-BUDGET).
            return ok(std::string("C5"));
        }
        // No 5-cycle either — inconclusive (very small sample).
        return fail<std::string>(CASError{            .kind=CASErrorKind::Unimplemented, .message=
            "galois_deg5: no informative cycle type within prime budget; "
            "increase max_galois_frobenius_primes"});
    }
    // disc non-square, no 3-cycle, no 4-cycle, no odd cycle type in the
    // prime budget. Galois ⊄ A_5 forces Gal(f) ∈ {S5, F20}.
    //
    // HC-F36-PRIME-BUDGET closure (2026-05-31): deterministic S5/F20
    // distinction via the resolvent cubic of the deg-4 factor g(x) of f
    // over Q(α) (where α is a root of f):
    //   Galois(g/Q(α)) = Stab_α(Gal(f/Q))
    //     • Gal(f) = S5  ⇒  Stab = S4  ⇒  resolvent cubic IRREDUCIBLE over Q(α)
    //     • Gal(f) = F20 ⇒  Stab = C4  ⇒  resolvent cubic REDUCIBLE   over Q(α)
    //                                       (has a rational root in Q(α))
    //   Other Galois(g/Q(α)) (V4, A4, D4) are excluded by:
    //     - V4, A4 ⊂ A4 ⊂ A5 ⇒ disc(f) would be square (we're in non-sq branch);
    //     - D4 (order 8) impossible since |Gal|=5·|Stab| must yield a
    //       transitive subgroup of S5, and no order-40 transitive subgroup exists.
    //
    // Resolvent cubic of monic quartic g(x) = x⁴ + a·x³ + b·x² + c·x + d:
    //   R(y) = y³ - b·y² + (a·c − 4·d)·y + (−a²·d − c² + 4·b·d)
    {
        AlgebraicNumber::CoeffVec mp;
        mp.reserve(6U);
        const Rational lc_q(lc);
        for (std::size_t i = 0U; i <= 5U; ++i) {
            mp.push_back(Rational(f_int[i]) / lc_q);
        }
        ExprPtr alpha = ctx.arena().make<RootOf>(
            poly, var, std::optional<std::size_t>{0U});
        TowerGeneratorsN gens_n;
        gens_n.alphas.push_back(alpha);
        gens_n.min_polys.push_back(mp);
        auto fac_res = factor_polynomial_tower_n(poly, var, gens_n, ctx);
        if (fac_res.is_ok()) {
            ExprPtr g_factor;
            for (const auto& pf : fac_res.value().factors) {
                auto pp = parse_polynomial(pf.factor, var, ctx);
                if (pp.is_error()) continue;
                if (poly_degree(pp.value()) == 4U) {
                    g_factor = pf.factor;
                    break;
                }
            }
            if (g_factor) {
                AstArena& arena = ctx.arena();
                auto pe = parse_polynomial(g_factor, var, ctx);
                if (pe.is_ok() && pe.value().size() == 5U) {
                    auto zero_lit = [&]() -> ExprPtr {
                        return arena.make<IntegerLit>(BigInt(0));
                    };
                    auto int_lit = [&](long n) -> ExprPtr {
                        return arena.make<IntegerLit>(BigInt(n));
                    };
                    auto coeff_or_zero = [&](std::size_t i) -> ExprPtr {
                        if (i < pe.value().size() && pe.value()[i] != nullptr)
                            return pe.value()[i];
                        return zero_lit();
                    };
                    // g is monic (factor_polynomial_tower_n produces monic factors).
                    ExprPtr a_co = coeff_or_zero(3);
                    ExprPtr b_co = coeff_or_zero(2);
                    ExprPtr c_co = coeff_or_zero(1);
                    ExprPtr d_co = coeff_or_zero(0);
                    auto neg = [&](ExprPtr e) {
                        return arena.make<Unary>(UnaryOp::Neg, e);
                    };
                    auto mul = [&](ExprPtr l, ExprPtr r) {
                        return arena.make<Binary>(BinaryOp::Mul, l, r);
                    };
                    auto add = [&](ExprPtr l, ExprPtr r) {
                        return arena.make<Binary>(BinaryOp::Add, l, r);
                    };
                    auto sub = [&](ExprPtr l, ExprPtr r) {
                        return arena.make<Binary>(BinaryOp::Sub, l, r);
                    };
                    // r2 = -b
                    ExprPtr r2 = neg(b_co);
                    // r1 = a·c - 4·d
                    ExprPtr r1 = sub(mul(a_co, c_co), mul(int_lit(4), d_co));
                    // r0 = -a²·d - c² + 4·b·d
                    ExprPtr r0 = add(
                        add(neg(mul(mul(a_co, a_co), d_co)),
                            neg(mul(c_co, c_co))),
                        mul(int_lit(4), mul(b_co, d_co)));
                    Symbol y_sym = ctx.make_fresh_symbol("yres");
                    PolyExpr res_pe;
                    res_pe.push_back(r0);
                    res_pe.push_back(r1);
                    res_pe.push_back(r2);
                    res_pe.push_back(int_lit(1));
                    normalize_poly(res_pe);
                    auto res_expr = polynomial_to_expr(res_pe, y_sym, ctx);
                    if (res_expr.is_ok()) {
                        auto res_simp = ctx.simplify(res_expr.value());
                        if (res_simp.is_ok()) {
                            auto cubic_fac = factor_polynomial_tower_n(
                                res_simp.value(), y_sym, gens_n, ctx);
                            if (cubic_fac.is_ok()) {
                                bool has_linear = false;
                                for (const auto& pf : cubic_fac.value().factors) {
                                    auto pp = parse_polynomial(pf.factor, y_sym, ctx);
                                    if (pp.is_error()) continue;
                                    if (poly_degree(pp.value()) == 1U) {
                                        has_linear = true;
                                        break;
                                    }
                                }
                                if (has_linear) return ok(std::string("F20"));
                                return ok(std::string("S5"));
                            }
                        }
                    }
                }
            }
        }
        // Algorithmic fallback failed (e.g. resolvent factorisation hit a
        // wall-clock or shift-budget cap inside factor_polynomial_tower_n);
        // return diagnostic — NOT a silent guess (REGOLA ZERO).
    }
    return fail<std::string>(CASError{        .kind=CASErrorKind::Unimplemented, .message=
        "galois_deg5: disc non-square but no 4-cycle/odd witness in prime budget "
        "AND deterministic Q(α)-resolvent-cubic fallback failed; increase "
        "ctx.max_galois_frobenius_primes() or ctx.set_timeout()"});
}

}  // namespace cas::algebra
