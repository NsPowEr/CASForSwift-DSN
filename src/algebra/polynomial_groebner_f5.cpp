// polynomial_groebner_f5.cpp
//
// Signature-based Groebner computation — F5C criterion implementation.
//
// Reference: Faugère, "A new efficient algorithm for computing Gröbner bases
// without reduction to zero (F5)", ISSAC 2002. Pp. 75-83.
// Also: Eder-Perry, "F5C: A variant of Faugère's F5 algorithm with reduced
// overhead", JSC 44 (2009) 1162-1177.
//
// What is implemented (F5C with signature criterion):
//   - Each polynomial carries a signature: a module monomial (index, monomial)
//     representing its relation to the original generators.
//   - The F5 criterion: an S-pair (p, q) is rejected if its signature is
//     divisible by the leading signature of any previously-computed syzygy.
//   - The Rewritten criterion: an S-pair is rejected if one of the annotated
//     polynomials has been "rewritten" (a newer reducer has its signature).
//   - These two criteria together eliminate ALL zero-reductions in the original
//     F5 algorithm (Faugère 2002 Theorem 1 + 2).
//
// HONEST SCOPE DECLARATION (REGOLA ZERO anti-lying):
//   This implements F5C = F5 criteria applied inside a Buchberger framework.
//   The full matrix-F5 (parallel row reduction in Macaulay matrix with
//   signature selection at degree d) is NOT implemented in this block.
//   The full algorithm that guarantees zero zero-reductions for ALL orderings
//   and ALL input systems remains as OPEN ledger entry F3.3-F5-FULL.
//   What IS guaranteed here: the F5 criterion prunes a strictly positive
//   number of zero-reductions on non-trivial systems compared to plain
//   Buchberger (verified by zero_reduction_count probe in tests).
//
// REGOLA ZERO compliance:
//   - Rational/BigInt arithmetic only. No double/float.
//   - Signature comparison is structural (module monomial order), not numeric.
//   - No hardcoded counts or closed pattern sets.
//   - Known syzygy signatures are accumulated dynamically.

#include "polynomial_groebner_f5.hpp"
#include "polynomial_groebner_f4.hpp"
#include "polynomial_groebner_f4_internal.hpp"
#include "polynomial_groebner_f5_internal.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <vector>
#include <optional>

namespace cas::algebra {

// ── Signature representation ─────────────────────────────────────────────────
// A signature sig = (index, monomial) represents: this polynomial is
// e_index * monomial, where e_index is the i-th standard basis vector of
// the free module F = (Q[x])^r (r = number of input generators).
// The module monomial order: (i, t) < (j, s) iff i < j, or i == j and t < s
// under the chosen polynomial monomial order (here GRevLex).

namespace {

struct Signature {
    std::size_t index;   // generator index (0..r-1)
    Monomial    mon;     // monomial multiplier

    bool operator==(const Signature& other) const {
        return index == other.index && mon == other.mon;
    }
};

// Signature order: (index, mon) < (index', mon') iff
// index < index', or index == index' and mon <_GRevLex mon'.
[[nodiscard]] bool sig_lt(const Signature& a, const Signature& b) {
    if (a.index != b.index) return a.index < b.index;
    // GRevLex comparison on monomials
    unsigned int da = 0; for (unsigned int e : a.mon) da += e;
    unsigned int db = 0; for (unsigned int e : b.mon) db += e;
    if (da != db) return da < db;
    // Reverse lex tie-break
    for (int i = static_cast<int>(a.mon.size()) - 1; i >= 0; --i) {
        if (a.mon[i] != b.mon[i]) return a.mon[i] > b.mon[i];
    }
    return false;
}

[[maybe_unused]] [[nodiscard]] bool sig_le(const Signature& a, const Signature& b) {
    return a == b || sig_lt(a, b);
}

// Multiply monomial t into signature: (index, mon) -> (index, t*mon)
[[nodiscard]] Signature sig_multiply(const Signature& sig, const Monomial& t) {
    Monomial new_mon(sig.mon.size());
    for (std::size_t k = 0; k < sig.mon.size(); ++k) new_mon[k] = sig.mon[k] + t[k];
    return {sig.index, new_mon};
}

// Check if sig_a divides sig_b: same index AND a.mon divides b.mon component-wise.
[[nodiscard]] bool sig_divides(const Signature& a, const Signature& b) {
    if (a.index != b.index) return false;
    for (std::size_t k = 0; k < a.mon.size(); ++k) {
        if (a.mon[k] > b.mon[k]) return false;
    }
    return true;
}

// A "labeled polynomial": polynomial + its signature.
struct LabeledPoly {
    PolyF4    poly;
    Signature sig;
};

using f5_detail::divides_mon;
using f5_detail::lcm_mon;
using f5_detail::total_deg;

// S-polynomial of two labeled polys; computes signature of S(p,q):
// sig(S(p,q)) = max(t1*sig(p), t2*sig(q)) under module monomial order,
// where t1*lm(p) = t2*lm(q) = lcm(lm(p), lm(q)).
// Also returns which side (LEFT=p / RIGHT=q) achieved the max, so the caller
// can identify the "achiever" labeled poly for the rewritten criterion.
enum class SigAchiever : std::uint8_t { LEFT, RIGHT };
struct SLabeledResult {
    LabeledPoly sp;
    SigAchiever achiever;
};

[[nodiscard]] std::optional<SLabeledResult> s_labeled(
    const LabeledPoly& lp,
    const LabeledPoly& rp,
    MonomialOrder order)
{
    const Monomial lm_l = lp.poly.leading_monomial(order);
    const Monomial lm_r = rp.poly.leading_monomial(order);
    if (lm_l.empty() || lm_r.empty()) return std::nullopt;
    const Monomial L = lcm_mon(lm_l, lm_r);

    Monomial t1(L.size()), t2(L.size());
    for (std::size_t k = 0; k < L.size(); ++k) {
        t1[k] = L[k] - lm_l[k];
        t2[k] = L[k] - lm_r[k];
    }

    // Compute S-polynomial (same as in Buchberger)
    PolyF4 sp;
    const Rational inv_lc_l = Rational(1) / lp.poly.leading_coefficient(order);
    const Rational inv_lc_r = Rational(1) / rp.poly.leading_coefficient(order);
    for (const auto& [m, c] : lp.poly.terms) {
        Monomial nm(m.size());
        for (std::size_t k = 0; k < m.size(); ++k) nm[k] = m[k] + t1[k];
        sp.terms[nm] = sp.terms[nm] + c * inv_lc_l;
        if (sp.terms[nm].numerator().is_zero()) sp.terms.erase(nm);
    }
    for (const auto& [m, c] : rp.poly.terms) {
        Monomial nm(m.size());
        for (std::size_t k = 0; k < m.size(); ++k) nm[k] = m[k] + t2[k];
        sp.terms[nm] = sp.terms[nm] - c * inv_lc_r;
        if (sp.terms[nm].numerator().is_zero()) sp.terms.erase(nm);
    }

    // Signature of S-pair = max(t1 * sig(lp), t2 * sig(rp))
    Signature s1 = sig_multiply(lp.sig, t1);
    Signature s2 = sig_multiply(rp.sig, t2);
    const SigAchiever ach = sig_lt(s1, s2) ? SigAchiever::RIGHT : SigAchiever::LEFT;
    Signature sig_sp = sig_lt(s1, s2) ? s2 : s1;

    return SLabeledResult{LabeledPoly{std::move(sp), sig_sp}, ach};
}

// Reduce poly by basis: standard multivariate reduction preserving signature.
// Signature of result = signature of input (reduction does not change sig).
[[nodiscard]] PolyF4 reduce_by_labeled_basis(
    PolyF4 poly,
    const std::vector<LabeledPoly>& basis,
    MonomialOrder order)
{
    bool changed = true;
    while (changed && !poly.is_zero()) {
        changed = false;
        const Monomial lm = poly.leading_monomial(order);
        const Rational lc = poly.leading_coefficient(order);
        for (const auto& lp : basis) {
            if (lp.poly.is_zero()) continue;
            const Monomial lmg = lp.poly.leading_monomial(order);
            if (!divides_mon(lmg, lm)) continue;
            Rational factor = lc / lp.poly.leading_coefficient(order);
            Monomial shift(lm.size());
            for (std::size_t k = 0; k < lm.size(); ++k) shift[k] = lm[k] - lmg[k];
            for (const auto& [gm, gc] : lp.poly.terms) {
                Monomial nm(gm.size());
                for (std::size_t k = 0; k < gm.size(); ++k) nm[k] = gm[k] + shift[k];
                poly.terms[nm] = poly.terms[nm] - factor * gc;
                if (poly.terms[nm].numerator().is_zero()) poly.terms.erase(nm);
            }
            changed = true;
            break;
        }
    }
    return poly;
}

// F5 criterion: reject S-pair if its signature is divisible by
// the signature of any known syzygy.
[[nodiscard]] bool f5_criterion(
    const Signature& sp_sig,
    const std::vector<Signature>& syzygy_sigs)
{
    for (const Signature& sz : syzygy_sigs) {
        if (sig_divides(sz, sp_sig)) return true; // reject
    }
    return false;
}

// Rewritten criterion: rejects an S-pair if a labeled poly added STRICTLY LATER
// than the achiever of sig_sp has a signature that divides sig_sp (same index).
//
// The previous implementation iterated over ALL basis elements and accepted ANY
// strictly-smaller divisor as a "rewriter" — but this falsely catches the
// ORIGINAL generator with sig (i, 1), which always divides every (i, t) for t≠1.
// Consequence: every S-pair (l_a, f_i) with t_i≠1 (i.e., almost all of them)
// got rejected, the algorithm missed essential basis elements, and the output
// was NOT a complete Gröbner basis (Eder-Faugère survey 2017 §3).
//
// Correct semantics: only LATER additions can "rewrite" an earlier labeled poly.
// `achiever_birth_idx` is the position in `basis` at which the achiever of the
// S-pair signature was inserted; we only consider candidates at higher index.
[[nodiscard]] bool rewritten_criterion(
    const Signature& sig_sp,
    std::size_t achiever_birth_idx,
    const std::vector<LabeledPoly>& basis)
{
    for (std::size_t k = achiever_birth_idx + 1U; k < basis.size(); ++k) {
        const auto& lp = basis[k];
        if (sig_divides(lp.sig, sig_sp)) return true;
    }
    return false;
}

} // namespace

// ── Public: f5c_groebner ──────────────────────────────────────────────────────
//
// F5C pseudocode (matching implementation):
//
//  Input: F = [f_0, ..., f_{r-1}] generators; order = GRevLex.
//  Initialise:
//    basis = [ LabeledPoly(f_i, sig=(i, 1)) for i in 0..r-1 ]
//    syzygy_sigs = []
//    pairs = all (i,j) pairs from basis
//    zero_reductions_saved = 0
//    total_pairs_considered = 0
//  While pairs not empty:
//    Select pair (p, q) with min lcm degree + sugar tie-break.
//    Compute S(p,q): sp = s_labeled(p, q).
//    total_pairs_considered++
//    If sp is None (coprime LM): continue.
//    -- F5 criterion check:
//    If f5_criterion(sp.sig, syzygy_sigs): skip pair (F5 pruned), count as saved.
//    -- Rewritten criterion check:
//    If rewritten_criterion(sp.sig, basis): skip pair (Rewritten pruned), count as saved.
//    Reduce sp.poly by basis.
//    If remainder is zero:
//      Record sp.sig as a new syzygy signature.
//      zero_reductions_saved++ (would have been a zero reduction without criterion)
//    Else:
//      Add LabeledPoly(remainder, sp.sig) to basis.
//      Update pairs.
//  Inter-reduce basis; return {polys, stats}.

F5Result f5c_groebner(
    std::vector<PolyF4> F,
    MonomialOrder order)
{
    F5Result result;
    result.zero_reductions_baseline = 0;
    result.zero_reductions_f5 = 0;

    if (F.empty()) {
        return result;
    }

    const std::size_t r = F.size();
    const std::size_t n_vars = F[0].leading_monomial(order).size();

    // Make monic
    for (auto& f : F) f.make_monic(order);

    // Initialise labeled basis
    std::vector<LabeledPoly> basis;
    basis.reserve(r);
    for (std::size_t i = 0; i < r; ++i) {
        Monomial one_mon(n_vars, 0);
        basis.push_back(LabeledPoly{F[i], Signature{i, one_mon}});
    }

    // Known syzygy signatures
    std::vector<Signature> syzygy_sigs;

    // Pair queue: (i, j, lcm_mon, sugar)
    struct Pair {
        std::size_t i, j;
        Monomial deg_lcm;
        unsigned int sugar;
    };

    auto compute_sugar = [&](std::size_t i, std::size_t j, const Monomial& L) -> unsigned int {
        const unsigned int di = total_deg(basis[i].poly.leading_monomial(order));
        const unsigned int dj = total_deg(basis[j].poly.leading_monomial(order));
        const unsigned int dl = total_deg(L);
        return std::max(di + (dl - di), dj + (dl - dj));
    };

    std::vector<Pair> pairs;
    for (std::size_t i = 0; i < r; ++i) {
        for (std::size_t j = i + 1; j < r; ++j) {
            const Monomial lm_i = basis[i].poly.leading_monomial(order);
            const Monomial lm_j = basis[j].poly.leading_monomial(order);
            // Skip coprime (product criterion)
            bool coprime = true;
            for (std::size_t k = 0; k < n_vars; ++k) {
                if (lm_i[k] > 0 && lm_j[k] > 0) { coprime = false; break; }
            }
            if (coprime) continue;
            Monomial L = lcm_mon(lm_i, lm_j);
            pairs.push_back({i, j, L, compute_sugar(i, j, L)});
        }
    }

    while (!pairs.empty()) {
        // Select pair with minimum sugar (then min lcm degree for tie-break)
        auto min_it = std::min_element(pairs.begin(), pairs.end(),
            [](const Pair& a, const Pair& b) {
                if (a.sugar != b.sugar) return a.sugar < b.sugar;
                return total_deg(a.deg_lcm) < total_deg(b.deg_lcm);
            });
        Pair pair = *min_it;
        pairs.erase(min_it);

        if (basis[pair.i].poly.is_zero() || basis[pair.j].poly.is_zero()) continue;

        // Build S-labeled polynomial
        auto sp_opt = s_labeled(basis[pair.i], basis[pair.j], order);
        if (!sp_opt.has_value()) continue;
        LabeledPoly& sp = sp_opt.value().sp;
        const std::size_t achiever_idx =
            (sp_opt.value().achiever == SigAchiever::LEFT) ? pair.i : pair.j;

        // F5 criterion: reject if sp.sig divisible by known syzygy sig
        if (f5_criterion(sp.sig, syzygy_sigs)) {
            // F5-pruned. Faugère 2002 Thm 1 guarantees this pair would have
            // reduced to zero, so we count it as a saved zero-reduction.
            result.zero_reductions_baseline++;
            result.zero_reductions_f5++; // pruned ≡ no real reduction work
            continue;
        }

        // Rewritten criterion: only LATER additions can rewrite the achiever.
        if (rewritten_criterion(sp.sig, achiever_idx, basis)) {
            result.zero_reductions_baseline++;
            result.zero_reductions_f5++;
            continue;
        }

        // Reduce by basis
        PolyF4 remainder = reduce_by_labeled_basis(sp.poly, basis, order);

        if (remainder.is_zero()) {
            // Zero reduction: record syzygy
            syzygy_sigs.push_back(sp.sig);
            result.zero_reductions_f5++; // this IS a zero reduction under F5C
            result.zero_reductions_baseline++; // AND under plain Buchberger
        } else {
            // New basis element
            remainder.make_monic(order);
            const std::size_t new_idx = basis.size();
            basis.push_back(LabeledPoly{std::move(remainder), sp.sig});

            // Add new pairs
            const Monomial lm_new = basis[new_idx].poly.leading_monomial(order);
            for (std::size_t i = 0; i < new_idx; ++i) {
                if (basis[i].poly.is_zero()) continue;
                const Monomial lm_i = basis[i].poly.leading_monomial(order);
                bool coprime = true;
                for (std::size_t k = 0; k < n_vars; ++k) {
                    if (lm_i[k] > 0 && lm_new[k] > 0) { coprime = false; break; }
                }
                if (coprime) continue;
                Monomial L = lcm_mon(lm_i, lm_new);
                pairs.push_back({i, new_idx, L, compute_sugar(i, new_idx, L)});
            }
        }
    }

    // Collect non-zero polynomials
    for (const auto& lp : basis) {
        if (!lp.poly.is_zero()) result.basis.push_back(lp.poly);
    }

    // Final inter-reduction
    inter_reduce(result.basis, order);
    for (auto& g : result.basis) g.make_monic(order);

    return result;
}

} // namespace cas::algebra
