#include "polynomial_groebner_f4_internal.hpp"

#include "algebra_internal.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace cas::algebra::detail {
namespace {

// Sugar strategy (Giovini-Mora-Niesi-Robbiano 1991): each S-polynomial
// carries a virtual homogenized degree ("sugar") which is monotonically
// non-decreasing along the algorithm execution. Selection by minimum
// sugar (tie-break: minimum lcm total degree) prevents the engine from
// jumping to high-degree pairs whose reduction would be obviated by
// lower-degree progress, and yields a basis of bounded cardinality.
//
//   sugar(input_i) = total_degree(lm(input_i))
//   sugar(S(f,g)) = max(
//       sugar(f) + deg(lcm/lm_f),
//       sugar(g) + deg(lcm/lm_g))
struct Pair {
    size_t i;
    size_t j;
    Monomial deg_lcm;
    unsigned int sugar;
};

[[nodiscard]] Monomial lcm_mon(const Monomial& a, const Monomial& b) {
    Monomial res(a.size());
    for (size_t k = 0; k < a.size(); ++k) res[k] = std::max(a[k], b[k]);
    return res;
}

[[nodiscard]] unsigned int total_degree(const Monomial& m) noexcept {
    unsigned int d = 0;
    for (unsigned int e : m) d += e;
    return d;
}

// Sugar of S-polynomial (GMNR 1991).
[[nodiscard]] unsigned int s_poly_sugar(unsigned int sugar_i, unsigned int sugar_j,
                                        const Monomial& lm_i, const Monomial& lm_j,
                                        const Monomial& lcm) noexcept {
    const unsigned int deg_i = total_degree(lm_i);
    const unsigned int deg_j = total_degree(lm_j);
    const unsigned int deg_lcm = total_degree(lcm);
    const unsigned int s_left = sugar_i + (deg_lcm - deg_i);
    const unsigned int s_right = sugar_j + (deg_lcm - deg_j);
    return std::max(s_left, s_right);
}

[[nodiscard]] bool divides(const Monomial& a, const Monomial& b) {
    for (size_t k = 0; k < a.size(); ++k) {
        if (a[k] > b[k]) return false;
    }
    return true;
}

// Product criterion: gcd(lm(i), lm(j)) == 1 → S-poly reduces to 0
[[nodiscard]] bool coprime(const Monomial& a, const Monomial& b) {
    for (size_t k = 0; k < a.size(); ++k) {
        if (a[k] > 0 && b[k] > 0) return false;
    }
    return true;
}

[[nodiscard]] PolyF4 multiply_by_monomial(const PolyF4& poly, const Monomial& shift, const Rational& factor) {
    PolyF4 result;
    for (const auto& [mon, coeff] : poly.terms) {
        Monomial shifted(mon.size(), 0);
        for (size_t k = 0; k < mon.size(); ++k) shifted[k] = mon[k] + shift[k];
        result.terms[shifted] = result.terms[shifted] + coeff * factor;
        if (result.terms[shifted].numerator().is_zero()) result.terms.erase(shifted);
    }
    return result;
}

void subtract_into(PolyF4& lhs, const PolyF4& rhs) {
    for (const auto& [mon, coeff] : rhs.terms) {
        lhs.terms[mon] = lhs.terms[mon] - coeff;
        if (lhs.terms[mon].numerator().is_zero()) lhs.terms.erase(mon);
    }
}

[[nodiscard]] PolyF4 reduce_by_basis(PolyF4 poly, const std::vector<PolyF4>& basis, MonomialOrder order) {
    while (!poly.is_zero()) {
        const Monomial lm = poly.leading_monomial(order);
        const Rational lc = poly.leading_coefficient(order);
        bool reduced = false;
        for (const PolyF4& g : basis) {
            const Monomial lm_g = g.leading_monomial(order);
            if (lm_g.empty() || !divides(lm_g, lm)) continue;

            Monomial shift(lm.size(), 0);
            for (size_t k = 0; k < lm.size(); ++k) shift[k] = lm[k] - lm_g[k];
            subtract_into(poly, multiply_by_monomial(g, shift, lc / g.leading_coefficient(order)));
            reduced = true;
            break;
        }
        if (!reduced) break;
    }
    return poly;
}

[[nodiscard]] PolyF4 s_polynomial(const PolyF4& lhs, const PolyF4& rhs, MonomialOrder order) {
    const Monomial lm_l = lhs.leading_monomial(order);
    const Monomial lm_r = rhs.leading_monomial(order);
    const Monomial common = lcm_mon(lm_l, lm_r);

    Monomial shift_l(common.size(), 0);
    Monomial shift_r(common.size(), 0);
    for (size_t k = 0; k < common.size(); ++k) {
        shift_l[k] = common[k] - lm_l[k];
        shift_r[k] = common[k] - lm_r[k];
    }

    PolyF4 result = multiply_by_monomial(lhs, shift_l, Rational(1) / lhs.leading_coefficient(order));
    subtract_into(result, multiply_by_monomial(rhs, shift_r, Rational(1) / rhs.leading_coefficient(order)));
    return result;
}

// Gebauer-Moeller update: add new_idx to basis, pruning pairs via GM criteria.
// After this call, `pairs` contains the updated pair set (old surviving + new).
// `basis_sugar[k]` carries the sugar of basis[k]; new pair sugars are
// computed via s_poly_sugar (GMNR 1991).
void gm_update(std::vector<Pair>& pairs, const std::vector<PolyF4>& basis,
               const std::vector<unsigned int>& basis_sugar,
               size_t new_idx, MonomialOrder order) {
    const Monomial lm_new = basis[new_idx].leading_monomial(order);
    const unsigned int sugar_new = basis_sugar[new_idx];

    // Candidate new pairs (i, new_idx) for all i < new_idx
    std::vector<Pair> cands;
    cands.reserve(new_idx);
    for (size_t i = 0; i < new_idx; ++i) {
        const Monomial lm_i = basis[i].leading_monomial(order);
        if (lm_i.empty()) continue;
        const Monomial lcm = lcm_mon(lm_i, lm_new);
        const unsigned int sugar = s_poly_sugar(basis_sugar[i], sugar_new, lm_i, lm_new, lcm);
        cands.push_back({i, new_idx, lcm, sugar});
    }

    // Product criterion: discard pairs where lm_i ⊥ lm_new (S-poly → 0)
    cands.erase(std::remove_if(cands.begin(), cands.end(),
        [&](const Pair& p) {
            return coprime(basis[p.i].leading_monomial(order), lm_new);
        }), cands.end());

    // Gebauer-Moeller criterion for new pairs:
    // Remove (i,new) if ∃(j,new) in cands with lm_j | lcm(i,new) and lcm(j,new)≠lcm(i,new)
    cands.erase(std::remove_if(cands.begin(), cands.end(),
        [&](const Pair& p) {
            for (const Pair& q : cands) {
                if (q.i == p.i) continue;
                const Monomial lm_q = basis[q.i].leading_monomial(order);
                if (divides(lm_q, p.deg_lcm) && q.deg_lcm != p.deg_lcm) return true;
            }
            return false;
        }), cands.end());

    // Chain criterion: remove existing pairs (i,j) where lm_new | lcm(i,j)
    // AND lcm(i,new)≠lcm(i,j) AND lcm(j,new)≠lcm(i,j)
    pairs.erase(std::remove_if(pairs.begin(), pairs.end(),
        [&](const Pair& p) {
            if (!divides(lm_new, p.deg_lcm)) return false;
            const Monomial lm_i = basis[p.i].leading_monomial(order);
            const Monomial lm_j = basis[p.j].leading_monomial(order);
            const Monomial lcm_in = lcm_mon(lm_i, lm_new);
            const Monomial lcm_jn = lcm_mon(lm_j, lm_new);
            return lcm_in != p.deg_lcm && lcm_jn != p.deg_lcm;
        }), pairs.end());

    for (auto& c : cands) pairs.push_back(std::move(c));
}

// Sugar selection strategy (Giovini-Mora-Niesi-Robbiano 1991):
//   primary key:   minimum sugar
//   tie-breaker:   minimum total lcm degree (normal selection within sugar class)
//
// The sugar field on `Pair` is non-decreasing in algorithm execution,
// so progress sweeps degree by degree. This avoids the brute-force
// "normal" strategy's tendency to jump to high-degree S-pairs that
// would be obviated by lower-degree progress.
Pair select_pair(std::vector<Pair>& pairs) {
    auto it = std::min_element(pairs.begin(), pairs.end(),
        [](const Pair& a, const Pair& b) {
            if (a.sugar != b.sugar) return a.sugar < b.sugar;
            const unsigned int da = total_degree(a.deg_lcm);
            const unsigned int db = total_degree(b.deg_lcm);
            return da < db;
        });
    Pair p = std::move(*it);
    pairs.erase(it);
    return p;
}

}  // namespace

Result<std::vector<PolyF4>> buchberger_groebner(std::vector<PolyF4> basis, MonomialOrder order) {
    for (auto& g : basis) g.make_monic(order);

    // Per-element sugar (GMNR 1991). For an input generator, sugar is the
    // total degree of its leading monomial (interpreted as a virtual
    // homogenisation height). Reduction can only lower the sugar of a
    // polynomial, so the sugar of a new basis element is bounded by the
    // sugar of the S-pair from which it was produced.
    std::vector<unsigned int> basis_sugar;
    basis_sugar.reserve(basis.size());
    for (const PolyF4& g : basis) {
        basis_sugar.push_back(total_degree(g.leading_monomial(order)));
    }

    // Build initial pair set via GM update for each initial generator
    std::vector<Pair> pairs;
    for (size_t idx = 1; idx < basis.size(); ++idx) {
        gm_update(pairs, basis, basis_sugar, idx, order);
    }

    // Termination: the pair queue empties in finite time (Buchberger
    // termination theorem; Hilbert basis theorem bounds |G|). With Sugar
    // selection + Gebauer-Moeller pruning + on-fly minimisation (below),
    // pair count stays close to the theoretical lower bound, so no
    // artificial guard on |pairs| or |basis| is required.
    while (!pairs.empty()) {
        Pair pair = select_pair(pairs);

        // Skip pairs referencing now-zero basis slots (Buchberger 1985
        // on-fly inter-reduction may have invalidated either side).
        if (basis[pair.i].is_zero() || basis[pair.j].is_zero()) continue;

        PolyF4 remainder = reduce_by_basis(s_polynomial(basis[pair.i], basis[pair.j], order), basis, order);
        if (remainder.is_zero()) continue;

        remainder.make_monic(order);
        const size_t new_idx = basis.size();
        const unsigned int new_sugar = pair.sugar;
        basis.push_back(std::move(remainder));
        basis_sugar.push_back(new_sugar);

        gm_update(pairs, basis, basis_sugar, new_idx, order);

        // Buchberger 1985 on-fly inter-reduction: zero out basis
        // elements whose LM is a strict multiple of lm(basis[new_idx]),
        // and re-reduce trailing terms of all live elements. Keeps |G|
        // at its theoretical minimum (Cox-Little-O'Shea Thm 2.7.4)
        // during the run rather than only at the end.
        std::vector<std::size_t> removed;
        inter_reduce_after_addition(basis, new_idx, order, removed);
        // Discard pending pairs that touch zeroed-out basis slots; the
        // index check inside the main loop also handles this lazily,
        // but pruning here keeps the queue small.
        if (!removed.empty()) {
            pairs.erase(std::remove_if(pairs.begin(), pairs.end(),
                [&](const Pair& p) {
                    return basis[p.i].is_zero() || basis[p.j].is_zero();
                }), pairs.end());
        }
    }

    // Compact zero placeholders before the final pass.
    basis.erase(std::remove_if(basis.begin(), basis.end(),
        [](const PolyF4& g) { return g.is_zero(); }), basis.end());
    inter_reduce(basis, order);
    for (auto& g : basis) g.make_monic(order);
    return ok(std::move(basis));
}

}  // namespace cas::algebra::detail
