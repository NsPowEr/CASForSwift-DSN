#include "polynomial_groebner_f4_internal.hpp"

#include "algebra_internal.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace cas::algebra::detail {
namespace {

struct Pair {
    size_t i;
    size_t j;
    Monomial deg_lcm;
};

[[nodiscard]] Monomial lcm_mon(const Monomial& a, const Monomial& b) {
    Monomial res(a.size());
    for (size_t k = 0; k < a.size(); ++k) res[k] = std::max(a[k], b[k]);
    return res;
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
void gm_update(std::vector<Pair>& pairs, const std::vector<PolyF4>& basis,
               size_t new_idx, MonomialOrder order) {
    const Monomial lm_new = basis[new_idx].leading_monomial(order);

    // Candidate new pairs (i, new_idx) for all i < new_idx
    std::vector<Pair> cands;
    cands.reserve(new_idx);
    for (size_t i = 0; i < new_idx; ++i) {
        const Monomial lm_i = basis[i].leading_monomial(order);
        if (lm_i.empty()) continue;
        cands.push_back({i, new_idx, lcm_mon(lm_i, lm_new)});
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

// Select pair with minimum total degree of lcm (normal selection strategy)
Pair select_pair(std::vector<Pair>& pairs) {
    auto it = std::min_element(pairs.begin(), pairs.end(),
        [](const Pair& a, const Pair& b) {
            size_t da = 0, db = 0;
            for (size_t k = 0; k < a.deg_lcm.size(); ++k) da += a.deg_lcm[k];
            for (size_t k = 0; k < b.deg_lcm.size(); ++k) db += b.deg_lcm[k];
            return da < db;
        });
    Pair p = std::move(*it);
    pairs.erase(it);
    return p;
}

}  // namespace

Result<std::vector<PolyF4>> buchberger_groebner(std::vector<PolyF4> basis, MonomialOrder order) {
    constexpr std::size_t kMaxBuchbergerPairs = 8192;
    constexpr std::size_t kMaxBasisSize = 256;

    for (auto& g : basis) g.make_monic(order);

    // Build initial pair set via GM update for each initial generator
    std::vector<Pair> pairs;
    for (size_t idx = 1; idx < basis.size(); ++idx) {
        gm_update(pairs, basis, idx, order);
    }

    std::size_t processed = 0;
    while (!pairs.empty()) {
        if (++processed > kMaxBuchbergerPairs || basis.size() > kMaxBasisSize) {
            return fail<std::vector<PolyF4>>(make_error(
                CASErrorKind::Timeout,
                "Buchberger fallback exceeded exact algebra resource guard"));
        }

        Pair pair = select_pair(pairs);

        PolyF4 remainder = reduce_by_basis(s_polynomial(basis[pair.i], basis[pair.j], order), basis, order);
        if (remainder.is_zero()) continue;

        remainder.make_monic(order);
        const size_t new_idx = basis.size();
        basis.push_back(std::move(remainder));

        gm_update(pairs, basis, new_idx, order);
    }

    inter_reduce(basis, order);
    for (auto& g : basis) g.make_monic(order);
    return ok(std::move(basis));
}

}  // namespace cas::algebra::detail
