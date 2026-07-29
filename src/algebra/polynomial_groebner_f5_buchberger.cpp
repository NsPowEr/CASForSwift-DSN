// polynomial_groebner_f5_buchberger.cpp
// Baseline plain Buchberger (zero-reduction counting), split out of
// polynomial_groebner_f5.cpp (T-049 anti-monolith).

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

using f5_detail::divides_mon;
using f5_detail::lcm_mon;
using f5_detail::total_deg;

// ── Plain Buchberger with zero-reduction counting ─────────────────────────────
// Used as baseline for F5 comparison probe.

BuchbergerCountResult buchberger_with_zero_count(
    std::vector<PolyF4> F,
    MonomialOrder order)
{
    BuchbergerCountResult result;
    result.zero_reductions = 0;

    if (F.empty()) return result;

    const std::size_t n_vars = F[0].leading_monomial(order).size();
    for (auto& f : F) f.make_monic(order);

    struct Pair {
        std::size_t i, j;
        Monomial deg_lcm;
        unsigned int sugar;
    };

    auto s_poly = [&](const PolyF4& lp, const PolyF4& rp) -> PolyF4 {
        const Monomial lm_l = lp.leading_monomial(order);
        const Monomial lm_r = rp.leading_monomial(order);
        const Monomial L = lcm_mon(lm_l, lm_r);
        Monomial t1(L.size()), t2(L.size());
        for (std::size_t k = 0; k < L.size(); ++k) {
            t1[k] = L[k] - lm_l[k];
            t2[k] = L[k] - lm_r[k];
        }
        PolyF4 sp;
        const Rational inv_l = Rational(1) / lp.leading_coefficient(order);
        const Rational inv_r = Rational(1) / rp.leading_coefficient(order);
        for (const auto& [m, c] : lp.terms) {
            Monomial nm(m.size());
            for (std::size_t k = 0; k < m.size(); ++k) nm[k] = m[k] + t1[k];
            sp.terms[nm] = sp.terms[nm] + c * inv_l;
            if (sp.terms[nm].numerator().is_zero()) sp.terms.erase(nm);
        }
        for (const auto& [m, c] : rp.terms) {
            Monomial nm(m.size());
            for (std::size_t k = 0; k < m.size(); ++k) nm[k] = m[k] + t2[k];
            sp.terms[nm] = sp.terms[nm] - c * inv_r;
            if (sp.terms[nm].numerator().is_zero()) sp.terms.erase(nm);
        }
        return sp;
    };

    auto reduce_fully = [&](PolyF4 poly, const std::vector<PolyF4>& basis) -> PolyF4 {
        bool changed = true;
        while (changed && !poly.is_zero()) {
            changed = false;
            const Monomial lm = poly.leading_monomial(order);
            const Rational lc = poly.leading_coefficient(order);
            for (const auto& g : basis) {
                if (g.is_zero()) continue;
                const Monomial lmg = g.leading_monomial(order);
                if (!divides_mon(lmg, lm)) continue;
                Rational factor = lc / g.leading_coefficient(order);
                Monomial shift(lm.size());
                for (std::size_t k = 0; k < lm.size(); ++k) shift[k] = lm[k] - lmg[k];
                for (const auto& [gm, gc] : g.terms) {
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
    };

    std::vector<PolyF4> basis = F;
    std::vector<Pair> pairs;
    auto make_pairs = [&](std::size_t new_idx) {
        const Monomial lm_new = basis[new_idx].leading_monomial(order);
        for (std::size_t i = 0; i < new_idx; ++i) {
            const Monomial lm_i = basis[i].leading_monomial(order);
            bool coprime = true;
            for (std::size_t k = 0; k < n_vars; ++k) {
                if (lm_i[k] > 0 && lm_new[k] > 0) { coprime = false; break; }
            }
            if (coprime) continue;
            Monomial L = lcm_mon(lm_i, lm_new);
            const unsigned int di = total_deg(lm_i);
            const unsigned int dn = total_deg(lm_new);
            const unsigned int dl = total_deg(L);
            unsigned int sugar = std::max(di + (dl - di), dn + (dl - dn));
            pairs.push_back({i, new_idx, L, sugar});
        }
    };

    for (std::size_t i = 1; i < basis.size(); ++i) make_pairs(i);

    while (!pairs.empty()) {
        auto min_it = std::min_element(pairs.begin(), pairs.end(),
            [](const Pair& a, const Pair& b) {
                if (a.sugar != b.sugar) return a.sugar < b.sugar;
                return total_deg(a.deg_lcm) < total_deg(b.deg_lcm);
            });
        Pair pair = *min_it;
        pairs.erase(min_it);

        if (basis[pair.i].is_zero() || basis[pair.j].is_zero()) continue;

        PolyF4 sp = s_poly(basis[pair.i], basis[pair.j]);
        PolyF4 rem = reduce_fully(std::move(sp), basis);

        if (rem.is_zero()) {
            result.zero_reductions++;
        } else {
            rem.make_monic(order);
            std::size_t new_idx = basis.size();
            basis.push_back(rem);
            make_pairs(new_idx);
        }
    }

    basis.erase(
        std::remove_if(basis.begin(), basis.end(), [](const PolyF4& g) { return g.is_zero(); }),
        basis.end());
    inter_reduce(basis, order);
    for (auto& g : basis) g.make_monic(order);
    result.basis = std::move(basis);
    return result;
}

} // namespace cas::algebra
