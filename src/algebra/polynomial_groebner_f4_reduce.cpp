#include "polynomial_groebner_f4.hpp"

#include <algorithm>
#include <utility>

namespace cas::algebra {
namespace {

[[nodiscard]] bool divides(const Monomial& a, const Monomial& b) {
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) return false;
    }
    return true;
}

// Reduce `poly` by `basis` (single-step multivariate reduction).
[[nodiscard]] PolyF4 reduce_fully(PolyF4 poly, const std::vector<PolyF4>& basis, MonomialOrder order) {
    bool any = true;
    while (any && !poly.is_zero()) {
        any = false;
        for (auto it = poly.terms.begin(); it != poly.terms.end(); ) {
            const Monomial& m = it->first;
            bool reduced = false;
            for (const PolyF4& g : basis) {
                const Monomial lmg = g.leading_monomial(order);
                if (!divides(lmg, m)) continue;
                Rational factor = it->second / g.leading_coefficient(order);
                Monomial shift(m.size());
                for (size_t k = 0; k < m.size(); ++k) shift[k] = m[k] - lmg[k];
                for (const auto& [gm, gc] : g.terms) {
                    Monomial nm(gm.size());
                    for (size_t k = 0; k < gm.size(); ++k) nm[k] = gm[k] + shift[k];
                    poly.terms[nm] = poly.terms[nm] - factor * gc;
                    if (poly.terms[nm].numerator().is_zero()) poly.terms.erase(nm);
                }
                any = true;
                reduced = true;
                it = poly.terms.begin();
                break;
            }
            if (!reduced) ++it;
        }
    }
    return poly;
}

// lcm of two monomials
[[nodiscard]] Monomial lcm_mon(const Monomial& a, const Monomial& b) {
    Monomial r(a.size());
    for (size_t k = 0; k < a.size(); ++k) r[k] = std::max(a[k], b[k]);
    return r;
}

// S-polynomial of two polynomials
[[nodiscard]] PolyF4 s_poly(const PolyF4& f, const PolyF4& g, MonomialOrder order) {
    const Monomial lmf = f.leading_monomial(order);
    const Monomial lmg = g.leading_monomial(order);
    const Monomial L   = lcm_mon(lmf, lmg);
    Monomial sf(lmf.size()), sg(lmg.size());
    for (size_t k = 0; k < L.size(); ++k) { sf[k] = L[k] - lmf[k]; sg[k] = L[k] - lmg[k]; }

    PolyF4 result;
    Rational inv_lc_f = Rational(1) / f.leading_coefficient(order);
    Rational inv_lc_g = Rational(1) / g.leading_coefficient(order);
    for (const auto& [m, c] : f.terms) {
        Monomial nm(m.size()); for (size_t k=0;k<m.size();++k) nm[k]=m[k]+sf[k];
        result.terms[nm] = result.terms[nm] + c * inv_lc_f;
        if (result.terms[nm].numerator().is_zero()) result.terms.erase(nm);
    }
    for (const auto& [m, c] : g.terms) {
        Monomial nm(m.size()); for (size_t k=0;k<m.size();++k) nm[k]=m[k]+sg[k];
        result.terms[nm] = result.terms[nm] - c * inv_lc_g;
        if (result.terms[nm].numerator().is_zero()) result.terms.erase(nm);
    }
    return result;
}

}  // namespace

void inter_reduce(std::vector<PolyF4>& G, MonomialOrder order) {  // NOLINT(bugprone-easily-swappable-parameters)
    if (G.empty()) return;
    size_t n_vars = G[0].leading_monomial(order).size();

    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < G.size(); ++i) {
            PolyF4 f = G[i];
            bool f_changed = false;
            auto it = f.terms.begin();
            while (it != f.terms.end()) {
                Monomial m = it->first;
                bool reduced = false;
                for (size_t j = 0; j < G.size(); ++j) {
                    if (i == j) continue;
                    Monomial lmj = G[j].leading_monomial(order);
                    if (divides(lmj, m)) {
                        Rational factor = it->second / G[j].leading_coefficient(order);
                        Monomial t(n_vars);
                        for (size_t k = 0; k < n_vars; ++k) t[k] = m[k] - lmj[k];

                        for (const auto& [mon_j, coeff_j] : G[j].terms) {
                            Monomial nm(n_vars);
                            for (size_t k = 0; k < n_vars; ++k) nm[k] = mon_j[k] + t[k];
                            f.terms[nm] = f.terms[nm] - factor * coeff_j;
                            if (f.terms[nm].numerator().is_zero()) f.terms.erase(nm);
                        }
                        f_changed = true;
                        reduced = true;
                        it = f.terms.begin();
                        break;
                    }
                }
                if (!reduced) ++it;
            }

            if (f_changed) {
                if (f.is_zero()) {
                    G.erase(G.begin() + i);
                    changed = true;
                    break;
                }
                f.make_monic(order);
                G[i] = std::move(f);
                changed = true;
            }
        }
    }
}

// Buchberger 1985 on-fly tail-reduction (Cox-Little-O'Shea Thm 2.7.4
// safe variant): after a new polynomial G[new_idx] is appended, re-
// reduce *trailing terms* (not leading) of every other live basis
// element by the rest of G. This keeps each basis element close to
// reduced form during the algorithm run, which feeds simpler S-polys
// into subsequent gm_update steps.
//
// Important: full Buchberger minimization (removing elements whose LM
// is divisible by another's LM) is NOT safe to do during the run —
// the pair queue still contains S-polynomials whose computation
// references these elements. Minimization is therefore deferred to
// the post-loop `inter_reduce` pass. On-fly we only tighten tails.
//
// `removed_indices` always returns empty here (no element is dropped
// during the run); the parameter is retained for API stability and to
// support a future full-minimization variant that synchronises the
// pair queue.
void inter_reduce_after_addition(
    std::vector<PolyF4>& G,
    std::size_t new_idx,
    MonomialOrder order,
    std::vector<std::size_t>& removed_indices) {
    removed_indices.clear();
    if (new_idx >= G.size()) return;
    const Monomial lm_new = G[new_idx].leading_monomial(order);
    if (lm_new.empty()) return;
    const std::size_t n_vars = lm_new.size();

    for (std::size_t i = 0; i < G.size(); ++i) {
        if (G[i].is_zero()) continue;
        PolyF4 f = G[i];
        const Monomial lm_f = f.leading_monomial(order);
        bool f_changed = false;
        auto it = f.terms.begin();
        while (it != f.terms.end()) {
            const Monomial m = it->first;
            // Preserve the leading monomial — reducing it would imply
            // dropping the element from the basis, which is unsafe
            // during the run (see header comment).
            if (m == lm_f) { ++it; continue; }
            bool reduced = false;
            for (std::size_t j = 0; j < G.size(); ++j) {
                if (i == j) continue;
                if (G[j].is_zero()) continue;
                const Monomial lmj = G[j].leading_monomial(order);
                if (lmj.empty()) continue;
                if (divides(lmj, m)) {
                    Rational factor = it->second / G[j].leading_coefficient(order);
                    Monomial t(n_vars);
                    for (std::size_t k = 0; k < n_vars; ++k) t[k] = m[k] - lmj[k];
                    for (const auto& [mon_j, coeff_j] : G[j].terms) {
                        Monomial nm(n_vars);
                        for (std::size_t k = 0; k < n_vars; ++k) nm[k] = mon_j[k] + t[k];
                        // Skip writes that would alter the leading term.
                        if (nm == lm_f) continue;
                        f.terms[nm] = f.terms[nm] - factor * coeff_j;
                        if (f.terms[nm].numerator().is_zero()) f.terms.erase(nm);
                    }
                    f_changed = true;
                    reduced = true;
                    it = f.terms.begin();
                    break;
                }
            }
            if (!reduced) ++it;
        }
        if (f_changed) {
            G[i] = std::move(f);
        }
    }
}

bool is_reduced_groebner_basis(const std::vector<PolyF4>& G, MonomialOrder order) {
    if (G.empty()) return true;

    for (size_t i = 0; i < G.size(); ++i) {
        if (G[i].is_zero()) return false;

        // (2) Monic check
        const Rational lc = G[i].leading_coefficient(order);
        if (lc != Rational(1)) return false;

        const Monomial lmi = G[i].leading_monomial(order);

        for (size_t j = 0; j < G.size(); ++j) {
            if (i == j) continue;
            const Monomial lmj = G[j].leading_monomial(order);

            // (3) No LM of G[j] divides any term of G[i]
            for (const auto& [m, _] : G[i].terms) {
                if (divides(lmj, m)) return false;
            }

            // (1) Groebner: S(i,j) reduces to 0
            PolyF4 sp = s_poly(G[i], G[j], order);
            // Build basis without G[i] and G[j] for reduction? No — reduce by full G.
            PolyF4 rem = reduce_fully(std::move(sp), G, order);
            if (!rem.is_zero()) return false;
        }
    }
    return true;
}

}  // namespace cas::algebra
