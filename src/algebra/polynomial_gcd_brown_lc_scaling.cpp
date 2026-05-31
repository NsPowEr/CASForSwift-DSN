// polynomial_gcd_brown_lc_scaling.cpp — Leading-coefficient pre-scaling helpers
// for Brown's modular multivariate GCD (Geddes-Czapor-Labahn §7.4.2, Algorithm 7.2).
//
// Resolves the residual case from F3.1-BROWN-MODULAR where the true gcd has a
// polynomial (non-scalar) leading coefficient in the main variable.  Without
// pre-scaling, CRT cross-prime reconstruction recovers a candidate whose
// lc-in-main-var is wrong, the divides_sparse_z certificate fails, and the
// dispatcher falls back to gcd_eval_interp_z (Z coefficient growth).
//
// Strategy (GCL §7.4.2 + Wang-style monic normalization):
//   1. Compute lc_bound_poly L = gcd(lc_main(P), lc_main(Q)) recursively in
//      Z[x_0,...,x_{n-2}].  The true gcd's lc-in-main-var divides L.
//   2. For each lucky prime p, reduce L mod p → Lp; reduce gcd-in-Fp gp_p as
//      before.  Find a scalar u ∈ Fp* such that lc_main(u * gp_p) = Lp (up to
//      monomial alignment).  Multiply gp_p by u.
//   3. After CRT lift, candidate has lc_main = L.  Either (a) certificate
//      divides_sparse_z(P, cand) ∧ divides_sparse_z(Q, cand) succeeds directly
//      (case L = lc_main(true_gcd)), or (b) candidate has a spurious factor
//      h = L / lc_main(true_gcd) in n-1 vars — remove via main-var content
//      extraction (recursive gcd of main-var-coefficient layer).
//
// Reference: Geddes-Czapor-Labahn, "Algorithms for Computer Algebra" §7.4.2
// Algorithm 7.2 (pages 313-317).

#include "cas/algebra.hpp"
#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_gcd_fp_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <utility>
#include <vector>

namespace cas::algebra {

using fp_helpers::BMonomial;
using fp_helpers::BSparsePoly;
using fp_helpers::pos_mod;
using fp_helpers::bigint_gcd;
using fp_helpers::reduce_sparse_mod_p;
using fp_helpers::deg_in_var;

// Extract the leading-coefficient layer of `sp` w.r.t. main_var.  Returned poly
// lives in the same variable index space but every monomial has main_var-exp = 0
// (i.e., it's a polynomial in the n-1 sub-variables).
BSparsePoly extract_lc_in_var(const BSparsePoly& sp, std::size_t main_var) {
    if (sp.empty()) return {};
    std::size_t d = deg_in_var(sp, main_var);
    BSparsePoly out;
    for (const auto& [m, c] : sp) {
        std::size_t k = (main_var < m.size()) ? m[main_var] : 0U;
        if (k != d) continue;
        BMonomial m2 = m;
        if (main_var < m2.size()) m2[main_var] = 0U;
        out[m2] = c;
    }
    return out;
}

// Sparse poly mult mod p.
BSparsePoly multiply_sparse_mod_p(const BSparsePoly& A, const BSparsePoly& B,
                                   const BigInt& p, std::size_t n_vars) {
    BSparsePoly out;
    if (A.empty() || B.empty()) return out;
    for (const auto& [ma, ca] : A) {
        for (const auto& [mb, cb] : B) {
            BMonomial nm(n_vars, 0U);
            for (std::size_t i = 0; i < n_vars; ++i) {
                unsigned int ea = (i < ma.size()) ? ma[i] : 0U;
                unsigned int eb = (i < mb.size()) ? mb[i] : 0U;
                nm[i] = ea + eb;
            }
            BigInt v = pos_mod(out[nm] + ca * cb, p);
            if (v.is_zero()) out.erase(nm);
            else out[nm] = v;
        }
    }
    return out;
}

// Scale gp in place by scalar u mod p.
void scale_by_lc(BSparsePoly& gp, const BigInt& u, const BigInt& p) {
    if (u == BigInt(1)) return;
    BSparsePoly out;
    for (auto& [m, c] : gp) {
        BigInt v = pos_mod(c * u, p);
        if (!v.is_zero()) out[m] = v;
    }
    gp = std::move(out);
}

// Compute scalar unit u ∈ Fp* such that u * lcg ≡ Lp (when they are associates,
// i.e., scalar multiples of one another).  Picks the lex-leading monomial of Lp,
// looks up the matching monomial in lcg, and solves u = Lp[m] / lcg[m] mod p.
// Returns nullopt if the polynomials are NOT scalar multiples (unlucky prime).
std::optional<BigInt> compute_lc_scalar_ratio(const BSparsePoly& Lp,
                                                const BSparsePoly& lcg,
                                                const BigInt& p) {
    if (Lp.empty() || lcg.empty()) return std::nullopt;
    // Lex-leading monomial of Lp (std::map gives ascending order; last is leading).
    auto it = std::prev(Lp.end());
    const BMonomial& m_lead = it->first;
    const BigInt& a = it->second;
    auto jt = lcg.find(m_lead);
    if (jt == lcg.end()) return std::nullopt;
    const BigInt& b = jt->second;
    if (b.is_zero()) return std::nullopt;
    auto inv_b = numtheory::modular_inverse(pos_mod(b, p), p);
    if (inv_b.is_error()) return std::nullopt;
    BigInt u = pos_mod(a * inv_b.value(), p);
    // Verify: u * lcg ≡ Lp on ALL monomials (scalar-multiple check).
    // If even one monomial mismatches, lcg and Lp differ in support → unlucky.
    if (lcg.size() != Lp.size()) return std::nullopt;
    for (const auto& [m, lc_c] : lcg) {
        auto kt = Lp.find(m);
        if (kt == Lp.end()) return std::nullopt;
        BigInt prod = pos_mod(lc_c * u, p);
        if (prod != pos_mod(kt->second, p)) return std::nullopt;
    }
    return u;
}

// Reduce lc_bound poly (integer-coefficient) mod p; if leading monomial coeff
// vanishes mod p the prime is unlucky (returns empty).
BSparsePoly reduce_lc_bound_mod_p(const BSparsePoly& L, const BigInt& p) {
    BSparsePoly Lp = reduce_sparse_mod_p(L, p);
    if (Lp.empty()) return Lp;
    // Check that leading monomial of L survives mod p (else degree drop → unlucky).
    auto it_orig = std::prev(L.end());
    auto it_red = std::prev(Lp.end());
    if (it_orig->first != it_red->first) return {};  // signal unlucky.
    return Lp;
}

// Convert a BSparsePoly (n-var indexed, main_var-exp = 0 throughout) into a
// MultivariatePolynomial over the n-1 sub-variables (vars without main_var).
MultivariatePolynomial sub_sparse_to_mv(const BSparsePoly& sp,
                                         const std::vector<Symbol>& vars,
                                         std::size_t main_var) {
    std::vector<MultivariateTerm> terms;
    for (const auto& [m, c] : sp) {
        if (c.is_zero()) continue;
        std::vector<std::pair<Symbol, unsigned int>> f;
        for (std::size_t i = 0; i < vars.size(); ++i) {
            if (i == main_var) continue;
            unsigned int e = (i < m.size()) ? m[i] : 0U;
            if (e > 0U) f.emplace_back(vars[i], e);
        }
        terms.push_back(MultivariateTerm{ .coefficient = c, .factors = std::move(f) });
    }
    return MultivariatePolynomial(std::move(terms));
}

// Inverse: MultivariatePolynomial over n-1 sub-vars → BSparsePoly in n-var space
// with main_var-exp = 0.
BSparsePoly mv_to_sub_sparse(const MultivariatePolynomial& p,
                              const std::vector<Symbol>& vars,
                              std::size_t main_var) {
    BSparsePoly out;
    const std::size_t n = vars.size();
    for (const auto& term : p.terms()) {
        BMonomial mono(n, 0U);
        for (const auto& [sym, exp] : term.factors)
            for (std::size_t i = 0; i < n; ++i)
                if (i != main_var && vars[i].name == sym.name) mono[i] = exp;
        BigInt& slot = out[mono];
        slot += term.coefficient;
        if (slot.is_zero()) out.erase(mono);
    }
    return out;
}

// Pseudo-division of `dividend` (n-var sparse poly) by `divisor` (n-var sparse
// poly).  Returns true with `quo` populated iff division is EXACT in Z.
bool exact_divide_sparse_z(const BSparsePoly& dividend,
                                   const BSparsePoly& divisor,
                                   std::size_t n_vars,
                                   BSparsePoly& quo) {
    quo.clear();
    if (divisor.empty()) return false;
    if (dividend.empty()) return true;
    BSparsePoly rem = dividend;
    auto [dlm, dlc] = *std::prev(divisor.end());
    const std::size_t budget = (rem.size() + 1U) * (divisor.size() + 1U) + 16U;
    std::size_t steps = 0;
    while (!rem.empty()) {
        if (++steps > budget) return false;
        auto [rlm, rlc] = *std::prev(rem.end());
        for (std::size_t i = 0; i < n_vars; ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            if (ri < di) return false;
        }
        if ((rlc % dlc) != BigInt(0)) return false;
        BigInt qc = rlc / dlc;
        BMonomial qm(n_vars, 0U);
        for (std::size_t i = 0; i < n_vars; ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            qm[i] = static_cast<unsigned int>(ri - di);
        }
        quo[qm] = qc;
        for (const auto& [dm, dc] : divisor) {
            BMonomial nm(n_vars, 0U);
            for (std::size_t i = 0; i < n_vars; ++i)
                nm[i] = qm[i] + ((i < dm.size()) ? dm[i] : 0U);
            rem[nm] -= qc * dc;
            if (rem[nm].is_zero()) rem.erase(nm);
        }
    }
    return true;
}

bool remove_spurious_main_var_factor(const BSparsePoly& cand,
                                      const std::vector<Symbol>& vars,
                                      std::size_t main_var,
                                      symbolic::CASContext& ctx,
                                      BSparsePoly& out) {
    out.clear();
    const std::size_t n = vars.size();
    // Decompose cand into layers indexed by main_var-degree.
    std::map<std::size_t, BSparsePoly> layers;
    for (const auto& [m, c] : cand) {
        std::size_t k = (main_var < m.size()) ? m[main_var] : 0U;
        BMonomial m2 = m;
        if (main_var < m2.size()) m2[main_var] = 0U;
        layers[k][m2] = c;
    }
    if (layers.empty()) return false;
    // Iteratively gcd all layers via recursive gcd_brown_modular in n-1 sub-vars.
    MultivariatePolynomial h_mv;
    bool first = true;
    for (auto& [k, layer] : layers) {
        MultivariatePolynomial layer_mv = sub_sparse_to_mv(layer, vars, main_var);
        if (first) { h_mv = std::move(layer_mv); first = false; continue; }
        auto rec = gcd_brown_modular(h_mv, layer_mv, ctx);
        if (rec.is_error()) return false;
        h_mv = rec.value();
        if (h_mv.is_zero()) return false;
    }
    BSparsePoly h_sp = mv_to_sub_sparse(h_mv, vars, main_var);
    if (h_sp.empty()) return false;
    // If h is a unit (scalar 1 or -1), no spurious factor — caller's certificate
    // already failed, so this attempt won't help.  Bail to avoid pointless work.
    bool h_is_unit = (h_sp.size() == 1U);
    if (h_is_unit) {
        auto it = h_sp.begin();
        bool all_zero_expo = std::all_of(it->first.begin(), it->first.end(),
                                          [](unsigned int e){ return e == 0U; });
        if (all_zero_expo && (it->second == BigInt(1) || it->second == BigInt(-1)))
            return false;
    }
    // Divide each layer exactly by h_sp, re-assemble into out at main-var-degree k.
    for (const auto& [k, layer] : layers) {
        BSparsePoly quo;
        if (!exact_divide_sparse_z(layer, h_sp, n, quo)) return false;
        for (const auto& [qm, qc] : quo) {
            BMonomial nm = qm;
            if (main_var >= nm.size()) nm.resize(main_var + 1U, 0U);
            nm[main_var] = static_cast<unsigned int>(k);
            out[nm] = qc;
        }
    }
    return !out.empty();
}

}  // namespace cas::algebra
