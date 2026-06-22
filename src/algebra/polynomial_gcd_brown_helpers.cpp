#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include "polynomial_gcd_fp_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <set>
#include <string>
#include <vector>

namespace cas::algebra::fp_helpers {

[[nodiscard]] BSparsePoly to_sparse(const MultivariatePolynomial& p,
                                    const std::vector<Symbol>& vars) {
    BSparsePoly sp;
    for (const auto& term : p.terms()) {
        BMonomial mono(vars.size(), 0U);
        for (const auto& [sym, exp] : term.factors)
            for (std::size_t i = 0; i < vars.size(); ++i)
                if (vars[i].name == sym.name) mono[i] = exp;
        sp[mono] += term.coefficient;
        if (sp[mono].is_zero()) sp.erase(mono);
    }
    return sp;
}

[[nodiscard]] MultivariatePolynomial from_sparse(const BSparsePoly& sp,
                                                 const std::vector<Symbol>& vars) {
    std::vector<MultivariateTerm> terms;
    for (const auto& [m, c] : sp) {
        if (c.is_zero()) continue;
        std::vector<std::pair<Symbol, unsigned int>> f;
        for (std::size_t i = 0; i < vars.size(); ++i)
            if (m[i] > 0U) f.emplace_back(vars[i], m[i]);
        terms.push_back(MultivariateTerm{ .coefficient = c, .factors = std::move(f) });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] std::vector<Symbol> collect_vars(const MultivariatePolynomial& p,
                                               const MultivariatePolynomial& q) {
    auto vp = p.variables();
    auto vq = q.variables();
    std::vector<Symbol> all = vp;
    for (const auto& s : vq)
        if (std::none_of(all.begin(), all.end(),
                         [&](const Symbol& a){ return a.name == s.name; }))
            all.push_back(s);
    std::sort(all.begin(), all.end(),
              [](const Symbol& a, const Symbol& b){ return a.name < b.name; });
    return all;
}

[[nodiscard]] BigInt multivar_mignotte_bound(const BSparsePoly& P, const BSparsePoly& Q) {
    std::size_t dP = 0, dQ = 0;
    for (const auto& [m, _] : P) {
        std::size_t td = 0; for (auto e : m) td += e;
        dP = std::max(dP, td);
    }
    for (const auto& [m, _] : Q) {
        std::size_t td = 0; for (auto e : m) td += e;
        dQ = std::max(dQ, td);
    }
    BigInt nP = sparse_inf_norm(P);
    BigInt nQ = sparse_inf_norm(Q);
    BigInt mn = (nP < nQ) ? nP : nQ;
    if (mn.is_zero()) mn = BigInt(1);
    return BigInt(1).shift_left_bits(std::min(dP, dQ)) * mn;
}

[[nodiscard]] bool divides_sparse_z(
    const BSparsePoly& dividend, const BSparsePoly& divisor,
    std::size_t n_vars) {
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

} // namespace cas::algebra::fp_helpers
