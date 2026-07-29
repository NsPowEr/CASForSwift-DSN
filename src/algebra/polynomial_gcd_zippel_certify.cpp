// polynomial_gcd_zippel_certify.cpp — certificates for the Zippel-Prony GCD (F3.1).
//
// Split out of polynomial_gcd_zippel_prony.cpp (anti-monolith, 500-line limit).
// Holds everything that validates a candidate produced by the sparse
// interpolation core, plus the sparse Z-polynomial bridge both the single-prime
// and the CRT path share:
//   * to_sparse_z / from_sparse_z — MultivariatePolynomial <-> exponent-vector map
//   * certify_divides             — exact division g | a in Z[x_1..x_n]
//   * is_maximal_gcd_candidate    — cofactor coprimality (A37/A40 soundness)
//   * finish_if_maximal           — sign normalization gated on maximality
//
// Ref: Zippel EUROSAM 1979; Geddes-Czapor-Labahn §7.6 (leading-coefficient
// problem in sparse GCD).

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include "polynomial_gcd_multivariate_helpers.hpp"
#include "polynomial_gcd_zippel_internal.hpp"

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra::zippel_detail {

ZSparsePoly to_sparse_z(const MultivariatePolynomial& p,
                         const std::vector<Symbol>& vars) {
    ZSparsePoly sp;
    for (const auto& term : p.terms()) {
        ZMonomial mono(vars.size(), 0U);
        for (const auto& [sym, exp] : term.factors)
            for (std::size_t i = 0; i < vars.size(); ++i)
                if (vars[i].name == sym.name) mono[i] = exp;
        sp[mono] += term.coefficient;
        if (sp[mono].is_zero()) sp.erase(mono);
    }
    return sp;
}

MultivariatePolynomial from_sparse_z(const ZSparsePoly& sp,
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

// Exact certificate: does b divide a in Z[x_1..x_n]? (sparse leading-term division)
bool certify_divides(const MultivariatePolynomial& a, const MultivariatePolynomial& b,
                     const std::vector<Symbol>& vars) {
    ZSparsePoly sa = to_sparse_z(a, vars);
    ZSparsePoly sb = to_sparse_z(b, vars);
    if (sb.empty()) return false;
    if (sa.empty()) return true;
    const std::size_t nv = vars.size();
    ZSparsePoly rem = sa;
    auto [dlm, dlc] = *std::prev(sb.end());
    const std::size_t budget = (rem.size() + 1U) * (sb.size() + 1U) + 16U;
    std::size_t steps = 0;
    while (!rem.empty()) {
        if (++steps > budget) return false;
        auto [rlm, rlc] = *std::prev(rem.end());
        for (std::size_t i = 0; i < nv; ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            if (ri < di) return false;
        }
        if ((rlc % dlc) != BigInt(0)) return false;
        BigInt qc = rlc / dlc;
        ZMonomial qm(nv, 0U);
        for (std::size_t i = 0; i < nv; ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            qm[i] = static_cast<unsigned int>(ri - di);
        }
        for (const auto& [dm, dc] : sb) {
            ZMonomial nm(nv, 0U);
            for (std::size_t i = 0; i < nv; ++i)
                nm[i] = qm[i] + ((i < dm.size()) ? dm[i] : 0U);
            rem[nm] -= qc * dc;
            if (rem[nm].is_zero()) rem.erase(nm);
        }
    }
    return true;
}

// A37: certify_divides only proves g_cand divides P and Q — a proper factor of
// the true GCD also divides both and would pass that check. The gap is
// structural: the per-sample Fp specialization in the Prony core forces every
// univariate gcd MONIC, which silently discards LC_{x1}(gcd(P,Q)) whenever that
// leading coefficient is a non-constant polynomial in the other variables.
// Phase 2's Prony interpolation then reconstructs a constant-1 leading
// coefficient instead of the true one, undercounting shared factors (e.g.
// gcd(x*y^2-x, x*y-x) = x*(y-1), but every Fp-monic sample collapses to "x"
// regardless of the evaluation point, so Prony faithfully — and wrongly —
// interpolates the constant polynomial 1). Maximality is certified
// independently via the cofactors: g is maximal iff gcd(P/g, Q/g) = 1, checked
// with gcd_multivariate_eval_interp — an exact Z evaluation/interpolation
// algorithm with no per-sample Fp-monic step, so it is not subject to the same
// bias.
[[nodiscard]] Result<bool> is_maximal_gcd_candidate(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    const MultivariatePolynomial& g_cand,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx) {
    auto cof_p = exact_quotient(P, g_cand, vars, ctx);
    if (cof_p.is_error()) {
        return fail<bool>(cof_p.error());
    }
    if (!cof_p.value().has_value()) {
        return fail<bool>(make_error(CASErrorKind::Unimplemented,
            "gcd_zippel_prony: maximality check — cofactor P/g not exact"));
    }
    auto cof_q = exact_quotient(Q, g_cand, vars, ctx);
    if (cof_q.is_error()) {
        return fail<bool>(cof_q.error());
    }
    if (!cof_q.value().has_value()) {
        return fail<bool>(make_error(CASErrorKind::Unimplemented,
            "gcd_zippel_prony: maximality check — cofactor Q/g not exact"));
    }
    auto cofactor_gcd = gcd_multivariate_eval_interp(cof_p.value().value(), cof_q.value().value(), ctx);
    if (cofactor_gcd.is_error()) {
        return fail<bool>(cofactor_gcd.error());
    }
    return ok(is_unit_polynomial(cofactor_gcd.value()));
}

// Sign-normalizes g_cand (lex-leading positive) and returns it only once
// maximality is certified (see is_maximal_gcd_candidate above); otherwise
// reports Unimplemented so the caller (gcd_zippel_sparse) falls through to
// gcd_brown_modular, which has no per-sample Fp-monic step and is therefore
// not exposed to the same bias.
[[nodiscard]] Result<MultivariatePolynomial> finish_if_maximal(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    MultivariatePolynomial g_cand,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx,
    std::size_t* out_samples_used,
    std::size_t samples_used,
    const char* stage) {
    auto maximal = is_maximal_gcd_candidate(P, Q, g_cand, vars, ctx);
    if (!maximal.is_ok() || !maximal.value()) {
        if (out_samples_used) *out_samples_used = samples_used;
        return make_unimplemented<MultivariatePolynomial>(
            "algebra", "gcd_zippel_prony", std::string("stage=") + stage,
            "ZIPPEL_PRONY_NOT_MAXIMAL",
            "Candidate divides P and Q but fails the cofactor-coprimality "
            "maximality certificate (A37, Fp-monic leading-coefficient bias) "
            "— fall back to gcd_brown_modular", "F3.1");
    }

    auto sp = to_sparse_z(g_cand, vars);
    if (!sp.empty()) {
        auto last = std::prev(sp.end());
        if (last->second.is_negative())
            for (auto& [_, c] : sp) c = -c;
        g_cand = from_sparse_z(sp, vars);
    }
    if (out_samples_used) *out_samples_used = samples_used;
    return ok(std::move(g_cand));
}

}  // namespace cas::algebra::zippel_detail
