// F3.2 — Conversion between AST expressions and the internal MPoly.

#include "factor_multivariate_internal.hpp"

#include <algorithm>
#include <optional>

namespace cas::algebra {

namespace {

[[nodiscard]] std::optional<std::size_t> var_index_of(
    const WangContext& wc, const std::string& name) {
    for (std::size_t i = 0; i < wc.vars.size(); ++i) {
        if (wc.vars[i].name == name) {
            return i;
        }
    }
    return std::nullopt;
}

}  // namespace

Result<MPoly> mpoly_from_expr(ExprPtr expr, const WangContext& wc, symbolic::CASContext& ctx) {
    auto parsed = parse_multivariate_polynomial(expr, ctx);
    if (parsed.is_error()) {
        return fail<MPoly>(parsed.error());
    }
    MPoly result;
    const std::size_t n = wc.nvars();
    for (const auto& term : parsed.value().terms()) {
        Monomial mono(n, 0U);
        for (const auto& [sym, exp] : term.factors) {
            auto idx = var_index_of(wc, sym.name);
            if (!idx.has_value()) {
                return fail<MPoly>(make_error(
                    CASErrorKind::Unimplemented,
                    "factor_multivariate: variable '" + sym.name +
                        "' not present in declared variable ordering"));
            }
            mono[*idx] += exp;
        }
        BigInt& slot = result.terms[mono];
        slot += term.coefficient;
        if (slot.is_zero()) {
            result.terms.erase(mono);
        }
    }
    return ok(std::move(result));
}

Result<ExprPtr> mpoly_to_expr(const MPoly& p, const WangContext& wc, symbolic::CASContext& ctx) {
    if (p.is_zero()) {
        return ok(make_integer(ctx.arena(), 0));
    }
    std::vector<MultivariateTerm> terms;
    terms.reserve(p.terms.size());
    for (const auto& [mono, coeff] : p.terms) {
        MultivariateTerm t;
        t.coefficient = coeff;
        for (std::size_t i = 0; i < mono.size(); ++i) {
            if (mono[i] > 0U) {
                t.factors.emplace_back(wc.vars[i], mono[i]);
            }
        }
        terms.push_back(std::move(t));
    }
    MultivariatePolynomial mv(std::move(terms));
    return multivariate_to_expr(mv, ctx);
}

}  // namespace cas::algebra
