#include "cas/symbolic.hpp"
#include "cas/ast.hpp"
#include "simplify_impl.hpp"
#include "cas/error_helpers.hpp"
#include <map>
#include <vector>

namespace cas::symbolic::detail {

[[nodiscard]] Result<ExprPtr> trunc_series_add(const SeriesExp& a, const SeriesExp& b, AstArena& arena, CASContext& ctx) {
    std::map<long long, ExprPtr> terms;
    for (const auto& [exp, coeff] : a.terms) terms[exp] = coeff;
    for (const auto& [exp, coeff] : b.terms) {
        if (terms.count(exp)) {
            auto sum = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{terms[exp], coeff}));
            if (sum.is_error()) return sum;
            terms[exp] = sum.value();
        } else {
            terms[exp] = coeff;
        }
    }
    
    std::vector<std::pair<long long, ExprPtr>> result_terms;
    long long order = std::min(a.order, b.order);
    for (auto& [exp, coeff] : terms) {
        if (exp < order) {
            result_terms.push_back({exp, coeff});
        }
    }
    
    return ok(arena.make<SeriesExp>(a.var, a.point, std::move(result_terms), order));
}

[[nodiscard]] Result<ExprPtr> trunc_series_mul(const SeriesExp& a, const SeriesExp& b, AstArena& arena, CASContext& ctx) {
    std::map<long long, ExprPtr> terms;
    long long order = std::min(a.order, b.order);
    
    for (const auto& [exp_a, coeff_a] : a.terms) {
        for (const auto& [exp_b, coeff_b] : b.terms) {
            long long exp = exp_a + exp_b;
            if (exp >= order) continue;
            
            auto prod = ctx.simplify(arena.make<Product>(std::vector<ExprPtr>{coeff_a, coeff_b}));
            if (prod.is_error()) return prod;
            
            if (terms.count(exp)) {
                auto sum = ctx.simplify(arena.make<Sum>(std::vector<ExprPtr>{terms[exp], prod.value()}));
                if (sum.is_error()) return sum;
                terms[exp] = sum.value();
            } else {
                terms[exp] = prod.value();
            }
        }
    }
    
    std::vector<std::pair<long long, ExprPtr>> result_terms;
    for (auto& [exp, coeff] : terms) {
        result_terms.push_back({exp, coeff});
    }
    
    return ok(arena.make<SeriesExp>(a.var, a.point, std::move(result_terms), order));
}

[[nodiscard]] Result<ExprPtr> trunc_series_inv(const SeriesExp& s, AstArena&, CASContext&) {
    if (s.terms.empty()) return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "Cannot invert empty series"));
    
    auto it = s.terms.begin();
    while (it != s.terms.end() && is_zero_expr(it->second)) ++it;
    if (it == s.terms.end()) return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "Cannot invert zero series"));
    
    // F0.8-MIGRATED
    return make_unimplemented<ExprPtr>(
        "symbolic", "trunc_series_inv",
        "truncated power series inversion (non-constant leading term)",
        error::reason_codes::SYMBOLIC_SERIES_INVERSION,
        "Implement full series inversion via Newton iteration on truncated series ring",
        "F1.x");
}

} // namespace cas::symbolic::detail
