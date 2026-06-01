#include "integrate_engine.hpp"
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include <set>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

struct ExprLess {
    bool operator()(ExprPtr lhs, ExprPtr rhs) const noexcept {
        return symbolic::canonical_compare(lhs, rhs) < 0;
    }
};

// Internal helper to collect potential u = g(x) candidates.
// Candidates are arguments of functions or bases of powers.
void collect_substitution_candidates(
    ExprPtr expr,
    const Symbol& var,
    std::set<ExprPtr, ExprLess>& candidates) {
    if (!expr || !integrate_detail::depends_on(expr, var)) {
        return;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        // The FuncCall itself is a u-substitution candidate: integrals of
        // the form ∫ f(g(x))·g'(x) dx benefit from u = g(x) where g may
        // itself be a function call (e.g. ∫ ln(x)/x dx with u = ln(x),
        // du = dx/x).  Add the call as candidate, plus its arguments
        // (for the "outer is f, inner is g" pattern), and recurse.
        candidates.insert(expr);
        for (const auto& arg : call->args) {
            if (integrate_detail::depends_on(arg, var) && !expr_is<Symbol>(arg)) {
                candidates.insert(arg);
            }
            collect_substitution_candidates(arg, var, candidates);
        }
    } else if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) {
             if (integrate_detail::depends_on(bin->left, var) && !expr_is<Symbol>(bin->left)) {
                 candidates.insert(bin->left);
             }
        }
        collect_substitution_candidates(bin->left, var, candidates);
        collect_substitution_candidates(bin->right, var, candidates);
    } else if (const auto* p = expr_cast<Product>(expr)) {
        for (const auto& f : p->factors) {
            collect_substitution_candidates(f, var, candidates);
        }
    } else if (const auto* s = expr_cast<Sum>(expr)) {
        for (const auto& t : s->terms) {
            collect_substitution_candidates(t, var, candidates);
        }
    } else if (const auto* u = expr_cast<Unary>(expr)) {
        collect_substitution_candidates(u->operand, var, candidates);
    }
}

// General expression substituter: replaces occurrences of 'pattern' with 'replacement'.
ExprPtr replace_expr(ExprPtr expr, ExprPtr pattern, ExprPtr replacement, AstArena& arena) {
    if (structural_equal(expr, pattern)) {
        return replacement;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        std::vector<ExprPtr> args;
        args.reserve(call->args.size());
        for (auto arg : call->args) args.push_back(replace_expr(arg, pattern, replacement, arena));
        return arena.make<FuncCall>(call->name, std::move(args));
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        return arena.make<Binary>(bin->op, 
            replace_expr(bin->left, pattern, replacement, arena),
            replace_expr(bin->right, pattern, replacement, arena));
    }
    if (const auto* p = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> factors;
        factors.reserve(p->factors.size());
        for (auto f : p->factors) factors.push_back(replace_expr(f, pattern, replacement, arena));
        return arena.make<Product>(std::move(factors));
    }
    if (const auto* s = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(s->terms.size());
        for (auto t : s->terms) terms.push_back(replace_expr(t, pattern, replacement, arena));
        return arena.make<Sum>(std::move(terms));
    }
    if (const auto* u = expr_cast<Unary>(expr)) {
        return arena.make<Unary>(u->op, replace_expr(u->operand, pattern, replacement, arena));
    }

    return expr;
}

} // anonymous namespace

Result<std::optional<ExprPtr>> integrate_by_substitution(
    const ExprPtr& integrand,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    
    // CAS-L2-16: Automated Variable Substitution (u-substitution recognition).
    // Strategy: find g(x) such that integrand = f(g(x)) * g'(x).
    
    std::set<ExprPtr, ExprLess> candidates;
    collect_substitution_candidates(integrand, var, candidates);
    
    // HC-004: Fresh symbol for u to avoid collisions in the symbolic workspace.
    const Symbol u_sym = ctx.make_fresh_symbol("u");
    const ExprPtr u_expr = ctx.arena().make<Symbol>(u_sym);
    
    for (const auto& g : candidates) {
        // 1. Compute du/dx = g'(x)
        auto dg_res = diff(g, var, 1U, ctx);
        if (dg_res.is_error()) continue;
        ExprPtr dg = dg_res.value();
        
        // Skip constants (handled by linearity/integrate_once)
        if (integrate_detail::is_rational_value(dg, 0, 1)) continue;
        
        // 2. Candidate integrand in terms of u: f(u) = integrand / g'(x)
        auto ratio = ctx.simplify(integrate_detail::make_binary(ctx.arena(), BinaryOp::Div, integrand, dg));
        if (ratio.is_error()) continue;
        
        // 3. Check if f(u) depends on x ONLY through g(x).
        ExprPtr f_u_raw = replace_expr(ratio.value(), g, u_expr, ctx.arena());
        auto f_u = ctx.simplify(f_u_raw);
        if (f_u.is_error()) continue;
        
        // If f_u no longer depends on x, then integrand = f(g(x)) * g'(x).
        if (!integrate_detail::depends_on(f_u.value(), var)) {
            // Success! Integrate f(u) du.
            auto primitive_u = integrate(f_u.value(), u_sym, ctx);
            if (primitive_u.is_ok()) {
                // Back-substitute u -> g(x) to get the result in terms of x.
                auto result = replace_expr(primitive_u.value(), u_expr, g, ctx.arena());
                
                // Mandatory Verification: diff(∫f dx, x) == f
                auto verification = diff(result, var, 1U, ctx);
                if (verification.is_ok()) {
                     auto s_diff = ctx.simplify(verification.value());
                     auto s_integrand = ctx.simplify(integrand);
                     if (s_diff.is_ok() && s_integrand.is_ok() && 
                         structural_equal(s_diff.value(), s_integrand.value())) {
                         return ok(std::make_optional(result));
                     }
                }
            }
        }
    }
    
    return ok(std::optional<ExprPtr>{std::nullopt});
}

namespace integrate_detail {

Result<ExprPtr> Integrator::try_u_substitution_for_product(const Product& product, const Symbol& var) {
    // Attempt general u-substitution for products.
    auto res = integrate_by_substitution(context_.arena().make<Product>(product), var, context_);
    if (res.is_ok() && res.value().has_value()) {
        return ok(res.value().value());
    }
    
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "No supported u-substitution pattern found"));
}

Result<ExprPtr> Integrator::integrate_via_partial_fractions(ExprPtr expr, const Symbol& var) {
    auto terms = algebra::partial_fractions(expr, var, context_);
    if (terms.is_error()) {
        return fail<ExprPtr>(terms.error());
    }

    if (terms.value().empty()) {
        return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Partial fraction decomposition produced no terms"));
    }

    if (terms.value().size() <= 1U) {
         return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Partial fractions did not decompose the expression"));
    }

    std::vector<ExprPtr> primitives;
    primitives.reserve(terms.value().size());
    for (ExprPtr term : terms.value()) {
        auto simplified_term = context_.simplify(term);
        if (simplified_term.is_error()) return fail<ExprPtr>(simplified_term.error());
        
        auto primitive = integrate_once(simplified_term.value(), var);
        if (primitive.is_error()) return primitive;
        primitives.push_back(primitive.value());
    }
    return ok(make_sum(arena_, std::move(primitives)));
}

}  // namespace integrate_detail

}  // namespace cas::calculus
