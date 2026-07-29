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

// Recognizes `arg` as k*w for a nonzero integer k (any sign), or w itself
// (k=1), via structural comparison against w (not name-based like the
// Weierstrass multi-angle helper — w here is an arbitrary expression, not
// necessarily a bare variable).
[[nodiscard]] std::optional<BigInt> integer_multiple_of_expr(ExprPtr arg, ExprPtr w) {
    if (!arg || !w) return std::nullopt;
    if (structural_equal(arg, w)) return BigInt(1);
    if (const auto* un = expr_cast<Unary>(arg)) {
        if (un->op != UnaryOp::Neg) return std::nullopt;
        auto inner = integer_multiple_of_expr(un->operand, w);
        if (!inner.has_value()) return std::nullopt;
        return -(*inner);
    }
    if (const auto* prod = expr_cast<Product>(arg)) {
        if (prod->factors.size() != 2U) return std::nullopt;
        const auto* int_a = expr_cast<IntegerLit>(prod->factors[0]);
        const auto* int_b = expr_cast<IntegerLit>(prod->factors[1]);
        if (int_a && structural_equal(prod->factors[1], w) && !int_a->value.is_zero()) return int_a->value;
        if (int_b && structural_equal(prod->factors[0], w) && !int_b->value.is_zero()) return int_b->value;
        return std::nullopt;
    }
    return std::nullopt;
}

// General expression substituter: replaces occurrences of 'pattern' with
// 'replacement'. When pattern is exp(w), also recognizes exp(k*w) for any
// nonzero integer k (including exp(-w)) as replacement^k — exp(-x) is
// algebraically 1/exp(x), but the two are different ExprPtr shapes, so a
// literal structural_equal match alone would miss it (e.g. candidate
// g=exp(x) in 1/(exp(x)+exp(-x)): without this, exp(-x) never becomes 1/u
// and the candidate looks like it still depends on x).
ExprPtr replace_expr(ExprPtr expr, ExprPtr pattern, ExprPtr replacement, AstArena& arena) {
    if (structural_equal(expr, pattern)) {
        return replacement;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::Exp && call->args.size() == 1U) {
            if (const auto* pattern_call = expr_cast<FuncCall>(pattern);
                pattern_call && pattern_call->func_id == BuiltinOp::Exp && pattern_call->args.size() == 1U) {
                if (auto k = integer_multiple_of_expr(call->args[0], pattern_call->args[0]); k.has_value()) {
                    return arena.make<Binary>(BinaryOp::Pow, replacement, arena.make<IntegerLit>(*k));
                }
            }
        }
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
        
        // 2. Candidate integrand in terms of u: f(u) = integrand / g'(x).
        //    Use the raw quotient (no simplifier pass) so trig-linearisation
        //    identities (sin(x)·cos(x) → sin(2x)/2, etc.) do not destroy the
        //    candidate sub-expression before replacement can lock onto it.
        ExprPtr ratio_raw = integrate_detail::make_binary(ctx.arena(), BinaryOp::Div, integrand, dg);
        // 3. Replace candidate g(x) by u in the raw ratio, then simplify.
        //    Simplifying AFTER the substitution lets `u` act as an opaque
        //    atom so cancellations like sin·u/sin → u proceed without the
        //    multi-angle expansion that fires on (sin·cos² / sin) directly.
        ExprPtr f_u_raw = replace_expr(ratio_raw, g, u_expr, ctx.arena());
        auto f_u = ctx.simplify(f_u_raw);
        if (f_u.is_error()) continue;

        // Fixpoint pass: algebraic simplification of the raw substitution
        // can surface *new* literal occurrences of g that were not visible
        // to the first replace_expr pass (e.g. for x^3*cos(x^2), g=x^2:
        // ratio_raw = x^3*cos(x^2)/(2x) only reduces x^3/(2x) -> x^2/2
        // *during* simplify, after replace_expr already ran, so the
        // emergent x^2 factor is missed and the candidate looks like it
        // still depends on x). Re-run replace+simplify while a candidate
        // still depends on var and each pass keeps changing the result;
        // bounded so a genuinely non-closing candidate cannot loop
        // unboundedly (a converging fixpoint stabilizes in 1-2 passes).
        for (unsigned int pass = 0; pass < ctx.max_substitution_fixpoint_passes()
             && integrate_detail::depends_on(f_u.value(), var); ++pass) {
            ExprPtr next_raw = replace_expr(f_u.value(), g, u_expr, ctx.arena());
            auto next = ctx.simplify(next_raw);
            if (next.is_error() || structural_equal(next.value(), f_u.value())) break;
            f_u = next;
        }

        // If f_u no longer depends on x, then integrand = f(g(x)) * g'(x).
        if (!integrate_detail::depends_on(f_u.value(), var)) {
            auto primitive_u = integrate(f_u.value(), u_sym, ctx);
            if (primitive_u.is_ok()) {
                // Back-substitute u -> g(x) to get the result in terms of x.
                auto result = replace_expr(primitive_u.value(), u_expr, g, ctx.arena());
                
                // Mandatory Verification: diff(∫f dx, x) == f. Compare via
                // together(D(result) - integrand) simplifying to the zero
                // literal, not raw structural_equal on two independently
                // simplified forms: simplify is not confluent on reciprocal
                // shapes (e.g. Pow(Product(x,log(x)),-1) vs
                // Product(Pow(x,-1),Pow(log(x),-1)) for the same value), so a
                // structurally-different-but-equal derivative would be
                // rejected and the correct substitution discarded. Same
                // zero-difference idiom as integrate_by_parts' verification.
                auto verification = diff(result, var, 1U, ctx);
                if (verification.is_ok()) {
                    ExprPtr delta = ctx.arena().make<Binary>(BinaryOp::Sub, verification.value(), integrand);
                    auto delta_tog = algebra::together(delta, ctx);
                    ExprPtr delta_for_simp = delta_tog.is_ok() ? delta_tog.value() : delta;
                    auto delta_simp = ctx.simplify(delta_for_simp);
                    bool is_zero = delta_simp.is_ok()
                        && expr_is<IntegerLit>(delta_simp.value())
                        && expr_ref<IntegerLit>(delta_simp.value()).value.is_zero();
                    if (is_zero) {
                        return ok(std::make_optional(result));
                    }
                }
            }
        }
    }
    
    return ok(std::optional<ExprPtr>{std::nullopt});
}

namespace integrate_detail {

Result<ExprPtr> Integrator::try_u_substitution_for_product_impl(const Product& product, const Symbol& var) {
    // Attempt general u-substitution for products.
    auto res = integrate_by_substitution(context_.arena().make<Product>(product), var, context_);
    if (res.is_ok() && res.value().has_value()) {
        return ok(res.value().value());
    }
    
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "No supported u-substitution pattern found"));
}

Result<ExprPtr> Integrator::integrate_via_partial_fractions_impl(ExprPtr expr, const Symbol& var) {
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
