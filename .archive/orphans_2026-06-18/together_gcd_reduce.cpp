// together_gcd_reduce.cpp — polynomial GCD content reduction helper with support for non-polynomial bases.
// Extracted from factorization_num_den.cpp to respect the 500-line anti-monolith limit.

#include "cas/algebra.hpp"
#include "cas/ast_nodes.hpp"
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <set>
#include <map>
#include <string>
#include <vector>

namespace cas::algebra {

namespace {

void collect_symbol_names(ExprPtr expr, std::set<std::string>& out) {
    if (!expr) return;
    if (const auto* s = expr_cast<Symbol>(expr)) {
        out.insert(s->name);
        return;
    }
    if (const auto* u = expr_cast<Unary>(expr)) {
        collect_symbol_names(u->operand, out);
        return;
    }
    if (const auto* b = expr_cast<Binary>(expr)) {
        collect_symbol_names(b->left, out);
        collect_symbol_names(b->right, out);
        return;
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr t : sum->terms) collect_symbol_names(t, out);
        return;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) collect_symbol_names(f, out);
        return;
    }
    if (const auto* fc = expr_cast<FuncCall>(expr)) {
        for (ExprPtr a : fc->args) collect_symbol_names(a, out);
        return;
    }
}

void collect_non_polynomial_bases(ExprPtr expr, std::vector<ExprPtr>& bases) {
    if (!expr) return;
    
    // Check if the node itself is in bases (using pointer equality since AST has hash-consing)
    if (std::find(bases.begin(), bases.end(), expr) != bases.end()) {
        return;
    }

    if (expr_cast<IntegerLit>(expr) || expr_cast<RationalLit>(expr) || expr_cast<DecimalLit>(expr) || expr_cast<Symbol>(expr)) {
        return;
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            collect_non_polynomial_bases(unary->operand, bases);
            return;
        }
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Add || binary->op == BinaryOp::Sub || binary->op == BinaryOp::Mul) {
            collect_non_polynomial_bases(binary->left, bases);
            collect_non_polynomial_bases(binary->right, bases);
            return;
        }
        if (binary->op == BinaryOp::Pow) {
            auto exp_res = poly_parse_nonnegative_integer_exponent(binary->right);
            if (exp_res.is_ok()) {
                collect_non_polynomial_bases(binary->left, bases);
                return;
            }
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr term : sum->terms) {
            collect_non_polynomial_bases(term, bases);
        }
        return;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            collect_non_polynomial_bases(factor, bases);
        }
        return;
    }

    // It's a non-polynomial base
    bases.push_back(expr);
}

Result<ExprPtr> substitute_multiple(ExprPtr expr, const std::map<ExprPtr, ExprPtr>& replacement_map, symbolic::CASContext& ctx) {
    if (!expr) return ok(expr);

    auto it = replacement_map.find(expr);
    if (it != replacement_map.end()) {
        return ok(it->second);
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        auto operand_res = substitute_multiple(unary->operand, replacement_map, ctx);
        if (operand_res.is_error()) return operand_res;
        if (operand_res.value() == unary->operand) return ok(expr);
        return ok(ctx.arena().make<Unary>(unary->op, operand_res.value()));
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        auto left_res = substitute_multiple(binary->left, replacement_map, ctx);
        if (left_res.is_error()) return left_res;
        auto right_res = substitute_multiple(binary->right, replacement_map, ctx);
        if (right_res.is_error()) return right_res;
        if (left_res.value() == binary->left && right_res.value() == binary->right) return ok(expr);
        return ok(ctx.arena().make<Binary>(binary->op, left_res.value(), right_res.value()));
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        bool changed = false;
        for (ExprPtr term : sum->terms) {
            auto term_res = substitute_multiple(term, replacement_map, ctx);
            if (term_res.is_error()) return term_res;
            terms.push_back(term_res.value());
            if (term_res.value() != term) changed = true;
        }
        if (!changed) return ok(expr);
        return ok(ctx.arena().make<Sum>(std::move(terms)));
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> factors;
        bool changed = false;
        for (ExprPtr factor : product->factors) {
            auto factor_res = substitute_multiple(factor, replacement_map, ctx);
            if (factor_res.is_error()) return factor_res;
            factors.push_back(factor_res.value());
            if (factor_res.value() != factor) changed = true;
        }
        if (!changed) return ok(expr);
        return ok(ctx.arena().make<Product>(std::move(factors)));
    }

    if (const auto* func = expr_cast<FuncCall>(expr)) {
        std::vector<ExprPtr> args;
        bool changed = false;
        for (ExprPtr arg : func->args) {
            auto arg_res = substitute_multiple(arg, replacement_map, ctx);
            if (arg_res.is_error()) return arg_res;
            args.push_back(arg_res.value());
            if (arg_res.value() != arg) changed = true;
        }
        if (!changed) return ok(expr);
        return ok(ctx.arena().make<FuncCall>(func->name, std::move(args)));
    }

    return ok(expr);
}

Result<ExprPtr> rebuild_rationalized(const PolyExpr& poly, ExprPtr S, ExprPtr U, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> terms;
    for (std::size_t k = 0; k < poly.size(); ++k) {
        ExprPtr coeff = poly[k];
        if (!coeff || is_zero_expr(coeff)) continue;
        
        ExprPtr power = nullptr;
        std::size_t m = k / 2;
        if (m > 0) {
            auto pow_res = pow_expr(U, m, ctx);
            if (pow_res.is_error()) return fail<ExprPtr>(pow_res.error());
            power = pow_res.value();
        }
        
        if (k % 2 == 1) {
            if (power) {
                auto mul_res = multiply_exprs(power, S, ctx);
                if (mul_res.is_error()) return fail<ExprPtr>(mul_res.error());
                power = mul_res.value();
            } else {
                power = S;
            }
        }
        
        ExprPtr term = coeff;
        if (power) {
            auto mul_res = multiply_exprs(coeff, power, ctx);
            if (mul_res.is_error()) return fail<ExprPtr>(mul_res.error());
            term = mul_res.value();
        }
        terms.push_back(term);
    }
    if (terms.empty()) return ok(make_integer(ctx.arena(), 0));
    if (terms.size() == 1) return ok(terms[0]);
    return ok(ctx.arena().make<Sum>(std::move(terms)));
}

Result<ReducedRational> try_rationalize_sqrt(ExprPtr N, ExprPtr D, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> bases;
    collect_non_polynomial_bases(D, bases);
    
    for (ExprPtr base : bases) {
        const auto* func = expr_cast<FuncCall>(base);
        if (!func || func->name != "sqrt" || func->args.size() != 1) {
            continue;
        }
        
        ExprPtr U = func->args[0];
        ExprPtr S = base;
        
        ExprPtr y_sym = ctx.arena().make<Symbol>("__temp_y_rationalize");
        std::map<ExprPtr, ExprPtr> fwd;
        fwd[S] = y_sym;
        
        auto N_y_res = substitute_multiple(N, fwd, ctx);
        if (N_y_res.is_error()) continue;
        auto D_y_res = substitute_multiple(D, fwd, ctx);
        if (D_y_res.is_error()) continue;
        
        Symbol y_var("__temp_y_rationalize");
        auto poly_N_res = parse_polynomial(N_y_res.value(), y_var, ctx);
        if (poly_N_res.is_error()) continue;
        auto poly_D_res = parse_polynomial(D_y_res.value(), y_var, ctx);
        if (poly_D_res.is_error()) continue;
        
        const auto& poly_N = poly_N_res.value();
        const auto& poly_D = poly_D_res.value();
        
        if (poly_D.degree() != 1U) {
            continue;
        }
        
        ExprPtr D_parts = poly_D[0];
        ExprPtr C = poly_D[1];
        
        std::vector<ExprPtr> conj_coeffs;
        conj_coeffs.push_back(D_parts);
        auto neg_C = negate_expr(C, ctx);
        if (neg_C.is_error()) continue;
        conj_coeffs.push_back(neg_C.value());
        PolyExpr poly_conj(std::move(conj_coeffs));
        
        auto poly_N_conj_res = poly_multiply(poly_N, poly_conj, ctx);
        if (poly_N_conj_res.is_error()) continue;
        auto poly_D_conj_res = poly_multiply(poly_D, poly_conj, ctx);
        if (poly_D_conj_res.is_error()) continue;
        
        auto N_new_res = rebuild_rationalized(poly_N_conj_res.value(), S, U, ctx);
        if (N_new_res.is_error()) continue;
        auto D_new_res = rebuild_rationalized(poly_D_conj_res.value(), S, U, ctx);
        if (D_new_res.is_error()) continue;
        
        auto N_simp = simplify_expr(N_new_res.value(), ctx);
        if (N_simp.is_error()) continue;
        auto D_simp = simplify_expr(D_new_res.value(), ctx);
        if (D_simp.is_error()) continue;
        
        return ok(ReducedRational{N_simp.value(), D_simp.value()});
    }
    
    return fail<ReducedRational>(make_error(CASErrorKind::Unimplemented, "No sqrt rationalization applied"));
}

} // namespace

Result<ReducedRational> reduce_rational_by_gcd(
    ExprPtr N, ExprPtr D, symbolic::CASContext& ctx)
{
    ReducedRational identity{N, D};
    if (!ctx.together_gcd_enabled()) return ok(identity);
    if (!N || !D) return ok(identity);
    if (is_zero_expr(N) || is_one_expr(D) || is_one_expr(N)) return ok(identity);

    // Try rationalizing sqrt first
    auto rationalized = try_rationalize_sqrt(N, D, ctx);
    if (rationalized.is_ok()) {
        N = rationalized.value().numerator;
        D = rationalized.value().denominator;
    }

    // Collect non-polynomial bases
    std::vector<ExprPtr> bases;
    collect_non_polynomial_bases(N, bases);
    collect_non_polynomial_bases(D, bases);

    // Create substitution maps
    std::map<ExprPtr, ExprPtr> forward_map;
    std::map<ExprPtr, ExprPtr> backward_map;
    for (std::size_t i = 0; i < bases.size(); ++i) {
        std::string var_name = "__alg_var_" + std::to_string(i);
        ExprPtr temp_sym = ctx.arena().make<Symbol>(var_name);
        forward_map[bases[i]] = temp_sym;
        backward_map[temp_sym] = bases[i];

        // If the base is sqrt(U), also map U -> temp_sym^2
        if (const auto* func = expr_cast<FuncCall>(bases[i])) {
            if (func->name == "sqrt" && func->args.size() == 1) {
                ExprPtr U = func->args[0];
                ExprPtr temp_sq = ctx.arena().make<Binary>(
                    BinaryOp::Pow,
                    temp_sym,
                    ctx.arena().make<IntegerLit>(BigInt(2))
                );
                forward_map[U] = temp_sq;
                backward_map[temp_sq] = U;
            }
        }
    }

    // Substitute non-polynomial bases
    auto N_sub_res = substitute_multiple(N, forward_map, ctx);
    if (N_sub_res.is_error()) return ok(identity);
    auto D_sub_res = substitute_multiple(D, forward_map, ctx);
    if (D_sub_res.is_error()) return ok(identity);

    ExprPtr N_sub = N_sub_res.value();
    ExprPtr D_sub = D_sub_res.value();

    // Guard: cap symbol count to avoid catastrophic GCD on wide expressions.
    std::set<std::string> shared;
    collect_symbol_names(N_sub, shared);
    std::set<std::string> dvars;
    collect_symbol_names(D_sub, dvars);
    
    // Intersect → symbols common to both (potential cancellation indeterminates).
    std::vector<std::string> common;
    std::set_intersection(shared.begin(), shared.end(), dvars.begin(), dvars.end(),
                          std::back_inserter(common));
    if (common.empty()) return ok(identity);
    if (common.size() > ctx.together_gcd_max_symbols()) return ok(identity);

    // Compute the multivariate GCD; soft-fail to identity on any error.
    auto g_res = polynomial_gcd_multivariate(N_sub, D_sub, ctx);
    if (g_res.is_error()) return ok(identity);
    ExprPtr g = g_res.value();
    if (is_zero_expr(g) || is_one_expr(g)) return ok(identity);

    // Pick the lexicographically-first shared symbol as the main variable for
    // univariate exact-divide. Any var present in g works because g divides
    // both N and D exactly; we deliberately choose deterministically.
    Symbol main_var(common.front());

    // Guard: skip reduction if degree in main_var blows up past the cap.
    auto deg_g = polynomial_degree(g, main_var, ctx);
    if (deg_g.is_error()) return ok(identity);
    if (deg_g.value() > ctx.together_gcd_max_degree()) return ok(identity);

    auto N_q_res = polynomial_exact_divide(N_sub, g, main_var, ctx);
    if (N_q_res.is_error()) return ok(identity);  // non-exact → fall back unreduced
    auto D_q_res = polynomial_exact_divide(D_sub, g, main_var, ctx);
    if (D_q_res.is_error()) return ok(identity);

    // Substitute back original bases
    auto N_q_final = substitute_multiple(N_q_res.value(), backward_map, ctx);
    if (N_q_final.is_error()) return ok(identity);
    auto D_q_final = substitute_multiple(D_q_res.value(), backward_map, ctx);
    if (D_q_final.is_error()) return ok(identity);

    return ok(ReducedRational{N_q_final.value(), D_q_final.value()});
}

} // namespace cas::algebra
