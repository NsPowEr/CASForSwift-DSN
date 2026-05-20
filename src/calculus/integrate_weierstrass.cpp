// CAS-L2-14 — Weierstrass substitution for ∫ R(sin x, cos x) dx.
//
// Substitution: t = tan(x/2), so
//   sin(x) = 2t / (1 + t²)
//   cos(x) = (1 - t²) / (1 + t²)
//   dx     = 2 / (1 + t²) dt
//
// Transforms any rational function of sin(x), cos(x) into a rational
// function of t, integrable via the existing rational pipeline.
// Back-substitution t = tan(x/2) recovers the closed form in x.
//
// Anti-furbizia: this is the GENERAL algorithm for R(sin, cos) — not a
// table of patterns. The dispatcher only checks that every trig
// sub-expression depending on var is sin(var) or cos(var), then
// delegates to the rational integrator.

#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"

#include <vector>

namespace cas::calculus {

namespace {

// Tiny structural check: x contains Symbol(var).
[[nodiscard]] bool contains_var(ExprPtr x, const Symbol& v) {
    if (!x) return false;
    if (auto* s = expr_cast<Symbol>(x)) return s->name == v.name;
    if (auto* un = expr_cast<Unary>(x)) return contains_var(un->operand, v);
    if (auto* bin = expr_cast<Binary>(x)) return contains_var(bin->left, v) || contains_var(bin->right, v);
    if (auto* sum = expr_cast<Sum>(x)) { for (auto t : sum->terms) if (contains_var(t, v)) return true; return false; }
    if (auto* prod = expr_cast<Product>(x)) { for (auto t : prod->factors) if (contains_var(t, v)) return true; return false; }
    if (auto* fc = expr_cast<FuncCall>(x)) { for (auto t : fc->args) if (contains_var(t, v)) return true; return false; }
    return false;
}

// Walk AST, return true if every FuncCall depending on var is sin(var)
// or cos(var). Rejects exp/ln/tan/asin/etc. on var. Constants and
// other-variable sub-expressions are fine.
[[nodiscard]] bool is_weierstrass_candidate_walk(ExprPtr expr, const Symbol& var) {
    if (!expr) return true;
    if (auto* call = expr_cast<FuncCall>(expr)) {
        // Trig call: must be sin or cos of exactly Symbol(var).
        const BuiltinOp f = call->func_id;
        bool depends_on_var = false;
        for (ExprPtr a : call->args)
            if (contains_var(a, var)) { depends_on_var = true; break; }
        if (depends_on_var) {
            // Sin/Cos with arg == var allowed; recurse into arg.
            const bool is_sin_or_cos = (f == BuiltinOp::Sin || f == BuiltinOp::Cos);
            if (!is_sin_or_cos) return false;
            if (call->args.size() != 1U) return false;
            const auto* sym = expr_cast<Symbol>(call->args[0]);
            if (!sym || sym->name != var.name) return false;
            return true;
        }
        // Function call independent of var — allowed (constant).
        return true;
    }
    if (auto* un = expr_cast<Unary>(expr)) {
        return is_weierstrass_candidate_walk(un->operand, var);
    }
    if (auto* bin = expr_cast<Binary>(expr)) {
        return is_weierstrass_candidate_walk(bin->left, var)
            && is_weierstrass_candidate_walk(bin->right, var);
    }
    if (auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr t : sum->terms)
            if (!is_weierstrass_candidate_walk(t, var)) return false;
        return true;
    }
    if (auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors)
            if (!is_weierstrass_candidate_walk(f, var)) return false;
        return true;
    }
    return true;  // literals, symbols, constants: OK
}

// Replace FuncCall(Sin, [var]) → sin_repl and FuncCall(Cos, [var]) → cos_repl
// recursively, preserving structural sharing on unchanged sub-trees.
[[nodiscard]] ExprPtr substitute_trig(
    ExprPtr expr, const Symbol& var, ExprPtr sin_repl, ExprPtr cos_repl,
    AstArena& arena) {
    if (!expr) return expr;
    if (auto* call = expr_cast<FuncCall>(expr)) {
        if (call->args.size() == 1U) {
            const auto* sym = expr_cast<Symbol>(call->args[0]);
            if (sym && sym->name == var.name) {
                if (call->func_id == BuiltinOp::Sin) return sin_repl;
                if (call->func_id == BuiltinOp::Cos) return cos_repl;
            }
        }
        // Recurse into args (in case nested).
        std::vector<ExprPtr> new_args;
        new_args.reserve(call->args.size());
        bool changed = false;
        for (ExprPtr a : call->args) {
            ExprPtr a2 = substitute_trig(a, var, sin_repl, cos_repl, arena);
            if (a2 != a) changed = true;
            new_args.push_back(a2);
        }
        if (!changed) return expr;
        return arena.make<FuncCall>(call->func_id, std::move(new_args));
    }
    if (auto* un = expr_cast<Unary>(expr)) {
        ExprPtr inner = substitute_trig(un->operand, var, sin_repl, cos_repl, arena);
        if (inner == un->operand) return expr;
        return arena.make<Unary>(un->op, inner);
    }
    if (auto* bin = expr_cast<Binary>(expr)) {
        ExprPtr l = substitute_trig(bin->left, var, sin_repl, cos_repl, arena);
        ExprPtr r = substitute_trig(bin->right, var, sin_repl, cos_repl, arena);
        if (l == bin->left && r == bin->right) return expr;
        return arena.make<Binary>(bin->op, l, r);
    }
    if (auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(sum->terms.size());
        bool changed = false;
        for (ExprPtr t : sum->terms) {
            ExprPtr t2 = substitute_trig(t, var, sin_repl, cos_repl, arena);
            if (t2 != t) changed = true;
            terms.push_back(t2);
        }
        if (!changed) return expr;
        return arena.make<Sum>(std::move(terms));
    }
    if (auto* prod = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> factors;
        factors.reserve(prod->factors.size());
        bool changed = false;
        for (ExprPtr f : prod->factors) {
            ExprPtr f2 = substitute_trig(f, var, sin_repl, cos_repl, arena);
            if (f2 != f) changed = true;
            factors.push_back(f2);
        }
        if (!changed) return expr;
        return arena.make<Product>(std::move(factors));
    }
    return expr;
}

}  // namespace

// Public entry. Attempts Weierstrass substitution for ∫ expr dx.
// Returns Unimplemented if integrand is not in the rational(sin,cos) class.
Result<ExprPtr> integrate_weierstrass(ExprPtr expr, const Symbol& var,
                                      symbolic::CASContext& ctx) {
    if (!is_weierstrass_candidate_walk(expr, var)) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::Unimplemented,
            "Weierstrass: integrand contains non-sin/cos trig of var",
            std::nullopt});
    }
    AstArena& arena = ctx.arena();

    Symbol t = ctx.make_fresh_symbol("w");
    ExprPtr t_sym = arena.make<Symbol>(t.name);
    ExprPtr two = arena.make<IntegerLit>(BigInt(2));
    ExprPtr one = arena.make<IntegerLit>(BigInt(1));
    ExprPtr t_sq = arena.make<Binary>(BinaryOp::Pow, t_sym, two);
    ExprPtr one_plus_t_sq = arena.make<Sum>(std::vector<ExprPtr>{one, t_sq});
    ExprPtr one_minus_t_sq = arena.make<Sum>(std::vector<ExprPtr>{
        arena.make<IntegerLit>(BigInt(1)),
        arena.make<Unary>(UnaryOp::Neg, t_sq)});
    ExprPtr sin_repl = arena.make<Binary>(BinaryOp::Div,
        arena.make<Product>(std::vector<ExprPtr>{
            arena.make<IntegerLit>(BigInt(2)), t_sym}),
        one_plus_t_sq);
    ExprPtr cos_repl = arena.make<Binary>(BinaryOp::Div,
        one_minus_t_sq, one_plus_t_sq);

    ExprPtr substituted = substitute_trig(expr, var, sin_repl, cos_repl, arena);
    // Jacobian dx/dt = 2/(1 + t²)
    ExprPtr jacobian = arena.make<Binary>(BinaryOp::Div,
        arena.make<IntegerLit>(BigInt(2)), one_plus_t_sq);
    ExprPtr in_t = arena.make<Product>(std::vector<ExprPtr>{substituted, jacobian});
    // Aggressive normalization: together() folds nested rational
    // (1+(1-t²)/(1+t²)) into single fraction so downstream integrator
    // sees a properly-formed rational(t) instead of nested mess.
    auto in_t_tog = algebra::together(in_t, ctx);
    if (in_t_tog.is_ok()) in_t = in_t_tog.value();
    auto in_t_simp = ctx.simplify(in_t);
    if (in_t_simp.is_ok()) in_t = in_t_simp.value();

    // Decompose into clean (numerator/denominator) so downstream
    // integrate_linear_over_quadratic / partial_fractions can pattern-
    // match against polynomials directly. Without this, in_t may stay
    // as a Product of mixed Pow factors that breaks dispatch.
    if (auto parts = algebra::apart_num_den(in_t, ctx); parts.is_ok()) {
        ExprPtr num = parts.value().numerator;
        ExprPtr den = parts.value().denominator;
        auto num_simp = ctx.simplify(num);
        auto den_simp = ctx.simplify(den);
        if (num_simp.is_ok()) num = num_simp.value();
        if (den_simp.is_ok()) den = den_simp.value();
        in_t = arena.make<Binary>(BinaryOp::Div, num, den);
    }

    auto integ = calculus::integrate(in_t, t, ctx);
    if (integ.is_error()) {
        return fail<ExprPtr>(CASError{
            CASErrorKind::Unimplemented,
            "Weierstrass: rational integration in t failed: "
                + integ.error().message,
            std::nullopt});
    }

    // Back-substitute t = tan(var/2)
    ExprPtr var_expr = arena.make<Symbol>(var);
    ExprPtr half_var = arena.make<Binary>(BinaryOp::Div, var_expr, two);
    ExprPtr tan_half = arena.make<FuncCall>(BuiltinOp::Tan,
        std::vector<ExprPtr>{half_var});
    auto back = symbolic::substitute(integ.value(), t, tan_half, ctx);
    if (back.is_error()) return back;
    return ctx.simplify(back.value());
}

}  // namespace cas::calculus
