#pragma once

// integrate_risch_internal.hpp — Shared helpers for the Risch integration
// pipeline.  This is an INTERNAL header: do NOT add to include/cas/.
//
// Included by:
//   integrate_risch.cpp      (dispatcher + integrate_log_polynomial_part)
//   integrate_risch_rde.cpp  (Q-base Risch DE solvers, Bronstein §5-6)
//
// Both translation units define anonymous namespaces that include these
// helpers; `inline` ensures a single linkage-compatible definition per TU.

#include "calculus_internal.hpp"

#include <cstddef>
#include <vector>

namespace cas::calculus {

namespace risch_helpers {

// Deep structural equality for two AST nodes.  Pointer identity O(1) fast
// path; recursive structural comparison for non-identical pointers.
// Used both by the Q-base RDE solver and by the top-level dispatcher.
[[nodiscard]] inline bool deep_struct_equal(ExprPtr a, ExprPtr b) {
    if (a == b) return true;
    if (!a || !b) return false;
    if (const auto* sa = expr_cast<Symbol>(a))
        if (const auto* sb = expr_cast<Symbol>(b)) return sa->name == sb->name;
    if (const auto* ia = expr_cast<IntegerLit>(a))
        if (const auto* ib = expr_cast<IntegerLit>(b)) return ia->value == ib->value;
    if (const auto* ra = expr_cast<RationalLit>(a))
        if (const auto* rb = expr_cast<RationalLit>(b))
            return ra->numerator == rb->numerator && ra->denominator == rb->denominator;
    if (const auto* ca = expr_cast<Constant>(a))
        if (const auto* cb = expr_cast<Constant>(b)) return ca->value == cb->value;
    if (const auto* fa = expr_cast<FuncCall>(a))
        if (const auto* fb = expr_cast<FuncCall>(b)) {
            if (fa->func_id != fb->func_id || fa->args.size() != fb->args.size()) return false;
            for (std::size_t i = 0; i < fa->args.size(); ++i)
                if (!deep_struct_equal(fa->args[i], fb->args[i])) return false;
            return true;
        }
    if (const auto* ba = expr_cast<Binary>(a))
        if (const auto* bb = expr_cast<Binary>(b))
            return ba->op == bb->op &&
                   deep_struct_equal(ba->left, bb->left) &&
                   deep_struct_equal(ba->right, bb->right);
    if (const auto* ua = expr_cast<Unary>(a))
        if (const auto* ub = expr_cast<Unary>(b))
            return ua->op == ub->op && deep_struct_equal(ua->operand, ub->operand);
    if (const auto* su_a = expr_cast<Sum>(a))
        if (const auto* su_b = expr_cast<Sum>(b)) {
            if (su_a->terms.size() != su_b->terms.size()) return false;
            for (std::size_t i = 0; i < su_a->terms.size(); ++i)
                if (!deep_struct_equal(su_a->terms[i], su_b->terms[i])) return false;
            return true;
        }
    if (const auto* pa = expr_cast<Product>(a))
        if (const auto* pb = expr_cast<Product>(b)) {
            if (pa->factors.size() != pb->factors.size()) return false;
            for (std::size_t i = 0; i < pa->factors.size(); ++i)
                if (!deep_struct_equal(pa->factors[i], pb->factors[i])) return false;
            return true;
        }
    return false;
}

// HPP-007 helper: deep tree walker that substitutes every occurrence of a
// target subexpression with a replacement.  Used in the log-derivative
// recognition pass to replace transcendental generators with fresh symbols
// so polynomial GCD cancellation becomes reachable.
[[nodiscard]] inline ExprPtr deep_replace_expr(
    ExprPtr e, ExprPtr target, ExprPtr replacement, AstArena& arena) {
    if (!e) return e;
    if (deep_struct_equal(e, target)) return replacement;
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        std::vector<ExprPtr> na;
        na.reserve(fc->args.size());
        for (ExprPtr a : fc->args)
            na.push_back(deep_replace_expr(a, target, replacement, arena));
        return arena.make<FuncCall>(fc->func_id, std::move(na));
    }
    if (const auto* bin = expr_cast<Binary>(e))
        return arena.make<Binary>(bin->op,
            deep_replace_expr(bin->left,  target, replacement, arena),
            deep_replace_expr(bin->right, target, replacement, arena));
    if (const auto* un = expr_cast<Unary>(e))
        return arena.make<Unary>(un->op,
            deep_replace_expr(un->operand, target, replacement, arena));
    if (const auto* sum = expr_cast<Sum>(e)) {
        std::vector<ExprPtr> nt;
        nt.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms)
            nt.push_back(deep_replace_expr(t, target, replacement, arena));
        return arena.make<Sum>(std::move(nt));
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        std::vector<ExprPtr> nf;
        nf.reserve(prod->factors.size());
        for (ExprPtr f : prod->factors)
            nf.push_back(deep_replace_expr(f, target, replacement, arena));
        return arena.make<Product>(std::move(nf));
    }
    return e;
}

} // namespace risch_helpers

// Forward declaration for the polynomial-quotient + Hermite/Trager assembly
// step of the Risch pipeline (Bronstein §5-6, steps 3 onwards).
// Defined in integrate_risch_hermite.cpp; called by integrate_risch().
//
// Parameters:
//   gen_expr           — integrand expressed in field generators
//   expr_original      — original integrand (reserved for round-trip checks)
//   var                — integration variable
//   field              — differential extension tower
//   poly_integral_part — antiderivative of the polynomial quotient (may be 0)
//   context            — CAS context
[[nodiscard]] Result<ExprPtr> integrate_risch_poly_and_rational_part(
    ExprPtr gen_expr,
    ExprPtr expr_original,
    const Symbol& var,
    const DifferentialField& field,
    ExprPtr poly_integral_part,
    symbolic::CASContext& context);

// Note: integrate_log_polynomial_part uses algebra::PolyExpr (defined in
// polynomial_internal.hpp, a private header).  Its declaration lives only
// in integrate_risch_hermite.cpp where that header is included.

} // namespace cas::calculus
