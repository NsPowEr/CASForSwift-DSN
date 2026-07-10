#pragma once
// Helpers shared (T-045 anti-monolith split) between summation_abramov.cpp
// (Gosper / polygamma / Abramov drivers) and summation_abramov_quadratic.cpp
// (RootOf-aware quadratic-atom antidifference). Not part of the public CAS API.

#include "cas/builtin_functions.hpp"
#include "cas/symbolic.hpp"
#include "summation_internal.hpp"  // rational_expr + ast/Rational/AstArena

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus::abramov_detail {

[[nodiscard]] inline bool extract_rational_const(ExprPtr e, Rational& out) {
    if (const auto* il = expr_cast<IntegerLit>(e)) {
        out = Rational(il->value);
        return true;
    }
    if (const auto* rl = expr_cast<RationalLit>(e)) {
        out = Rational(rl->numerator, rl->denominator);
        return true;
    }
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (extract_rational_const(un->operand, out)) {
            out = -out;
            return true;
        }
    }
    return false;
}

// Build the antidifference of A/(k+a)^m as A·(−1)^(m−1)/(m−1)!·ψ^(m−1)(k+a).
// All arguments are ExprPtrs — symbolic or algebraic coefficients pass through.
[[nodiscard]] inline ExprPtr polygamma_antidiff(
    ExprPtr A, ExprPtr k_plus_a, unsigned int m, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    BigInt sign = (m % 2U == 1U) ? BigInt(1) : BigInt(-1);
    BigInt fact(1);
    for (unsigned int i = 2U; i < m; ++i)
        fact *= BigInt(static_cast<long long>(i));
    ExprPtr scale_expr = rational_expr(arena, Rational(sign, fact));
    ExprPtr coefficient = arena.make<Binary>(BinaryOp::Mul, A, scale_expr);

    ExprPtr ant;
    if (m == 1U) {
        std::vector<ExprPtr> args{k_plus_a};
        ant = arena.make<FuncCall>(BuiltinOp::Digamma, std::move(args));
    } else {
        std::vector<ExprPtr> args{
            arena.make<IntegerLit>(BigInt(static_cast<long long>(m - 1U))),
            k_plus_a,
        };
        ant = arena.make<FuncCall>(BuiltinOp::Polygamma, std::move(args));
    }
    return arena.make<Binary>(BinaryOp::Mul, coefficient, ant);
}

// Fast structural check: does `expr` contain a Div node or Pow with negative
// exponent?  If not, the polygamma/Abramov path is skipped.
[[nodiscard]] inline bool has_rational_dependency(
    ExprPtr expr, const std::string& k_name) {
    (void)k_name;
    if (!expr) return false;
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Div) return true;
        if (bin->op == BinaryOp::Pow) {
            if (const auto* exp_lit = expr_cast<IntegerLit>(bin->right);
                exp_lit && exp_lit->value.is_negative()) {
                return true;
            }
        }
        return has_rational_dependency(bin->left, k_name) ||
               has_rational_dependency(bin->right, k_name);
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr t : sum->terms) {
            if (has_rational_dependency(t, k_name)) return true;
        }
        return false;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr t : prod->factors) {
            if (has_rational_dependency(t, k_name)) return true;
        }
        return false;
    }
    if (const auto* un = expr_cast<Unary>(expr)) {
        return has_rational_dependency(un->operand, k_name);
    }
    return false;
}

[[nodiscard]] std::optional<ExprPtr> try_polygamma_antidiff(
    ExprPtr term, const Symbol& k, symbolic::CASContext& ctx);

// RootOf-aware digamma/polygamma antidifference for Q-irreducible quadratic
// atoms. Defined in summation_abramov_quadratic.cpp.
[[nodiscard]] std::optional<ExprPtr> try_quadratic_atom_antidiff(
    ExprPtr term, const Symbol& k, symbolic::CASContext& ctx);

// Partial-fraction coefficients of (A1·k + A0)/((k−α)^m (k−β)^m) over Q(α,β):
//   (A1·k + A0)/((k−α)^m(k−β)^m) = Σ_{j=1..m} [ C_j/(k−α)^j + D_j/(k−β)^j ].
// Returned in j order (index 0 ↔ j=1). Pure AST construction, no simplify:
// α and β may be RootOf nodes or plain rationals — the identity is the same,
// which is what makes the formula directly testable with rational roots.
// Defined in summation_abramov_quadratic.cpp (F5.7-B6BIS).
[[nodiscard]] std::vector<std::pair<ExprPtr, ExprPtr>> quadratic_pf_coeffs(
    ExprPtr A0, ExprPtr A1, ExprPtr alpha, ExprPtr beta, unsigned int m,
    AstArena& arena);

}  // namespace cas::calculus::abramov_detail
