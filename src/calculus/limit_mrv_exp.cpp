// limit_mrv_exp.cpp — Exponential-product limit analysis for the Gruntz MRV
// algorithm (Gruntz 1996 §3.5).
//
// Provides:
//   - append_scaled_exponential_term     (file-local)
//   - collect_exponential_product_terms  (file-local)
//   - try_exponential_product_limit      (file-local; wrapper: mrv_try_exponential_product_limit)
//
// Recognises expressions of the form  c · exp(f₁)^a₁ · exp(f₂)^a₂ · …
// and resolves the limit x→+∞ by finding the dominant exponential exponent
// via compare_growth, then reading off the sign of its integer coefficient.

#include "limit_mrv_internal.hpp"

#include <algorithm>
#include <optional>
#include <vector>

namespace cas::calculus {

namespace {

// ---------------------------------------------------------------------------
// append_scaled_exponential_term
//
// Recursion depth bound (1024) is derived from AST nesting depth, not
// arbitrary: each recursive call strictly descends one AST level.  1024 is a
// conservative ceiling that returns false (Unimplemented path) on
// pathologically deep inputs rather than overflowing the call stack.
// See SOURCE NOTE in test/unit/symbolic/test_limit_tower_adaptive.cpp.
// ---------------------------------------------------------------------------
[[nodiscard]] bool append_scaled_exponential_term(
    std::vector<ExponentialTerm>& terms,
    ExprPtr exponent,
    long long coefficient,
    AstArena& arena,
    unsigned int max_depth,
    unsigned int depth = 0U) {
    if (depth >= max_depth) return false;
    if (coefficient == 0) return true;

    if (const auto* unary = expr_cast<Unary>(exponent)) {
        if (unary->op == UnaryOp::Neg) {
            if (coefficient == std::numeric_limits<long long>::min()) return false;
            return append_scaled_exponential_term(
                terms, unary->operand, -coefficient, arena, max_depth, depth + 1U);
        }
    }

    if (const auto* binary = expr_cast<Binary>(exponent)) {
        if (binary->op == BinaryOp::Mul) {
            if (auto lhs_coeff = integer_value(binary->left)) {
                long long scaled{};
                if (!safe_mul_i64(coefficient, *lhs_coeff, scaled)) return false;
                return append_scaled_exponential_term(
                    terms, binary->right, scaled, arena, max_depth, depth + 1U);
            }
            if (auto rhs_coeff = integer_value(binary->right)) {
                long long scaled{};
                if (!safe_mul_i64(coefficient, *rhs_coeff, scaled)) return false;
                return append_scaled_exponential_term(
                    terms, binary->left, scaled, arena, max_depth, depth + 1U);
            }
        }
    }

    if (const auto* product = expr_cast<Product>(exponent)) {
        long long scaled = coefficient;
        std::vector<ExprPtr> symbolic_factors;
        for (ExprPtr factor : product->factors) {
            if (auto factor_coeff = integer_value(factor)) {
                if (!safe_mul_i64(scaled, *factor_coeff, scaled)) return false;
            } else {
                symbolic_factors.push_back(factor);
            }
        }
        if (symbolic_factors.empty()) return true;
        ExprPtr symbolic_exponent = symbolic_factors.size() == 1U
            ? symbolic_factors.front()
            : arena.make<Product>(std::move(symbolic_factors));
        // Guard: if Product reconstruction yields the same pointer, we would
        // recurse forever.  The depth bound provides a second safety net.
        if (symbolic_exponent == exponent) {
            terms.push_back(ExponentialTerm{.exponent = exponent, .coefficient = coefficient});
            return true;
        }
        return append_scaled_exponential_term(
            terms, symbolic_exponent, scaled, arena, max_depth, depth + 1U);
    }

    terms.push_back(ExponentialTerm{.exponent = exponent, .coefficient = coefficient});
    return true;
}

[[nodiscard]] bool collect_exponential_product_terms(
    ExprPtr expr,
    const Symbol& var,
    long long multiplier,
    std::vector<ExponentialTerm>& terms,
    AstArena& arena,
    unsigned int max_depth) {
    if (multiplier == 0) return true;
    if (!expr) return false;

    if (!depends_on(expr, var)) {
        return integer_value(expr).has_value();
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->func_id == BuiltinOp::Exp && call->args.size() == 1U) {
            return append_scaled_exponential_term(
                terms, call->args.front(), multiplier, arena, max_depth);
        }
        return false;
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op != UnaryOp::Neg) return false;
        return collect_exponential_product_terms(
            unary->operand, var, multiplier, terms, arena, max_depth);
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Mul) {
            return collect_exponential_product_terms(
                       binary->left, var, multiplier, terms, arena, max_depth) &&
                   collect_exponential_product_terms(
                       binary->right, var, multiplier, terms, arena, max_depth);
        }
        if (binary->op == BinaryOp::Div) {
            long long denominator_multiplier{};
            if (!safe_mul_i64(multiplier, -1LL, denominator_multiplier)) return false;
            return collect_exponential_product_terms(
                       binary->left, var, multiplier, terms, arena, max_depth) &&
                   collect_exponential_product_terms(
                       binary->right, var, denominator_multiplier, terms, arena, max_depth);
        }
        if (binary->op == BinaryOp::Pow) {
            const auto exponent = integer_value(binary->right);
            if (!exponent.has_value()) return false;
            long long scaled{};
            if (!safe_mul_i64(multiplier, *exponent, scaled)) return false;
            return collect_exponential_product_terms(
                binary->left, var, scaled, terms, arena, max_depth);
        }
        return false;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            if (!collect_exponential_product_terms(factor, var, multiplier, terms, arena, max_depth))
                return false;
        }
        return true;
    }

    return integer_value(expr).has_value();
}

[[nodiscard]] std::optional<Result<ExprPtr>> try_exponential_product_limit(
    ExprPtr expr,
    const Symbol& var,
    AstArena& arena,
    symbolic::CASContext& ctx) {
    std::vector<ExponentialTerm> raw_terms;
    if (!collect_exponential_product_terms(expr, var, 1LL, raw_terms, arena, ctx.mrv_max_append_depth())
        || raw_terms.empty()) {
        return std::nullopt;
    }

    // Merge terms with structurally identical exponents.
    std::vector<ExponentialTerm> terms;
    for (const auto& raw : raw_terms) {
        if (!depends_on(raw.exponent, var)) continue;
        bool merged = false;
        for (auto& term : terms) {
            if (structural_equal(term.exponent, raw.exponent)) {
                term.coefficient += raw.coefficient;
                merged = true;
                break;
            }
        }
        if (!merged) terms.push_back(raw);
    }

    terms.erase(
        std::remove_if(terms.begin(), terms.end(),
            [](const ExponentialTerm& t) { return t.coefficient == 0; }),
        terms.end());
    if (terms.empty()) return std::nullopt;

    // Repeatedly identify and eliminate the dominant exponential exponent.
    while (!terms.empty()) {
        std::size_t dominant = 0U;
        bool undecidable_tie = false;
        for (std::size_t i = 1U; i < terms.size(); ++i) {
            const int cmp = compare_growth(
                terms[i].exponent, terms[dominant].exponent, var, ctx);
            if (cmp > 0) {
                dominant = i;
                undecidable_tie = false;
            } else if (cmp == 0
                       && !structural_equal(terms[i].exponent, terms[dominant].exponent)) {
                undecidable_tie = true;
            }
        }

        long long dominant_coefficient = terms[dominant].coefficient;
        for (std::size_t i = 0U; i < terms.size(); ++i) {
            if (i == dominant) continue;
            if (structural_equal(terms[i].exponent, terms[dominant].exponent)) {
                dominant_coefficient += terms[i].coefficient;
            }
        }

        if (undecidable_tie) return std::nullopt;
        if (dominant_coefficient > 0) return ok(arena.make<Constant>(MathConstant::Infinity));
        if (dominant_coefficient < 0) return ok(limit_make_integer(arena, 0));

        ExprPtr dominant_expr = terms[dominant].exponent;
        terms.erase(
            std::remove_if(terms.begin(), terms.end(),
                [&](const ExponentialTerm& t) {
                    return structural_equal(t.exponent, dominant_expr);
                }),
            terms.end());
    }

    return std::nullopt;
}

} // namespace

// ---------------------------------------------------------------------------
// Public wrapper
// ---------------------------------------------------------------------------
std::optional<Result<ExprPtr>> mrv_try_exponential_product_limit(
    ExprPtr expr, const Symbol& var, AstArena& arena, symbolic::CASContext& ctx) {
    return try_exponential_product_limit(expr, var, arena, ctx);
}

} // namespace cas::calculus
