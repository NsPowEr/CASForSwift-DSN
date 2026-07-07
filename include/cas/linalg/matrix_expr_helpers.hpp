#ifndef CAS_LINALG_MATRIX_EXPR_HELPERS_HPP
#define CAS_LINALG_MATRIX_EXPR_HELPERS_HPP

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/result.hpp"

#include <tuple>
#include <vector>

namespace cas::linalg {

/**
 * @brief Helper to create an integer literal expression.
 */
[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, long long value);
[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, const BigInt& value);

/**
 * @brief Helper to simplify an expression using the provided context.
 */
[[nodiscard]] Result<ExprPtr> simplify(symbolic::CASContext& ctx, ExprPtr expr);

/**
 * @brief Helper to add two expressions and simplify the result.
 */
[[nodiscard]] Result<ExprPtr> add_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs);

/**
 * @brief Helper to negate an expression and simplify the result.
 */
[[nodiscard]] Result<ExprPtr> negate_expr(symbolic::CASContext& ctx, ExprPtr expr);

/**
 * @brief Helper to subtract two expressions and simplify the result.
 */
[[nodiscard]] Result<ExprPtr> sub_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs);

/**
 * @brief Helper to multiply two expressions and simplify the result.
 */
[[nodiscard]] Result<ExprPtr> mul_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs);

/**
 * @brief Helper to divide two expressions and simplify the result.
 */
[[nodiscard]] Result<ExprPtr> div_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs);

/**
 * @brief Checks if an expression is a numeric literal zero.
 */
[[nodiscard]] bool is_zero_expr(ExprPtr expr);

/**
 * @brief Checks if an expression is a numeric literal one.
 */
[[nodiscard]] bool is_one_expr(ExprPtr expr);

/**
 * @brief Tries to extract a BigInt from an expression.
 */
[[nodiscard]] std::optional<BigInt> try_get_bigint(ExprPtr expr);

/**
 * @brief Estimates the complexity (AST size/node count) of an expression.
 */
[[nodiscard]] std::size_t estimate_complexity(ExprPtr expr);

/**
 * @brief Checks if an expression is structurally non-zero.
 */
[[nodiscard]] bool is_structurally_nonzero(ExprPtr expr);

/**
 * @brief Checks if an expression is known to be non-zero using assumptions in the context.
 */
[[nodiscard]] bool is_known_nonzero(ExprPtr expr, symbolic::CASContext& ctx);

/**
 * @brief Estimates the total symbolic degree of an expression.
 */
[[nodiscard]] std::size_t total_degree(ExprPtr expr);

/**
 * @brief Lexicographic pivot quality score for selection algorithms.
 */
struct PivotScore {
    int certainty;
    int neg_complexity;
    int neg_total_degree;

    [[nodiscard]] auto as_tuple() const noexcept {
        return std::tie(certainty, neg_complexity, neg_total_degree);
    }
    [[nodiscard]] bool operator<(const PivotScore& other) const noexcept {
        return as_tuple() < other.as_tuple();
    }
    [[nodiscard]] bool operator==(const PivotScore& other) const noexcept {
        return as_tuple() == other.as_tuple();
    }
    [[nodiscard]] bool operator>(const PivotScore& other) const noexcept {
        return other < *this;
    }
};

/**
 * @brief Creates a PivotScore for the given expression based on the context's assumptions.
 */
[[nodiscard]] PivotScore make_pivot_score(ExprPtr val, symbolic::CASContext& ctx);

/**
 * @brief Sum of squares Σ v_i² (no sqrt). Used as squared-norm for symbolic vectors.
 */
[[nodiscard]] Result<ExprPtr> sym_norm_sq(const std::vector<ExprPtr>& v, symbolic::CASContext& ctx);

/**
 * @brief Euclidean norm sqrt(Σ v_i²), simplified.
 */
[[nodiscard]] Result<ExprPtr> sym_norm(const std::vector<ExprPtr>& v, symbolic::CASContext& ctx);

} // namespace cas::linalg

#endif // CAS_LINALG_MATRIX_EXPR_HELPERS_HPP
