// builtin_rewrite_internal.hpp — Private declarations for builtin rewrite.
#pragma once

#include "cas/symbolic.hpp"
#include "cas/rational.hpp"

#include <vector>
#include <optional>
#include <string>

namespace cas::symbolic {

// Helper utilities

[[nodiscard]] bool exact_expr_is_positive(ExprPtr expr);
[[nodiscard]] bool exact_expr_is_nonnegative(ExprPtr expr);
[[nodiscard]] bool exact_expr_is_negative(ExprPtr expr);
[[nodiscard]] bool expr_is_positive_under_assumptions(ExprPtr expr, const Assumptions* assumptions);
[[nodiscard]] bool expr_is_nonnegative_under_assumptions(ExprPtr expr, const Assumptions* assumptions);

// Trig rules and helpers
void add_trig_rules(std::vector<RewriteRule>& rules, AstArena& arena);
[[nodiscard]] Result<ExprPtr> try_rewrite_trig(ExprPtr expr, AstArena& arena);

// Log/Exp rules and helpers
void add_log_exp_rules(std::vector<RewriteRule>& rules, AstArena& arena);
[[nodiscard]] Result<ExprPtr> try_rewrite_log_exp(ExprPtr expr, AstArena& arena, const Assumptions* assumptions);

// Algebraic rules and helpers
[[nodiscard]] Result<ExprPtr> try_rewrite_algebraic(
    ExprPtr expr,
    AstArena& arena,
    const Assumptions* assumptions,
    CASContext* context);

// Shared trig/parity checkers
[[nodiscard]] bool is_odd_parity_function(BuiltinOp func_id);
[[nodiscard]] bool is_even_parity_function(BuiltinOp func_id);
[[nodiscard]] bool is_parity_rewrite_function(BuiltinOp func_id);

enum class SquareFunctionKind : std::uint8_t {
    None,
    Sin,
    Cos,
};
[[nodiscard]] SquareFunctionKind square_function_kind(ExprPtr expr);

[[nodiscard]] bool may_match_builtin_rewrite(ExprPtr expr);

}  // namespace cas::symbolic
