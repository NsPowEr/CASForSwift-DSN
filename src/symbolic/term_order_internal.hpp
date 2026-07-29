// term_order_internal.hpp — Private declarations for term ordering.
#pragma once

#include "cas/symbolic.hpp"
#include <vector>
#include <string>

namespace cas::symbolic {

// Helpers declared in term_order.cpp
[[nodiscard]] int term_kind_rank(ExprKind kind) noexcept;
[[nodiscard]] int compare_kind_precedence(ExprKind lhs, ExprKind rhs) noexcept;
[[nodiscard]] int compare_string_precedence(const std::string& lhs, const std::string& rhs) noexcept;
[[nodiscard]] int get_builtin_precedence(BuiltinOp op) noexcept;
[[nodiscard]] int get_binary_op_precedence(BinaryOp op) noexcept;
[[nodiscard]] int compare_bigint(const BigInt& lhs, const BigInt& rhs) noexcept;
[[nodiscard]] std::vector<ExprPtr> term_order_children(ExprPtr expr);

[[nodiscard]] int compare_head_precedence(ExprPtr lhs, ExprPtr rhs) noexcept;

// Functions implemented in term_order_lpo.cpp
[[nodiscard]] bool term_order_gt(ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] bool term_order_ge(ExprPtr lhs, ExprPtr rhs);

// Functions implemented in term_order_kbw.cpp
[[nodiscard]] TermOrderRelation relation_from_compare(int cmp) noexcept;
[[nodiscard]] TermOrderRelation compare_knuth_bendix_weight_order(ExprPtr lhs, ExprPtr rhs);

}  // namespace cas::symbolic
