#pragma once

#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include <string>
#include <vector>
#include <optional>
#include <unordered_map>

namespace cas::symbolic {
namespace detail {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message);
[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr);
[[nodiscard]] ExprPtr negate_expr(ExprPtr expr, AstArena& arena);
void append_difference_terms(ExprPtr expr, bool negate, std::vector<ExprPtr>& terms, AstArena& arena);
[[nodiscard]] int compare_exact_scalars(const Rational& lhs, const Rational& rhs);
[[nodiscard]] bool range_is_exact_zero(ExprPtr lower, ExprPtr upper);
[[nodiscard]] bool exact_range_excludes_zero(ExprPtr lower, ExprPtr upper);
[[nodiscard]] bool is_wildcard_name(const std::string& name);
[[nodiscard]] std::size_t expr_weight(ExprPtr expr);
[[nodiscard]] int term_kind_rank(ExprKind kind) noexcept;
[[nodiscard]] int compare_kind_precedence(ExprKind lhs, ExprKind rhs) noexcept;
[[nodiscard]] int compare_string_precedence(const std::string& lhs, const std::string& rhs) noexcept;
[[nodiscard]] int compare_head_precedence(ExprPtr lhs, ExprPtr rhs) noexcept;
[[nodiscard]] std::vector<ExprPtr> term_order_children(ExprPtr expr);
[[nodiscard]] TermOrderRelation relation_from_compare(int cmp) noexcept;
[[nodiscard]] bool term_order_ge(ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] bool term_order_gt(ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] TermOrderRelation compare_knuth_bendix_weight_order(ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] TermOrderRelation compare_rewrite_terms_impl(ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] bool match_pattern_impl(ExprPtr expr, ExprPtr pattern, MatchMap& matches);
[[nodiscard]] bool match_sequence_ac_internal(const std::vector<ExprPtr>& exprs, const std::vector<ExprPtr>& patterns, MatchMap& matches, std::vector<bool>* matched_exprs);

} // namespace detail
} // namespace cas::symbolic
