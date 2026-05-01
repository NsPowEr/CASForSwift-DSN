#pragma once

#include "cas/symbolic.hpp"
#include <chrono>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace cas::symbolic {

// Internal matching utilities
[[nodiscard]] bool is_wildcard_name(const std::string& name);
[[nodiscard]] bool match_pattern_impl(ExprPtr expr, ExprPtr pattern, MatchMap& matches);
[[nodiscard]] bool match_ac_pattern_impl(ExprPtr expr, ExprPtr pattern, MatchMap& matches);
[[nodiscard]] bool match_sequence_ac_internal(
    const std::vector<ExprPtr>& exprs,
    const std::vector<ExprPtr>& patterns,
    MatchMap& matches,
    std::vector<bool>* matched_exprs);
[[nodiscard]] bool match_sequence_ac_exact(
    const std::vector<ExprPtr>& exprs,
    const std::vector<ExprPtr>& patterns,
    MatchMap& matches);

// Internal rewrite utilities
[[nodiscard]] bool is_strict_rewrite_reduction(ExprPtr before, ExprPtr after);
[[nodiscard]] bool rewrite_rule_is_oriented_impl(const RewriteRule& rule);
[[nodiscard]] ExprPtr instantiate_pattern(ExprPtr pattern, const MatchMap& matches, AstArena& arena);
[[nodiscard]] Result<ExprPtr> try_apply_rule_here(ExprPtr expr, const RewriteRule& rule, AstArena& arena);
[[nodiscard]] Result<ExprPtr> apply_rule_impl(ExprPtr expr, const RewriteRule& rule, TraversalStrategy strategy, AstArena& arena);
[[nodiscard]] Result<ExprPtr> apply_rule_set_impl(ExprPtr expr, const std::vector<RewriteRule>& rules, AstArena& arena);

// Internal scalar and expression helpers
[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr);
[[nodiscard]] int compare_exact_scalars(const Rational& lhs, const Rational& rhs);
[[nodiscard]] ExprPtr negate_expr(ExprPtr expr, AstArena& arena);
[[nodiscard]] bool range_is_exact_zero(ExprPtr lower, ExprPtr upper);
[[nodiscard]] bool exact_range_excludes_zero(ExprPtr lower, ExprPtr upper);
[[nodiscard]] std::size_t expr_weight(ExprPtr expr);
[[nodiscard]] int compare_head_precedence(ExprPtr lhs, ExprPtr rhs) noexcept;
[[nodiscard]] TermOrderRelation compare_rewrite_terms_impl(ExprPtr lhs, ExprPtr rhs);

// Internal CASContext helpers
[[nodiscard]] const ComputationTrace& empty_trace() noexcept;
[[nodiscard]] Result<ExprPtr> materialize_expr_impl(ExprPtr expr, AstArena& arena);

constexpr std::uint64_t kTimeoutCheckInterval = 1024U;

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message);

} // namespace cas::symbolic
