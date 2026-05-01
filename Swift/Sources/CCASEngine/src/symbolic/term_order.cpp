#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace cas::symbolic {
namespace {

[[nodiscard]] int term_kind_rank(ExprKind kind) noexcept {
    switch (kind) {
    case ExprKind::Null:
        return 0;
    case ExprKind::IntegerLit:
    case ExprKind::RationalLit:
    case ExprKind::DecimalLit:
    case ExprKind::Constant:
        return 1;
    case ExprKind::Symbol:
        return 2;
    case ExprKind::Unary:
        return 3;
    case ExprKind::Binary:
        return 4;
    case ExprKind::Sum:
        return 5;
    case ExprKind::Product:
        return 6;
    case ExprKind::FuncCall:
        return 7;
    case ExprKind::Integral:
    case ExprKind::Derivative:
    case ExprKind::Limit:
    case ExprKind::RootOf:
    case ExprKind::Matrix:
        return 8;
    }

    return 9;
}

[[nodiscard]] int compare_kind_precedence(ExprKind lhs, ExprKind rhs) noexcept {
    const int lhs_rank = term_kind_rank(lhs);
    const int rhs_rank = term_kind_rank(rhs);
    if (lhs_rank < rhs_rank) {
        return -1;
    }
    if (lhs_rank > rhs_rank) {
        return 1;
    }
    return 0;
}

[[nodiscard]] int compare_string_precedence(const std::string& lhs, const std::string& rhs) noexcept {
    if (lhs < rhs) {
        return -1;
    }
    if (rhs < lhs) {
        return 1;
    }
    return 0;
}

[[nodiscard]] std::vector<ExprPtr> term_order_children(ExprPtr expr) {
    return visit_expr(
        expr,
        [](const auto& node) -> std::vector<ExprPtr> {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<Node, Unary>) {
                return {node.operand};
            } else if constexpr (std::is_same_v<Node, Binary>) {
                return {node.left, node.right};
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                return node.args;
            } else if constexpr (std::is_same_v<Node, Sum>) {
                return node.terms;
            } else if constexpr (std::is_same_v<Node, Product>) {
                return node.factors;
            } else if constexpr (std::is_same_v<Node, Integral>) {
                std::vector<ExprPtr> children{node.integrand};
                if (node.lower.has_value()) {
                    children.push_back(*node.lower);
                }
                if (node.upper.has_value()) {
                    children.push_back(*node.upper);
                }
                return children;
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                return {node.expression};
            } else if constexpr (std::is_same_v<Node, Limit>) {
                return {node.expression, node.point};
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                return {node.polynomial};
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                return node.elements;
            } else {
                return {};
            }
        });
}

[[nodiscard]] bool term_order_ge(ExprPtr lhs, ExprPtr rhs);

[[nodiscard]] bool term_order_gt(ExprPtr lhs, ExprPtr rhs) {
    if (!lhs || !rhs || structural_equal(lhs, rhs)) {
        return false;
    }

    const std::vector<ExprPtr> lhs_children = term_order_children(lhs);
    for (ExprPtr child : lhs_children) {
        if (term_order_ge(child, rhs)) {
            return true;
        }
    }

    const std::vector<ExprPtr> rhs_children = term_order_children(rhs);
    const auto dominates_rhs_children = [&]() {
        return std::all_of(
            rhs_children.begin(),
            rhs_children.end(),
            [&](ExprPtr child) {
                return term_order_gt(lhs, child);
            });
    };

    const int head_cmp = compare_head_precedence(lhs, rhs);
    if (head_cmp > 0 && dominates_rhs_children()) {
        return true;
    }

    if (head_cmp == 0 && lhs_children.size() == rhs_children.size() && dominates_rhs_children()) {
        for (std::size_t index = 0; index < lhs_children.size(); ++index) {
            if (structural_equal(lhs_children[index], rhs_children[index])) {
                continue;
            }
            return term_order_gt(lhs_children[index], rhs_children[index]);
        }
    }

    return false;
}

[[nodiscard]] bool term_order_ge(ExprPtr lhs, ExprPtr rhs) {
    return structural_equal(lhs, rhs) || term_order_gt(lhs, rhs);
}

[[nodiscard]] TermOrderRelation relation_from_compare(int cmp) noexcept {
    if (cmp < 0) {
        return TermOrderRelation::Less;
    }
    if (cmp > 0) {
        return TermOrderRelation::Greater;
    }
    return TermOrderRelation::Equivalent;
}

[[nodiscard]] TermOrderRelation compare_knuth_bendix_weight_order(ExprPtr lhs, ExprPtr rhs) {
    if (!lhs || !rhs) {
        return lhs == rhs ? TermOrderRelation::Equivalent : TermOrderRelation::Incomparable;
    }
    if (structural_equal(lhs, rhs)) {
        return TermOrderRelation::Equivalent;
    }

    const std::size_t lhs_weight = expr_weight(lhs);
    const std::size_t rhs_weight = expr_weight(rhs);
    if (lhs_weight < rhs_weight) {
        return TermOrderRelation::Less;
    }
    if (lhs_weight > rhs_weight) {
        return TermOrderRelation::Greater;
    }

    const int head_cmp = compare_head_precedence(lhs, rhs);
    if (head_cmp != 0) {
        return relation_from_compare(head_cmp);
    }

    const std::vector<ExprPtr> lhs_children = term_order_children(lhs);
    const std::vector<ExprPtr> rhs_children = term_order_children(rhs);
    const std::size_t shared_size = std::min(lhs_children.size(), rhs_children.size());
    for (std::size_t index = 0; index < shared_size; ++index) {
        const TermOrderRelation child_cmp = compare_knuth_bendix_weight_order(lhs_children[index], rhs_children[index]);
        if (child_cmp != TermOrderRelation::Equivalent) {
            return child_cmp;
        }
    }
    if (lhs_children.size() < rhs_children.size()) {
        return TermOrderRelation::Less;
    }
    if (lhs_children.size() > rhs_children.size()) {
        return TermOrderRelation::Greater;
    }

    return relation_from_compare(canonical_compare(lhs, rhs));
}

} // namespace

int compare_head_precedence(ExprPtr lhs, ExprPtr rhs) noexcept {
    const ExprKind lhs_kind = expr_kind(lhs);
    const ExprKind rhs_kind = expr_kind(rhs);
    const int kind_cmp = compare_kind_precedence(lhs_kind, rhs_kind);
    if (kind_cmp != 0) {
        return kind_cmp;
    }

    if (const auto* lhs_unary = expr_cast<Unary>(lhs)) {
        const auto* rhs_unary = expr_cast<Unary>(rhs);
        if (lhs_unary->op < rhs_unary->op) {
            return -1;
        }
        if (rhs_unary->op < lhs_unary->op) {
            return 1;
        }
    }

    if (const auto* lhs_binary = expr_cast<Binary>(lhs)) {
        const auto* rhs_binary = expr_cast<Binary>(rhs);
        if (lhs_binary->op < rhs_binary->op) {
            return -1;
        }
        if (rhs_binary->op < lhs_binary->op) {
            return 1;
        }
    }

    if (const auto* lhs_call = expr_cast<FuncCall>(lhs)) {
        const auto* rhs_call = expr_cast<FuncCall>(rhs);
        return compare_string_precedence(lhs_call->name, rhs_call->name);
    }

    if (const auto* lhs_symbol = expr_cast<Symbol>(lhs)) {
        const auto* rhs_symbol = expr_cast<Symbol>(rhs);
        return compare_string_precedence(lhs_symbol->name, rhs_symbol->name);
    }

    return 0;
}

TermOrderRelation compare_rewrite_terms_impl(ExprPtr lhs, ExprPtr rhs) {
    if (!lhs || !rhs) {
        return lhs == rhs ? TermOrderRelation::Equivalent : TermOrderRelation::Incomparable;
    }
    if (structural_equal(lhs, rhs)) {
        return TermOrderRelation::Equivalent;
    }

    const bool lhs_gt_rhs = term_order_gt(lhs, rhs);
    const bool rhs_gt_lhs = term_order_gt(rhs, lhs);
    if (lhs_gt_rhs != rhs_gt_lhs) {
        return lhs_gt_rhs ? TermOrderRelation::Greater : TermOrderRelation::Less;
    }

    return compare_knuth_bendix_weight_order(lhs, rhs);
}

bool is_strict_rewrite_reduction(ExprPtr before, ExprPtr after) {
    return compare_rewrite_terms_impl(after, before) == TermOrderRelation::Less;
}

bool rewrite_rule_is_oriented_impl(const RewriteRule& rule) {
    return static_cast<bool>(rule.pattern) &&
        static_cast<bool>(rule.replacement) &&
        is_strict_rewrite_reduction(rule.pattern, rule.replacement);
}

TermOrderRelation compare_rewrite_terms(ExprPtr lhs, ExprPtr rhs) {
    return compare_rewrite_terms_impl(lhs, rhs);
}

bool rewrite_rule_is_oriented(const RewriteRule& rule) {
    return rewrite_rule_is_oriented_impl(rule);
}

bool is_strongly_normalizing(const std::vector<RewriteRule>& rules) {
    return std::all_of(
        rules.begin(),
        rules.end(),
        [](const RewriteRule& rule) {
            return rewrite_rule_is_oriented_impl(rule);
        });
}

} // namespace cas::symbolic
