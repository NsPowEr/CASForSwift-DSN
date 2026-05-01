#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <functional>
#include <vector>

namespace cas::symbolic {

[[nodiscard]] bool rule_condition_satisfied(const RewriteRule& rule, const MatchMap& matches) {
    return !rule.condition || rule.condition(matches);
}

[[nodiscard]] ExprPtr build_ac_expr(ExprKind kind, std::vector<ExprPtr> operands, AstArena& arena) {
    if (operands.empty()) {
        return kind == ExprKind::Sum
            ? arena.make<IntegerLit>(BigInt(0))
            : arena.make<IntegerLit>(BigInt(1));
    }
    if (operands.size() == 1U) {
        return operands.front();
    }
    if (kind == ExprKind::Sum) {
        return arena.make<Sum>(std::move(operands));
    }
    return arena.make<Product>(std::move(operands));
}

[[nodiscard]] Result<ExprPtr> try_apply_rule_here(ExprPtr expr, const RewriteRule& rule, AstArena& arena) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot rewrite null expression"));
    }
    MatchMap matches;
    if (match_pattern_impl(expr, rule.pattern, matches) && rule_condition_satisfied(rule, matches)) {
        ExprPtr rewritten = instantiate_pattern(rule.replacement, matches, arena);
        if (is_strict_rewrite_reduction(expr, rewritten)) {
            return ok(rewritten);
        }
    }

    const ExprKind expr_kind_value = expr_kind(expr);
    if (expr_kind_value == ExprKind::Sum && expr_kind(rule.pattern) == ExprKind::Sum) {
        const auto& expr_terms = expr_ref<Sum>(expr).terms;
        const auto& pattern_terms = expr_ref<Sum>(rule.pattern).terms;
        std::vector<bool> matched_terms;
        MatchMap ac_matches;
        if (expr_terms.size() > pattern_terms.size() &&
            match_sequence_ac_internal(expr_terms, pattern_terms, ac_matches, &matched_terms) &&
            rule_condition_satisfied(rule, ac_matches)) {
            std::vector<ExprPtr> rewritten_terms;
            rewritten_terms.reserve(expr_terms.size() - pattern_terms.size() + 1U);
            for (std::size_t index = 0; index < expr_terms.size(); ++index) {
                if (!matched_terms[index]) {
                    rewritten_terms.push_back(expr_terms[index]);
                }
            }
            rewritten_terms.push_back(instantiate_pattern(rule.replacement, ac_matches, arena));
            ExprPtr rewritten = build_ac_expr(ExprKind::Sum, std::move(rewritten_terms), arena);
            if (is_strict_rewrite_reduction(expr, rewritten)) {
                return ok(rewritten);
            }
        }
    }

    if (expr_kind_value == ExprKind::Product && expr_kind(rule.pattern) == ExprKind::Product) {
        const auto& expr_factors = expr_ref<Product>(expr).factors;
        const auto& pattern_factors = expr_ref<Product>(rule.pattern).factors;
        std::vector<bool> matched_factors;
        MatchMap ac_matches;
        if (expr_factors.size() > pattern_factors.size() &&
            match_sequence_ac_internal(expr_factors, pattern_factors, ac_matches, &matched_factors) &&
            rule_condition_satisfied(rule, ac_matches)) {
            std::vector<ExprPtr> rewritten_factors;
            rewritten_factors.reserve(expr_factors.size() - pattern_factors.size() + 1U);
            for (std::size_t index = 0; index < expr_factors.size(); ++index) {
                if (!matched_factors[index]) {
                    rewritten_factors.push_back(expr_factors[index]);
                }
            }
            rewritten_factors.push_back(instantiate_pattern(rule.replacement, ac_matches, arena));
            ExprPtr rewritten = build_ac_expr(ExprKind::Product, std::move(rewritten_factors), arena);
            if (is_strict_rewrite_reduction(expr, rewritten)) {
                return ok(rewritten);
            }
        }
    }

    return ok(expr);
}

using RewriteStep = std::function<Result<ExprPtr>(ExprPtr)>;

[[nodiscard]] Result<ExprPtr> rewrite_children(ExprPtr expr, const RewriteStep& recurse, AstArena& arena) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot rewrite null expression"));
    }

    return visit_expr(
        expr,
        [&](const auto& node) -> Result<ExprPtr> {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, Symbol> ||
                std::is_same_v<Node, Constant>) {
                return ok(expr);
            } else if constexpr (std::is_same_v<Node, Unary>) {
                auto operand = recurse(node.operand);
                if (operand.is_error()) {
                    return operand;
                }
                if (operand.value() == node.operand) {
                    return ok(expr);
                }
                return ok(arena.make<Unary>(node.op, operand.value()));
            } else if constexpr (std::is_same_v<Node, Binary>) {
                auto left = recurse(node.left);
                if (left.is_error()) {
                    return left;
                }
                auto right = recurse(node.right);
                if (right.is_error()) {
                    return right;
                }
                if (left.value() == node.left && right.value() == node.right) {
                    return ok(expr);
                }
                return ok(arena.make<Binary>(node.op, left.value(), right.value()));
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> args;
                args.reserve(node.args.size());
                bool changed = false;
                for (ExprPtr arg : node.args) {
                    auto rewritten = recurse(arg);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || rewritten.value() != arg;
                    args.push_back(rewritten.value());
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(arena.make<FuncCall>(node.name, std::move(args)));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> terms;
                terms.reserve(node.terms.size());
                bool changed = false;
                for (ExprPtr term : node.terms) {
                    auto rewritten = recurse(term);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || rewritten.value() != term;
                    terms.push_back(rewritten.value());
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(arena.make<Sum>(std::move(terms)));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> factors;
                factors.reserve(node.factors.size());
                bool changed = false;
                for (ExprPtr factor : node.factors) {
                    auto rewritten = recurse(factor);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || rewritten.value() != factor;
                    factors.push_back(rewritten.value());
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(arena.make<Product>(std::move(factors)));
            } else if constexpr (std::is_same_v<Node, Integral>) {
                auto integrand = recurse(node.integrand);
                if (integrand.is_error()) {
                    return integrand;
                }

                std::optional<ExprPtr> lower = node.lower;
                std::optional<ExprPtr> upper = node.upper;
                bool changed = integrand.value() != node.integrand;
                if (node.lower.has_value()) {
                    auto rewritten_lower = recurse(*node.lower);
                    if (rewritten_lower.is_error()) {
                        return rewritten_lower;
                    }
                    changed = changed || rewritten_lower.value() != *node.lower;
                    lower = rewritten_lower.value();
                }
                if (node.upper.has_value()) {
                    auto rewritten_upper = recurse(*node.upper);
                    if (rewritten_upper.is_error()) {
                        return rewritten_upper;
                    }
                    changed = changed || rewritten_upper.value() != *node.upper;
                    upper = rewritten_upper.value();
                }

                if (!changed) {
                    return ok(expr);
                }
                return ok(arena.make<Integral>(integrand.value(), node.variable, lower, upper));
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                auto expression = recurse(node.expression);
                if (expression.is_error()) {
                    return expression;
                }
                if (expression.value() == node.expression) {
                    return ok(expr);
                }
                return ok(arena.make<Derivative>(expression.value(), node.variable, node.order));
            } else if constexpr (std::is_same_v<Node, Limit>) {
                auto expression = recurse(node.expression);
                if (expression.is_error()) {
                    return expression;
                }
                auto point = recurse(node.point);
                if (point.is_error()) {
                    return point;
                }
                if (expression.value() == node.expression && point.value() == node.point) {
                    return ok(expr);
                }
                return ok(arena.make<Limit>(expression.value(), node.variable, point.value(), node.direction));
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                auto polynomial = recurse(node.polynomial);
                if (polynomial.is_error()) {
                    return polynomial;
                }
                if (polynomial.value() == node.polynomial) {
                    return ok(expr);
                }
                return ok(arena.make<RootOf>(polynomial.value(), node.variable, node.root_index));
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::vector<ExprPtr> elements;
                elements.reserve(node.elements.size());
                bool changed = false;
                for (ExprPtr element : node.elements) {
                    auto rewritten = recurse(element);
                    if (rewritten.is_error()) {
                        return rewritten;
                    }
                    changed = changed || rewritten.value() != element;
                    elements.push_back(rewritten.value());
                }
                if (!changed) {
                    return ok(expr);
                }
                return ok(arena.make<Matrix>(node.rows, node.cols, std::move(elements)));
            } else {
                return ok(expr);
            }
        });
}

[[nodiscard]] Result<ExprPtr> apply_rule_bottom_up(ExprPtr expr, const RewriteRule& rule, AstArena& arena) {
    auto rebuilt = rewrite_children(
        expr,
        [&](ExprPtr child) {
            return apply_rule_bottom_up(child, rule, arena);
        },
        arena);
    if (rebuilt.is_error()) {
        return rebuilt;
    }
    return try_apply_rule_here(rebuilt.value(), rule, arena);
}

[[nodiscard]] Result<ExprPtr> apply_rule_top_down(ExprPtr expr, const RewriteRule& rule, AstArena& arena) {
    auto rewritten = try_apply_rule_here(expr, rule, arena);
    if (rewritten.is_error()) {
        return rewritten;
    }
    return rewrite_children(
        rewritten.value(),
        [&](ExprPtr child) {
            return apply_rule_top_down(child, rule, arena);
        },
        arena);
}

[[nodiscard]] Result<ExprPtr> apply_rule_impl(ExprPtr expr, const RewriteRule& rule, TraversalStrategy strategy, AstArena& arena) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot rewrite null expression"));
    }

    if (!rule.pattern || !rule.replacement) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Rewrite rule must define pattern and replacement"));
    }
    if (!rewrite_rule_is_oriented_impl(rule)) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "Rewrite rule is not strictly decreasing under the term order"));
    }

    switch (strategy) {
        case TraversalStrategy::BottomUp:
            return apply_rule_bottom_up(expr, rule, arena);
        case TraversalStrategy::TopDown:
            return apply_rule_top_down(expr, rule, arena);
        case TraversalStrategy::FixPoint: {
            ExprPtr current = expr;
            while (true) {
                auto rewritten = apply_rule_bottom_up(current, rule, arena);
                if (rewritten.is_error()) {
                    return rewritten;
                }
                if (rewritten.value() == current) {
                    return rewritten;
                }
                current = rewritten.value();
            }
        }
    }

    return ok(expr);
}

[[nodiscard]] Result<ExprPtr> apply_rule_set_impl(ExprPtr expr, const std::vector<RewriteRule>& rules, AstArena& arena) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot rewrite null expression"));
    }
    if (!is_strongly_normalizing(rules)) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::InvalidArgument,
            "Rewrite rule set contains a rule that is not strictly decreasing under the term order"));
    }

    ExprPtr current = expr;
    while (true) {
        bool changed = false;
        for (const RewriteRule& rule : rules) {
            auto rewritten = apply_rule_impl(current, rule, TraversalStrategy::BottomUp, arena);
            if (rewritten.is_error()) {
                return rewritten;
            }
            if (rewritten.value() != current) {
                current = rewritten.value();
                changed = true;
            }
        }
        if (!changed) {
            return ok(current);
        }
    }
}

Result<ExprPtr> apply_rule(ExprPtr expr, const RewriteRule& rule, TraversalStrategy strategy, AstArena& arena) {
    return apply_rule_impl(expr, rule, strategy, arena);
}

Result<ExprPtr> apply_rule_set(ExprPtr expr, const std::vector<RewriteRule>& rules, AstArena& arena) {
    return apply_rule_set_impl(expr, rules, arena);
}

} // namespace cas::symbolic
