#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <string>
#include <vector>

namespace cas::symbolic {

[[nodiscard]] bool is_wildcard_name(const std::string& name) {
    return name.size() > 1U && name.back() == '_';
}

[[nodiscard]] bool match_pattern_impl(ExprPtr expr, ExprPtr pattern, MatchMap& matches) {
    if (!expr || !pattern) {
        return expr == pattern;
    }

    if (const auto* wildcard = expr_cast<Symbol>(pattern)) {
        if (is_wildcard_name(wildcard->name)) {
            const auto existing = matches.find(wildcard->name);
            if (existing == matches.end()) {
                matches.emplace(wildcard->name, expr);
                return true;
            }
            return structural_equal(existing->second, expr);
        }
    }

    if (expr_kind(expr) != expr_kind(pattern)) {
        return false;
    }

    if (const auto* lhs_unary = expr_cast<Unary>(expr)) {
        const auto* rhs_unary = expr_cast<Unary>(pattern);
        return lhs_unary->op == rhs_unary->op && match_pattern_impl(lhs_unary->operand, rhs_unary->operand, matches);
    }

    if (const auto* lhs_binary = expr_cast<Binary>(expr)) {
        const auto* rhs_binary = expr_cast<Binary>(pattern);
        return lhs_binary->op == rhs_binary->op &&
            match_pattern_impl(lhs_binary->left, rhs_binary->left, matches) &&
            match_pattern_impl(lhs_binary->right, rhs_binary->right, matches);
    }

    if (const auto* lhs_call = expr_cast<FuncCall>(expr)) {
        const auto* rhs_call = expr_cast<FuncCall>(pattern);
        if (lhs_call->name != rhs_call->name || lhs_call->args.size() != rhs_call->args.size()) {
            return false;
        }

        for (std::size_t index = 0; index < lhs_call->args.size(); ++index) {
            if (!match_pattern_impl(lhs_call->args[index], rhs_call->args[index], matches)) {
                return false;
            }
        }
        return true;
    }

    if (const auto* lhs_sum = expr_cast<Sum>(expr)) {
        const auto* rhs_sum = expr_cast<Sum>(pattern);
        return match_sequence_ac_exact(lhs_sum->terms, rhs_sum->terms, matches);
    }

    if (const auto* lhs_product = expr_cast<Product>(expr)) {
        const auto* rhs_product = expr_cast<Product>(pattern);
        return match_sequence_ac_exact(lhs_product->factors, rhs_product->factors, matches);
    }

    if (const auto* lhs_integral = expr_cast<Integral>(expr)) {
        const auto* rhs_integral = expr_cast<Integral>(pattern);
        const bool same_bounds =
            lhs_integral->lower.has_value() == rhs_integral->lower.has_value() &&
            lhs_integral->upper.has_value() == rhs_integral->upper.has_value();
        if (!same_bounds || lhs_integral->variable.name != rhs_integral->variable.name) {
            return false;
        }
        if (!match_pattern_impl(lhs_integral->integrand, rhs_integral->integrand, matches)) {
            return false;
        }
        if (lhs_integral->lower.has_value() &&
            (!match_pattern_impl(*lhs_integral->lower, *rhs_integral->lower, matches) ||
             !match_pattern_impl(*lhs_integral->upper, *rhs_integral->upper, matches))) {
            return false;
        }
        return true;
    }

    if (const auto* lhs_derivative = expr_cast<Derivative>(expr)) {
        const auto* rhs_derivative = expr_cast<Derivative>(pattern);
        return lhs_derivative->order == rhs_derivative->order &&
            lhs_derivative->variable.name == rhs_derivative->variable.name &&
            match_pattern_impl(lhs_derivative->expression, rhs_derivative->expression, matches);
    }

    if (const auto* lhs_limit = expr_cast<Limit>(expr)) {
        const auto* rhs_limit = expr_cast<Limit>(pattern);
        return lhs_limit->direction == rhs_limit->direction &&
            lhs_limit->variable.name == rhs_limit->variable.name &&
            match_pattern_impl(lhs_limit->expression, rhs_limit->expression, matches) &&
            match_pattern_impl(lhs_limit->point, rhs_limit->point, matches);
    }

    if (const auto* lhs_root = expr_cast<RootOf>(expr)) {
        const auto* rhs_root = expr_cast<RootOf>(pattern);
        return lhs_root->variable.name == rhs_root->variable.name &&
            lhs_root->root_index == rhs_root->root_index &&
            match_pattern_impl(lhs_root->polynomial, rhs_root->polynomial, matches);
    }

    if (const auto* lhs_matrix = expr_cast<Matrix>(expr)) {
        const auto* rhs_matrix = expr_cast<Matrix>(pattern);
        if (lhs_matrix->rows != rhs_matrix->rows ||
            lhs_matrix->cols != rhs_matrix->cols ||
            lhs_matrix->elements.size() != rhs_matrix->elements.size()) {
            return false;
        }

        for (std::size_t index = 0; index < lhs_matrix->elements.size(); ++index) {
            if (!match_pattern_impl(lhs_matrix->elements[index], rhs_matrix->elements[index], matches)) {
                return false;
            }
        }
        return true;
    }

    return structural_equal(expr, pattern);
}

[[nodiscard]] bool match_sequence_ac_internal(
    const std::vector<ExprPtr>& exprs,
    const std::vector<ExprPtr>& patterns,
    MatchMap& matches,
    std::vector<bool>* matched_exprs) {
    if (exprs.size() < patterns.size()) {
        return false;
    }

    std::vector<std::size_t> ordered_patterns(patterns.size());
    for (std::size_t index = 0; index < patterns.size(); ++index) {
        ordered_patterns[index] = index;
    }
    std::sort(
        ordered_patterns.begin(),
        ordered_patterns.end(),
        [&](std::size_t lhs, std::size_t rhs) {
            return expr_weight(patterns[lhs]) > expr_weight(patterns[rhs]);
        });

    std::vector<bool> used(exprs.size(), false);
    std::function<bool(std::size_t, MatchMap&)> backtrack =
        [&](std::size_t pattern_index, MatchMap& current_matches) -> bool {
            if (pattern_index == ordered_patterns.size()) {
                return true;
            }

            const ExprPtr current_pattern = patterns[ordered_patterns[pattern_index]];
            for (std::size_t expr_index = 0; expr_index < exprs.size(); ++expr_index) {
                if (used[expr_index]) {
                    continue;
                }

                MatchMap trial_matches = current_matches;
                if (!match_pattern_impl(exprs[expr_index], current_pattern, trial_matches)) {
                    continue;
                }

                used[expr_index] = true;
                if (backtrack(pattern_index + 1U, trial_matches)) {
                    current_matches = std::move(trial_matches);
                    return true;
                }
                used[expr_index] = false;
            }

            return false;
        };

    const bool matched = backtrack(0U, matches);
    if (matched && matched_exprs != nullptr) {
        *matched_exprs = used;
    }
    return matched;
}

[[nodiscard]] bool match_sequence_ac_exact(
    const std::vector<ExprPtr>& exprs,
    const std::vector<ExprPtr>& patterns,
    MatchMap& matches) {
    if (exprs.size() != patterns.size()) {
        return false;
    }
    return match_sequence_ac_internal(exprs, patterns, matches, nullptr);
}

[[nodiscard]] bool match_ac_pattern_impl(ExprPtr expr, ExprPtr pattern, MatchMap& matches) {
    if (!expr || !pattern) {
        return expr == pattern;
    }

    if (expr_kind(expr) == ExprKind::Sum && expr_kind(pattern) == ExprKind::Sum) {
        return match_sequence_ac_exact(expr_ref<Sum>(expr).terms, expr_ref<Sum>(pattern).terms, matches);
    }
    if (expr_kind(expr) == ExprKind::Product && expr_kind(pattern) == ExprKind::Product) {
        return match_sequence_ac_exact(expr_ref<Product>(expr).factors, expr_ref<Product>(pattern).factors, matches);
    }
    return match_pattern_impl(expr, pattern, matches);
}

bool match_pattern(ExprPtr expr, ExprPtr pattern, MatchMap& out_matches) {
    return match_pattern_impl(expr, pattern, out_matches);
}

bool match_ac_pattern(ExprPtr expr, ExprPtr pattern, MatchMap& out_matches) {
    return match_ac_pattern_impl(expr, pattern, out_matches);
}

} // namespace cas::symbolic
