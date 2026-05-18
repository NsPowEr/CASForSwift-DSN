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
    case ExprKind::SeriesExp: return 85;
    case ExprKind::Quantity: return 90;
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

[[nodiscard]] int get_builtin_precedence(BuiltinOp op) noexcept {
    switch (op) {
    case BuiltinOp::Exp: return 100;
    case BuiltinOp::Ln:
    case BuiltinOp::Log: return 90;
    case BuiltinOp::Sqrt: return 85;
    case BuiltinOp::Sin:
    case BuiltinOp::Cos:
    case BuiltinOp::Tan: return 80;
    case BuiltinOp::Asin:
    case BuiltinOp::Acos:
    case BuiltinOp::Atan: return 75;
    case BuiltinOp::Sinh:
    case BuiltinOp::Cosh:
    case BuiltinOp::Tanh: return 70;
    default: return 50;
    }
}

[[nodiscard]] int get_binary_op_precedence(BinaryOp op) noexcept {
    switch (op) {
    case BinaryOp::Pow: return 60;
    case BinaryOp::Mul: return 40;
    case BinaryOp::Div: return 40;
    case BinaryOp::Add: return 30;
    case BinaryOp::Sub: return 30;
    case BinaryOp::Mod: return 20;
    case BinaryOp::Equal:
    case BinaryOp::Less:
    case BinaryOp::Greater:
    case BinaryOp::LessEqual:
    case BinaryOp::GreaterEqual:
        return 10;
    }
    return 0;
}

[[nodiscard]] int compare_bigint(const BigInt& lhs, const BigInt& rhs) noexcept {
    if (lhs < rhs) return -1;
    if (rhs < lhs) return 1;
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
            } else if constexpr (std::is_same_v<Node, SeriesExp>) {
                std::vector<ExprPtr> children{node.point};
                for (const auto& [exp, coeff] : node.terms) children.push_back(coeff);
                return children;
            } else {
                return {};
            }
        });
}

} // namespace

int canonical_compare(ExprPtr lhs, ExprPtr rhs) noexcept {
    if (lhs == rhs) return 0;
    if (!lhs) return -1;
    if (!rhs) return 1;

    const ExprKind lhs_kind = expr_kind(lhs);
    const ExprKind rhs_kind = expr_kind(rhs);
    const int lhs_rank = term_kind_rank(lhs_kind);
    const int rhs_rank = term_kind_rank(rhs_kind);
    if (lhs_rank != rhs_rank) return lhs_rank < rhs_rank ? -1 : 1;

    // Se hanno lo stesso rango, confrontiamo il tipo specifico prima di castare
    if (lhs_kind != rhs_kind) return static_cast<int>(lhs_kind) < static_cast<int>(rhs_kind) ? -1 : 1;

    if (const auto* lhs_integer = expr_cast<IntegerLit>(lhs)) {
        return compare_bigint(lhs_integer->value, expr_cast<IntegerLit>(rhs)->value);
    }

    if (const auto* lhs_rational = expr_cast<RationalLit>(lhs)) {
        const auto* rhs_rational = expr_cast<RationalLit>(rhs);
        const int num_cmp = compare_bigint(lhs_rational->numerator, rhs_rational->numerator);
        if (num_cmp != 0) return num_cmp;
        return compare_bigint(lhs_rational->denominator, rhs_rational->denominator);
    }

    if (const auto* lhs_constant = expr_cast<Constant>(lhs)) {
        const auto* rhs_constant = expr_cast<Constant>(rhs);
        if (lhs_constant->value < rhs_constant->value) return -1;
        if (rhs_constant->value < lhs_constant->value) return 1;
        return 0;
    }

    if (const auto* lhs_symbol = expr_cast<Symbol>(lhs)) {
        const auto* rhs_symbol = expr_cast<Symbol>(rhs);
        if (lhs_symbol->name < rhs_symbol->name) return -1;
        if (rhs_symbol->name < lhs_symbol->name) return 1;
        return 0;
    }

    if (const auto* lhs_derivative = expr_cast<Derivative>(lhs)) {
        const auto* rhs_derivative = expr_cast<Derivative>(rhs);
        if (lhs_derivative->order < rhs_derivative->order) return -1;
        if (lhs_derivative->order > rhs_derivative->order) return 1;
        int var_cmp = compare_string_precedence(lhs_derivative->variable.name, rhs_derivative->variable.name);
        if (var_cmp != 0) return var_cmp;
        return canonical_compare(lhs_derivative->expression, rhs_derivative->expression);
    }

    if (const auto* lhs_integral = expr_cast<Integral>(lhs)) {
        const auto* rhs_integral = expr_cast<Integral>(rhs);
        int var_cmp = compare_string_precedence(lhs_integral->variable.name, rhs_integral->variable.name);
        if (var_cmp != 0) return var_cmp;
        int int_cmp = canonical_compare(lhs_integral->integrand, rhs_integral->integrand);
        if (int_cmp != 0) return int_cmp;
        // Compare optional limits
        auto compare_opt = [](const std::optional<ExprPtr>& l, const std::optional<ExprPtr>& r) -> int {
            if (!l.has_value() && !r.has_value()) return 0;
            if (!l.has_value()) return -1;
            if (!r.has_value()) return 1;
            return canonical_compare(*l, *r);
        };
        int low_cmp = compare_opt(lhs_integral->lower, rhs_integral->lower);
        if (low_cmp != 0) return low_cmp;
        return compare_opt(lhs_integral->upper, rhs_integral->upper);
    }

    if (const auto* lhs_limit = expr_cast<Limit>(lhs)) {
        const auto* rhs_limit = expr_cast<Limit>(rhs);
        int var_cmp = compare_string_precedence(lhs_limit->variable.name, rhs_limit->variable.name);
        if (var_cmp != 0) return var_cmp;
        if (lhs_limit->direction < rhs_limit->direction) return -1;
        if (lhs_limit->direction > rhs_limit->direction) return 1;
        int expr_cmp = canonical_compare(lhs_limit->expression, rhs_limit->expression);
        if (expr_cmp != 0) return expr_cmp;
        return canonical_compare(lhs_limit->point, rhs_limit->point);
    }

    if (const auto* lhs_unary = expr_cast<Unary>(lhs)) {
        const auto* rhs_unary = expr_cast<Unary>(rhs);
        if (lhs_unary->op < rhs_unary->op) return -1;
        if (lhs_unary->op > rhs_unary->op) return 1;
        return canonical_compare(lhs_unary->operand, rhs_unary->operand);
    }

    if (const auto* lhs_binary = expr_cast<Binary>(lhs)) {
        const auto* rhs_binary = expr_cast<Binary>(rhs);
        if (lhs_binary->op < rhs_binary->op) return -1;
        if (lhs_binary->op > rhs_binary->op) return 1;
        int left_cmp = canonical_compare(lhs_binary->left, rhs_binary->left);
        if (left_cmp != 0) return left_cmp;
        return canonical_compare(lhs_binary->right, rhs_binary->right);
    }

    if (const auto* lhs_call = expr_cast<FuncCall>(lhs)) {
        const auto* rhs_call = expr_cast<FuncCall>(rhs);
        if (lhs_call->func_id != rhs_call->func_id) {
            const int lp = get_builtin_precedence(lhs_call->func_id);
            const int rp = get_builtin_precedence(rhs_call->func_id);
            if (lp != rp) return lp < rp ? -1 : 1;
            return lhs_call->func_id < rhs_call->func_id ? -1 : 1;
        }
        // Same func_id — for Unknown, compare names
        if (lhs_call->func_id == BuiltinOp::Unknown) {
            const int name_cmp = compare_string_precedence(lhs_call->name, rhs_call->name);
            if (name_cmp != 0) return name_cmp;
        }
        // Same function: compare arguments
        const auto& la = lhs_call->args;
        const auto& ra = rhs_call->args;
        const std::size_t shared = std::min(la.size(), ra.size());
        for (std::size_t i = 0; i < shared; ++i) {
            const int c = canonical_compare(la[i], ra[i]);
            if (c != 0) return c;
        }
        if (la.size() < ra.size()) return -1;
        if (ra.size() < la.size()) return 1;
        return 0;
    }

    const auto lhs_children = term_order_children(lhs);
    const auto rhs_children = term_order_children(rhs);
    const auto shared_size = std::min(lhs_children.size(), rhs_children.size());
    for (std::size_t i = 0; i < shared_size; ++i) {
        const int cmp = canonical_compare(lhs_children[i], rhs_children[i]);
        if (cmp != 0) return cmp;
    }
    if (lhs_children.size() < rhs_children.size()) return -1;
    if (rhs_children.size() < lhs_children.size()) return 1;
    return 0;
}

int compare_head_precedence(ExprPtr lhs, ExprPtr rhs) noexcept {
    if (lhs == rhs) return 0;
    const ExprKind lhs_kind = expr_kind(lhs);
    const ExprKind rhs_kind = expr_kind(rhs);
    const int kind_cmp = compare_kind_precedence(lhs_kind, rhs_kind);
    if (kind_cmp != 0) {
        return kind_cmp;
    }

    if (const auto* lhs_unary = expr_cast<Unary>(lhs)) {
        const auto* rhs_unary = expr_cast<Unary>(rhs);
        if (lhs_unary->op != rhs_unary->op) {
            return lhs_unary->op < rhs_unary->op ? -1 : 1;
        }
    }

    if (const auto* lhs_binary = expr_cast<Binary>(lhs)) {
        const auto* rhs_binary = expr_cast<Binary>(rhs);
        if (lhs_binary->op != rhs_binary->op) {
            const int lp = get_binary_op_precedence(lhs_binary->op);
            const int rp = get_binary_op_precedence(rhs_binary->op);
            if (lp != rp) return lp < rp ? -1 : 1;
            return lhs_binary->op < rhs_binary->op ? -1 : 1;
        }
    }

    if (const auto* lhs_call = expr_cast<FuncCall>(lhs)) {
        const auto* rhs_call = expr_cast<FuncCall>(rhs);
        if (lhs_call->func_id != rhs_call->func_id) {
            const int lp = get_builtin_precedence(lhs_call->func_id);
            const int rp = get_builtin_precedence(rhs_call->func_id);
            if (lp != rp) return lp < rp ? -1 : 1;
            return lhs_call->func_id < rhs_call->func_id ? -1 : 1;
        }
        return compare_string_precedence(lhs_call->name, rhs_call->name);
    }

    if (const auto* lhs_symbol = expr_cast<Symbol>(lhs)) {
        const auto* rhs_symbol = expr_cast<Symbol>(rhs);
        return compare_string_precedence(lhs_symbol->name, rhs_symbol->name);
    }

    return 0;
}

namespace {

[[nodiscard]] bool term_order_ge(ExprPtr lhs, ExprPtr rhs);

[[nodiscard]] bool term_order_gt(ExprPtr lhs, ExprPtr rhs) {
    if (!lhs || !rhs || lhs == rhs) {
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
    return lhs == rhs || term_order_gt(lhs, rhs);
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
    if (lhs == rhs) {
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

TermOrderRelation compare_rewrite_terms_impl(ExprPtr lhs, ExprPtr rhs) {
    if (!lhs || !rhs) {
        return lhs == rhs ? TermOrderRelation::Equivalent : TermOrderRelation::Incomparable;
    }
    if (lhs == rhs) {
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
