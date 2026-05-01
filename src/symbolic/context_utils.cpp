#include "context_impl.hpp"
#include <algorithm>
#include <type_traits>

namespace cas::symbolic::detail {

[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr) {
    if (!expr) return std::nullopt;
    if (const auto* integer = expr_cast<IntegerLit>(expr)) return Rational(integer->value);
    if (const auto* rational = expr_cast<RationalLit>(expr)) return Rational(rational->numerator, rational->denominator);
    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) {
        if (auto inner = exact_scalar_from_expr(unary->operand); inner.has_value()) return -(*inner);
    }
    return std::nullopt;
}

[[nodiscard]] ExprPtr negate_expr(ExprPtr expr, AstArena& arena) {
    if (!expr) return expr;
    if (const auto* integer = expr_cast<IntegerLit>(expr)) return arena.make<IntegerLit>(-integer->value);
    if (const auto* rational = expr_cast<RationalLit>(expr)) return arena.make<RationalLit>(-rational->numerator, rational->denominator);
    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) return unary->operand;
    return arena.make<Unary>(UnaryOp::Neg, expr);
}

void append_difference_terms(ExprPtr expr, bool negate, std::vector<ExprPtr>& terms, AstArena& arena) {
    if (!expr) return;
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr term : sum->terms) append_difference_terms(term, negate, terms, arena);
        return;
    }
    terms.push_back(negate ? negate_expr(expr, arena) : expr);
}

[[nodiscard]] int compare_exact_scalars(const Rational& lhs, const Rational& rhs) {
    const BigInt left_cross = lhs.numerator() * rhs.denominator();
    const BigInt right_cross = rhs.numerator() * lhs.denominator();
    return (left_cross < right_cross) ? -1 : (left_cross > right_cross ? 1 : 0);
}


[[nodiscard]] bool is_wildcard_name(const std::string& name) {
    return name.size() > 1U && name.back() == '_';
}

[[nodiscard]] std::size_t expr_weight(ExprPtr expr) {
    if (!expr) return 0U;
    return 1U + visit_expr(expr, [](const auto& node) -> std::size_t {
        using Node = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<Node, IntegerLit> || std::is_same_v<Node, RationalLit> || std::is_same_v<Node, DecimalLit> || std::is_same_v<Node, Symbol> || std::is_same_v<Node, Constant>) return 0U;
        else if constexpr (std::is_same_v<Node, Unary>) return expr_weight(node.operand);
        else if constexpr (std::is_same_v<Node, Binary>) return expr_weight(node.left) + expr_weight(node.right);
        else if constexpr (std::is_same_v<Node, FuncCall>) { std::size_t w = 0; for (ExprPtr a : node.args) w += expr_weight(a); return w; }
        else if constexpr (std::is_same_v<Node, Sum>) { std::size_t w = 0; for (ExprPtr t : node.terms) w += expr_weight(t); return w; }
        else if constexpr (std::is_same_v<Node, Product>) { std::size_t w = 0; for (ExprPtr f : node.factors) w += expr_weight(f); return w; }
        else if constexpr (std::is_same_v<Node, Integral>) return expr_weight(node.integrand) + (node.lower ? expr_weight(*node.lower) : 0) + (node.upper ? expr_weight(*node.upper) : 0);
        else if constexpr (std::is_same_v<Node, Derivative>) return expr_weight(node.expression);
        else if constexpr (std::is_same_v<Node, Limit>) return expr_weight(node.expression) + expr_weight(node.point);
        else if constexpr (std::is_same_v<Node, RootOf>) return expr_weight(node.polynomial);
        else if constexpr (std::is_same_v<Node, Matrix>) { std::size_t w = 0; for (ExprPtr e : node.elements) w += expr_weight(e); return w; }
        return 0U;
    });
}

} // namespace cas::symbolic::detail

namespace cas::symbolic {

[[nodiscard]] bool range_is_exact_zero(ExprPtr lower, ExprPtr upper) {
    auto exact_lower = detail::exact_scalar_from_expr(lower);
    auto exact_upper = detail::exact_scalar_from_expr(upper);
    return exact_lower && exact_upper && exact_lower->numerator().is_zero() && exact_upper->numerator().is_zero();
}

[[nodiscard]] bool exact_range_excludes_zero(ExprPtr lower, ExprPtr upper) {
    auto exact_lower = detail::exact_scalar_from_expr(lower);
    auto exact_upper = detail::exact_scalar_from_expr(upper);
    if (!exact_lower || !exact_upper) return false;
    return detail::compare_exact_scalars(*exact_lower, Rational(BigInt(0))) > 0
        || detail::compare_exact_scalars(*exact_upper, Rational(BigInt(0))) < 0;
}

} // namespace cas::symbolic
