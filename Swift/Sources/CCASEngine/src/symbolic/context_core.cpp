#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <chrono>
#include <functional>
#include <utility>
#include <vector>

namespace cas::symbolic {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr) {
    if (!expr) {
        return std::nullopt;
    }

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }

    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            auto inner = exact_scalar_from_expr(unary->operand);
            if (inner.has_value()) {
                return -(*inner);
            }
        }
    }

    return std::nullopt;
}

[[nodiscard]] ExprPtr negate_expr(ExprPtr expr, AstArena& arena) {
    if (!expr) {
        return expr;
    }

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return arena.make<IntegerLit>(-integer->value);
    }

    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return arena.make<RationalLit>(-rational->numerator, rational->denominator);
    }

    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) {
        return unary->operand;
    }

    return arena.make<Unary>(UnaryOp::Neg, expr);
}

[[nodiscard]] int compare_exact_scalars(const Rational& lhs, const Rational& rhs) {
    const BigInt left_cross = lhs.numerator() * rhs.denominator();
    const BigInt right_cross = rhs.numerator() * lhs.denominator();
    if (left_cross < right_cross) {
        return -1;
    }
    if (left_cross > right_cross) {
        return 1;
    }
    return 0;
}

[[nodiscard]] std::size_t expr_weight(ExprPtr expr) {
    if (!expr) {
        return 0U;
    }

    return 1U + visit_expr(
        expr,
        [](const auto& node) -> std::size_t {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, Symbol> ||
                std::is_same_v<Node, Constant>) {
                return 0U;
            } else if constexpr (std::is_same_v<Node, Unary>) {
                return expr_weight(node.operand);
            } else if constexpr (std::is_same_v<Node, Binary>) {
                return expr_weight(node.left) + expr_weight(node.right);
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::size_t weight = 0U;
                for (ExprPtr arg : node.args) {
                    weight += expr_weight(arg);
                }
                return weight;
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::size_t weight = 0U;
                for (ExprPtr term : node.terms) {
                    weight += expr_weight(term);
                }
                return weight;
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::size_t weight = 0U;
                for (ExprPtr factor : node.factors) {
                    weight += expr_weight(factor);
                }
                return weight;
            } else if constexpr (std::is_same_v<Node, Integral>) {
                return expr_weight(node.integrand) +
                    (node.lower.has_value() ? expr_weight(*node.lower) : 0U) +
                    (node.upper.has_value() ? expr_weight(*node.upper) : 0U);
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                return expr_weight(node.expression);
            } else if constexpr (std::is_same_v<Node, Limit>) {
                return expr_weight(node.expression) + expr_weight(node.point);
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                return expr_weight(node.polynomial);
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::size_t weight = 0U;
                for (ExprPtr element : node.elements) {
                    weight += expr_weight(element);
                }
                return weight;
            } else {
                return 0U;
            }
        });
}

[[nodiscard]] const ComputationTrace& empty_trace() noexcept {
    static const ComputationTrace trace;
    return trace;
}

[[nodiscard]] Result<ExprPtr> materialize_expr_impl(ExprPtr expr, AstArena& arena) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot clone null expression"));
    }

    return ok(visit_expr(
        expr,
        [&](const auto& node) -> ExprPtr {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, Symbol> ||
                std::is_same_v<Node, Constant>) {
                return arena.make<Node>(node);
            } else if constexpr (std::is_same_v<Node, Unary>) {
                return arena.make<Unary>(node.op, materialize_expr_impl(node.operand, arena).value());
            } else if constexpr (std::is_same_v<Node, Binary>) {
                return arena.make<Binary>(
                    node.op,
                    materialize_expr_impl(node.left, arena).value(),
                    materialize_expr_impl(node.right, arena).value());
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> args;
                args.reserve(node.args.size());
                for (ExprPtr arg : node.args) {
                    args.push_back(materialize_expr_impl(arg, arena).value());
                }
                return arena.make<FuncCall>(node.name, std::move(args));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> terms;
                terms.reserve(node.terms.size());
                for (ExprPtr term : node.terms) {
                    terms.push_back(materialize_expr_impl(term, arena).value());
                }
                return arena.make<Sum>(std::move(terms));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> factors;
                factors.reserve(node.factors.size());
                for (ExprPtr factor : node.factors) {
                    factors.push_back(materialize_expr_impl(factor, arena).value());
                }
                return arena.make<Product>(std::move(factors));
            } else if constexpr (std::is_same_v<Node, Integral>) {
                return arena.make<Integral>(
                    materialize_expr_impl(node.integrand, arena).value(),
                    node.variable,
                    node.lower.has_value()
                        ? std::optional<ExprPtr>(materialize_expr_impl(*node.lower, arena).value())
                        : std::nullopt,
                    node.upper.has_value()
                        ? std::optional<ExprPtr>(materialize_expr_impl(*node.upper, arena).value())
                        : std::nullopt);
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                return arena.make<Derivative>(
                    materialize_expr_impl(node.expression, arena).value(),
                    node.variable,
                    node.order);
            } else if constexpr (std::is_same_v<Node, Limit>) {
                return arena.make<Limit>(
                    materialize_expr_impl(node.expression, arena).value(),
                    node.variable,
                    materialize_expr_impl(node.point, arena).value(),
                    node.direction);
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                return arena.make<RootOf>(
                    materialize_expr_impl(node.polynomial, arena).value(),
                    node.variable,
                    node.root_index);
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::vector<ExprPtr> elements;
                elements.reserve(node.elements.size());
                for (ExprPtr element : node.elements) {
                    elements.push_back(materialize_expr_impl(element, arena).value());
                }
                return arena.make<Matrix>(node.rows, node.cols, std::move(elements));
            } else {
                return ExprPtr{};
            }
        }));
}

[[nodiscard]] ExprPtr instantiate_pattern(ExprPtr pattern, const MatchMap& matches, AstArena& arena) {
    if (!pattern) {
        return ExprPtr{};
    }

    if (const auto* wildcard = expr_cast<Symbol>(pattern)) {
        if (is_wildcard_name(wildcard->name)) {
            const auto found = matches.find(wildcard->name);
            if (found != matches.end()) {
                return found->second;
            }
        }
    }

    return visit_expr(
        pattern,
        [&](const auto& node) -> ExprPtr {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (
                std::is_same_v<Node, IntegerLit> ||
                std::is_same_v<Node, RationalLit> ||
                std::is_same_v<Node, DecimalLit> ||
                std::is_same_v<Node, Symbol> ||
                std::is_same_v<Node, Constant>) {
                return arena.make<Node>(node);
            } else if constexpr (std::is_same_v<Node, Unary>) {
                return arena.make<Unary>(node.op, instantiate_pattern(node.operand, matches, arena));
            } else if constexpr (std::is_same_v<Node, Binary>) {
                return arena.make<Binary>(
                    node.op,
                    instantiate_pattern(node.left, matches, arena),
                    instantiate_pattern(node.right, matches, arena));
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> args;
                args.reserve(node.args.size());
                for (ExprPtr arg : node.args) {
                    args.push_back(instantiate_pattern(arg, matches, arena));
                }
                return arena.make<FuncCall>(node.name, std::move(args));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> terms;
                terms.reserve(node.terms.size());
                for (ExprPtr term : node.terms) {
                    terms.push_back(instantiate_pattern(term, matches, arena));
                }
                return arena.make<Sum>(std::move(terms));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> factors;
                factors.reserve(node.factors.size());
                for (ExprPtr factor : node.factors) {
                    factors.push_back(instantiate_pattern(factor, matches, arena));
                }
                return arena.make<Product>(std::move(factors));
            } else if constexpr (std::is_same_v<Node, Integral>) {
                return arena.make<Integral>(
                    instantiate_pattern(node.integrand, matches, arena),
                    node.variable,
                    node.lower.has_value()
                        ? std::optional<ExprPtr>(instantiate_pattern(*node.lower, matches, arena))
                        : std::nullopt,
                    node.upper.has_value()
                        ? std::optional<ExprPtr>(instantiate_pattern(*node.upper, matches, arena))
                        : std::nullopt);
            } else if constexpr (std::is_same_v<Node, Derivative>) {
                return arena.make<Derivative>(
                    instantiate_pattern(node.expression, matches, arena),
                    node.variable,
                    node.order);
            } else if constexpr (std::is_same_v<Node, Limit>) {
                return arena.make<Limit>(
                    instantiate_pattern(node.expression, matches, arena),
                    node.variable,
                    instantiate_pattern(node.point, matches, arena),
                    node.direction);
            } else if constexpr (std::is_same_v<Node, RootOf>) {
                return arena.make<RootOf>(
                    instantiate_pattern(node.polynomial, matches, arena),
                    node.variable,
                    node.root_index);
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                std::vector<ExprPtr> elements;
                elements.reserve(node.elements.size());
                for (ExprPtr element : node.elements) {
                    elements.push_back(instantiate_pattern(element, matches, arena));
                }
                return arena.make<Matrix>(node.rows, node.cols, std::move(elements));
            } else {
                return ExprPtr{};
            }
        });
}

CASContext::CASContext() : rewrite_provider_(&default_rewrite_provider()) {}

void CASContext::define(const Symbol& symbol, ExprPtr value) {
    variables_[symbol.name] = value;
}

std::optional<ExprPtr> CASContext::lookup(const Symbol& symbol) const {
    const auto found = variables_.find(symbol.name);
    if (found == variables_.end()) {
        return std::nullopt;
    }
    return found->second;
}

Assumptions& CASContext::assumptions() noexcept {
    return assumptions_;
}

const Assumptions& CASContext::assumptions() const noexcept {
    return assumptions_;
}

AstArena& CASContext::arena() noexcept {
    return arena_;
}

const AstArena& CASContext::arena() const noexcept {
    return arena_;
}

void CASContext::set_rewrite_provider(const RewriteProvider* provider) noexcept {
    rewrite_provider_ = provider;
}

const RewriteProvider* CASContext::rewrite_provider() const noexcept {
    return rewrite_provider_;
}

void CASContext::enable_trace(bool enabled) noexcept {
    trace_enabled_ = enabled;
    if (!enabled) {
        trace_.clear();
    }
}

const ComputationTrace& CASContext::get_trace() const noexcept {
    return trace_enabled_ ? trace_ : empty_trace();
}

void CASContext::set_timeout(std::chrono::milliseconds timeout) noexcept {
    timeout_ = timeout;
}

Result<ExprPtr> CASContext::simplify(ExprPtr expr) {
    const bool owns_operation = !operation_active_;
    if (owns_operation) {
        operation_active_ = true;
        trace_capture_active_ = trace_enabled_;
        trace_.clear();
        ops_count_ = 0;
        operation_started_at_ = std::chrono::steady_clock::now();
    }
    auto result = symbolic::simplify(expr, *this);
    if (owns_operation) {
        operation_active_ = false;
        trace_capture_active_ = false;
        ops_count_ = 0;
    }
    return result;
}

Result<ExprPtr> materialize_expr(ExprPtr expr, AstArena& arena) {
    return materialize_expr_impl(expr, arena);
}

Result<bool> mathematically_equal(ExprPtr lhs, ExprPtr rhs, CASContext& context) {
    const bool owns_operation = !context.operation_active_;
    if (owns_operation) {
        context.operation_active_ = true;
        context.trace_capture_active_ = false;
        context.trace_.clear();
    }
    
    auto lhs_simplified = context.simplify(lhs);
    if (lhs_simplified.is_error()) return fail<bool>(lhs_simplified.error());
    
    auto rhs_simplified = context.simplify(rhs);
    if (rhs_simplified.is_error()) return fail<bool>(rhs_simplified.error());
    
    bool equal = structural_equal(lhs_simplified.value(), rhs_simplified.value());
    
    if (owns_operation) {
        context.operation_active_ = false;
    }
    
    return ok(equal);
}

} // namespace cas::symbolic
