#include "cas/numeric.hpp"
#include <cmath>
#include <numbers>

namespace cas::numeric {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

} // namespace

Result<double> NumericEvaluator::evaluate(ExprPtr expr) {
    if (!expr) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument, "Cannot evaluate null expression"));
    }

    return visit_expr(
        expr,
        [this](const auto& node) -> Result<double> {
            using NodeT = std::decay_t<decltype(node)>;

            if constexpr (std::is_same_v<NodeT, IntegerLit>) {
                return ok(node.value.to_double());
            } else if constexpr (std::is_same_v<NodeT, RationalLit>) {
                auto num = node.numerator.to_double();
                auto den = node.denominator.to_double();
                return ok(num / den);
            } else if constexpr (std::is_same_v<NodeT, DecimalLit>) {
                return ok(node.to_double());
            } else if constexpr (std::is_same_v<NodeT, Symbol>) {
                auto it = env_.find(node.name);
                if (it != env_.end()) {
                    return ok(it->second);
                }
                return fail<double>(make_error(CASErrorKind::InvalidArgument, "Symbol '" + node.name + "' not found in numeric environment"));
            } else if constexpr (std::is_same_v<NodeT, Constant>) {
                switch (node.value) {
                    case MathConstant::Pi: return ok(std::numbers::pi);
                    case MathConstant::E: return ok(std::numbers::e);
                    case MathConstant::Infinity: return ok(static_cast<double>(INFINITY));
                    default: return fail<double>(make_error(CASErrorKind::Unimplemented, "Constant not supported in numeric evaluator"));
                }
            } else if constexpr (std::is_same_v<NodeT, Unary>) {
                auto op_res = evaluate(node.operand);
                if (op_res.is_error()) return op_res;
                switch (node.op) {
                    case UnaryOp::Neg: return ok(-op_res.value());
                    default: return fail<double>(make_error(CASErrorKind::Unimplemented, "Unary op not supported in numeric evaluator"));
                }
            } else if constexpr (std::is_same_v<NodeT, Binary>) {
                auto left = evaluate(node.left);
                if (left.is_error()) return left;
                auto right = evaluate(node.right);
                if (right.is_error()) return right;

                switch (node.op) {
                    case BinaryOp::Add: return ok(left.value() + right.value());
                    case BinaryOp::Sub: return ok(left.value() - right.value());
                    case BinaryOp::Mul: return ok(left.value() * right.value());
                    case BinaryOp::Div: {
                        if (right.value() == 0.0) return fail<double>(make_error(CASErrorKind::Undefined, "Division by zero in numeric evaluator"));
                        return ok(left.value() / right.value());
                    }
                    case BinaryOp::Pow: {
                        double val = std::pow(left.value(), right.value());
                        if (std::isnan(val)) return fail<double>(make_error(CASErrorKind::Undefined, "Math domain error in power operation"));
                        return ok(val);
                    }
                    default: return fail<double>(make_error(CASErrorKind::Unimplemented, "Binary op not supported in numeric evaluator"));
                }
            } else if constexpr (std::is_same_v<NodeT, Sum>) {
                double total = 0.0;
                for (const auto& term : node.terms) {
                    auto res = evaluate(term);
                    if (res.is_error()) return res;
                    total += res.value();
                }
                return ok(total);
            } else if constexpr (std::is_same_v<NodeT, Product>) {
                double total = 1.0;
                for (const auto& factor : node.factors) {
                    auto res = evaluate(factor);
                    if (res.is_error()) return res;
                    total *= res.value();
                }
                return ok(total);
            } else if constexpr (std::is_same_v<NodeT, FuncCall>) {
                std::vector<double> args;
                for (auto arg : node.args) {
                    auto res = evaluate(arg);
                    if (res.is_error()) return res;
                    args.push_back(res.value());
                }

                if (node.func_id == BuiltinOp::Sin) return ok(std::sin(args[0]));
                if (node.func_id == BuiltinOp::Cos) return ok(std::cos(args[0]));
                if (node.func_id == BuiltinOp::Tan) {
                    double val = std::tan(args[0]);
                    if (std::isinf(val)) return fail<double>(make_error(CASErrorKind::Undefined, "Tan: Singular point"));
                    return ok(val);
                }
                if (node.func_id == BuiltinOp::Exp) return ok(std::exp(args[0]));
                if (node.func_id == BuiltinOp::Ln || node.func_id == BuiltinOp::Log) {
                    if (args[0] <= 0) return fail<double>(make_error(CASErrorKind::Undefined, "Logarithm of non-positive number"));
                    return ok(std::log(args[0]));
                }
                if (node.func_id == BuiltinOp::Sqrt) {
                    if (args[0] < 0) return fail<double>(make_error(CASErrorKind::Undefined, "Square root of negative number"));
                    return ok(std::sqrt(args[0]));
                }
                if (node.func_id == BuiltinOp::Abs) return ok(std::abs(args[0]));
                if (node.name == "atan2") return ok(std::atan2(args[0], args[1]));
                
                return fail<double>(make_error(CASErrorKind::Unimplemented, "Function '" + node.name + "' not supported in numeric evaluator"));
            } else {
                return fail<double>(make_error(CASErrorKind::Unimplemented, "Node type not supported in numeric evaluator"));
            }
        });
}

Result<double> eval(ExprPtr expr, const NumericEnv& env) {
    return NumericEvaluator(env).evaluate(expr);
}

} // namespace cas::numeric
