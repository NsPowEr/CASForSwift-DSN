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
                    case MathConstant::EulerGamma: return ok(0.5772156649015328606);
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
            } else if constexpr (std::is_same_v<NodeT, RootOf>) {
                // Newton-Raphson per trovare una radice numerica del polinomio P(x) = 0
                const std::string& var_name = node.variable.name;
                
                // Salvataggio ambiente
                auto old_it = env_.find(var_name);
                std::optional<double> old_val;
                if (old_it != env_.end()) old_val = old_it->second;

                auto f = [&](double val) -> Result<double> {
                    env_[var_name] = val;
                    return evaluate(node.polynomial);
                };

                auto df = [&](double val) -> Result<double> {
                    const double h = 1e-7;
                    auto f1 = f(val + h);
                    if (f1.is_error()) return f1;
                    auto f2 = f(val - h);
                    if (f2.is_error()) return f2;
                    return ok((f1.value() - f2.value()) / (2.0 * h));
                };

                // Guess iniziale: usiamo l'indice per diversificare i punti di partenza
                double x = 0.0;
                if (node.root_index.has_value()) {
                    // Mappa 0 -> 1.0, 1 -> -1.0, 2 -> 2.0, 3 -> -2.0, etc.
                    // Questo aiuta a trovare radici diverse per indici diversi
                    double idx = static_cast<double>(*node.root_index);
                    if (*node.root_index % 2 == 0) x = (idx / 2.0) + 1.0;
                    else x = -((idx + 1.0) / 2.0);
                }
                
                bool found = false;
                for (int attempt = 0; attempt < 10 && !found; ++attempt) {
                    if (attempt > 0) {
                        // Se il primo tentativo fallisce, proviamo una spirale
                        x = (attempt % 2 == 0) ? static_cast<double>(attempt) : -static_cast<double>(attempt);
                    }
                    
                    for (int i = 0; i < 50; ++i) {
                        auto fx = f(x);
                        if (fx.is_error()) break;
                        if (std::abs(fx.value()) < 1e-10) { found = true; break; }

                        auto dfx = df(x);
                        if (dfx.is_error() || std::abs(dfx.value()) < 1e-15) break;

                        double next_x = x - fx.value() / dfx.value();
                        if (std::abs(next_x - x) < 1e-12) { found = true; x = next_x; break; }
                        x = next_x;
                    }
                }

                // Ripristino ambiente
                if (old_val) env_[var_name] = *old_val;
                else env_.erase(var_name);

                if (found) return ok(x);
                return fail<double>(make_error(CASErrorKind::Undefined, "Newton-Raphson failed to converge for RootOf"));
            } else {
                return fail<double>(make_error(CASErrorKind::Unimplemented, "Node type not supported in numeric evaluator"));
            }
        });
}

Result<double> eval(ExprPtr expr, const NumericEnv& env) {
    return NumericEvaluator(env).evaluate(expr);
}

} // namespace cas::numeric
