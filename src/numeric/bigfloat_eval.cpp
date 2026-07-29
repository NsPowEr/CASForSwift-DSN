#include "cas/bigfloat.hpp"
#include "cas/numeric.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/error_helpers.hpp"
#include <cmath>
#include <string>
#include <unordered_set>

namespace cas::numeric {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}

} // namespace

// Recursive BigFloat evaluator.
// Follows the same structure as NumericEvaluator but uses BigFloat throughout.
class BigFloatEvaluator {
public:
    explicit BigFloatEvaluator(mpfr_prec_t prec, const NumericEnv& env, std::size_t max_recursion_depth = 256U)
        : prec_(prec), env_(env), max_recursion_depth_(max_recursion_depth) {}

    [[nodiscard]] Result<BigFloat> evaluate(ExprPtr expr) {
        if (!expr)
            return fail<BigFloat>(make_error(CASErrorKind::InvalidArgument,
                "Cannot evaluate null expression"));

        if (current_depth_ >= max_recursion_depth_) {
            return make_unimplemented<BigFloat>(
                "numeric", "BigFloatEvaluator::evaluate",
                "expression recursion depth",
                error::reason_codes::RECURSION_DEPTH_EXCEEDED,
                "Increase max_recursion_depth",
                "A20",
                "BigFloat evaluation recursion depth limit exceeded");
        }

        auto [it, inserted] = active_nodes_.insert(expr);
        if (!inserted) {
            return make_unimplemented<BigFloat>(
                "numeric", "BigFloatEvaluator::evaluate",
                "expression cycle",
                error::reason_codes::CYCLE_DETECTED,
                "Ensure expression AST has no cyclic references",
                "A20",
                "Cyclic evaluation detected during BigFloat evaluation");
        }

        struct ScopeGuard {
            std::size_t& depth;
            std::unordered_set<ExprPtr, ExprHash>& nodes;
            ExprPtr expr;
            ScopeGuard(std::size_t& d, std::unordered_set<ExprPtr, ExprHash>& n, ExprPtr e)
                : depth(d), nodes(n), expr(e) { ++depth; }
            ~ScopeGuard() { --depth; nodes.erase(expr); }
        } scope_guard(current_depth_, active_nodes_, expr);

        return visit_expr(expr,
            [this](const auto& node) -> Result<BigFloat> {
                using NodeT = std::decay_t<decltype(node)>;

                if constexpr (std::is_same_v<NodeT, IntegerLit>) {
                    return ok(BigFloat::from_integer_string(
                        node.value.decimal(), prec_));
                }
                else if constexpr (std::is_same_v<NodeT, RationalLit>) {
                    return ok(BigFloat::from_rational_parts(
                        node.numerator.decimal(),
                        node.denominator.decimal(),
                        prec_));
                }
                else if constexpr (std::is_same_v<NodeT, DecimalLit>) {
                    // Parse decimal string directly for full precision
                    return ok(BigFloat::from_integer_string(node.text, prec_));
                }
                else if constexpr (std::is_same_v<NodeT, Symbol>) {
                    auto it = env_.find(node.name);
                    if (it != env_.end())
                        return ok(BigFloat::from_double(it->second, prec_));
                    return fail<BigFloat>(make_error(CASErrorKind::InvalidArgument,
                        "Symbol '" + node.name + "' not in numeric environment"));
                }
                else if constexpr (std::is_same_v<NodeT, Constant>) {
                    switch (node.value) {
                    case MathConstant::Pi:
                        return ok(BigFloat::pi(prec_));
                    case MathConstant::E:
                        return ok(BigFloat::e(prec_));
                    case MathConstant::EulerGamma:
                        return ok(BigFloat::euler_gamma(prec_));
                    case MathConstant::Infinity:
                        return fail<BigFloat>(make_error(CASErrorKind::Undefined,
                            "Infinity has no finite BigFloat representation"));
                    case MathConstant::NegInfinity:
                        return fail<BigFloat>(make_error(CASErrorKind::Undefined,
                            "NegInfinity has no finite BigFloat representation"));
                    case MathConstant::ComplexInfinity:
                        return fail<BigFloat>(make_error(CASErrorKind::Undefined,
                            "ComplexInfinity has no finite BigFloat representation"));
                    case MathConstant::Indeterminate:
                        return fail<BigFloat>(make_error(CASErrorKind::Undefined,
                            "Indeterminate form has no numeric value"));
                    default:
                        // F0.8-MIGRATED
                        return make_unimplemented<BigFloat>(
                            "numeric", "BigFloatEvaluator::evaluate",
                            "math constant not yet mapped to BigFloat (MPFR) function",
                            error::reason_codes::NUMERIC_UNSUPPORTED_CONSTANT,
                            "Add MPFR constant via BigFloat::from_mpfr_const() for the missing MathConstant",
                            "F1.x");
                    }
                }
                else if constexpr (std::is_same_v<NodeT, Unary>) {
                    auto op = evaluate(node.operand);
                    if (op.is_error()) return op;
                    switch (node.op) {
                    case UnaryOp::Neg: return ok(-op.value());
                    default:
                        // F0.8-MIGRATED
                        return make_unimplemented<BigFloat>(
                            "numeric", "BigFloatEvaluator::evaluate",
                            "unary operator not dispatched in BigFloat evaluator",
                            error::reason_codes::NUMERIC_UNSUPPORTED_UNARY_OP,
                            "Add BigFloat dispatch for the missing UnaryOp in bigfloat_eval.cpp",
                            "F1.x");
                    }
                }
                else if constexpr (std::is_same_v<NodeT, Binary>) {
                    auto lv = evaluate(node.left);
                    if (lv.is_error()) return lv;
                    auto rv = evaluate(node.right);
                    if (rv.is_error()) return rv;
                    switch (node.op) {
                    case BinaryOp::Add: return ok(lv.value() + rv.value());
                    case BinaryOp::Sub: return ok(lv.value() - rv.value());
                    case BinaryOp::Mul: return ok(lv.value() * rv.value());
                    case BinaryOp::Div:
                        if (rv.value().is_zero())
                            return fail<BigFloat>(make_error(CASErrorKind::Undefined,
                                "BigFloat division by zero"));
                        return ok(lv.value() / rv.value());
                    case BinaryOp::Pow:
                        return ok(BigFloat::pow(lv.value(), rv.value()));
                    default:
                        // F0.8-MIGRATED
                        return make_unimplemented<BigFloat>(
                            "numeric", "BigFloatEvaluator::evaluate",
                            "binary operator not dispatched in BigFloat evaluator",
                            error::reason_codes::NUMERIC_UNSUPPORTED_BINARY_OP,
                            "Add BigFloat dispatch for the missing BinaryOp in bigfloat_eval.cpp",
                            "F1.x");
                    }
                }
                else if constexpr (std::is_same_v<NodeT, Sum>) {
                    BigFloat total(prec_);
                    for (auto term : node.terms) {
                        auto r = evaluate(term);
                        if (r.is_error()) return r;
                        total = total + r.value();
                    }
                    return ok(std::move(total));
                }
                else if constexpr (std::is_same_v<NodeT, Product>) {
                    BigFloat total = BigFloat::from_double(1.0, prec_);
                    for (auto factor : node.factors) {
                        auto r = evaluate(factor);
                        if (r.is_error()) return r;
                        total = total * r.value();
                    }
                    return ok(std::move(total));
                }
                else if constexpr (std::is_same_v<NodeT, FuncCall>) {
                    // Evaluate arguments
                    std::vector<BigFloat> args;
                    args.reserve(node.args.size());
                    for (auto a : node.args) {
                        auto r = evaluate(a);
                        if (r.is_error()) return r;
                        args.push_back(std::move(r.value()));
                    }
                    if (args.empty())
                        return fail<BigFloat>(make_error(CASErrorKind::InvalidArgument,
                            "Function with no arguments"));

                    const auto op = node.func_id;
                    if (op == BuiltinOp::Sin)  return ok(BigFloat::sin(args[0]));
                    if (op == BuiltinOp::Cos)  return ok(BigFloat::cos(args[0]));
                    if (op == BuiltinOp::Tan)  {
                        auto r = BigFloat::tan(args[0]);
                        if (r.is_inf())
                            return fail<BigFloat>(make_error(CASErrorKind::Undefined,
                                "tan: singular point"));
                        return ok(std::move(r));
                    }
                    if (op == BuiltinOp::Asin) return ok(BigFloat::asin(args[0]));
                    if (op == BuiltinOp::Acos) return ok(BigFloat::acos(args[0]));
                    if (op == BuiltinOp::Atan) return ok(BigFloat::atan(args[0]));
                    if (op == BuiltinOp::Exp)  return ok(BigFloat::exp(args[0]));
                    if (op == BuiltinOp::Ln || op == BuiltinOp::Log) {
                        if (args[0].is_negative() || args[0].is_zero())
                            return fail<BigFloat>(make_error(CASErrorKind::Undefined,
                                "ln: argument must be positive"));
                        return ok(BigFloat::ln(args[0]));
                    }
                    if (op == BuiltinOp::Sqrt) {
                        if (args[0].is_negative())
                            return fail<BigFloat>(make_error(CASErrorKind::Undefined,
                                "sqrt: argument must be non-negative"));
                        return ok(BigFloat::sqrt(args[0]));
                    }
                    if (op == BuiltinOp::Abs)     return ok(BigFloat::abs(args[0]));
                    if (op == BuiltinOp::Gamma)   return ok(BigFloat::gamma(args[0]));
                    if (op == BuiltinOp::Erf)     return ok(BigFloat::erf(args[0]));
                    if (op == BuiltinOp::Sinh)    return ok(BigFloat::sinh(args[0]));
                    if (op == BuiltinOp::Cosh)    return ok(BigFloat::cosh(args[0]));
                    if (op == BuiltinOp::Tanh)    return ok(BigFloat::tanh(args[0]));

                    // F0.8-MIGRATED
                    return make_unimplemented<BigFloat>(
                        "numeric", "BigFloatEvaluator::evaluate",
                        "function '" + node.name + "' not dispatched in BigFloat evaluator",
                        error::reason_codes::NUMERIC_UNSUPPORTED_FUNCTION,
                        "Add BigFloat::fn() wrapper for '" + node.name + "' via mpfr_fn() in bigfloat.cpp",
                        "F1.x");
                }
                else {
                    // F0.8-MIGRATED
                    return make_unimplemented<BigFloat>(
                        "numeric", "BigFloatEvaluator::evaluate",
                        "AST node type not handled in BigFloat evaluator dispatch",
                        error::reason_codes::NUMERIC_UNSUPPORTED_NODE_TYPE,
                        "Add a constexpr branch for the missing node type in bigfloat_eval.cpp",
                        "F1.x");
                }
            });
    }

private:
    mpfr_prec_t prec_;
    const NumericEnv& env_;
    std::size_t max_recursion_depth_{256U};
    std::size_t current_depth_{0U};
    std::unordered_set<ExprPtr, ExprHash> active_nodes_;
};

// ── Public API ────────────────────────────────────────────────────────────────

Result<BigFloat> eval_bigfloat(ExprPtr expr, mpfr_prec_t prec_bits,
    const NumericEnv& env, std::size_t max_recursion_depth = 256U) {
    return BigFloatEvaluator(prec_bits, env, max_recursion_depth).evaluate(expr);
}

Result<std::string> eval_mpfr(ExprPtr expr, unsigned int decimal_digits,
    const NumericEnv& env) {
    if (decimal_digits == 0) decimal_digits = 15;
    const mpfr_prec_t bits = decimal_digits_to_bits(decimal_digits);
    auto result = eval_bigfloat(expr, bits, env, 256U);
    if (result.is_error()) return fail<std::string>(result.error());
    return ok(result.value().to_string(static_cast<int>(decimal_digits)));
}

// L3-03 overload context-aware: applica precision dal CASContext.
Result<std::string> eval_mpfr(ExprPtr expr, symbolic::CASContext& ctx,
    const NumericEnv& env) {
    unsigned int decimal_digits = ctx.numeric_precision_digits();
    if (decimal_digits == 0) decimal_digits = 15;
    const mpfr_prec_t bits = decimal_digits_to_bits(decimal_digits);
    auto result = eval_bigfloat(expr, bits, env, ctx.max_recursion_depth());
    if (result.is_error()) return fail<std::string>(result.error());
    return ok(result.value().to_string(static_cast<int>(decimal_digits)));
}

} // namespace cas::numeric
