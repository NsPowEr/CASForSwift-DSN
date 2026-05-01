#include "simplify_impl.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/numeric.hpp"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const FuncCall& node) {
    std::vector<ExprPtr> args;
    args.reserve(node.args.size());
    for (std::size_t i = 0; i < node.args.size(); ++i) {
        auto simplify_arg = [&]() -> Result<ExprPtr> {
            if (trace_enabled_) {
                std::vector<ExprPtr> current_args = node.args;
                for (std::size_t j = 0; j < args.size(); ++j) current_args[j] = args[j];
                ScopedFrame frame(*this, [this, current_args = std::move(current_args), i, name = node.name](ExprPtr value) mutable {
                    current_args[i] = value;
                    return arena_.make<FuncCall>(name, std::move(current_args));
                });
                return simplify_expr(node.args[i]);
            }
            return simplify_expr(node.args[i]);
        };
        auto res = simplify_arg();
        if (res.is_error()) return res;
        args.push_back(res.value());
    }

    const ExprPtr target_before = expr_ptr_sequence_identical(args, node.args) ? original : (trace_enabled_ ? arena_.make<FuncCall>(node.name, args) : ExprPtr{});

    if (rewrite_provider_ != nullptr && may_rewrite_function_call(node.name, args)) {
        ExprPtr rewrite_target = expr_ptr_sequence_identical(args, node.args) ? original : arena_.make<FuncCall>(node.name, args);
        auto rewritten = rewrite_provider_->try_rewrite(rewrite_target, arena_, assumptions_);
        if (rewritten.is_ok() && rewritten.value() != rewrite_target) {
            append_trace(RuleId::RewriteProviderApplied, rewrite_target, rewritten.value());
            return simplify_expr(rewritten.value());
        }
    }

    if (node.name == "sin" && args.size() == 1U) {
        if (is_zero_expr(args.front())) return traced_result(RuleId::SimplifySinZero, target_before, make_integer(arena_, BigInt(0)));
        if (is_constant_expr(args.front(), MathConstant::Pi)) return traced_result(RuleId::SimplifySinPi, target_before, make_integer(arena_, BigInt(0)));
    }
    if (node.name == "cos" && args.size() == 1U) {
        if (is_zero_expr(args.front())) return traced_result(RuleId::SimplifyCosZero, target_before, make_integer(arena_, BigInt(1)));
        if (is_constant_expr(args.front(), MathConstant::Pi)) return traced_result(RuleId::SimplifyCosPi, target_before, make_integer(arena_, BigInt(-1)));
    }
    if (node.name == "exp" && args.size() == 1U) {
        if (is_zero_expr(args.front())) return traced_result(RuleId::SimplifyExpZero, target_before, make_integer(arena_, BigInt(1)));
        if (is_one_expr(args.front())) return traced_result(RuleId::SimplifyExpOne, target_before, make_constant(arena_, MathConstant::E));
        if (const auto* sum = expr_cast<Sum>(args.front())) {
            std::vector<ExprPtr> factors;
            for (ExprPtr term : sum->terms) factors.push_back(arena_.make<FuncCall>("exp", std::vector<ExprPtr>{term}));
            auto rewritten = simplify_product_factors(factors, arena_.make<Product>(factors));
            if (rewritten.is_ok()) { append_trace(RuleId::SimplifyExpSum, target_before, rewritten.value()); return rewritten; }
        }
    }
    if (node.name == "ln" && args.size() == 1U) {
        if (is_one_expr(args.front())) return traced_result(RuleId::SimplifyLnOne, target_before, make_integer(arena_, BigInt(0)));
        if (is_constant_expr(args.front(), MathConstant::E)) return traced_result(RuleId::SimplifyLnE, target_before, make_integer(arena_, BigInt(1)));
        if (const auto* power = expr_cast<Binary>(args.front()); power != nullptr && power->op == BinaryOp::Pow && is_constant_expr(power->left, MathConstant::E)) {
            return traced_result(RuleId::SimplifyLnExp, target_before, power->right);
        }
    }
    if (node.name == "sqrt" && args.size() == 1U) {
        if (const auto* power = expr_cast<Binary>(args.front()); power != nullptr && power->op == BinaryOp::Pow) {
            if (auto exp = try_get_integer_exponent(power->right); exp.has_value() && *exp == BigInt(2)) {
                if (is_known_nonnegative(power->left)) { append_assumption(target_before); return traced_result(RuleId::SimplifySqrtSquare, target_before, power->left); }
                return traced_result(RuleId::SimplifySqrtSquare, target_before, arena_.make<FuncCall>("abs", std::vector<ExprPtr>{power->left}));
            }
        }
    }

    if (context_ != nullptr && args.size() == 1U && expr_is<Matrix>(args.front())) {
        const auto& m_node = expr_ref<Matrix>(args.front());
        cas::linalg::MatrixExpr m_expr(m_node.rows, m_node.cols, m_node.elements);

        if (node.name == "det") {
            auto res = cas::linalg::determinant(m_expr, *context_);
            if (res.is_ok()) return ok(res.value());
        }
        if (node.name == "rank") {
            auto res = cas::linalg::rank(m_expr, *context_);
            if (res.is_ok()) return ok(make_integer(arena_, static_cast<long long>(res.value())));
        }
        if (node.name == "trace") {
            auto res = cas::linalg::trace(m_expr, *context_);
            if (res.is_ok()) return ok(res.value());
        }
        if (node.name == "inv") {
            auto res = cas::linalg::inverse(m_expr, *context_);
            if (res.is_ok()) return ok(arena_.make<Matrix>(res.value().rows(), res.value().cols(), res.value().elements()));
        }
        if (node.name == "transpose") {
            auto res = cas::linalg::transpose(m_expr);
            if (res.is_ok()) return ok(arena_.make<Matrix>(res.value().rows(), res.value().cols(), res.value().elements()));
        }
    }

    if (node.name == "N" && args.size() == 1U) {
        if (const auto* matrix = expr_cast<Matrix>(args.front())) {
            std::vector<ExprPtr> numeric_elements;
            for (auto elem : matrix->elements) {
                auto val = cas::numeric::eval(elem);
                if (val.is_ok()) {
                    numeric_elements.push_back(arena_.make<DecimalLit>(val.value()));
                } else {
                    numeric_elements.push_back(elem);
                }
            }
            return ok(arena_.make<Matrix>(matrix->rows, matrix->cols, std::move(numeric_elements)));
        }

        auto val = cas::numeric::eval(args.front());
        if (val.is_ok()) {
            return ok(arena_.make<DecimalLit>(val.value()));
        }
    }

    if (expr_ptr_sequence_identical(args, node.args)) return ok(original);
    return ok(arena_.make<FuncCall>(node.name, std::move(args)));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Integral& node) { return simplify_passthrough(original, node); }
Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Derivative& node) { return simplify_passthrough(original, node); }
Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Limit& node) { return simplify_passthrough(original, node); }
Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const RootOf& node) { return simplify_passthrough(original, node); }

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Matrix& node) {
    std::vector<ExprPtr> elements;
    elements.reserve(node.elements.size());
    for (std::size_t i = 0; i < node.elements.size(); ++i) {
        auto simplify_el = [&]() -> Result<ExprPtr> {
            if (trace_enabled_) {
                std::vector<ExprPtr> current = node.elements;
                for (std::size_t j = 0; j < elements.size(); ++j) current[j] = elements[j];
                ScopedFrame frame(*this, [this, current = std::move(current), i, r = node.rows, c = node.cols](ExprPtr value) mutable {
                    current[i] = value; return arena_.make<Matrix>(r, c, std::move(current));
                });
                return simplify_expr(node.elements[i]);
            }
            return simplify_expr(node.elements[i]);
        };
        auto res = simplify_el();
        if (res.is_error()) return res;
        elements.push_back(res.value());
    }
    if (expr_ptr_sequence_identical(elements, node.elements)) return ok(original);
    return ok(arena_.make<Matrix>(node.rows, node.cols, std::move(elements)));
}

template <typename Node>
Result<ExprPtr> Simplifier::simplify_passthrough(ExprPtr original, const Node&) { return ok(original); }

bool Simplifier::may_rewrite_function_call(const std::string& name, const std::vector<ExprPtr>& args) const {
    if (args.size() != 1U) return false;
    if (name == "tan") return true;
    if (is_parity_rewrite_function(name)) {
        if (is_zero_expr(args.front())) return name == "sin" || name == "cos" || name == "tan" || name == "sinh" || name == "tanh";
        const auto* unary = expr_cast<Unary>(args.front());
        return unary != nullptr && unary->op == UnaryOp::Neg;
    }
    if (name == "exp") return is_zero_expr(args.front()) || is_one_expr(args.front());
    if (name == "ln") {
        if (is_one_expr(args.front()) || is_constant_expr(args.front(), MathConstant::E)) return true;
        if (const auto* quot = expr_cast<Binary>(args.front()); quot != nullptr && quot->op == BinaryOp::Div) return is_known_positive(quot->left) && is_known_positive(quot->right);
        if (const auto* pow = expr_cast<Binary>(args.front()); pow != nullptr && pow->op == BinaryOp::Pow) return is_constant_expr(pow->left, MathConstant::E) || is_known_positive(pow->left);
        if (const auto* prod = expr_cast<Product>(args.front())) return std::all_of(prod->factors.begin(), prod->factors.end(), [this](ExprPtr f) { return is_known_positive(f); });
        const auto* sqrt = expr_cast<FuncCall>(args.front());
        return sqrt != nullptr && sqrt->name == "sqrt" && sqrt->args.size() == 1U && is_known_positive(sqrt->args.front());
    }
    if (name != "sqrt") {
        if (name == "det" || name == "rank" || name == "trace" || name == "inv" || name == "transpose" || name == "N") {
            return !args.empty() && (expr_is<Matrix>(args.front()) || name == "N");
        }
        return false;
    }
    if (const auto* quot = expr_cast<Binary>(args.front()); quot != nullptr && quot->op == BinaryOp::Div) return is_known_nonnegative(quot->left) && is_known_positive(quot->right);
    const auto* prod = expr_cast<Product>(args.front());
    return prod != nullptr && std::all_of(prod->factors.begin(), prod->factors.end(), [this](ExprPtr f) { return is_known_nonnegative(f); });
}

bool Simplifier::may_rewrite_sum_terms(const std::vector<ExprPtr>& terms) const {
    ExprPtr s_arg, c_arg;
    for (ExprPtr term : terms) {
        const auto* pow = expr_cast<Binary>(term);
        if (pow == nullptr || pow->op != BinaryOp::Pow) continue;
        if (auto exp = try_get_integer_exponent(pow->right); !exp.has_value() || *exp != BigInt(2)) continue;
        const auto* call = expr_cast<FuncCall>(pow->left);
        if (call == nullptr || call->args.size() != 1U) continue;
        if (call->name == "sin") {
            if (c_arg && structural_equal(call->args.front(), c_arg)) return true;
            if (!s_arg) s_arg = call->args.front();
        } else if (call->name == "cos") {
            if (s_arg && structural_equal(call->args.front(), s_arg)) return true;
            if (!c_arg) c_arg = call->args.front();
        }
    }
    return false;
}

bool Simplifier::may_rewrite_power(ExprPtr base, ExprPtr exponent) const {
    if (!is_constant_expr(base, MathConstant::E)) return false;
    const auto* ln = expr_cast<FuncCall>(exponent);
    return ln != nullptr && ln->name == "ln" && ln->args.size() == 1U && is_known_positive(ln->args.front());
}

bool Simplifier::is_known_positive(ExprPtr expr) const {
    if (!expr) return false;
    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value()) return !rat.value.numerator().is_zero() && !rat.value.numerator().is_negative();
    if (const auto* symbol = expr_cast<Symbol>(expr)) return assumptions_ != nullptr && assumptions_->is_positive(*symbol);
    if (const auto* constant = expr_cast<Constant>(expr)) return is_known_positive_constant(constant->value);
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) return is_known_positive(bin->left);
        if (bin->op == BinaryOp::Div) return is_known_positive(bin->left) && is_known_positive(bin->right);
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) if (!is_known_positive(f)) return false;
        return !prod->factors.empty();
    }
    return false;
}

bool Simplifier::is_known_nonnegative(ExprPtr expr) const {
    if (!expr) return false;
    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value()) return !rat.value.numerator().is_negative();
    if (const auto* symbol = expr_cast<Symbol>(expr)) return assumptions_ != nullptr && assumptions_->is_positive(*symbol);
    if (const auto* constant = expr_cast<Constant>(expr)) return is_known_nonnegative_constant(constant->value);
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) return is_known_positive(bin->left);
        if (bin->op == BinaryOp::Div) return is_known_nonnegative(bin->left) && is_known_positive(bin->right);
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) if (!is_known_nonnegative(f)) return false;
        return !prod->factors.empty();
    }
    return false;
}

} // namespace cas::symbolic::detail
