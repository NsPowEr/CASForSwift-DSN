// Sign predicates for the Simplifier — extracted from simplify_functions.cpp
// (2026-06-19 anti-monolith split). is_known_positive / _nonnegative / _negative
// and is_assumed_nonzero: conservative, assumption-aware sign analysis over the
// AST. Each returns true only when the sign is provable (never guesses), so the
// engine stays correct for symbolic inputs of unknown sign. Dependencies via
// simplify_impl.hpp.

#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

bool Simplifier::is_known_positive(ExprPtr expr) const {
    if (!expr) return false;
    if (assumptions_ != nullptr && assumptions_->is_positive(expr)) return true;
    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value())
        return !rat.value.numerator().is_zero() && !rat.value.numerator().is_negative();
    if (const auto* constant = expr_cast<Constant>(expr))
        return is_known_positive_constant(constant->value);
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) return is_known_positive(bin->left);
        if (bin->op == BinaryOp::Div)
            return is_known_positive(bin->left) && is_known_positive(bin->right);
    }
    // √x > 0 when x > 0 (principal real branch). Mirrors the Pow(x,1/2) case
    // above for the FuncCall(Sqrt) surface form, so e.g. is_known_positive(√3)
    // holds regardless of which canonical form √3 currently has.
    if (const auto* fc = expr_cast<FuncCall>(expr);
        fc != nullptr && fc->func_id == BuiltinOp::Sqrt && fc->args.size() == 1U)
        return is_known_positive(fc->args.front());
    if (const auto* prod = expr_cast<Product>(expr)) {
        if (prod->factors.empty()) return false;
        for (ExprPtr f : prod->factors) if (!is_known_positive(f)) return false;
        return true;
    }
    return false;
}

bool Simplifier::is_known_nonnegative(ExprPtr expr) const {
    if (!expr) return false;
    if (assumptions_ != nullptr && assumptions_->is_nonnegative(expr)) return true;
    if (is_known_positive(expr)) return true;
    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value())
        return !rat.value.numerator().is_negative();
    if (const auto* constant = expr_cast<Constant>(expr))
        return is_known_nonnegative_constant(constant->value);
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Div)
            return is_known_nonnegative(bin->left) && is_known_positive(bin->right);
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        if (prod->factors.empty()) return false;
        for (ExprPtr f : prod->factors) if (!is_known_nonnegative(f)) return false;
        return true;
    }
    return false;
}

bool Simplifier::is_known_negative(ExprPtr expr) const {
    if (!expr) return false;
    if (assumptions_ != nullptr && assumptions_->is_negative(expr)) return true;
    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value())
        return rat.value.numerator().is_negative();
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Div)
            return (is_known_negative(bin->left) && is_known_positive(bin->right))
                || (is_known_positive(bin->left) && is_known_negative(bin->right));
        // base^n < 0 ⟺ base < 0 and n is an odd integer (e.g. (−√3)^(−1) < 0).
        // Even/non-integer powers of a negative base are positive or non-real,
        // so we only claim negativity for the odd-integer case.
        if (bin->op == BinaryOp::Pow && is_known_negative(bin->left)) {
            if (auto n = try_get_integer_exponent(bin->right); n.has_value())
                return (*n % BigInt(2)) != BigInt(0);
        }
    }
    if (const auto* fc = expr_cast<FuncCall>(expr);
        fc != nullptr && fc->func_id == BuiltinOp::Sqrt && fc->args.size() == 1U)
        return false;  // principal √ is never negative (≥ 0 on its real domain)
    if (const auto* prod = expr_cast<Product>(expr)) {
        int negative_count = 0;
        for (ExprPtr f : prod->factors) {
            if (is_known_negative(f)) negative_count++;
            else if (!is_known_positive(f)) return false;
        }
        return (negative_count % 2 != 0);
    }
    if (const auto* unary = expr_cast<Unary>(expr);
        unary != nullptr && unary->op == UnaryOp::Neg)
        return is_known_positive(unary->operand);
    return false;
}

bool Simplifier::is_assumed_nonzero(ExprPtr expr) const {
    if (!expr) return false;
    if (assumptions_ != nullptr && assumptions_->is_nonzero(expr)) return true;
    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value())
        return !rat.value.numerator().is_zero();
    return is_known_positive(expr) || is_known_negative(expr);
}

} // namespace cas::symbolic::detail
