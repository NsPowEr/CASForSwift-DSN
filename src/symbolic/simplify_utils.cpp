#include "simplify_impl.hpp"
#include <algorithm>
#include <limits>

namespace cas::symbolic {

// Forward declaration if needed, but it's defined in simplify_impl.hpp or similar headers.
int canonical_compare(ExprPtr lhs, ExprPtr rhs) noexcept;

namespace detail {

thread_local int simplification_depth = 0;

DepthGuard::DepthGuard() { ++simplification_depth; }
DepthGuard::~DepthGuard() { --simplification_depth; }
bool DepthGuard::exceeded() const { return simplification_depth > MAX_SIMPLIFICATION_DEPTH; }

[[nodiscard]] bool is_odd_parity_function(BuiltinOp op) {
    return op == BuiltinOp::Sin || op == BuiltinOp::Tan || op == BuiltinOp::Cot || op == BuiltinOp::Csc ||
           op == BuiltinOp::Sinh || op == BuiltinOp::Tanh || op == BuiltinOp::Coth;
}

[[nodiscard]] bool is_even_parity_function(BuiltinOp op) {
    return op == BuiltinOp::Cos || op == BuiltinOp::Sec || op == BuiltinOp::Cosh;
}

[[nodiscard]] bool is_parity_rewrite_function(BuiltinOp op) {
    return is_odd_parity_function(op) || is_even_parity_function(op);
}

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] ExprPtr make_integer(AstArena& arena, BigInt value) {
    return arena.make<IntegerLit>(std::move(value));
}

[[nodiscard]] ExprPtr make_rational(AstArena& arena, Rational value) {
    if (value.is_integer()) {
        return make_integer(arena, value.numerator());
    }
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] Result<ExprPtr> make_rational_result(AstArena& arena, Rational value) {
    return ok(make_rational(arena, std::move(value)));
}

[[nodiscard]] ExprPtr make_constant(AstArena& arena, MathConstant value) {
    return arena.make<Constant>(value);
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* integer = expr_cast<IntegerLit>(expr)) return integer->value.is_zero();
    if (const auto* rational = expr_cast<RationalLit>(expr)) return rational->numerator.is_zero();
    return false;
}

[[nodiscard]] bool is_one_expr(ExprPtr expr) {
    if (!expr) return false;
    static const BigInt kOne(1);
    if (const auto* integer = expr_cast<IntegerLit>(expr)) return integer->value == kOne;
    if (const auto* rational = expr_cast<RationalLit>(expr)) return rational->numerator == kOne && rational->denominator == kOne;
    return false;
}

[[nodiscard]] bool is_constant_expr(ExprPtr expr, MathConstant constant) {
    const auto* value = expr_cast<Constant>(expr);
    return value != nullptr && value->value == constant;
}

[[nodiscard]] Result<bool> try_get_exact_rational(ExprPtr expr, LiteralRational& out) {
    if (!expr) return ok(false);
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        out = LiteralRational{.value = Rational(integer->value), .exact = true};
        return ok(true);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        auto exact = Rational::make(rational->numerator, rational->denominator);
        if (exact.is_error()) return fail<bool>(exact.error());
        out = LiteralRational{.value = std::move(exact.value()), .exact = true};
        return ok(true);
    }
    return ok(false);
}

[[nodiscard]] std::optional<BigInt> try_get_integer_exponent(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) return integer->value;
    return std::nullopt;
}

[[nodiscard]] Rational decimal_to_rational(const DecimalLit& node) {
    std::string text = node.text;
    std::size_t dot_pos = text.find('.');
    if (dot_pos == std::string::npos) {
        return Rational(BigInt::parse(text).value());
    }
    
    std::size_t e_pos = text.find('e');
    if (e_pos == std::string::npos) e_pos = text.find('E');
    
    std::string mantissa;
    if (e_pos != std::string::npos) {
        mantissa = text.substr(0, e_pos);
    } else {
        mantissa = text;
    }
    
    std::string fraction_part = mantissa.substr(dot_pos + 1);
    mantissa.erase(dot_pos, 1);
    
    BigInt numerator = BigInt::parse(mantissa).value();
    BigInt denominator = BigInt::parse(std::string("1") + std::string(fraction_part.length(), '0')).value();
    
    if (e_pos != std::string::npos) {
        int exp = std::stoi(text.substr(e_pos + 1));
        if (exp > 0) {
            numerator *= BigInt::parse(std::string("1") + std::string(exp, '0')).value();
        } else if (exp < 0) {
            denominator *= BigInt::parse(std::string("1") + std::string(-exp, '0')).value();
        }
    }
    
    return Rational(numerator, denominator);
}

[[nodiscard]] Rational pow_rational_nonnegative(Rational base, BigInt exponent) {
    Rational result(BigInt(1));
    BigInt remaining = std::move(exponent);
    while (!remaining.is_zero()) {
        if ((remaining % BigInt(2)) == BigInt(1)) result *= base;
        remaining /= BigInt(2);
        if (!remaining.is_zero()) base *= base;
    }
    return result;
}

[[nodiscard]] bool is_known_positive_constant(MathConstant value) noexcept {
    switch (value) {
    case MathConstant::Pi:
    case MathConstant::E:
    case MathConstant::Infinity: return true;
    default: return false;
    }
}

[[nodiscard]] bool is_known_nonnegative_constant(MathConstant value) noexcept {
    return is_known_positive_constant(value);
}

[[nodiscard]] bool expr_ptr_sequence_identical(const std::vector<ExprPtr>& lhs, const std::vector<ExprPtr>& rhs) noexcept {
    if (lhs.size() != rhs.size()) return false;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i] != rhs[i]) return false;
    }
    return true;
}

[[nodiscard]] int saturating_add_degree(int lhs, int rhs) noexcept {
    if (lhs < 0 || rhs < 0) return -1;
    if (lhs > std::numeric_limits<int>::max() - rhs) return std::numeric_limits<int>::max();
    return lhs + rhs;
}

template <typename UInt>
[[nodiscard]] std::optional<UInt> parse_bounded_unsigned_decimal(const std::string& decimal) {
    UInt value = 0;
    for (char ch : decimal) {
        const unsigned int digit = static_cast<unsigned int>(ch - '0');
        if (value > (std::numeric_limits<UInt>::max() - static_cast<UInt>(digit)) / static_cast<UInt>(10)) return std::nullopt;
        value = static_cast<UInt>(value * static_cast<UInt>(10) + static_cast<UInt>(digit));
    }
    return value;
}

[[nodiscard]] int polynomial_degree(ExprPtr expr) noexcept {
    if (!expr) return -1;
    if (expr_is<IntegerLit>(expr) || expr_is<RationalLit>(expr) || expr_is<DecimalLit>(expr) || expr_is<Constant>(expr)) return 0;
    if (expr_is<Symbol>(expr)) return 1;
    if (const auto* product = expr_cast<Product>(expr)) {
        int degree = 0;
        for (ExprPtr factor : product->factors) degree = saturating_add_degree(degree, polynomial_degree(factor));
        return degree;
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Pow) {
            const auto exponent = try_get_integer_exponent(binary->right);
            if (exponent.has_value() && !exponent->is_negative()) {
                if (expr_is<Symbol>(binary->left)) {
                    const auto parsed = parse_bounded_unsigned_decimal<unsigned int>(exponent->decimal());
                    if (!parsed.has_value() || parsed.value() > static_cast<unsigned int>(std::numeric_limits<int>::max())) return std::numeric_limits<int>::max();
                    return static_cast<int>(parsed.value());
                }
            }
        }
    }
    return 1;
}

} // namespace detail

// canonical_compare implementation has been moved to term_order.cpp for unification.

} // namespace cas::symbolic
