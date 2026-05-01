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

[[nodiscard]] bool is_odd_parity_function(const std::string& name) {
    return name == "sin" || name == "tan" || name == "cot" || name == "csc" ||
           name == "sinh" || name == "tanh" || name == "coth";
}

[[nodiscard]] bool is_even_parity_function(const std::string& name) {
    return name == "cos" || name == "sec" || name == "cosh";
}

[[nodiscard]] bool is_parity_rewrite_function(const std::string& name) {
    return is_odd_parity_function(name) || is_even_parity_function(name);
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

[[nodiscard]] int kind_rank(ExprKind kind) noexcept {
    switch (kind) {
    case ExprKind::Null: return 0;
    case ExprKind::IntegerLit:
    case ExprKind::RationalLit:
    case ExprKind::DecimalLit:
    case ExprKind::Constant: return 1;
    case ExprKind::Symbol: return 2;
    case ExprKind::Unary: return 3;
    case ExprKind::Binary: return 4;
    case ExprKind::FuncCall: return 5;
    case ExprKind::Product: return 6;
    case ExprKind::Sum: return 7;
    default: return 8;
    }
}

[[nodiscard]] int compare_bigint(const BigInt& lhs, const BigInt& rhs) noexcept {
    if (lhs < rhs) return -1;
    if (rhs < lhs) return 1;
    return 0;
}

[[nodiscard]] int compare_expr_vectors(const std::vector<ExprPtr>& lhs, const std::vector<ExprPtr>& rhs) noexcept {
    const auto shared_size = std::min(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < shared_size; ++i) {
        const int cmp = canonical_compare(lhs[i], rhs[i]);
        if (cmp != 0) return cmp;
    }
    if (lhs.size() < rhs.size()) return -1;
    if (rhs.size() < lhs.size()) return 1;
    return 0;
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

int canonical_compare(ExprPtr lhs, ExprPtr rhs) noexcept {
    if (lhs == rhs) return 0;
    if (!lhs) return -1;
    if (!rhs) return 1;

    const ExprKind lhs_kind = expr_kind(lhs);
    const ExprKind rhs_kind = expr_kind(rhs);
    const int lhs_rank = detail::kind_rank(lhs_kind);
    const int rhs_rank = detail::kind_rank(rhs_kind);
    if (lhs_rank != rhs_rank) return lhs_rank < rhs_rank ? -1 : 1;
    if (lhs_kind != rhs_kind) return static_cast<int>(lhs_kind) < static_cast<int>(rhs_kind) ? -1 : 1;

    if (const auto* lhs_integer = expr_cast<IntegerLit>(lhs)) {
        return detail::compare_bigint(lhs_integer->value, expr_cast<IntegerLit>(rhs)->value);
    }

    if (const auto* lhs_rational = expr_cast<RationalLit>(lhs)) {
        const auto* rhs_rational = expr_cast<RationalLit>(rhs);
        const int num_cmp = detail::compare_bigint(lhs_rational->numerator, rhs_rational->numerator);
        if (num_cmp != 0) return num_cmp;
        return detail::compare_bigint(lhs_rational->denominator, rhs_rational->denominator);
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
        if (lhs_call->name < rhs_call->name) return -1;
        if (lhs_call->name > rhs_call->name) return 1;
        if (lhs_call->args.size() < rhs_call->args.size()) return -1;
        if (lhs_call->args.size() > rhs_call->args.size()) return 1;
        for (std::size_t i = 0; i < lhs_call->args.size(); ++i) {
            int cmp = canonical_compare(lhs_call->args[i], rhs_call->args[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

    if (const auto* lhs_product = expr_cast<Product>(lhs)) return detail::compare_expr_vectors(lhs_product->factors, expr_cast<Product>(rhs)->factors);
    if (const auto* lhs_sum = expr_cast<Sum>(lhs)) return detail::compare_expr_vectors(lhs_sum->terms, expr_cast<Sum>(rhs)->terms);

    return structural_equal(lhs, rhs) ? 0 : (expr_kind_name(lhs_kind) < expr_kind_name(rhs_kind) ? -1 : 1);
}

} // namespace cas::symbolic
