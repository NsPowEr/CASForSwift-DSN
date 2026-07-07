#include "simplify_impl.hpp"
#include "cas/complex_rational.hpp"
#include <algorithm>
#include <limits>

namespace cas::symbolic {

// Forward declaration if needed, but it's defined in simplify_impl.hpp or similar headers.
int canonical_compare(ExprPtr lhs, ExprPtr rhs) noexcept;

namespace detail {

thread_local int simplification_depth = 0;
thread_local std::vector<ExprPtr> simplify_ancestor_path;

DepthGuard::DepthGuard(int max_depth) : max_depth_(max_depth) { ++simplification_depth; }
DepthGuard::~DepthGuard() { --simplification_depth; }
bool DepthGuard::exceeded() const { return simplification_depth > max_depth_; }

// F7.0-A3.2 async depth propagation.
int current_simplify_depth() noexcept { return simplification_depth; }

AsyncDepthScope::AsyncDepthScope(int inherited_depth) noexcept
    : prev_(simplification_depth) {
    simplification_depth = inherited_depth;
}
AsyncDepthScope::~AsyncDepthScope() noexcept {
    simplification_depth = prev_;
}

}  // namespace detail

// F7.0-A4.2 — strict-canonical-form invariant check.
namespace {

[[nodiscard]] bool is_exact_zero_lit(ExprPtr e) noexcept {
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
    return false;
}
[[nodiscard]] bool is_exact_one_lit(ExprPtr e) noexcept {
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value == BigInt(1);
    if (const auto* rl = expr_cast<RationalLit>(e)) {
        return rl->numerator == BigInt(1) && rl->denominator == BigInt(1);
    }
    return false;
}

}  // namespace

bool is_strictly_canonical(ExprPtr expr) noexcept {
    if (!expr) return true;
    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() < 2U) return false;        // singleton not collapsed
        for (std::size_t i = 0; i < sum->terms.size(); ++i) {
            const ExprPtr t = sum->terms[i];
            if (expr_is<Sum>(t)) return false;            // nested Sum
            if (is_exact_zero_lit(t)) return false;       // zero summand
            if (!is_strictly_canonical(t)) return false;
            if (i > 0) {
                int prev_deg = detail::polynomial_degree(sum->terms[i - 1]);
                int curr_deg = detail::polynomial_degree(t);
                if (prev_deg < curr_deg) return false;
                if (prev_deg == curr_deg && canonical_compare(sum->terms[i - 1], t) > 0) {
                    return false;
                }
            }
        }
        return true;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        if (prod->factors.size() < 2U) return false;
        for (std::size_t i = 0; i < prod->factors.size(); ++i) {
            const ExprPtr f = prod->factors[i];
            if (expr_is<Product>(f)) return false;        // nested Product
            if (is_exact_one_lit(f)) return false;        // one factor
            if (!is_strictly_canonical(f)) return false;
            if (i > 0 && canonical_compare(prod->factors[i - 1], f) > 0) {
                return false;
            }
        }
        return true;
    }
    // Recurse into other node kinds via children. Conservative default true
    // for leaves and structurally-trivial nodes.
    if (const auto* bin = expr_cast<Binary>(expr)) {
        return is_strictly_canonical(bin->left) && is_strictly_canonical(bin->right);
    }
    if (const auto* un = expr_cast<Unary>(expr)) {
        return is_strictly_canonical(un->operand);
    }
    if (const auto* fc = expr_cast<FuncCall>(expr)) {
        for (const auto& a : fc->args) {
            if (!is_strictly_canonical(a)) return false;
        }
        return true;
    }
    return true;
}

namespace detail {

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

[[nodiscard]] ExprPtr make_complex(AstArena& arena, const ComplexRational& value) {
    if (value.is_real()) {
        return make_rational(arena, value.real());
    }
    return arena.make<ComplexLit>(
        value.real().numerator(), value.real().denominator(),
        value.imag().numerator(), value.imag().denominator());
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (!expr) return false;
    if (const auto* integer = expr_cast<IntegerLit>(expr)) return integer->value.is_zero();
    if (const auto* rational = expr_cast<RationalLit>(expr)) return rational->numerator.is_zero();
    if (const auto* complex = expr_cast<ComplexLit>(expr)) return complex->re_num.is_zero() && complex->im_num.is_zero();
    return false;
}

[[nodiscard]] bool is_one_expr(ExprPtr expr) {
    if (!expr) return false;
    static const BigInt kOne(1);
    if (const auto* integer = expr_cast<IntegerLit>(expr)) return integer->value == kOne;
    if (const auto* rational = expr_cast<RationalLit>(expr)) return rational->numerator == kOne && rational->denominator == kOne;
    if (const auto* complex = expr_cast<ComplexLit>(expr))
        return complex->re_num == kOne && complex->re_den == kOne && complex->im_num.is_zero();
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

[[nodiscard]] Result<bool> try_get_exact_complex(ExprPtr expr, LiteralComplex& out) {
    if (!expr) return ok(false);
    if (const auto* complex = expr_cast<ComplexLit>(expr)) {
        auto re = Rational::make(complex->re_num, complex->re_den);
        auto im = Rational::make(complex->im_num, complex->im_den);
        if (re.is_error()) return fail<bool>(re.error());
        if (im.is_error()) return fail<bool>(im.error());
        out = LiteralComplex{.value = ComplexRational(std::move(re.value()), std::move(im.value())), .exact = true};
        return ok(true);
    }
    if (const auto* constant_node = expr_cast<Constant>(expr);
        constant_node != nullptr && constant_node->value == MathConstant::I) {
        out = LiteralComplex{.value = ComplexRational::imag_unit(), .exact = true};
        return ok(true);
    }
    if (const auto* u = expr_cast<Unary>(expr);
        u != nullptr && u->op == UnaryOp::Neg)
    {
        LiteralComplex inner;
        auto rec = try_get_exact_complex(u->operand, inner);
        if (rec.is_ok() && rec.value()) {
            out = LiteralComplex{.value = -inner.value, .exact = true};
            return ok(true);
        }
    }
    LiteralRational lr;
    auto res = try_get_exact_rational(expr, lr);
    if (res.is_ok() && res.value()) {
        out = LiteralComplex{.value = ComplexRational(std::move(lr.value)), .exact = true};
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
    case MathConstant::EulerGamma:
    case MathConstant::Infinity: return true;
    case MathConstant::NegInfinity:
    case MathConstant::ComplexInfinity:
    case MathConstant::Indeterminate:
    case MathConstant::I:
    case MathConstant::NaN: return false;
    }
    return false;
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
    if (expr_is<IntegerLit>(expr) || expr_is<RationalLit>(expr) || expr_is<DecimalLit>(expr) || expr_is<Constant>(expr) || expr_is<ComplexLit>(expr)) return 0;
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
                    const auto exponent_val = exponent->to_u64();
                    if (exponent_val > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) return std::numeric_limits<int>::max();
                    return static_cast<int>(exponent_val);
                }
            }
        }
    }
    return 1;
}

} // namespace detail

// canonical_compare implementation has been moved to term_order.cpp for unification.

} // namespace cas::symbolic
