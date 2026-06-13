#pragma once

#include "cas/ast.hpp"
#include "cas/algebraic_tower.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"

#include <chrono>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas {
namespace symbolic {
class CASContext;
}

namespace algebra {

template <typename Coeff>
class UnivariatePolynomial {
public:
    UnivariatePolynomial() = default;
    explicit UnivariatePolynomial(std::vector<Coeff> coefficients)
        : coefficients_(std::move(coefficients)) {}

    [[nodiscard]] bool empty() const noexcept { return coefficients_.empty(); }
    [[nodiscard]] bool is_zero() const noexcept { return coefficients_.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return coefficients_.size(); }
    [[nodiscard]] std::size_t degree() const noexcept { 
        if (coefficients_.empty()) return 0U; 
        std::size_t d = coefficients_.size() - 1U;
        // In case of hidden zeroes at the end, though normalize() should prevent this.
        return d;
    }
    [[nodiscard]] Coeff constant_term() const { return coefficients_.empty() ? Coeff{} : coefficients_.front(); }
    [[nodiscard]] const Coeff& leading_coeff() const { 
        static const Coeff kZero{};
        return coefficients_.empty() ? kZero : coefficients_.back(); 
    }
    [[nodiscard]] const std::vector<Coeff>& coefficients() const noexcept { return coefficients_; }
    [[nodiscard]] std::vector<Coeff>& coefficients() noexcept { return coefficients_; }

    void reserve(std::size_t size) { coefficients_.reserve(size); }
    void resize(std::size_t size, const Coeff& value = Coeff{}) { coefficients_.resize(size, value); }
    void push_back(Coeff value) { coefficients_.push_back(std::move(value)); }

    Coeff& operator[](std::size_t index) { return coefficients_[index]; }
    const Coeff& operator[](std::size_t index) const { 
        if (index >= coefficients_.size()) {
            static const Coeff kZero{};
            return kZero;
        }
        return coefficients_[index]; 
    }

    template <typename ZeroPredicate>
    void normalize(ZeroPredicate&& is_zero) {
        while (!coefficients_.empty() && is_zero(coefficients_.back())) {
            coefficients_.pop_back();
        }
    }

private:
    std::vector<Coeff> coefficients_;
};

using IntPoly = UnivariatePolynomial<BigInt>;
using RatPoly = UnivariatePolynomial<Rational>;
using PolyExpr = UnivariatePolynomial<ExprPtr>;

struct IntegerExponent {
    std::size_t magnitude{0U};
    bool negative{false};
};

struct PolyDivisionResult {
    PolyExpr quotient;
    PolyExpr remainder;
};

struct RationalRootCandidate {
    BigInt numerator;
    BigInt denominator;
};

enum class IntegerGcdPath {
    Subresultant,
    PrimitiveFallback,
    PrimitiveFallbackPsi,
    PrimitiveFallbackBeta,
    HalfGcd,    // A4: Half-GCD (Knuth-Schönhage O(M(n)log n)) dispatch path
    ModularCrt, // B2.1: Multi-prime modular GCD via CRT (F2.1)
};

struct IntegerGcdResult {
    IntPoly gcd;
    IntegerGcdPath path{IntegerGcdPath::PrimitiveFallback};
};

[[nodiscard]] BigInt bigint_pow_nonnegative(BigInt base, std::size_t exponent);
[[nodiscard]] ExprPtr poly_make_integer(AstArena& arena, long long value);

struct PolyXGCDResult {
    PolyExpr gcd;
    PolyExpr s; // s*a + t*b = gcd
    PolyExpr t;
};
[[nodiscard]] Result<PolyXGCDResult> poly_extended_gcd(
    const PolyExpr& a, const PolyExpr& b, symbolic::CASContext& ctx);
[[nodiscard]] bool poly_is_zero_expr(ExprPtr expr);
[[nodiscard]] bool poly_is_one_expr(ExprPtr expr);
[[nodiscard]] bool poly_is_minus_one_expr(ExprPtr expr);
[[nodiscard]] bool poly_depends_on(ExprPtr expr, const std::string& variable_name);
[[nodiscard]] bool poly_contains_decimal_literal(ExprPtr expr);
[[nodiscard]] Result<ExprPtr> poly_simplify_expr(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> poly_clone_into_context(ExprPtr expr, symbolic::CASContext& ctx);
[[nodiscard]] PolyExpr poly_make_monomial(ExprPtr coefficient, std::size_t degree);
[[nodiscard]] Result<PolyExpr> poly_make_constant_poly(ExprPtr coefficient, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_add(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_negate(const PolyExpr& poly, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_subtract(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_multiply(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_divide_by_scalar(const PolyExpr& poly, ExprPtr scalar, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> poly_pow(PolyExpr base, std::size_t exponent, symbolic::CASContext& ctx);
[[nodiscard]] Result<std::size_t> poly_parse_nonnegative_integer_exponent(ExprPtr expr);

template <typename UInt>
[[nodiscard]] std::optional<UInt> parse_bounded_unsigned_decimal(const std::string& decimal) {
    static_assert(std::numeric_limits<UInt>::is_integer && !std::numeric_limits<UInt>::is_signed);

    UInt value = 0;
    for (char ch : decimal) {
        const unsigned int digit = static_cast<unsigned int>(ch - '0');
        if (value > (std::numeric_limits<UInt>::max() - static_cast<UInt>(digit)) / static_cast<UInt>(10)) {
            return std::nullopt;
        }
        value = static_cast<UInt>(value * static_cast<UInt>(10) + static_cast<UInt>(digit));
    }
    return value;
}

void normalize_poly(PolyExpr& poly);
[[nodiscard]] bool is_zero_poly(const PolyExpr& poly);
[[nodiscard]] std::size_t poly_degree(const PolyExpr& poly);
[[nodiscard]] ExprPtr leading_coefficient(const PolyExpr& poly);
[[nodiscard]] Result<PolyDivisionResult> divide_poly_with_remainder(
    const PolyExpr& dividend,
    const PolyExpr& divisor,
    symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> normalize_poly_monic(const PolyExpr& poly, symbolic::CASContext& ctx);
[[nodiscard]] Result<PolyExpr> parse_polynomial(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> polynomial_to_expr(const PolyExpr& poly, const Symbol& var, symbolic::CASContext& ctx);

[[nodiscard]] Result<IntPoly> poly_to_integer_poly(const PolyExpr& poly);
[[nodiscard]] Result<std::vector<BigInt>> poly_to_integer_coefficients(const PolyExpr& poly);
void normalize_integer_poly(IntPoly& coefficients);
[[nodiscard]] BigInt integer_content(const IntPoly& coefficients);
void divide_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar);
void multiply_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar);
[[nodiscard]] IntPoly primitive_integer_poly(IntPoly coefficients);
[[nodiscard]] IntPoly pseudo_remainder_integer_poly(const IntPoly& dividend, const IntPoly& divisor);
[[nodiscard]] bool try_divide_integer_coefficients_by_scalar(IntPoly& coefficients, const BigInt& scalar);
[[nodiscard]] Result<ExprPtr> integer_coefficients_to_expr(const IntPoly& poly, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] bool is_zero_integer_poly(const IntPoly& coefficients);
[[nodiscard]] IntPoly gcd_integer_poly_primitive(IntPoly lhs, IntPoly rhs);
[[nodiscard]] IntegerGcdResult gcd_integer_poly_with_subresultant(IntPoly lhs, IntPoly rhs);
// Dispatches to Half-GCD when min(deg) > ctx.half_gcd_degree_threshold(), else subresultant.
[[nodiscard]] IntegerGcdResult gcd_integer_poly_dispatch(
    IntPoly lhs, IntPoly rhs, const symbolic::CASContext& ctx);
[[nodiscard]] std::optional<RationalRootCandidate> find_rational_root_candidate(const IntPoly& coefficients);
[[nodiscard]] std::vector<BigInt> positive_divisors_or_one(const BigInt& value);
[[nodiscard]] Result<RatPoly> poly_to_rational_poly(const PolyExpr& poly);
[[nodiscard]] Result<std::vector<Rational>> poly_to_rational_coefficients(const PolyExpr& poly);
[[nodiscard]] Rational evaluate_integer_polynomial_at(const IntPoly& coefficients, const Rational& value);
[[nodiscard]] Rational evaluate_rational_polynomial_at(const RatPoly& coefficients, const Rational& value);

template <typename Coeff>
[[nodiscard]] Result<std::vector<Coeff>> pseudo_remainder_generic(
    std::vector<Coeff> a,
    const std::vector<Coeff>& b,
    symbolic::CASContext* ctx = nullptr);

// F3.5-DEBT-01 fix: optional wall-clock deadline for interruptible Res.
// When the subresultant chain exceeds the deadline, returns Unimplemented.
using ResultantDeadline = std::optional<std::chrono::steady_clock::time_point>;

template <typename Coeff>
[[nodiscard]] Result<Coeff> resultant_generic(
    std::vector<Coeff> a,
    std::vector<Coeff> b,
    symbolic::CASContext* ctx = nullptr,
    const ResultantDeadline& deadline = std::nullopt);

void normalize_rational_coefficients(RatPoly& coefficients);
[[nodiscard]] RatPoly add_rational_poly(const RatPoly& a, const RatPoly& b);
[[nodiscard]] RatPoly sub_rational_poly(const RatPoly& a, const RatPoly& b);
[[nodiscard]] RatPoly mul_rational_poly(const RatPoly& a, const RatPoly& b);
[[nodiscard]] std::pair<RatPoly, RatPoly> div_rem_rational_poly(const RatPoly& a, const RatPoly& b);
[[nodiscard]] std::tuple<RatPoly, RatPoly, RatPoly> extended_gcd_rational_poly(const RatPoly& a, const RatPoly& b);
[[nodiscard]] inline std::tuple<RatPoly, RatPoly, RatPoly> bezout_polynomials(const RatPoly& a, const RatPoly& b) {
    return extended_gcd_rational_poly(a, b);
}

// Hensel Lifting
[[nodiscard]] Result<std::pair<IntPoly, IntPoly>> hensel_lift(
    const IntPoly& f, 
    const IntPoly& g, 
    const IntPoly& h, 
    const BigInt& p, 
    std::size_t k);

[[nodiscard]] Result<std::vector<IntPoly>> hensel_lift_multi(
    const IntPoly& f,
    const std::vector<IntPoly>& factors,
    const BigInt& p,
    std::size_t k);

[[nodiscard]] Result<IntPoly> exact_divide_integer_poly(
    const IntPoly& dividend,
    const IntPoly& divisor,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<IntPoly>> factorize_univariate_hensel_or_kronecker(
    const IntPoly& f,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<IntPoly>> factorize_kronecker(
    const IntPoly& f,
    symbolic::CASContext& ctx);

// LLL and Lattice
using LatticeVector = std::vector<Rational>;
using LatticeMatrix = std::vector<LatticeVector>;
void lll_reduction(LatticeMatrix& b, double delta = 0.75);

[[nodiscard]] std::optional<IntPoly> find_factor_lll(
    const IntPoly& f,
    const IntPoly& g,
    const BigInt& pk,
    std::size_t max_deg,
    double delta = 0.75);

[[nodiscard]] std::optional<IntPoly> find_factor_by_hensel_recombination(
    const IntPoly& f,
    const std::vector<IntPoly>& modular_factors,
    const BigInt& prime,
    std::size_t lift_steps,
    std::size_t max_degree);

// Van Hoeij knapsack-lattice recombination (van Hoeij 2002 §2-4, Belabas 2004 §4).
// Finds one factor of f ∈ Z[x] from r lifted modular factors mod p^a = pk.
// Uses Newton power sum lattice (additive in selection): finding short vectors
// in the LLL-reduced lattice identifies {0,1}^r selection vectors for true factors.
// Bound: Mignotte (TAOCP §4.6.2 Thm F): ||p_k(h)||_∞ ≤ n·2^{n/2}·||f||_∞ < p^a/2.
// Returns nullopt if no factor found within the lattice (caller should report
// Unimplemented if all paths exhausted, not fall back silently).
// Primary path for |modular_factors| ≥ van_hoeij_threshold (CASContext default 8,
// where 2^8 = 256 enumeration budget is the crossover point).
// lll_threshold: r above which LLL knapsack is used instead of enumeration.
//   Default 10 (C(10,5)=252 with Mignotte pruning always fast).
//   Set to 0 to force LLL for all r (useful in tests to verify LLL path).
[[nodiscard]] std::optional<IntPoly> van_hoeij_knapsack_factor(
    const IntPoly& f,
    const std::vector<IntPoly>& modular_factors,
    const BigInt& pk,
    double delta = 0.75,
    std::size_t lll_threshold = 10U,
    symbolic::CASContext* ctx = nullptr);

// Modular Factoring
[[nodiscard]] Result<std::vector<IntPoly>> factor_polynomial_mod_p(IntPoly f, const BigInt& p, symbolic::CASContext* ctx = nullptr);

// Half-GCD (Knuth-Schönhage O(M(n)log n)) for integer polynomials (A4, F2 Block A).
// Dispatcher calls this when min(deg(a), deg(b)) > ctx.half_gcd_degree_threshold().
[[nodiscard]] IntPoly half_gcd_integer_poly(const IntPoly& a, const IntPoly& b);

// Berlekamp factorization in Fp[x] (A3, F2 Block A).
// Alternative to Cantor-Zassenhaus; exact for small p where deg(f)*p ≤ max_matrix_size.
// Returns Unimplemented (reason BERLEKAMP_MATRIX_TOO_LARGE) for larger inputs;
// caller should fall back to factor_polynomial_mod_p (Cantor-Zassenhaus).
// Pass ctx.max_berlekamp_matrix_size() as max_matrix_size; default 1024 matches
// the CASContext default so callers without a context can use the bare signature.
[[nodiscard]] Result<std::vector<IntPoly>> berlekamp_factor_mod_p(
    IntPoly f, const BigInt& p, std::size_t max_matrix_size = 1024U);

// Cyclotomic Polynomials
// max_n: upper bound on n to check; -1 = derive from degree (φ(n)=d → n ≤ 2*(d+1))
[[nodiscard]] std::optional<int> is_cyclotomic(const IntPoly& poly, int max_n = -1);
[[nodiscard]] std::vector<ExprPtr> cyclotomic_roots(int n, const Symbol& var, AstArena& arena);
// Φ_n(x): n-th cyclotomic polynomial over Z.  Returns empty IntPoly for n > 2^20.
[[nodiscard]] IntPoly compute_cyclotomic(int n);

// B2.1 — Modular GCD via multi-prime CRT (F2.1).
// Returns ok(gcd) certified by divisibility check; Unimplemented if prime budget exhausted.
// Bound: Mignotte (vzGG §6.7).  Bad-prime detection: deg-stability + lc-divisibility.
[[nodiscard]] Result<IntPoly> gcd_integer_poly_crt(
    const IntPoly& f, const IntPoly& g, const symbolic::CASContext& ctx);

// B2.2 — Modular resultant via multi-prime CRT (F2.2, Collins 1971).
// Returns ok(res ∈ Z) or Unimplemented if budget exhausted.
// Bound: Hadamard (vzGG §6.8).
[[nodiscard]] Result<BigInt> resultant_integer_poly_crt(
    const IntPoly& f, const IntPoly& g, const symbolic::CASContext& ctx);

// B2.3 — Bivariate resultant via evaluation-interpolation (vzGG §6.3).
// f_as_y_poly[i] = coefficient of y^i in f(x,y), as a polynomial in x (IntPoly).
// Returns res_y(f,g) ∈ Z[x] or Unimplemented.
// Degree bound: deg_x(res) ≤ deg_x(f)*deg_y(g) + deg_x(g)*deg_y(f) → exact interpolation.
[[nodiscard]] Result<IntPoly> resultant_bivariate_eval_interp(
    const std::vector<IntPoly>& f_as_y_poly,
    const std::vector<IntPoly>& g_as_y_poly,
    const symbolic::CASContext& ctx);

}  // namespace algebra
}  // namespace cas

#include "polynomial_resultant_generic.hpp"
