#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <cstdint>
#include <optional>
#include <vector>

namespace cas::algebra {

[[nodiscard]] static std::optional<std::vector<long long>> try_small_divisors(const BigInt& n) {
    if (n.is_zero()) return std::nullopt;
    const BigInt abs_n = n.abs();
    if (abs_n.bit_length() > 30U) return std::nullopt;
    const long long val = static_cast<long long>(abs_n.to_u64());
    if (val <= 0) return std::nullopt;
    std::vector<long long> divs;
    for (long long d = 1; d * d <= val; ++d) {
        if (val % d == 0) {
            divs.push_back(d);
            if (d != val / d) divs.push_back(val / d);
        }
    }
    return divs;
}

[[nodiscard]] static BigInt eval_int_poly_at_int(const IntPoly& p, long long x) {
    if (p.empty()) return BigInt(0);
    BigInt result(0);
    BigInt x_pow(1);
    const BigInt x_bi(x);
    for (const BigInt& coeff : p.coefficients()) {
        result += coeff * x_pow;
        x_pow *= x_bi;
    }
    return result;
}

[[nodiscard]] static std::optional<IntPoly> try_kronecker_quadratic_factor(
    const IntPoly& f) {
    if (f.size() < 5U) return std::nullopt;

    const BigInt v0 = eval_int_poly_at_int(f, 0);
    const BigInt v1 = eval_int_poly_at_int(f, 1);
    const BigInt v2 = eval_int_poly_at_int(f, -1);

    auto divs0 = try_small_divisors(v0);
    auto divs1 = try_small_divisors(v1);
    auto divs2 = try_small_divisors(v2);
    if (!divs0 || !divs1 || !divs2) return std::nullopt;

    for (long long d0 : *divs0) {
        for (int s0 : {1, -1}) {
            const long long sd0 = s0 * d0;
            for (long long d1 : *divs1) {
                for (int s1 : {1, -1}) {
                    const long long sd1 = s1 * d1;
                    for (long long d2 : *divs2) {
                        for (int s2 : {1, -1}) {
                            const long long sd2 = s2 * d2;
                            if (((sd1 % 2) + 2) % 2 != ((sd2 % 2) + 2) % 2) continue;
                            long long a_v = (sd1 + sd2) / 2 - sd0;
                            long long b_v = (sd1 - sd2) / 2;
                            long long c_v = sd0;
                            if (a_v == 0) continue;
                            if (a_v < 0) { a_v = -a_v; b_v = -b_v; c_v = -c_v; }
                            IntPoly g(std::vector<BigInt>{BigInt(c_v), BigInt(b_v), BigInt(a_v)});
                            normalize_integer_poly(g);
                            if (g.size() != 3U) continue;

                            // Divisibility test in pure BigInt: prem(f, g) = 0
                            // iff g | f in Q[x], the same predicate the old
                            // PolyExpr division computed — but that one interned
                            // an arena IntegerLit for every coefficient of f on
                            // every one of the ~17.6k candidates, at ~8 ms each.
                            IntPoly rem = pseudo_remainder_integer_poly(f, g);
                            normalize_integer_poly(rem);
                            if (rem.is_zero()) {
                                return g;
                            }
                        }
                    }
                }
            }
        }
    }
    return std::nullopt;
}

[[nodiscard]] Result<IntPoly> divide_integer_poly_by_linear_factor(
    const IntPoly& coefficients,
    const RationalRootCandidate& root,
    symbolic::CASContext& ctx) {
    
    std::vector<ExprPtr> f_coeffs, l_coeffs;
    for(auto c : coefficients.coefficients()) f_coeffs.push_back(ctx.arena().make<IntegerLit>(c));
    l_coeffs.push_back(ctx.arena().make<IntegerLit>(-root.numerator));
    l_coeffs.push_back(ctx.arena().make<IntegerLit>(root.denominator));
    
    PolyExpr dividend(f_coeffs);
    PolyExpr divisor(l_coeffs);

    auto division = divide_poly_with_remainder(dividend, divisor, ctx);
    if (division.is_error()) {
        return fail<IntPoly>(division.error());
    }
    if (!is_zero_poly(division.value().remainder)) {
        return fail<IntPoly>(make_error(
            CASErrorKind::Unimplemented,
            "La deflazione lineare non ha prodotto resto nullo"));
    }

    return poly_to_integer_poly(division.value().quotient);
}

void append_factor_with_multiplicity(
    std::vector<PolynomialFactor>& factors,
    ExprPtr factor,
    unsigned int multiplicity) {
    if (!factors.empty() && structural_equal(factors.back().factor, factor)) {
        factors.back().multiplicity += multiplicity;
        return;
    }
    factors.push_back(PolynomialFactor{
        .factor = factor,
        .multiplicity = multiplicity,
    });
}

[[nodiscard]] bool is_unit_integer_poly(const IntPoly& coefficients) {
    return coefficients.size() == 1U;
}

[[nodiscard]] IntPoly formal_derivative_integer_poly(const IntPoly& coefficients) {
    if (coefficients.size() <= 1U) {
        return IntPoly{};
    }

    IntPoly derivative;
    derivative.resize(coefficients.size() - 1U, BigInt(0));
    for (std::size_t degree = 1U; degree < coefficients.size(); ++degree) {
        derivative[degree - 1U] = coefficients[degree] * BigInt(static_cast<long long>(degree));
    }
    normalize_integer_poly(derivative);
    return derivative;
}

[[nodiscard]] Result<IntPoly> exact_divide_integer_poly(
    const IntPoly& dividend,
    const IntPoly& divisor,
    symbolic::CASContext& ctx) {
    if (divisor.empty()) {
        return fail<IntPoly>(make_error(
            CASErrorKind::Undefined,
            "La divisione esatta polinomiale ha ricevuto un divisore nullo"));
    }

    std::vector<ExprPtr> f_coeffs, d_coeffs;
    for(auto c : dividend.coefficients()) f_coeffs.push_back(ctx.arena().make<IntegerLit>(c));
    for(auto c : divisor.coefficients()) d_coeffs.push_back(ctx.arena().make<IntegerLit>(c));

    PolyExpr dividend_expr(f_coeffs);
    PolyExpr divisor_expr(d_coeffs);
    
    auto division = divide_poly_with_remainder(dividend_expr, divisor_expr, ctx);
    if (division.is_error()) {
        return fail<IntPoly>(division.error());
    }
    if (!is_zero_poly(division.value().remainder)) {
        return fail<IntPoly>(make_error(
            CASErrorKind::Unimplemented,
            "La divisione polinomiale interna non e' risultata esatta"));
    }
    return poly_to_integer_poly(division.value().quotient);
}

[[nodiscard]] Result<std::vector<IntegerSquareFreeFactor>> square_free_factorize_integer_poly(
    const IntPoly& primitive,
    symbolic::CASContext& ctx) {
    std::vector<IntegerSquareFreeFactor> factors;
    if (primitive.size() <= 1U) {
        return ok(std::move(factors));
    }

    IntPoly derivative = formal_derivative_integer_poly(primitive);
    if (derivative.empty()) {
        factors.push_back(IntegerSquareFreeFactor{
            .factor = primitive,
            .multiplicity = 1U,
        });
        return ok(std::move(factors));
    }

    IntPoly g = gcd_integer_poly_dispatch(primitive, derivative, ctx).gcd;
    auto w_result = exact_divide_integer_poly(primitive, g, ctx);
    if (w_result.is_error()) {
        return fail<std::vector<IntegerSquareFreeFactor>>(w_result.error());
    }

    IntPoly w = primitive_integer_poly(std::move(w_result.value()));
    unsigned int multiplicity = 1U;
    while (!is_unit_integer_poly(w)) {
        IntPoly y = gcd_integer_poly_dispatch(w, g, ctx).gcd;
        auto z_result = exact_divide_integer_poly(w, y, ctx);
        if (z_result.is_error()) {
            return fail<std::vector<IntegerSquareFreeFactor>>(z_result.error());
        }

        IntPoly z = primitive_integer_poly(std::move(z_result.value()));
        if (!is_unit_integer_poly(z)) {
            factors.push_back(IntegerSquareFreeFactor{
                .factor = std::move(z),
                .multiplicity = multiplicity,
            });
        }

        auto next_g = exact_divide_integer_poly(g, y, ctx);
        if (next_g.is_error()) {
            return fail<std::vector<IntegerSquareFreeFactor>>(next_g.error());
        }

        w = primitive_integer_poly(std::move(y));
        g = primitive_integer_poly(std::move(next_g.value()));
        ++multiplicity;
    }

    return ok(std::move(factors));
}

[[nodiscard]] static std::size_t poly_hash_fnv1a(const IntPoly& f) {
    constexpr std::size_t kFnvPrime = 0x00000100000001B3ULL;
    constexpr std::size_t kFnvBasis = 0xcbf29ce484222325ULL;
    std::size_t h = kFnvBasis;
    for (std::size_t i = 0; i < f.size(); ++i) {
        const std::string& dec = f[i].decimal();
        for (unsigned char c : dec) {
            h ^= static_cast<std::size_t>(c);
            h *= kFnvPrime;
        }
        h ^= 0xABU;
        h *= kFnvPrime;
    }
    return h;
}

[[nodiscard]] static bool is_small_prime(std::uint64_t n) {
    if (n < 2U) return false;
    if (n == 2U) return true;
    if ((n % 2U) == 0U) return false;
    for (std::uint64_t d = 3U; d <= n / d; d += 2U) {
        if ((n % d) == 0U) return false;
    }
    return true;
}

[[nodiscard]] BigInt select_factorization_prime(const IntPoly& f) {
    static constexpr int kPool[] = {
        13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
        53, 59, 61, 67, 71, 73, 79, 83, 89, 97
    };
    static constexpr std::size_t kPoolSize = sizeof(kPool) / sizeof(kPool[0]);
    BigInt lc = f.leading_coeff();
    const std::size_t start = poly_hash_fnv1a(f) % kPoolSize;
    for (std::size_t i = 0; i < kPoolSize; ++i) {
        const int p_val = kPool[(start + i) % kPoolSize];
        BigInt bp(p_val);
        if (!(lc % bp).is_zero()) return bp;
    }

    std::uint64_t candidate = 101U + 2U * (poly_hash_fnv1a(f) % 128U);
    if ((candidate % 2U) == 0U) ++candidate;
    for (;;) {
        if (is_small_prime(candidate)) {
            BigInt p(static_cast<long long>(candidate));
            if (!(lc % p).is_zero()) return p;
        }
        candidate += 2U;
    }
}

[[nodiscard]] Result<void> append_integer_factor_component(
    Factorization& factorization,
    const IntPoly& component,
    unsigned int multiplicity,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    if (component.size() <= 1U) {
        return ok();
    }

    IntPoly remaining = component;
    while (remaining.size() > 1U) {
        const auto root = find_rational_root_candidate(remaining);
        if (!root.has_value()) {
            break;
        }

        auto factor_expr = integer_coefficients_to_expr(
            IntPoly(std::vector<BigInt>{-root->numerator, root->denominator}),
            var,
            ctx);
        if (factor_expr.is_error()) {
            return fail<void>(factor_expr.error());
        }
        append_factor_with_multiplicity(factorization.factors, factor_expr.value(), multiplicity);

        auto quotient = divide_integer_poly_by_linear_factor(remaining, *root, ctx);
        if (quotient.is_error()) {
            return fail<void>(quotient.error());
        }
        remaining = primitive_integer_poly(std::move(quotient.value()));
    }

    while (remaining.size() >= 5U && remaining.size() <= 7U) {
        auto quad = try_kronecker_quadratic_factor(remaining);
        if (!quad.has_value()) break;
        auto quad_expr = integer_coefficients_to_expr(*quad, var, ctx);
        if (quad_expr.is_error()) return fail<void>(quad_expr.error());
        append_factor_with_multiplicity(factorization.factors, quad_expr.value(), multiplicity);
        auto quotient = exact_divide_integer_poly(remaining, *quad, ctx);
        if (quotient.is_error()) break;
        remaining = primitive_integer_poly(std::move(quotient.value()));
    }

    if (remaining.size() > 7U) {
        auto factors_res = factorize_univariate_hensel_or_kronecker(remaining, ctx);
        if (factors_res.is_error()) {
            return fail<void>(factors_res.error());
        }
        for (const auto& h : factors_res.value()) {
            auto e = integer_coefficients_to_expr(h, var, ctx);
            if (!e.is_ok()) return fail<void>(e.error());
            append_factor_with_multiplicity(factorization.factors, e.value(), multiplicity);
        }
        remaining = IntPoly{};
    }

    if (remaining.size() > 1U) {
        auto residual = integer_coefficients_to_expr(remaining, var, ctx);
        if (residual.is_error()) {
            return fail<void>(residual.error());
        }
        append_factor_with_multiplicity(factorization.factors, residual.value(), multiplicity);
    }

    return ok();
}

Result<Factorization> factor_over_integers(ExprPtr poly, const Symbol& var, symbolic::CASContext& ctx) {
    if (!poly) {
        return fail<Factorization>(make_error(
            CASErrorKind::InvalidArgument,
            "factor_over_integers richiede un polinomio non nullo"));
    }

    auto parsed = parse_polynomial(poly, var, ctx);
    if (parsed.is_error()) {
        return fail<Factorization>(parsed.error());
    }

    auto integer_coefficients = poly_to_integer_poly(parsed.value());
    if (integer_coefficients.is_error()) {
        return fail<Factorization>(integer_coefficients.error());
    }
    if (integer_coefficients.value().empty()) {
        return fail<Factorization>(make_error(
            CASErrorKind::InvalidArgument,
            "Il polinomio nullo non ha una fattorizzazione intera canonica"));
    }

    IntPoly primitive = std::move(integer_coefficients.value());
    BigInt content = integer_content(primitive);
    divide_integer_coefficients_by_scalar(primitive, content);

    if (!primitive.empty() && primitive.leading_coeff().is_negative()) {
        content = -content;
        for (BigInt& coefficient : primitive.coefficients()) {
            coefficient = -coefficient;
        }
    }

    Factorization factorization{
        .content = ctx.arena().make<IntegerLit>(content),
        .factors = {},
    };

    if (primitive.size() <= 1U) {
        return ok(std::move(factorization));
    }

    auto square_free = square_free_factorize_integer_poly(primitive, ctx);
    if (square_free.is_error()) {
        return fail<Factorization>(square_free.error());
    }

    for (const IntegerSquareFreeFactor& component : square_free.value()) {
        auto appended = append_integer_factor_component(
            factorization,
            component.factor,
            component.multiplicity,
            var,
            ctx);
        if (appended.is_error()) {
            return fail<Factorization>(appended.error());
        }
    }

    return ok(std::move(factorization));
}

} // namespace cas::algebra
