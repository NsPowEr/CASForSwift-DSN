#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <optional>
#include <vector>

namespace cas::algebra {

[[nodiscard]] static BigInt bigint_isqrt(const BigInt& n) {
    if (n.is_zero() || n.is_negative()) return BigInt(0);
    const std::size_t bits = n.bit_length();
    BigInt x = n.shift_right_bits(bits / 2U);
    if (x.is_zero()) x = BigInt(1);
    while (true) {
        BigInt x1 = (x + n / x) / BigInt(2);
        if (x1 >= x) return x;
        x = std::move(x1);
    }
}

[[nodiscard]] static std::optional<BigInt> perfect_square_root(const BigInt& n) {
    if (n.is_negative()) return std::nullopt;
    if (n.is_zero()) return BigInt(0);
    BigInt r = bigint_isqrt(n);
    if (r * r == n) return r;
    return std::nullopt;
}

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_degree_one(
    const IntPoly& coeffs,
    symbolic::CASContext& ctx) {
    const BigInt& b = coeffs[0];
    const BigInt& a = coeffs[1];
    Rational root(-b, a);
    return ok(std::vector<ExprPtr>{make_rational_expr(ctx.arena(), root)});
}

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_degree_two(
    const IntPoly& coeffs,
    const Symbol& var,
    symbolic::CASContext& ctx,
    ExprPtr original_poly) {
    const BigInt& c_coeff = coeffs[0];
    const BigInt& b_coeff = coeffs[1];
    const BigInt& a_coeff = coeffs[2];

    const BigInt disc = b_coeff * b_coeff - BigInt(4) * a_coeff * c_coeff;

    if (disc.is_zero()) {
        Rational root(-b_coeff, BigInt(2) * a_coeff);
        ExprPtr r = make_rational_expr(ctx.arena(), root);
        return ok(std::vector<ExprPtr>{r, r});
    }

    if (disc.is_negative()) {
        return ok(std::vector<ExprPtr>{
            ctx.arena().make<RootOf>(original_poly, var, std::optional<std::size_t>{0U}),
            ctx.arena().make<RootOf>(original_poly, var, std::optional<std::size_t>{1U}),
        });
    }

    auto maybe_r = perfect_square_root(disc);
    if (maybe_r.has_value()) {
        Rational root1(-b_coeff + maybe_r.value(), BigInt(2) * a_coeff);
        Rational root2(-b_coeff - maybe_r.value(), BigInt(2) * a_coeff);
        return ok(std::vector<ExprPtr>{
            make_rational_expr(ctx.arena(), root1),
            make_rational_expr(ctx.arena(), root2),
        });
    }

    ExprPtr neg_b = ctx.arena().make<IntegerLit>(-b_coeff);
    ExprPtr disc_expr = ctx.arena().make<IntegerLit>(disc);
    ExprPtr sqrt_disc = ctx.arena().make<FuncCall>(
        "sqrt", std::vector<ExprPtr>{disc_expr});
    ExprPtr two_a = ctx.arena().make<IntegerLit>(BigInt(2) * a_coeff);

    auto r1 = simplify_expr(
        ctx.arena().make<Binary>(BinaryOp::Div,
            ctx.arena().make<Sum>(std::vector<ExprPtr>{neg_b, sqrt_disc}),
            two_a),
        ctx);
    if (r1.is_error()) return fail<std::vector<ExprPtr>>(r1.error());

    auto r2 = simplify_expr(
        ctx.arena().make<Binary>(BinaryOp::Div,
            ctx.arena().make<Binary>(BinaryOp::Sub, neg_b, sqrt_disc),
            two_a),
        ctx);
    if (r2.is_error()) return fail<std::vector<ExprPtr>>(r2.error());

    return ok(std::vector<ExprPtr>{r1.value(), r2.value()});
}

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_factor(
    ExprPtr factor_expr,
    const Symbol& var,
    unsigned int multiplicity,
    symbolic::CASContext& ctx);

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_by_factoring(
    ExprPtr poly,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    auto factorization = factor_over_integers(poly, var, ctx);
    if (factorization.is_error()) {
        auto parsed = parse_polynomial(poly, var, ctx);
        if (parsed.is_error()) return fail<std::vector<ExprPtr>>(parsed.error());
        const std::size_t deg = poly_degree(parsed.value());
        std::vector<ExprPtr> roots;
        roots.reserve(deg);
        for (std::size_t k = 0; k < deg; ++k) {
            roots.push_back(ctx.arena().make<RootOf>(poly, var, std::optional<std::size_t>{k}));
        }
        return ok(std::move(roots));
    }

    std::vector<ExprPtr> all_roots;
    for (const PolynomialFactor& pf : factorization.value().factors) {
        auto factor_roots = solve_factor(pf.factor, var, pf.multiplicity, ctx);
        if (factor_roots.is_error()) return fail<std::vector<ExprPtr>>(factor_roots.error());
        all_roots.insert(all_roots.end(),
            factor_roots.value().begin(), factor_roots.value().end());
    }
    return ok(std::move(all_roots));
}

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_factor(
    ExprPtr factor_expr,
    const Symbol& var,
    unsigned int multiplicity,
    symbolic::CASContext& ctx) {
    auto factor_poly = parse_polynomial(factor_expr, var, ctx);
    if (factor_poly.is_error()) return fail<std::vector<ExprPtr>>(factor_poly.error());

    const std::size_t deg = poly_degree(factor_poly.value());
    std::vector<ExprPtr> roots;

    if (deg == 0U) return ok(std::move(roots));

    auto int_coeffs = poly_to_integer_poly(factor_poly.value());

    if (int_coeffs.is_ok() && deg == 1U) {
        auto r = solve_degree_one(int_coeffs.value(), ctx);
        if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
        for (unsigned int m = 0; m < multiplicity; ++m) {
            roots.insert(roots.end(), r.value().begin(), r.value().end());
        }
        return ok(std::move(roots));
    }

    if (int_coeffs.is_ok() && deg == 2U) {
        auto r = solve_degree_two(int_coeffs.value(), var, ctx, factor_expr);
        if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
        for (unsigned int m = 0; m < multiplicity; ++m) {
            roots.insert(roots.end(), r.value().begin(), r.value().end());
        }
        return ok(std::move(roots));
    }

    for (std::size_t k = 0; k < deg; ++k) {
        for (unsigned int m = 0; m < multiplicity; ++m) {
            roots.push_back(ctx.arena().make<RootOf>(
                factor_expr, var, std::optional<std::size_t>{k}));
        }
    }
    return ok(std::move(roots));
}

Result<std::vector<ExprPtr>> solve_polynomial(
    ExprPtr poly,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    if (!poly) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::InvalidArgument,
            "solve_polynomial richiede un polinomio non nullo"));
    }

    auto parsed = parse_polynomial(poly, var, ctx);
    if (parsed.is_error()) return fail<std::vector<ExprPtr>>(parsed.error());

    if (is_zero_poly(parsed.value())) {
        return fail<std::vector<ExprPtr>>(make_error(
            CASErrorKind::InvalidArgument,
            "solve_polynomial ha ricevuto il polinomio zero (soluzioni infinite)"));
    }

    const std::size_t deg = poly_degree(parsed.value());

    if (deg == 0U) {
        return ok(std::vector<ExprPtr>{});
    }

    auto int_coeffs = poly_to_integer_poly(parsed.value());

    if (int_coeffs.is_ok()) {
        switch (deg) {
        case 1U:
            return solve_degree_one(int_coeffs.value(), ctx);
        case 2U:
            return solve_degree_two(int_coeffs.value(), var, ctx, poly);
        default:
            return solve_by_factoring(poly, var, ctx);
        }
    }

    if (deg <= 4U) {
        return solve_by_factoring(poly, var, ctx);
    }

    std::vector<ExprPtr> roots;
    roots.reserve(deg);
    for (std::size_t k = 0; k < deg; ++k) {
        roots.push_back(ctx.arena().make<RootOf>(poly, var, std::optional<std::size_t>{k}));
    }
    return ok(std::move(roots));
}

} // namespace cas::algebra
