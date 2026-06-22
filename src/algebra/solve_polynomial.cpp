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

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_degree_one_expr(
    const PolyExpr& poly,
    symbolic::CASContext& ctx) {
    if (poly.degree() != 1U) return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument, "solve_degree_one_expr richiede un polinomio di grado 1"));
    AstArena& arena = ctx.arena();
    auto mk_int = [&](long long v) { return arena.make<IntegerLit>(BigInt(v)); };
    ExprPtr b = poly[0] ? poly[0] : mk_int(0);
    ExprPtr a = poly[1] ? poly[1] : mk_int(0);
    ExprPtr root = arena.make<Binary>(BinaryOp::Div, 
        arena.make<Binary>(BinaryOp::Sub, arena.make<IntegerLit>(BigInt(0)), b), a);
    auto r = simplify_expr(root, ctx);
    if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
    return ok(std::vector<ExprPtr>{r.value()});
}

[[nodiscard]] Result<std::vector<ExprPtr>> solve_degree_two_expr(
    const PolyExpr& poly,
    symbolic::CASContext& ctx) {
    if (poly.degree() != 2U) return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument, "solve_degree_two_expr richiede un polinomio di grado 2"));
    AstArena& arena = ctx.arena();
    auto mk_int = [&](long long v) { return arena.make<IntegerLit>(BigInt(v)); };
    ExprPtr c = poly[0] ? poly[0] : mk_int(0);
    ExprPtr b = poly[1] ? poly[1] : mk_int(0);
    ExprPtr a = poly[2] ? poly[2] : mk_int(0);
    auto mk_add = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Add, x, y); };
    auto mk_sub = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Sub, x, y); };
    auto mk_mul = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Mul, x, y); };
    auto mk_div = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Div, x, y); };
    auto mk_pow = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Pow, x, y); };
    
    ExprPtr four = mk_int(4);
    ExprPtr two = mk_int(2);
    ExprPtr b2 = mk_pow(b, two);
    ExprPtr four_ac = mk_mul(four, mk_mul(a, c));
    ExprPtr disc = mk_sub(b2, four_ac);
    ExprPtr sqrt_disc = arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{disc});
    ExprPtr neg_b = mk_sub(mk_int(0), b);
    ExprPtr two_a = mk_mul(two, a);
    
    auto s1 = simplify_expr(mk_div(mk_add(neg_b, sqrt_disc), two_a), ctx);
    auto s2 = simplify_expr(mk_div(mk_sub(neg_b, sqrt_disc), two_a), ctx);
    if (s1.is_error()) return fail<std::vector<ExprPtr>>(s1.error());
    if (s2.is_error()) return fail<std::vector<ExprPtr>>(s2.error());
    return ok(std::vector<ExprPtr>{s1.value(), s2.value()});
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
    symbolic::CASContext& ctx) {
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
        // Roots: (-b ± i*sqrt(-disc)) / (2a)
        const BigInt neg_disc = -disc;
        Rational real_rat(-b_coeff, BigInt(2) * a_coeff);
        ExprPtr real_part = make_rational_expr(ctx.arena(), real_rat);
        auto maybe_sqrt = perfect_square_root(neg_disc);
        ExprPtr I_const = ctx.arena().make<Constant>(MathConstant::I);
        ExprPtr imag_part;
        if (maybe_sqrt.has_value()) {
            Rational imag_rat(maybe_sqrt.value(), BigInt(2) * a_coeff);
            imag_part = make_rational_expr(ctx.arena(), imag_rat);
        } else {
            ExprPtr sqrt_neg_disc = ctx.arena().make<FuncCall>(
                "sqrt", std::vector<ExprPtr>{ctx.arena().make<IntegerLit>(neg_disc)});
            ExprPtr two_a = ctx.arena().make<IntegerLit>(BigInt(2) * a_coeff);
            auto r = simplify_expr(
                ctx.arena().make<Binary>(BinaryOp::Div, sqrt_neg_disc, two_a), ctx);
            if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
            imag_part = r.value();
        }
        ExprPtr i_imag = ctx.arena().make<Binary>(BinaryOp::Mul, I_const, imag_part);
        auto r1 = simplify_expr(
            ctx.arena().make<Binary>(BinaryOp::Add, real_part, i_imag), ctx);
        auto r2 = simplify_expr(
            ctx.arena().make<Binary>(BinaryOp::Sub, real_part, i_imag), ctx);
        if (r1.is_error()) return fail<std::vector<ExprPtr>>(r1.error());
        if (r2.is_error()) return fail<std::vector<ExprPtr>>(r2.error());
        return ok(std::vector<ExprPtr>{r1.value(), r2.value()});
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

        switch (deg) {
        case 1U: return solve_degree_one_expr(parsed.value(), ctx);
        case 2U: return solve_degree_two_expr(parsed.value(), ctx);
        case 3U: return solve_degree_three_expr(parsed.value(), ctx);
        case 4U: return solve_degree_four_expr(parsed.value(), ctx);
        default: break;
        }

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

    if (deg == 1U) {
        auto r = solve_degree_one_expr(factor_poly.value(), ctx);
        if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
        for (unsigned int m = 0; m < multiplicity; ++m) {
            roots.insert(roots.end(), r.value().begin(), r.value().end());
        }
        return ok(std::move(roots));
    }

    if (int_coeffs.is_ok() && deg == 2U) {
        auto r = solve_degree_two(int_coeffs.value(), ctx);
        if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
        for (unsigned int m = 0; m < multiplicity; ++m) {
            roots.insert(roots.end(), r.value().begin(), r.value().end());
        }
        return ok(std::move(roots));
    }

    if (deg == 2U) {
        auto r = solve_degree_two_expr(factor_poly.value(), ctx);
        if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
        for (unsigned int m = 0; m < multiplicity; ++m) {
            roots.insert(roots.end(), r.value().begin(), r.value().end());
        }
        return ok(std::move(roots));
    }

    if (deg == 3U) {
        auto r = solve_degree_three_expr(factor_poly.value(), ctx);
        if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
        for (unsigned int m = 0; m < multiplicity; ++m) {
            roots.insert(roots.end(), r.value().begin(), r.value().end());
        }
        return ok(std::move(roots));
    }

    /*
    if (deg == 4U) {
        auto r = solve_degree_four_expr(factor_poly.value(), ctx);
        if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
        for (unsigned int m = 0; m < multiplicity; ++m) {
            roots.insert(roots.end(), r.value().begin(), r.value().end());
        }
        return ok(std::move(roots));
    }
    */

    if (int_coeffs.is_ok()) {
        if (auto n = is_cyclotomic(int_coeffs.value(), ctx.max_cyclotomic_n())) {
            auto r = cyclotomic_roots(*n, var, ctx.arena());
            for (unsigned int m = 0; m < multiplicity; ++m) {
                roots.insert(roots.end(), r.begin(), r.end());
            }
            return ok(std::move(roots));
        }
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
            return solve_degree_two(int_coeffs.value(), ctx);
        default:
            break;
        }
    } else {
        // Exact solvers for non-integer coefficients (e.g., from back-substitution)
        switch (deg) {
        case 1U: return solve_degree_one_expr(parsed.value(), ctx);
        case 2U: return solve_degree_two_expr(parsed.value(), ctx);
        case 3U: return solve_degree_three_expr(parsed.value(), ctx);
        case 4U: return solve_degree_four_expr(parsed.value(), ctx);
        default: break;
        }
    }

    // Attempt factorization for higher degrees
    if (deg >= 3U) {
        auto factored = solve_by_factoring(poly, var, ctx);
        if (factored.is_ok() && !factored.value().empty()) {
            return factored;
        }
    }

    std::vector<ExprPtr> roots;
    roots.reserve(deg);
    for (std::size_t k = 0; k < deg; ++k) {
        roots.push_back(ctx.arena().make<RootOf>(poly, var, std::optional<std::size_t>{k}));
    }
    return ok(std::move(roots));
}

} // namespace cas::algebra
