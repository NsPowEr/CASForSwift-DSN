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

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_degree_three_expr(
    const PolyExpr& poly,
    symbolic::CASContext& ctx) {
    if (poly.degree() != 3U) return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument, "solve_degree_three_expr richiede un polinomio di grado 3"));

    AstArena& arena = ctx.arena();
    auto mk_int = [&](long long v) { return arena.make<IntegerLit>(BigInt(v)); };

    ExprPtr d = poly[0] ? poly[0] : mk_int(0);
    ExprPtr c = poly[1] ? poly[1] : mk_int(0);
    ExprPtr b = poly[2] ? poly[2] : mk_int(0);
    ExprPtr a = poly[3] ? poly[3] : mk_int(0);
    auto mk_add = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Add, x, y); };
    auto mk_sub = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Sub, x, y); };
    auto mk_mul = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Mul, x, y); };
    auto mk_div = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Div, x, y); };
    auto mk_pow = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Pow, x, y); };
    
    ExprPtr zero = mk_int(0);
    ExprPtr one = mk_int(1);
    ExprPtr two = mk_int(2);
    ExprPtr three = mk_int(3);
    ExprPtr four = mk_int(4);
    ExprPtr twenty_seven = mk_int(27);

    // p = b/a, q = c/a, r = d/a
    ExprPtr p = mk_div(b, a);
    ExprPtr q = mk_div(c, a);
    ExprPtr r = mk_div(d, a);

    // P = q - p^2 / 3
    ExprPtr p2 = mk_pow(p, two);
    ExprPtr p2_3 = mk_div(p2, three);
    ExprPtr bigP = mk_sub(q, p2_3);

    // Q = 2p^3 / 27 - pq / 3 + r
    ExprPtr p3 = mk_pow(p, three);
    ExprPtr two_p3 = mk_mul(two, p3);
    ExprPtr two_p3_27 = mk_div(two_p3, twenty_seven);
    ExprPtr pq = mk_mul(p, q);
    ExprPtr pq_3 = mk_div(pq, three);
    ExprPtr bigQ_part1 = mk_sub(two_p3_27, pq_3);
    ExprPtr bigQ = mk_add(bigQ_part1, r);

    // Delta = Q^2 / 4 + P^3 / 27
    ExprPtr bigQ2 = mk_pow(bigQ, two);
    ExprPtr bigQ2_4 = mk_div(bigQ2, four);
    ExprPtr bigP3 = mk_pow(bigP, three);
    ExprPtr bigP3_27 = mk_div(bigP3, twenty_seven);
    ExprPtr delta = mk_add(bigQ2_4, bigP3_27);

    // Delta check for multiple roots
    auto delta_simp = simplify_expr(delta, ctx);
    if (delta_simp.is_ok() && is_zero_expr(delta_simp.value())) {
        auto p_simp = simplify_expr(bigP, ctx);
        auto q_simp = simplify_expr(bigQ, ctx);
        if (p_simp.is_ok() && is_zero_expr(p_simp.value()) &&
            q_simp.is_ok() && is_zero_expr(q_simp.value())) {
            // Triple root case: x^3 = 0 -> x = -p/3
            auto root_res = simplify_expr(mk_sub(zero, mk_div(p, three)), ctx);
            if (root_res.is_error()) return fail<std::vector<ExprPtr>>(root_res.error());
            ExprPtr root = root_res.value();
            return ok(std::vector<ExprPtr>{root, root, root});
        }
        // Double root case: Delta = 0, P, Q != 0
        // roots are 2u - p/3, -u - p/3, -u - p/3 where u = cbrt(-Q/2)
        ExprPtr neg_Q_2 = mk_div(mk_sub(zero, bigQ), two);
        ExprPtr u = mk_pow(neg_Q_2, mk_div(one, three));
        ExprPtr p_3 = mk_div(p, three);
        auto s1 = simplify_expr(mk_sub(mk_mul(two, u), p_3), ctx);
        auto s2 = simplify_expr(mk_sub(mk_sub(zero, u), p_3), ctx);
        if (s1.is_error()) return fail<std::vector<ExprPtr>>(s1.error());
        if (s2.is_error()) return fail<std::vector<ExprPtr>>(s2.error());
        return ok(std::vector<ExprPtr>{s1.value(), s2.value(), s2.value()});
    }

    // u = cbrt(-Q/2 + sqrt(Delta)), v = cbrt(-Q/2 - sqrt(Delta))
    ExprPtr neg_Q = mk_sub(zero, bigQ);
    ExprPtr neg_Q_2 = mk_div(neg_Q, two);
    ExprPtr sqrt_delta = arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{delta});
    
    ExprPtr inner_u = mk_add(neg_Q_2, sqrt_delta);
    ExprPtr inner_v = mk_sub(neg_Q_2, sqrt_delta);
    
    ExprPtr one_third = mk_div(one, three);
    ExprPtr u = mk_pow(inner_u, one_third);
    ExprPtr v = mk_pow(inner_v, one_third);

    // omega = (-1 + I*sqrt(3)) / 2
    ExprPtr i_const = arena.make<Constant>(MathConstant::I);
    ExprPtr sqrt3 = arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{three});
    ExprPtr i_sqrt3 = mk_mul(i_const, sqrt3);
    ExprPtr neg1 = mk_sub(zero, one);
    ExprPtr w1_num = mk_add(neg1, i_sqrt3);
    ExprPtr w1 = mk_div(w1_num, two);
    ExprPtr w2_num = mk_sub(neg1, i_sqrt3);
    ExprPtr w2 = mk_div(w2_num, two);

    // t0 = u + v
    // t1 = w1 u + w2 v
    // t2 = w2 u + w1 v
    ExprPtr t0 = mk_add(u, v);
    ExprPtr t1 = mk_add(mk_mul(w1, u), mk_mul(w2, v));
    ExprPtr t2 = mk_add(mk_mul(w2, u), mk_mul(w1, v));

    ExprPtr p_3 = mk_div(p, three);
    
    auto s0 = simplify_expr(mk_sub(t0, p_3), ctx);
    auto s1 = simplify_expr(mk_sub(t1, p_3), ctx);
    auto s2 = simplify_expr(mk_sub(t2, p_3), ctx);

    if (s0.is_error()) return fail<std::vector<ExprPtr>>(s0.error());
    if (s1.is_error()) return fail<std::vector<ExprPtr>>(s1.error());
    if (s2.is_error()) return fail<std::vector<ExprPtr>>(s2.error());

    return ok(std::vector<ExprPtr>{s0.value(), s1.value(), s2.value()});
}

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_degree_two_expr(
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

[[nodiscard]] static Result<std::vector<ExprPtr>> solve_degree_four_expr(
    const PolyExpr& poly,
    symbolic::CASContext& ctx) {
    if (poly.degree() != 4U) return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument, "solve_degree_four_expr richiede un polinomio di grado 4"));
    AstArena& arena = ctx.arena();
    auto mk_int = [&](long long v) { return arena.make<IntegerLit>(BigInt(v)); };
    ExprPtr e = poly[0] ? poly[0] : mk_int(0);
    ExprPtr d = poly[1] ? poly[1] : mk_int(0);
    ExprPtr c = poly[2] ? poly[2] : mk_int(0);
    ExprPtr b = poly[3] ? poly[3] : mk_int(0);
    ExprPtr a = poly[4] ? poly[4] : mk_int(0);
    auto mk_add = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Add, x, y); };
    auto mk_sub = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Sub, x, y); };
    auto mk_mul = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Mul, x, y); };
    auto mk_div = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Div, x, y); };
    auto mk_pow = [&](ExprPtr x, ExprPtr y) { return arena.make<Binary>(BinaryOp::Pow, x, y); };

    ExprPtr zero = mk_int(0);
    ExprPtr two = mk_int(2);
    ExprPtr three = mk_int(3);
    ExprPtr four = mk_int(4);
    ExprPtr eight = mk_int(8);
    ExprPtr sixteen = mk_int(16);
    ExprPtr two_fifty_six = mk_int(256);

    // p = b/a, q = c/a, r = d/a, s = e/a
    ExprPtr p = mk_div(b, a);
    ExprPtr q = mk_div(c, a);
    ExprPtr r = mk_div(d, a);
    ExprPtr s = mk_div(e, a);

    // P = q - 3p^2/8
    ExprPtr bigP = mk_sub(q, mk_div(mk_mul(three, mk_pow(p, two)), eight));
    // Q = r - pq/2 + p^3/8
    ExprPtr bigQ = mk_add(mk_sub(r, mk_div(mk_mul(p, q), two)), mk_div(mk_pow(p, three), eight));
    // R = s - pr/4 + p^2q/16 - 3p^4/256
    ExprPtr bigR = mk_sub(mk_add(mk_sub(s, mk_div(mk_mul(p, r), four)), 
                         mk_div(mk_mul(mk_pow(p, two), q), sixteen)),
                         mk_div(mk_mul(three, mk_pow(p, four)), two_fifty_six));

    // Resolvent cubic: y^3 + 2Py^2 + (P^2 - 4R)y - Q^2 = 0
    PolyExpr resolvent;
    resolvent.push_back(mk_sub(zero, mk_pow(bigQ, two)));
    resolvent.push_back(mk_sub(mk_pow(bigP, two), mk_mul(four, bigR)));
    resolvent.push_back(mk_mul(two, bigP));
    resolvent.push_back(mk_int(1));
    
    auto cubic_roots = solve_degree_three_expr(resolvent, ctx);
    if (cubic_roots.is_error()) return cubic_roots;
    
    // Pick first root y0
    ExprPtr y0 = cubic_roots.value()[0];
    
    // If Q is zero, it's biquadratic. Check symbolically:
    auto q_simplified = simplify_expr(bigQ, ctx);
    if (q_simplified.is_ok() && is_zero_expr(q_simplified.value())) {
        // u^4 + Pu^2 + R = 0
        PolyExpr biquad;
        biquad.push_back(bigR); biquad.push_back(bigP); biquad.push_back(mk_int(1));
        auto z_roots = solve_degree_two_expr(biquad, ctx);
        if (z_roots.is_error()) return z_roots;
        std::vector<ExprPtr> roots;
        ExprPtr p_4 = mk_div(p, four);
        for (ExprPtr z : z_roots.value()) {
            ExprPtr sqrt_z = arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{z});
            roots.push_back(simplify_expr(mk_sub(sqrt_z, p_4), ctx).value());
            roots.push_back(simplify_expr(mk_sub(mk_sub(zero, sqrt_z), p_4), ctx).value());
        }
        return ok(roots);
    }

    // General Ferrari
    ExprPtr A = arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{y0});
    // B = 1/2 * (P + y0 - Q/A)
    ExprPtr B = mk_div(mk_sub(mk_add(bigP, y0), mk_div(bigQ, A)), two);
    // C = 1/2 * (P + y0 + Q/A)
    ExprPtr C = mk_div(mk_add(mk_add(bigP, y0), mk_div(bigQ, A)), two);

    // Solve u^2 + Au + B = 0 and u^2 - Au + C = 0
    PolyExpr q1; q1.push_back(B); q1.push_back(A); q1.push_back(mk_int(1));
    PolyExpr q2; q2.push_back(C); q2.push_back(mk_sub(zero, A)); q2.push_back(mk_int(1));
    
    auto r1 = solve_degree_two_expr(q1, ctx);
    auto r2 = solve_degree_two_expr(q2, ctx);
    if (r1.is_error()) return r1;
    if (r2.is_error()) return r2;
    
    std::vector<ExprPtr> final_roots;
    ExprPtr p_4 = mk_div(p, four);
    for (ExprPtr u : r1.value()) final_roots.push_back(simplify_expr(mk_sub(u, p_4), ctx).value());
    for (ExprPtr u : r2.value()) final_roots.push_back(simplify_expr(mk_sub(u, p_4), ctx).value());
    
    return ok(final_roots);
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
