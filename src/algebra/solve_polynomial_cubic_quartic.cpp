#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <optional>
#include <vector>

namespace cas::algebra {

[[nodiscard]] Result<std::vector<ExprPtr>> solve_degree_three_expr(
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

[[nodiscard]] Result<std::vector<ExprPtr>> solve_degree_four_expr(
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

} // namespace cas::algebra
