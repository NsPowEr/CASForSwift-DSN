#include "residue_theorem_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/algebraic_number.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/calculus.hpp"
#include <string>
#include <vector>

namespace cas::calculus {

namespace {
// Build the RootOf node for an irreducible quadratic m(y) = y^2 + b*y + c,
// using the supplied generator variable name `gen_name`.  The minimal
// polynomial is the same m, expressed in `gen_var`.
[[nodiscard]] ExprPtr build_min_poly_expr(
    const Rational& b, const Rational& c, const Symbol& gen_var, AstArena& arena) {
    ExprPtr y = arena.make<Symbol>(gen_var);
    ExprPtr y2 = arena.make<Binary>(BinaryOp::Pow, y, make_int(arena, 2));
    ExprPtr by = arena.make<Binary>(BinaryOp::Mul, make_rational_expr(arena, b), y);
    ExprPtr c_e = make_rational_expr(arena, c);
    ExprPtr sum1 = arena.make<Binary>(BinaryOp::Add, y2, by);
    return arena.make<Binary>(BinaryOp::Add, sum1, c_e);
}

// Build  y⁴ + b·y² + c  as ExprPtr.
[[nodiscard]] ExprPtr build_biquadratic_min_poly_expr(
    const Rational& b, const Rational& c, const Symbol& gen_var, AstArena& arena) {
    ExprPtr y = arena.make<Symbol>(gen_var);
    ExprPtr y2 = arena.make<Binary>(BinaryOp::Pow, y, make_int(arena, 2));
    ExprPtr y4 = arena.make<Binary>(BinaryOp::Pow, y, make_int(arena, 4));
    ExprPtr by2 = arena.make<Binary>(BinaryOp::Mul, make_rational_expr(arena, b), y2);
    ExprPtr c_e = make_rational_expr(arena, c);
    ExprPtr sum1 = arena.make<Binary>(BinaryOp::Add, y4, by2);
    return arena.make<Binary>(BinaryOp::Add, sum1, c_e);
}
} // namespace

Result<ExprPtr> contribution_from_quadratic(
    const Rational& b,
    const Rational& c,
    const Rational& discriminant,
    ExprPtr N_over_D,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    Symbol root_var = var;
    ExprPtr min_poly_in_x_raw = build_min_poly_expr(b, c, root_var, arena);
    auto min_poly_simp = simplify_or_fail(min_poly_in_x_raw, ctx);
    if (min_poly_simp.is_error()) return fail<ExprPtr>(min_poly_simp.error());
    ExprPtr alpha_expr = arena.make<RootOf>(min_poly_simp.value(), root_var, std::nullopt);

    auto res = residue(N_over_D, var, alpha_expr, ctx);
    if (res.is_error()) return fail<ExprPtr>(res.error());

    auto reduced = algebra::simplify_in_q_alpha(res.value(), ctx);
    if (reduced.is_error()) return fail<ExprPtr>(reduced.error());

    algebra::AlgebraicNumber::CoeffVec min_poly_coeffs;
    min_poly_coeffs.push_back(c);
    min_poly_coeffs.push_back(b);
    min_poly_coeffs.push_back(Rational(BigInt(1)));

    auto expressed = algebra::try_express_in_q_alpha(
        reduced.value(), alpha_expr, min_poly_coeffs, ctx);
    if (expressed.is_error()) return fail<ExprPtr>(expressed.error());
    if (!expressed.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Residue theorem: residue not expressible in Q(α) for this quadratic factor"});
    }

    const algebra::AlgebraicNumber& an = *expressed.value();
    const auto& value = an.value();
    Rational e_part(BigInt(0));
    Rational f_part(BigInt(0));
    if (value.size() >= 1U) e_part = value[0];
    if (value.size() >= 2U) f_part = value[1];

    Rational neg_disc = -discriminant;
    ExprPtr neg_disc_expr = make_rational_expr(arena, neg_disc);
    ExprPtr sqrt_neg_disc = arena.make<FuncCall>(
        std::string("sqrt"),
        std::vector<ExprPtr>{neg_disc_expr});

    ExprPtr f_expr = make_rational_expr(arena, f_part);
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr neg_pi = arena.make<Binary>(BinaryOp::Mul, make_int(arena, -1), pi);
    ExprPtr term = arena.make<Binary>(BinaryOp::Mul, neg_pi, f_expr);
    ExprPtr full = arena.make<Binary>(BinaryOp::Mul, term, sqrt_neg_disc);

    (void)e_part;

    return simplify_or_fail(full, ctx);
}

Result<ExprPtr> contribution_from_irreducible_biquadratic(
    const Rational& b,
    const Rational& c,
    ExprPtr N_over_D,
    const Symbol& var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    Symbol gen = var;
    ExprPtr min_poly_raw = build_biquadratic_min_poly_expr(b, c, gen, arena);
    auto min_poly_simp = simplify_or_fail(min_poly_raw, ctx);
    if (min_poly_simp.is_error()) return fail<ExprPtr>(min_poly_simp.error());
    ExprPtr alpha_expr = arena.make<RootOf>(min_poly_simp.value(), gen, std::nullopt);

    auto res = residue(N_over_D, var, alpha_expr, ctx);
    if (res.is_error()) return fail<ExprPtr>(res.error());

    auto reduced = algebra::simplify_in_q_alpha(res.value(), ctx);
    if (reduced.is_error()) return fail<ExprPtr>(reduced.error());

    algebra::AlgebraicNumber::CoeffVec min_poly_coeffs;
    min_poly_coeffs.push_back(c);
    min_poly_coeffs.push_back(Rational(BigInt(0)));
    min_poly_coeffs.push_back(b);
    min_poly_coeffs.push_back(Rational(BigInt(0)));
    min_poly_coeffs.push_back(Rational(BigInt(1)));

    auto expressed = algebra::try_express_in_q_alpha(
        reduced.value(), alpha_expr, min_poly_coeffs, ctx);
    if (expressed.is_error()) return fail<ExprPtr>(expressed.error());
    if (!expressed.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Residue theorem: residue not expressible in Q(α) for this biquadratic factor"});
    }

    const algebra::AlgebraicNumber& an = *expressed.value();
    const auto& value = an.value();
    Rational c0(BigInt(0)), c1(BigInt(0)), c2(BigInt(0)), c3(BigInt(0));
    if (value.size() >= 1U) c0 = value[0];
    if (value.size() >= 2U) c1 = value[1];
    if (value.size() >= 3U) c2 = value[2];
    if (value.size() >= 4U) c3 = value[3];

    ExprPtr c_expr = make_rational_expr(arena, c);
    ExprPtr sqrt_c = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{c_expr});

    ExprPtr two_sqrt_c = arena.make<Binary>(BinaryOp::Mul, make_int(arena, 2), sqrt_c);
    ExprPtr b_expr = make_rational_expr(arena, b);
    ExprPtr radicand = arena.make<Binary>(BinaryOp::Add, two_sqrt_c, b_expr);
    ExprPtr sqrt_radicand = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{radicand});

    ExprPtr sqrt_c_minus_b = arena.make<Binary>(BinaryOp::Sub, sqrt_c, b_expr);
    ExprPtr c3_expr = make_rational_expr(arena, c3);
    ExprPtr c3_times = arena.make<Binary>(BinaryOp::Mul, c3_expr, sqrt_c_minus_b);
    ExprPtr c1_expr = make_rational_expr(arena, c1);
    ExprPtr inner_sum = arena.make<Binary>(BinaryOp::Add, c1_expr, c3_times);

    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr neg_two_pi = arena.make<Binary>(BinaryOp::Mul, make_int(arena, -2), pi);
    ExprPtr part1 = arena.make<Binary>(BinaryOp::Mul, neg_two_pi, sqrt_radicand);
    ExprPtr full = arena.make<Binary>(BinaryOp::Mul, part1, inner_sum);

    (void)c0;
    (void)c2;
    return simplify_or_fail(full, ctx);
}

}  // namespace cas::calculus
