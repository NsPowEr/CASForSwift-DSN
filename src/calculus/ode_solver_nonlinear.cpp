// ode_solver_nonlinear.cpp — F5.3 nonlinear 1st-order ODE solvers.
//
// Currently handles two Riccati sub-families with closed-form solutions:
//
//   1.  q_0(x) ≡ 0                          (homogeneous Bernoulli reduction)
//          y' = q_1(x)·y + q_2(x)·y²
//       transforms via v = 1/y into the linear ODE
//          v' + q_1(x)·v = -q_2(x)
//       which is solved by `solve_ode_1st_order(Linear1stOrder)`.
//
//   2.  q_0, q_1, q_2 all constant in x   (autonomous Riccati)
//          y' = a + b·y + c·y²,    a, b, c ∈ Q.
//       Solved in closed form via the partial-fraction
//          ∫ dy / (a + b·y + c·y²) = x − C.
//       Two sub-cases on the discriminant Δ = b² − 4ac:
//          Δ > 0  →  y = (1/(2c))·( r_1 + r_2·tanh-style expansion ) using
//                    the two real roots r_1, r_2 of c·y² + b·y + a = 0.
//          Δ = 0  →  y = -1/(c·(x − C)) − b/(2c).
//          Δ < 0  →  y = (1/(2c))·( -b + √(-Δ)·tan( √(-Δ)·(x − C)/2 ) ).
//
// Other Riccati instances (variable coefficients with no known particular
// solution, Airy/Bessel reductions, etc.) return Unimplemented with an
// explicit diagnostic; closing those requires either a particular-solution
// oracle (Risch-DE, Bronstein ch.8) or the full 2nd-order variable-coefficient
// linear ODE solver (currently absent — tracked in F5.3 / B2 continuation).
//
// Clairaut and d'Alembert solvers will be appended here in a follow-up commit.

#include "cas/ode.hpp"
#include "cas/bigint.hpp"
#include "cas/rational.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] bool is_zero_expr(ExprPtr expr, symbolic::CASContext& ctx) {
    auto res = ctx.simplify(expr);
    if (res.is_error()) return false;
    if (const auto* il = expr_cast<IntegerLit>(res.value())) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(res.value())) return rl->numerator.is_zero();
    return false;
}

// Decide whether `e` is a pure rational/integer literal (constant in x).
[[nodiscard]] std::optional<Rational> as_rational(ExprPtr e, symbolic::CASContext& ctx) {
    auto s = ctx.simplify(e);
    if (s.is_error()) return std::nullopt;
    ExprPtr v = s.value();
    if (const auto* il = expr_cast<IntegerLit>(v)) {
        return Rational(il->value);
    }
    if (const auto* rl = expr_cast<RationalLit>(v)) {
        return Rational(rl->numerator, rl->denominator);
    }
    if (const auto* un = expr_cast<Unary>(v); un && un->op == UnaryOp::Neg) {
        auto inner = as_rational(un->operand, ctx);
        if (inner) return -*inner;
    }
    return std::nullopt;
}

[[nodiscard]] ExprPtr make_rational(AstArena& arena, const Rational& r) {
    if (r.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(r.numerator());
    }
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

// Build  y = -1/(c·(x − C)) − b/(2c)  for the Δ = 0 case.
[[nodiscard]] ExprPtr build_double_root_solution(
    const Rational& b, const Rational& c, const Symbol& x,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    Symbol C_sym = ctx.make_fresh_symbol("C");
    ExprPtr C = arena.make<Symbol>(C_sym.name);
    ExprPtr xs = arena.make<Symbol>(x.name);
    ExprPtr x_minus_C = arena.make<Binary>(BinaryOp::Sub, xs, C);
    ExprPtr c_expr = make_rational(arena, c);
    ExprPtr denom = arena.make<Binary>(BinaryOp::Mul, c_expr, x_minus_C);
    ExprPtr first = arena.make<Unary>(UnaryOp::Neg,
        arena.make<Binary>(BinaryOp::Div, arena.make<IntegerLit>(BigInt(1)), denom));
    Rational shift = -b / (Rational(BigInt(2)) * c);
    ExprPtr second = make_rational(arena, shift);
    return arena.make<Binary>(BinaryOp::Add, first, second);
}

// Build  y = (1/(2c))·( -b + √(-Δ)·tan( √(-Δ)·(x − C)/2 ) )  for Δ < 0.
[[nodiscard]] ExprPtr build_negative_disc_solution(
    const Rational& b, const Rational& c, const Rational& neg_disc,
    const Symbol& x, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    Symbol C_sym = ctx.make_fresh_symbol("C");
    ExprPtr C = arena.make<Symbol>(C_sym.name);
    ExprPtr xs = arena.make<Symbol>(x.name);
    ExprPtr sqrt_nd = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_rational(arena, neg_disc)});
    ExprPtr half = arena.make<RationalLit>(BigInt(1), BigInt(2));
    ExprPtr x_minus_C = arena.make<Binary>(BinaryOp::Sub, xs, C);
    ExprPtr tan_arg = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Mul, sqrt_nd, half), x_minus_C);
    ExprPtr tan_expr = arena.make<FuncCall>(BuiltinOp::Tan, std::vector<ExprPtr>{tan_arg});
    ExprPtr inner = arena.make<Binary>(BinaryOp::Add,
        arena.make<Unary>(UnaryOp::Neg, make_rational(arena, b)),
        arena.make<Binary>(BinaryOp::Mul, sqrt_nd, tan_expr));
    Rational two_c = Rational(BigInt(2)) * c;
    ExprPtr scale = make_rational(arena, Rational(BigInt(1)) / two_c);
    return arena.make<Binary>(BinaryOp::Mul, scale, inner);
}

// Build  y = (1/(2c))·( -b + √Δ · tanh( √Δ·(x − C)/2 ) ) for Δ > 0.
// Derivation: ∫ dy/(c·y² + b·y + a) = ∫ dx.  For Δ > 0 the antiderivative
// is (2/√Δ)·atanh((2c·y + b)/√Δ); inverting and using atanh ↔ tanh closes
// the form.
[[nodiscard]] ExprPtr build_positive_disc_solution(
    const Rational& b, const Rational& c, const Rational& disc,
    const Symbol& x, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    Symbol C_sym = ctx.make_fresh_symbol("C");
    ExprPtr C = arena.make<Symbol>(C_sym.name);
    ExprPtr xs = arena.make<Symbol>(x.name);
    ExprPtr sqrt_d = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_rational(arena, disc)});
    ExprPtr half = arena.make<RationalLit>(BigInt(1), BigInt(2));
    ExprPtr x_minus_C = arena.make<Binary>(BinaryOp::Sub, xs, C);
    ExprPtr tanh_arg = arena.make<Binary>(BinaryOp::Mul,
        arena.make<Binary>(BinaryOp::Mul, sqrt_d, half), x_minus_C);
    ExprPtr tanh_expr = arena.make<FuncCall>(BuiltinOp::Tanh,
        std::vector<ExprPtr>{tanh_arg});
    ExprPtr inner = arena.make<Binary>(BinaryOp::Add,
        arena.make<Unary>(UnaryOp::Neg, make_rational(arena, b)),
        arena.make<Binary>(BinaryOp::Mul, sqrt_d, tanh_expr));
    Rational two_c = Rational(BigInt(2)) * c;
    ExprPtr scale = make_rational(arena, Rational(BigInt(1)) / two_c);
    return arena.make<Binary>(BinaryOp::Mul, scale, inner);
}

// Bernoulli reduction:  y' = q_1·y + q_2·y² (q_0 ≡ 0).  Substitute v = 1/y:
//   v' = -y'/y² = -(q_1/y + q_2) = -q_1·v - q_2,
// a linear 1st-order ODE in v.  After v(x) is known, y = 1/v.
[[nodiscard]] Result<ExprPtr> solve_bernoulli_reduction(
    ExprPtr q1, ExprPtr q2, const Symbol& y, const Symbol& x,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    // Synthesize the linear ODE  v' + q_1·v + q_2 = 0  ⇒  P=q_1, Q=-q_2.
    OdeClassification fake(OdeType::Linear1stOrder,
        /*equation*/arena.make<IntegerLit>(BigInt(0)), y, x);
    fake.components.push_back(q1);
    fake.components.push_back(arena.make<Unary>(UnaryOp::Neg, q2));
    auto v_sol = solve_ode_1st_order(fake, ctx);
    if (v_sol.is_error()) return v_sol;
    // solve_ode_1st_order returns v(x); reciprocate.
    ExprPtr y_sol = arena.make<Binary>(BinaryOp::Div,
        arena.make<IntegerLit>(BigInt(1)), v_sol.value());
    auto simp = ctx.simplify(y_sol);
    if (simp.is_error()) return simp;
    return ok(arena.make<Binary>(BinaryOp::Equal,
        arena.make<Symbol>(y.name), simp.value()));
}

}  // namespace

Result<ExprPtr> solve_ode_nonlinear(
    const OdeClassification& cls,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    if (cls.type == OdeType::Riccati) {
        if (cls.components.size() != 3U) {
            return fail<ExprPtr>(make_error(CASErrorKind::InternalError,
                "Riccati classification missing q_0/q_1/q_2 components"));
        }
        ExprPtr q0 = cls.components[0];
        ExprPtr q1 = cls.components[1];
        ExprPtr q2 = cls.components[2];

        // Sub-family 1: q_0 ≡ 0 → reduce to linear via v = 1/y.
        if (is_zero_expr(q0, ctx)) {
            return solve_bernoulli_reduction(q1, q2, cls.y, cls.x, ctx);
        }

        // Sub-family 2: all coefficients constant in x.
        auto a = as_rational(q0, ctx);
        auto b = as_rational(q1, ctx);
        auto c = as_rational(q2, ctx);
        if (a.has_value() && b.has_value() && c.has_value()) {
            if (c->numerator().is_zero()) {
                return fail<ExprPtr>(make_error(CASErrorKind::InternalError,
                    "Riccati solver invoked with q_2 = 0 (should have been "
                    "classified as Linear1stOrder)"));
            }
            Rational disc = *b * *b - Rational(BigInt(4)) * (*a) * (*c);
            ExprPtr y_sol;
            if (disc.numerator().is_zero()) {
                y_sol = build_double_root_solution(*b, *c, cls.x, ctx);
            } else if (disc.numerator().is_negative()) {
                Rational neg_disc = -disc;
                y_sol = build_negative_disc_solution(*b, *c, neg_disc, cls.x, ctx);
            } else {
                y_sol = build_positive_disc_solution(*b, *c, disc, cls.x, ctx);
            }
            auto simp = ctx.simplify(y_sol);
            if (simp.is_error()) return simp;
            return ok(arena.make<Binary>(BinaryOp::Equal,
                arena.make<Symbol>(cls.y.name), simp.value()));
        }

        // General Riccati with variable coefficients and no known particular
        // solution: explicit Unimplemented with diagnostic (no silent failure,
        // no hardcode-of-passage).
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Riccati with variable coefficients and no particular "
                       "solution: requires 2nd-order linear variable-coefficient "
                       "solver (tracked in F5.3 / B2 continuation) or a "
                       "particular-solution oracle (Risch DE, Bronstein ch.8)",
            .hint = std::string("Provide a particular solution y_p via assumption "
                                "or specialise the equation (q_0 = 0, constant "
                                "coefficients, or known Airy/Bessel reduction).")});
    }

    // Clairaut / d'Alembert solvers will hook in here in a follow-up commit.
    return fail<ExprPtr>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Nonlinear ODE family not yet supported by solve_ode_nonlinear",
        .hint = std::nullopt});
}

}  // namespace cas::calculus
