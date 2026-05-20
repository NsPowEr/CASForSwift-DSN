// CAS-L3-12 — Symbolic finite-difference numeric derivatives.
//
// Builds symbolic approximations to derivatives via finite-difference
// formulas. Useful when symbolic diff fails or produces unwieldy form,
// or when the user wants explicit O(h^k) approximation order.
//
// Formulas (Abramowitz-Stegun 25.3):
//   Forward 1st order:    f'(x) ≈ (f(x+h) - f(x)) / h          O(h)
//   Central 2nd order:    f'(x) ≈ (f(x+h) - f(x-h)) / (2h)     O(h²)
//   Central 4th order:    f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)  O(h⁴)
//
// All formulas return symbolic expressions in x and h; user substitutes
// numeric h or symbolic step.

#include "calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"

namespace cas::calculus {

namespace {

[[nodiscard]] ExprPtr subst_x(ExprPtr expr, const Symbol& var, ExprPtr value,
                              symbolic::CASContext& ctx) {
    auto r = symbolic::substitute(expr, var, value, ctx);
    if (r.is_ok()) return r.value();
    return expr;
}

}  // namespace

[[nodiscard]] Result<ExprPtr> numeric_diff(
    ExprPtr expr, const Symbol& var, ExprPtr h,
    FiniteDiffOrder order, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr x = arena.make<Symbol>(var);

    auto eval_at = [&](ExprPtr offset) -> ExprPtr {
        if (!offset) return expr;
        ExprPtr point = arena.make<Sum>(std::vector<ExprPtr>{x, offset});
        return subst_x(expr, var, point, ctx);
    };

    if (order == FiniteDiffOrder::Forward1) {
        // (f(x+h) - f(x)) / h
        ExprPtr fxh = eval_at(h);
        ExprPtr fx = expr;
        ExprPtr num = arena.make<Binary>(BinaryOp::Sub, fxh, fx);
        ExprPtr result = arena.make<Binary>(BinaryOp::Div, num, h);
        return ctx.simplify(result);
    }
    if (order == FiniteDiffOrder::Central2) {
        // (f(x+h) - f(x-h)) / (2h)
        ExprPtr neg_h = arena.make<Unary>(UnaryOp::Neg, h);
        ExprPtr fxh = eval_at(h);
        ExprPtr fxmh = eval_at(neg_h);
        ExprPtr num = arena.make<Binary>(BinaryOp::Sub, fxh, fxmh);
        ExprPtr two = arena.make<IntegerLit>(BigInt(2));
        ExprPtr two_h = arena.make<Product>(std::vector<ExprPtr>{two, h});
        ExprPtr result = arena.make<Binary>(BinaryOp::Div, num, two_h);
        return ctx.simplify(result);
    }
    if (order == FiniteDiffOrder::Central4) {
        // (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
        ExprPtr two = arena.make<IntegerLit>(BigInt(2));
        ExprPtr two_h = arena.make<Product>(std::vector<ExprPtr>{two, h});
        ExprPtr neg_two_h = arena.make<Unary>(UnaryOp::Neg, two_h);
        ExprPtr neg_h = arena.make<Unary>(UnaryOp::Neg, h);
        ExprPtr f_plus_2h = eval_at(two_h);
        ExprPtr f_plus_h = eval_at(h);
        ExprPtr f_minus_h = eval_at(neg_h);
        ExprPtr f_minus_2h = eval_at(neg_two_h);
        ExprPtr eight = arena.make<IntegerLit>(BigInt(8));
        ExprPtr twelve = arena.make<IntegerLit>(BigInt(12));
        ExprPtr term_a = arena.make<Unary>(UnaryOp::Neg, f_plus_2h);
        ExprPtr term_b = arena.make<Product>(std::vector<ExprPtr>{eight, f_plus_h});
        ExprPtr term_c = arena.make<Unary>(UnaryOp::Neg,
            arena.make<Product>(std::vector<ExprPtr>{eight, f_minus_h}));
        ExprPtr term_d = f_minus_2h;
        ExprPtr num = arena.make<Sum>(std::vector<ExprPtr>{term_a, term_b, term_c, term_d});
        ExprPtr den = arena.make<Product>(std::vector<ExprPtr>{twelve, h});
        ExprPtr result = arena.make<Binary>(BinaryOp::Div, num, den);
        return ctx.simplify(result);
    }
    return fail<ExprPtr>(CASError{
        CASErrorKind::Unimplemented, "Unknown FiniteDiffOrder", std::nullopt});
}

}  // namespace cas::calculus
