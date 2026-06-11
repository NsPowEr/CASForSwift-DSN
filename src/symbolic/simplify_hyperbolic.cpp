// F7.5.A4 follow-up — exact-value identities for the hyperbolic builtins.
//
// Closes the limit-corpus regressions where lim x→0 of cosh(x), sinh(x),
// tanh(x), etc. left the canonical-value substitution unevaluated
// (CAS returned `cosh(0)` instead of `1`, `sinh(0)` instead of `0`).
// The simplifier dispatcher previously had no hyperbolic branch — the
// switch in `simplify_funcall` fell through to default, leaving the
// FuncCall intact. This file adds:
//
//   sinh(0) = 0       sinh(-x) = -sinh(x)
//   cosh(0) = 1       cosh(-x) = cosh(x)     (even)
//   tanh(0) = 0       tanh(-x) = -tanh(x)
//   coth(-x) = -coth(x)  (coth(0) is undefined — leave unchanged)
//
// All rules are exact algebraic identities (Bronstein §4.3), not
// numerical approximations.

#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_funcall_hyperbolic(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr /*target_before*/)
{
    if (args.size() != 1U) {
        const auto& orig_args = expr_ref<FuncCall>(original).args;
        if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
        return ok(arena_.make<FuncCall>(op, std::move(args)));
    }

    ExprPtr a = args.front();

    // Zero-value identities.
    if (is_zero_expr(a)) {
        if (op == BuiltinOp::Sinh || op == BuiltinOp::Tanh) {
            return ok(make_integer(arena_, BigInt(0)));
        }
        if (op == BuiltinOp::Cosh) {
            return ok(make_integer(arena_, BigInt(1)));
        }
        // coth(0) is undefined (singularity) — leave the FuncCall intact
        // so downstream code can decide whether to treat it as
        // ComplexInfinity in an extended-real context.
    }

    // Parity reduction: sinh, tanh, coth are odd; cosh is even.
    if (const auto* u = expr_cast<Unary>(a); u != nullptr && u->op == UnaryOp::Neg) {
        ExprPtr inner_call = arena_.make<FuncCall>(
            op, std::vector<ExprPtr>{u->operand});
        if (op == BuiltinOp::Cosh) {
            return simplify_expr(inner_call);
        }
        if (op == BuiltinOp::Sinh || op == BuiltinOp::Tanh ||
            op == BuiltinOp::Coth) {
            return simplify_expr(arena_.make<Unary>(UnaryOp::Neg, inner_call));
        }
    }

    // No rule fired: preserve identity when args are pointer-equal to
    // the original FuncCall args.
    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

}  // namespace cas::symbolic::detail
