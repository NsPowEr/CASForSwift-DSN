#include "simplify_impl.hpp"

// F1-DEBT1 split: arc-trig inverse simplifications extracted from simplify_trig.cpp
// to keep simplify_trig.cpp ≤500 LOC.
//
// Implements Simplifier::simplify_funcall_arc_trig for:
//   Asin, Acos, Atan
// All rules are standard branch-cut / identity reductions.

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_funcall_arc_trig(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr /*target_before*/)
{
    if (op == BuiltinOp::Asin && args.size() == 1U) {
        // asin(0) = 0
        if (is_zero_expr(args.front())) return ok(make_integer(arena_, BigInt(0)));
        // asin(1) = π/2, asin(-1) = -π/2
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            if (il->value == BigInt(1)) {
                return ok(arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi),
                    make_integer(arena_, BigInt(2))));
            }
            if (il->value == BigInt(-1)) {
                return ok(arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Binary>(BinaryOp::Div,
                        arena_.make<Constant>(MathConstant::Pi),
                        make_integer(arena_, BigInt(2)))));
            }
        }
        // asin(sin(x)) -> x if -π/2 ≤ x ≤ π/2
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Sin) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi_2 = arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi), make_integer(arena_, BigInt(2)));
                ExprPtr neg_pi_2 = arena_.make<Unary>(UnaryOp::Neg, pi_2);
                if (assumptions_->is_greater_equal(x, neg_pi_2)
                    && assumptions_->is_greater_equal(pi_2, x))
                    return ok(x);
            }
        }
    }
    if (op == BuiltinOp::Acos && args.size() == 1U) {
        // acos(0) = π/2, acos(1) = 0, acos(-1) = π
        if (is_zero_expr(args.front())) {
            return ok(arena_.make<Binary>(BinaryOp::Div,
                arena_.make<Constant>(MathConstant::Pi),
                make_integer(arena_, BigInt(2))));
        }
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            if (il->value == BigInt(1)) return ok(make_integer(arena_, BigInt(0)));
            if (il->value == BigInt(-1))
                return ok(arena_.make<Constant>(MathConstant::Pi));
        }
        // acos(cos(x)) -> x if 0 ≤ x ≤ π
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Cos) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi = arena_.make<Constant>(MathConstant::Pi);
                ExprPtr zero = make_integer(arena_, BigInt(0));
                if (assumptions_->is_greater_equal(x, zero)
                    && assumptions_->is_greater_equal(pi, x))
                    return ok(x);
            }
        }
    }
    if (op == BuiltinOp::Atan && args.size() == 1U) {
        // atan(+∞) = π/2, atan(-∞) = -π/2 (canonical real-axis limits).
        if (const auto* c = expr_cast<Constant>(args.front());
            c != nullptr && c->value == MathConstant::Infinity) {
            return ok(arena_.make<Binary>(BinaryOp::Div,
                arena_.make<Constant>(MathConstant::Pi),
                make_integer(arena_, BigInt(2))));
        }
        if (const auto* u = expr_cast<Unary>(args.front());
            u != nullptr && u->op == UnaryOp::Neg) {
            if (const auto* c = expr_cast<Constant>(u->operand);
                c != nullptr && c->value == MathConstant::Infinity) {
                return ok(arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Binary>(BinaryOp::Div,
                        arena_.make<Constant>(MathConstant::Pi),
                        make_integer(arena_, BigInt(2)))));
            }
        }
        // atan(tan(x)) -> x if -π/2 < x < π/2
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Tan) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi_2 = arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi), make_integer(arena_, BigInt(2)));
                ExprPtr neg_pi_2 = arena_.make<Unary>(UnaryOp::Neg, pi_2);
                if (assumptions_->is_greater(x, neg_pi_2)
                    && assumptions_->is_greater(pi_2, x))
                    return ok(x);
            }
        }
        if (is_zero_expr(args.front())) return ok(make_integer(arena_, BigInt(0)));
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            if (il->value == BigInt(1))
                return simplify_expr(arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi),
                    make_integer(arena_, BigInt(4))));
            if (il->value == BigInt(-1))
                return simplify_expr(arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Binary>(BinaryOp::Div,
                        arena_.make<Constant>(MathConstant::Pi),
                        make_integer(arena_, BigInt(4)))));
        }
        // atan special values for the constructible √3 family, completing the
        // exact table {0, π/4, π/2} with {π/6, π/3}. Detection is form-independent:
        // we test tan² = arg², keyed by sign. For x>0:
        //   atan(x)=π/6 ⟺ x=1/√3 ⟺ x²=1/3 ;  atan(x)=π/3 ⟺ x=√3 ⟺ x²=3.
        // The positive square root is unique, so this is exact (no false match).
        // arg arrives as the ratio b/a from arg(a+bi); a bare Pow(Div,2) does NOT
        // collapse √3 in the denominator, so for a Div we square numerator and
        // denominator separately (each atomic radical reduces: √3²→3) and divide.
        // Negative arguments (Q2/Q3 ratios) take the mirrored angle, independent of
        // whether the form is Unary::Neg or a negative Product (the odd-symmetry
        // rule below only catches the former). Anything else falls through.
        {
            const bool pos = is_known_positive(args.front());
            const bool neg = is_known_negative(args.front());
            if (pos || neg) {
                // Robustly compute tan² = arg², independent of the (many) surface
                // forms a √3 ratio can take. (−x)²=x², and atomic radicals collapse
                // (√3²→3) where a composite Pow(Neg,·)/Pow(Div,2) may not, so:
                //   • strip a leading Neg;
                //   • for base^exp, square as (sign-stripped base)^(2·exp);
                //   • for a Div, square numerator and denominator separately.
                auto strip_neg = [](ExprPtr e) -> ExprPtr {
                    if (const auto* u = expr_cast<Unary>(e); u != nullptr && u->op == UnaryOp::Neg)
                        return u->operand;
                    return e;
                };
                auto sq_atom = [&](ExprPtr e) -> std::optional<ExprPtr> {
                    ExprPtr b = strip_neg(e);
                    if (const auto* p = expr_cast<Binary>(b); p != nullptr && p->op == BinaryOp::Pow) {
                        auto two_exp = simplify_expr(arena_.make<Binary>(BinaryOp::Mul,
                            p->right, make_integer(arena_, BigInt(2))));
                        if (two_exp.is_ok()) {
                            auto q = simplify_expr(arena_.make<Binary>(BinaryOp::Pow,
                                strip_neg(p->left), two_exp.value()));
                            if (q.is_ok()) return q.value();
                        }
                        return std::nullopt;
                    }
                    auto q = simplify_expr(arena_.make<Binary>(BinaryOp::Pow,
                        b, make_integer(arena_, BigInt(2))));
                    if (q.is_ok()) return q.value();
                    return std::nullopt;
                };
                std::optional<ExprPtr> tan2;
                ExprPtr base_arg = strip_neg(args.front());
                if (const auto* d = expr_cast<Binary>(base_arg); d != nullptr && d->op == BinaryOp::Div) {
                    auto n2 = sq_atom(d->left);
                    auto m2 = sq_atom(d->right);
                    if (n2 && m2) {
                        auto q = simplify_expr(arena_.make<Binary>(BinaryOp::Div, *n2, *m2));
                        if (q.is_ok()) tan2 = q.value();
                    }
                } else {
                    tan2 = sq_atom(base_arg);
                }
                if (tan2) {
                    int denom = 0;  // emit π/denom, 0 = no match
                    if (const auto* il = expr_cast<IntegerLit>(*tan2)) {
                        if (il->value == BigInt(3)) denom = 3;          // tan²=3 → π/3
                    } else if (const auto* rl = expr_cast<RationalLit>(*tan2)) {
                        if (rl->numerator == BigInt(1) && rl->denominator == BigInt(3))
                            denom = 6;                                   // tan²=1/3 → π/6
                    }
                    if (denom != 0) {
                        ExprPtr angle = arena_.make<Binary>(BinaryOp::Div,
                            arena_.make<Constant>(MathConstant::Pi),
                            make_integer(arena_, BigInt(denom)));
                        if (neg) angle = arena_.make<Unary>(UnaryOp::Neg, angle);
                        return simplify_expr(angle);
                    }
                }
            }
        }
        // atan odd: atan(-x) = -atan(x)
        if (const auto* un = expr_cast<Unary>(args.front()); un && un->op == UnaryOp::Neg) {
            ExprPtr inner_atan = arena_.make<FuncCall>(
                BuiltinOp::Atan, std::vector<ExprPtr>{un->operand});
            return simplify_expr(arena_.make<Unary>(UnaryOp::Neg, inner_atan));
        }
    }

    // No rule fired: rebuild unchanged.
    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

} // namespace cas::symbolic::detail
