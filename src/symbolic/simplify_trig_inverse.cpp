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
