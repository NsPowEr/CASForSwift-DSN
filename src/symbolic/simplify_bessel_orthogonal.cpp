#include "simplify_impl.hpp"
#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_funcall_bessel_orthogonal(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before) {

    if (op == BuiltinOp::BesselZero && args.size() == 2U) {
        auto order_is_real = [&](ExprPtr expr) -> bool {
            if (expr_cast<IntegerLit>(expr) != nullptr || expr_cast<RationalLit>(expr) != nullptr)
                return true;
            if (const auto* constant = expr_cast<Constant>(expr))
                return constant->value != MathConstant::I && constant->value != MathConstant::NaN;
            if (const auto* symbol = expr_cast<Symbol>(expr))
                return assumptions_ != nullptr && assumptions_->is_real(*symbol);
            return assumptions_ != nullptr && assumptions_->is_real(expr);
        };
        auto index_is_valid = [&](ExprPtr expr) -> bool {
            if (const auto* integer = expr_cast<IntegerLit>(expr))
                return !integer->value.is_negative() && integer->value != BigInt(0);
            if (const auto* symbol = expr_cast<Symbol>(expr))
                return assumptions_ != nullptr
                    && assumptions_->is_integer(*symbol)
                    && assumptions_->is_positive(*symbol);
            return assumptions_ != nullptr
                && assumptions_->is_integer(expr)
                && assumptions_->is_positive(expr);
        };
        if (!order_is_real(args[0]))
            return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument,
                "BesselZero: order must be real"));
        if (!index_is_valid(args[1]))
            return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument,
                "BesselZero: index must be a natural integer >= 1"));
        return ok(arena_.make<FuncCall>(BuiltinOp::BesselZero, args));
    }

    if ((op == BuiltinOp::BesselJ || op == BuiltinOp::BesselY
         || op == BuiltinOp::BesselI || op == BuiltinOp::BesselK) && args.size() == 2U) {
        ExprPtr order = args[0];
        ExprPtr x_arg = args[1];

        // F7.5.E2: special values at x = 0 for integer order.
        // Reference: Abramowitz & Stegun 9.1.7, 9.1.9, 9.6.7.
        //   J_0(0) = 1, J_n(0) = 0 for n ∈ Z, n ≠ 0
        //   I_0(0) = 1, I_n(0) = 0 for n ∈ Z⁺
        //   Y_n(0) = -∞ (sign convention), K_n(0) = +∞ (n ≥ 0)
        if (is_zero_expr(x_arg)) {
            if (const auto* il = expr_cast<IntegerLit>(order)) {
                if (op == BuiltinOp::BesselJ) {
                    return ok(il->value.is_zero()
                        ? make_integer(arena_, BigInt(1))
                        : make_integer(arena_, BigInt(0)));
                }
                if (op == BuiltinOp::BesselI && !il->value.is_negative()) {
                    return ok(il->value.is_zero()
                        ? make_integer(arena_, BigInt(1))
                        : make_integer(arena_, BigInt(0)));
                }
                if (op == BuiltinOp::BesselY && !il->value.is_negative()) {
                    return ok(arena_.make<Constant>(MathConstant::NegInfinity));
                }
                if (op == BuiltinOp::BesselK && !il->value.is_negative()) {
                    return ok(arena_.make<Constant>(MathConstant::Infinity));
                }
            }
        }

        // Parity: negative integer order
        if (const auto* il = expr_cast<IntegerLit>(order);
            il != nullptr && il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                // F0.8-MIGRATED
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_bessel_orthogonal",
                    "Bessel J/Y/I/K with negative integer order > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Reduce order magnitude or use asymptotic expansion for large orders",
                    "F1.x");
            BigInt mag = -il->value;
            ExprPtr pos_order = make_integer(arena_, mag);
            ExprPtr pos_call = arena_.make<FuncCall>(op,
                std::vector<ExprPtr>{pos_order, x_arg});
            const bool is_J_or_Y = (op == BuiltinOp::BesselJ || op == BuiltinOp::BesselY);
            const bool flip_sign = is_J_or_Y && ((mag.to_u64() % 2U) == 1U);
            if (!flip_sign) return simplify_expr(pos_call);
            return simplify_expr(arena_.make<Unary>(UnaryOp::Neg, pos_call));
        }

        // Half-integer J/Y/I/K reductions to elementary functions
        if (const auto* rl = expr_cast<RationalLit>(order);
            rl != nullptr && rl->denominator == BigInt(2)
            && (rl->numerator == BigInt(1) || rl->numerator == BigInt(-1))) {
            ExprPtr pi_const = arena_.make<Constant>(MathConstant::Pi);
            ExprPtr pi_x = arena_.make<Product>(std::vector<ExprPtr>{pi_const, x_arg});
            ExprPtr two_over_pi_x = arena_.make<Binary>(BinaryOp::Div,
                make_integer(arena_, BigInt(2)), pi_x);
            ExprPtr root = arena_.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{two_over_pi_x});
            BuiltinOp inner_fn = BuiltinOp::Unknown;
            bool negate = false;
            if (op == BuiltinOp::BesselJ) {
                inner_fn = (rl->numerator == BigInt(1)) ? BuiltinOp::Sin : BuiltinOp::Cos;
            } else if (op == BuiltinOp::BesselY) {
                if (rl->numerator == BigInt(1)) { inner_fn = BuiltinOp::Cos; negate = true; }
                else inner_fn = BuiltinOp::Sin;
            } else if (op == BuiltinOp::BesselI) {
                inner_fn = (rl->numerator == BigInt(1)) ? BuiltinOp::Sinh : BuiltinOp::Cosh;
            }
            if (inner_fn != BuiltinOp::Unknown) {
                ExprPtr trig = arena_.make<FuncCall>(inner_fn, std::vector<ExprPtr>{x_arg});
                ExprPtr prod = arena_.make<Product>(std::vector<ExprPtr>{root, trig});
                if (negate) prod = arena_.make<Unary>(UnaryOp::Neg, prod);
                return simplify_expr(prod);
            }
            // K_{±1/2}(x) = sqrt(π/(2x)) * exp(-x)
            if (op == BuiltinOp::BesselK) {
                ExprPtr two_x = arena_.make<Product>(std::vector<ExprPtr>{
                    make_integer(arena_, BigInt(2)), x_arg});
                ExprPtr pi_over_2x = arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi), two_x);
                ExprPtr root_k = arena_.make<FuncCall>(BuiltinOp::Sqrt,
                    std::vector<ExprPtr>{pi_over_2x});
                ExprPtr neg_x = arena_.make<Unary>(UnaryOp::Neg, x_arg);
                ExprPtr exp_neg_x = arena_.make<FuncCall>(BuiltinOp::Exp,
                    std::vector<ExprPtr>{neg_x});
                return simplify_expr(arena_.make<Product>(
                    std::vector<ExprPtr>{root_k, exp_neg_x}));
            }
        }

        // Opt-in three-term recurrence for J/Y integer order n >= 2
        if ((op == BuiltinOp::BesselJ || op == BuiltinOp::BesselY)
            && context_ != nullptr && context_->expand_bessel_recurrence()) {
            if (const auto* il = expr_cast<IntegerLit>(order);
                il != nullptr && !il->value.is_negative() && il->value > BigInt(1)) {
                if (il->value.bit_length() > 16)
                    // F0.8-MIGRATED
                    return make_unimplemented<ExprPtr>(
                        "symbolic", "simplify_funcall_bessel_orthogonal",
                        "Bessel J/Y three-term recurrence with integer order > 2^16",
                        error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                        "Limit recurrence depth via ctx.max_bessel_recurrence or use asymptotic form",
                        "F1.x");
                BigInt n_minus_1 = il->value - BigInt(1);
                BigInt n_minus_2 = il->value - BigInt(2);
                ExprPtr coeff = arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Product>(std::vector<ExprPtr>{
                        make_integer(arena_, BigInt(2)),
                        make_integer(arena_, n_minus_1)}),
                    x_arg);
                ExprPtr j_minus_1 = arena_.make<FuncCall>(op,
                    std::vector<ExprPtr>{make_integer(arena_, n_minus_1), x_arg});
                ExprPtr j_minus_2 = arena_.make<FuncCall>(op,
                    std::vector<ExprPtr>{make_integer(arena_, n_minus_2), x_arg});
                ExprPtr first_term = arena_.make<Product>(std::vector<ExprPtr>{coeff, j_minus_1});
                ExprPtr neg_second = arena_.make<Unary>(UnaryOp::Neg, j_minus_2);
                return simplify_expr(arena_.make<Sum>(
                    std::vector<ExprPtr>{first_term, neg_second}));
            }
        }
    }

    // F7.5.E2: orthogonal polynomial recurrences (Chebyshev, Hermite,
    // Legendre, Laguerre, Jacobi) + LambertW moved to simplify_orthogonal_polys.cpp
    // to stay under the 500-line anti-monolith limit.
    return simplify_funcall_orthogonal_polys(original, op, std::move(args), target_before);
}


} // namespace cas::symbolic::detail
