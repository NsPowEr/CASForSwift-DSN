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

        // Half-integer J/Y/I/K reductions to elementary closed form.
        // RationalLit invariant: denominator > 0 and gcd(numerator, denom) = 1.
        // With denom = 2 the numerator is therefore odd, hence ν = num/2 is
        // a half-integer. Reduction:
        //   |num| == 1  →  direct base-case build (J=sin/cos, Y=∓cos/sin,
        //                                          I=sinh/cosh, K=√π/(2x)·e^{-x})
        //   |num| > 1   →  iterate three-term recurrence from base ±1/2
        // Bounded by ctx.max_bessel_half_integer_order() — no silent truncation.
        if (const auto* rl = expr_cast<RationalLit>(order);
            rl != nullptr && rl->denominator == BigInt(2)) {
            auto reduced = bessel_half_integer_reduce(op, rl->numerator, x_arg);
            if (reduced.is_ok()) return reduced;
            // On Unimplemented (bound exceeded) propagate; on any other error too.
            return reduced;
        }

        // Opt-in three-term recurrence for J/Y/I/K with integer order n >= 2.
        // References (DLMF 10.6.1, 10.29.1):
        //   J_n = (2(n-1)/x)·J_{n-1} - J_{n-2}        (same for Y)
        //   I_n = I_{n-2} - (2(n-1)/x)·I_{n-1}
        //   K_n = K_{n-2} + (2(n-1)/x)·K_{n-1}
        if (context_ != nullptr && context_->expand_bessel_recurrence()) {
            if (const auto* il = expr_cast<IntegerLit>(order);
                il != nullptr && !il->value.is_negative() && il->value > BigInt(1)) {
                if (il->value.bit_length() > 16)
                    return make_unimplemented<ExprPtr>(
                        "symbolic", "simplify_funcall_bessel_orthogonal",
                        "Bessel three-term recurrence with integer order > 2^16",
                        error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                        "Use asymptotic form for large integer orders",
                        "F1.x");
                BigInt n_minus_1 = il->value - BigInt(1);
                BigInt n_minus_2 = il->value - BigInt(2);
                ExprPtr coeff = arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Product>(std::vector<ExprPtr>{
                        make_integer(arena_, BigInt(2)),
                        make_integer(arena_, n_minus_1)}),
                    x_arg);
                ExprPtr b_minus_1 = arena_.make<FuncCall>(op,
                    std::vector<ExprPtr>{make_integer(arena_, n_minus_1), x_arg});
                ExprPtr b_minus_2 = arena_.make<FuncCall>(op,
                    std::vector<ExprPtr>{make_integer(arena_, n_minus_2), x_arg});
                ExprPtr coeff_times = arena_.make<Product>(std::vector<ExprPtr>{coeff, b_minus_1});

                ExprPtr expansion;
                switch (op) {
                    case BuiltinOp::BesselJ:
                    case BuiltinOp::BesselY:
                        expansion = arena_.make<Sum>(std::vector<ExprPtr>{
                            coeff_times,
                            arena_.make<Unary>(UnaryOp::Neg, b_minus_2)});
                        break;
                    case BuiltinOp::BesselI:
                        expansion = arena_.make<Sum>(std::vector<ExprPtr>{
                            b_minus_2,
                            arena_.make<Unary>(UnaryOp::Neg, coeff_times)});
                        break;
                    case BuiltinOp::BesselK:
                        expansion = arena_.make<Sum>(std::vector<ExprPtr>{
                            b_minus_2, coeff_times});
                        break;
                    default: break;
                }
                if (expansion) return simplify_expr(expansion);
            }
        }
    }

    // F7.5.E2: orthogonal polynomial recurrences (Chebyshev, Hermite,
    // Legendre, Laguerre, Jacobi) + LambertW moved to simplify_orthogonal_polys.cpp
    // to stay under the 500-line anti-monolith limit.
    return simplify_funcall_orthogonal_polys(original, op, std::move(args), target_before);
}

// ─────────────────────────────────────────────────────────────────────────────
// F7.5.E2: half-integer Bessel reduction.
//
// For ν = p_num/2 with p_num odd and |p_num| ≥ 1 the value Bν(x) admits an
// elementary closed form. Base cases (DLMF 10.16.1, 10.39.1):
//
//   J_{ 1/2}(x) =  √(2/(πx))·sin(x)        J_{-1/2}(x) =  √(2/(πx))·cos(x)
//   Y_{ 1/2}(x) = -√(2/(πx))·cos(x)        Y_{-1/2}(x) =  √(2/(πx))·sin(x)
//   I_{ 1/2}(x) =  √(2/(πx))·sinh(x)       I_{-1/2}(x) =  √(2/(πx))·cosh(x)
//   K_{±1/2}(x) =  √(π/(2x))·exp(-x)
//
// Higher half-integer orders are generated by the standard three-term
// recurrence (DLMF 10.6.1, 10.29.1) applied iteratively. With ν = m/2 (m odd),
// the coefficient 2ν is the integer m and the recurrence becomes:
//
//   J_{ν+1} =  (m/x)·J_ν - J_{ν-1}         (same shape for Y)
//   I_{ν+1} =  I_{ν-1} - (m/x)·I_ν
//   K_{ν+1} =  K_{ν-1} + (m/x)·K_ν
//
// Backward (ν → ν-1) inverts each formula. Iteration count is (|p_num|-1)/2.
// Bound: ctx.max_bessel_half_integer_order() guards against pathological input;
// exceeding the bound produces an explicit Unimplemented diagnostic rather than
// silent truncation.
Result<ExprPtr> Simplifier::bessel_half_integer_reduce(
    BuiltinOp op, const BigInt& p_num, ExprPtr x_arg) {

    // Build closed-form base for ν = sign·1/2 (sign = +1 or -1).
    auto build_base = [&](int sign) -> ExprPtr {
        // Common scalar √(2/(πx)) for J/Y/I; √(π/(2x)) for K.
        ExprPtr pi_const = arena_.make<Constant>(MathConstant::Pi);
        if (op == BuiltinOp::BesselK) {
            ExprPtr two_x = arena_.make<Product>(std::vector<ExprPtr>{
                make_integer(arena_, BigInt(2)), x_arg});
            ExprPtr pi_over_2x = arena_.make<Binary>(BinaryOp::Div, pi_const, two_x);
            ExprPtr root = arena_.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{pi_over_2x});
            ExprPtr neg_x = arena_.make<Unary>(UnaryOp::Neg, x_arg);
            ExprPtr exp_neg_x = arena_.make<FuncCall>(BuiltinOp::Exp,
                std::vector<ExprPtr>{neg_x});
            return arena_.make<Product>(std::vector<ExprPtr>{root, exp_neg_x});
        }
        ExprPtr pi_x = arena_.make<Product>(std::vector<ExprPtr>{pi_const, x_arg});
        ExprPtr two_over_pi_x = arena_.make<Binary>(BinaryOp::Div,
            make_integer(arena_, BigInt(2)), pi_x);
        ExprPtr root = arena_.make<FuncCall>(BuiltinOp::Sqrt,
            std::vector<ExprPtr>{two_over_pi_x});
        BuiltinOp inner = BuiltinOp::Unknown;
        bool negate = false;
        switch (op) {
            case BuiltinOp::BesselJ:
                inner = (sign > 0) ? BuiltinOp::Sin : BuiltinOp::Cos; break;
            case BuiltinOp::BesselY:
                if (sign > 0) { inner = BuiltinOp::Cos; negate = true; }
                else inner = BuiltinOp::Sin;
                break;
            case BuiltinOp::BesselI:
                inner = (sign > 0) ? BuiltinOp::Sinh : BuiltinOp::Cosh; break;
            default: return ExprPtr{};
        }
        ExprPtr trig = arena_.make<FuncCall>(inner, std::vector<ExprPtr>{x_arg});
        ExprPtr prod = arena_.make<Product>(std::vector<ExprPtr>{root, trig});
        if (negate) prod = arena_.make<Unary>(UnaryOp::Neg, prod);
        return prod;
    };

    // Base cases.
    if (p_num == BigInt(1))  return simplify_expr(build_base(+1));
    if (p_num == BigInt(-1)) return simplify_expr(build_base(-1));

    // Bound check (configurable, no silent truncation).
    const unsigned int bound = (context_ != nullptr)
        ? context_->max_bessel_half_integer_order()
        : 64U;
    BigInt abs_num = p_num.is_negative() ? -p_num : p_num;
    if (abs_num.bit_length() > 16U
        || abs_num > BigInt(static_cast<long long>(bound))) {
        return make_unimplemented<ExprPtr>(
            "symbolic", "bessel_half_integer_reduce",
            "Bessel half-integer order exceeds configured bound",
            error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
            "Increase ctx.max_bessel_half_integer_order() or reduce input magnitude",
            "F7.5.E2");
    }

    const long long target = p_num.is_negative()
        ? -static_cast<long long>(abs_num.to_u64())
        :  static_cast<long long>(abs_num.to_u64());
    const bool forward = target > 0;

    // Initial state:
    //   forward  : b_curr = B_{+1/2}, b_prev = B_{-1/2}, nu = +1/2  (encoded as 1)
    //   backward : b_curr = B_{-1/2}, b_prev = B_{+1/2}, nu = -1/2  (encoded as -1)
    // In both cases "b_prev" refers to B at ν shifted opposite to the step
    // direction, so the same shift `(b_prev, b_curr) ← (b_curr, b_next)`
    // applies after each iteration.
    ExprPtr b_prev = forward ? build_base(-1) : build_base(+1);
    ExprPtr b_curr = forward ? build_base(+1) : build_base(-1);
    long long nu_num = forward ? 1 : -1;  // numerator of current ν (over 2)

    const long long step = forward ? +2 : -2;
    while (nu_num != target) {
        // 2ν as integer = nu_num.
        ExprPtr coeff = arena_.make<Binary>(BinaryOp::Div,
            make_integer(arena_, BigInt(nu_num)), x_arg);
        ExprPtr coeff_curr = arena_.make<Product>(std::vector<ExprPtr>{coeff, b_curr});

        ExprPtr b_next;
        switch (op) {
            case BuiltinOp::BesselJ:
            case BuiltinOp::BesselY:
                // Forward : b_{ν+1} = (2ν/x)·b_ν - b_{ν-1}
                // Backward: b_{ν-1} = (2ν/x)·b_ν - b_{ν+1}  (same shape on stored prev)
                b_next = arena_.make<Sum>(std::vector<ExprPtr>{
                    coeff_curr,
                    arena_.make<Unary>(UnaryOp::Neg, b_prev)});
                break;
            case BuiltinOp::BesselI:
                if (forward) {
                    // b_{ν+1} = b_{ν-1} - (2ν/x)·b_ν
                    b_next = arena_.make<Sum>(std::vector<ExprPtr>{
                        b_prev,
                        arena_.make<Unary>(UnaryOp::Neg, coeff_curr)});
                } else {
                    // b_{ν-1} = b_{ν+1} + (2ν/x)·b_ν
                    b_next = arena_.make<Sum>(std::vector<ExprPtr>{b_prev, coeff_curr});
                }
                break;
            case BuiltinOp::BesselK:
                if (forward) {
                    // b_{ν+1} = b_{ν-1} + (2ν/x)·b_ν
                    b_next = arena_.make<Sum>(std::vector<ExprPtr>{b_prev, coeff_curr});
                } else {
                    // b_{ν-1} = b_{ν+1} - (2ν/x)·b_ν
                    b_next = arena_.make<Sum>(std::vector<ExprPtr>{
                        b_prev,
                        arena_.make<Unary>(UnaryOp::Neg, coeff_curr)});
                }
                break;
            default:
                return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument,
                    "bessel_half_integer_reduce: unsupported Bessel op"));
        }

        b_prev = b_curr;
        b_curr = b_next;
        nu_num += step;
    }

    return simplify_expr(b_curr);
}

} // namespace cas::symbolic::detail
