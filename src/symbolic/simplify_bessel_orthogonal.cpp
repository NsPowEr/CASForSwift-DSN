#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_funcall_bessel_orthogonal(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr /*target_before*/) {

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

        // Parity: negative integer order
        if (const auto* il = expr_cast<IntegerLit>(order);
            il != nullptr && il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "Bessel: negative integer order too large"));
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
                    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                        "Bessel recurrence: integer order too large"));
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

    // Chebyshev T_n — three-term recurrence, scales to any n
    if (op == BuiltinOp::ChebyshevT && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && !il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "ChebyshevT: degree too large"));
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            if (n == 1U) return ok(x);
            ExprPtr t_prev = make_integer(arena_, BigInt(1));
            ExprPtr t_curr = x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr two_x_tk = arena_.make<Product>(std::vector<ExprPtr>{
                    make_integer(arena_, BigInt(2)), x, t_curr});
                ExprPtr neg_prev = arena_.make<Unary>(UnaryOp::Neg, t_prev);
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{two_x_tk, neg_prev});
                auto s = simplify_expr(next);
                if (s.is_error()) return s;
                t_prev = t_curr;
                t_curr = s.value();
            }
            return ok(t_curr);
        }
    }

    // Chebyshev U_n — three-term recurrence
    if (op == BuiltinOp::ChebyshevU && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && !il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "ChebyshevU: degree too large"));
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            ExprPtr two_x = arena_.make<Product>(std::vector<ExprPtr>{
                make_integer(arena_, BigInt(2)), x});
            if (n == 1U) return simplify_expr(two_x);
            ExprPtr u_prev = make_integer(arena_, BigInt(1));
            ExprPtr u_curr = two_x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr two_x_uk = arena_.make<Product>(std::vector<ExprPtr>{
                    make_integer(arena_, BigInt(2)), x, u_curr});
                ExprPtr neg_prev = arena_.make<Unary>(UnaryOp::Neg, u_prev);
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{two_x_uk, neg_prev});
                auto s = simplify_expr(next);
                if (s.is_error()) return s;
                u_prev = u_curr;
                u_curr = s.value();
            }
            return ok(u_curr);
        }
    }

    // Hermite H_n (physicist) — three-term recurrence
    if (op == BuiltinOp::HermiteH && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && !il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "HermiteH: degree too large"));
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            ExprPtr two_x = arena_.make<Product>(std::vector<ExprPtr>{
                make_integer(arena_, BigInt(2)), x});
            if (n == 1U) return simplify_expr(two_x);
            ExprPtr h_prev = make_integer(arena_, BigInt(1));
            ExprPtr h_curr = two_x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr two_x_hk = arena_.make<Product>(std::vector<ExprPtr>{
                    make_integer(arena_, BigInt(2)), x, h_curr});
                ExprPtr two_k = make_integer(arena_, BigInt(static_cast<std::int64_t>(2U * k)));
                ExprPtr neg = arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Product>(std::vector<ExprPtr>{two_k, h_prev}));
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{two_x_hk, neg});
                auto s = simplify_expr(next);
                if (s.is_error()) return s;
                h_prev = h_curr;
                h_curr = s.value();
            }
            return ok(h_curr);
        }
    }

    // Hermite He_n (probabilist) — three-term recurrence
    if (op == BuiltinOp::HermiteHe && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && !il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "HermiteHe: degree too large"));
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            if (n == 1U) return ok(x);
            ExprPtr h_prev = make_integer(arena_, BigInt(1));
            ExprPtr h_curr = x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr x_hk = arena_.make<Product>(std::vector<ExprPtr>{x, h_curr});
                ExprPtr k_const = make_integer(arena_, BigInt(static_cast<std::int64_t>(k)));
                ExprPtr neg = arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Product>(std::vector<ExprPtr>{k_const, h_prev}));
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{x_hk, neg});
                auto s = simplify_expr(next);
                if (s.is_error()) return s;
                h_prev = h_curr;
                h_curr = s.value();
            }
            return ok(h_curr);
        }
    }

    // Legendre P_n — Bonnet recurrence, pure algorithm, no polynomial table
    if (op == BuiltinOp::LegendreP && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && il->value >= BigInt(0)) {
            if (il->value.bit_length() > 16)
                return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
                    "LegendreP: degree too large for symbolic expansion"));
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            if (n == 1U) return ok(x);
            ExprPtr p_prev = make_integer(arena_, BigInt(1));
            ExprPtr p_curr = x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                Rational a(BigInt(static_cast<std::int64_t>(2U * k + 1U)),
                           BigInt(static_cast<std::int64_t>(k + 1U)));
                Rational b(BigInt(static_cast<std::int64_t>(k)),
                           BigInt(static_cast<std::int64_t>(k + 1U)));
                ExprPtr term1 = arena_.make<Product>(std::vector<ExprPtr>{
                    make_rational(arena_, a), x, p_curr});
                ExprPtr neg_term2 = arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Product>(std::vector<ExprPtr>{make_rational(arena_, b), p_prev}));
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{term1, neg_term2});
                auto simp = simplify_expr(next);
                if (simp.is_error()) return simp;
                p_prev = p_curr;
                p_curr = simp.value();
            }
            return ok(p_curr);
        }
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

} // namespace cas::symbolic::detail
