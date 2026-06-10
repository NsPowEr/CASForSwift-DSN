// F7.5.E1: split from simplify_special_fn.cpp.
// Hosts hypergeometric (0F1, 1F1, 2F1) and elliptic (K, E, F, Π) identity
// rewrites — see Abramowitz-Stegun §15.1 and §17.4–§17.7.

#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_funcall_hyper_elliptic(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before) {
    (void)target_before;

    // Hypergeometric functions (F5.9 / Task #18).
    //
    //   0F1(;b;z),  1F1(a;b;z),  2F1(a,b;c;z).
    //
    // Identità riconosciute (Abramowitz-Stegun §15.1):
    //   pFq(...; 0) = 1                          (valore base in z = 0)
    //   1F1(a; a; z) = e^z                       (Kummer degenerate)
    //   2F1(a, b; b; z) = (1 − z)^(−a)           (parametri degeneri)
    //   2F1(1, 1; 2; -z) = ln(1 + z) / z         (log closed-form)
    //
    // Non riconosciute: lascia FuncCall opaca per consumi successivi.
    if (op == BuiltinOp::Hypergeometric0F1 && args.size() == 2U) {
        if (is_zero_expr(args[1])) return ok(make_integer(arena_, BigInt(1)));
    }
    auto symbols_equal = [](ExprPtr x, ExprPtr y) -> bool {
        if (x == y) return true;
        if (!x || !y) return false;
        const auto* sx = expr_cast<Symbol>(x);
        const auto* sy = expr_cast<Symbol>(y);
        if (sx && sy) return sx->name == sy->name;
        const auto* ix = expr_cast<IntegerLit>(x);
        const auto* iy = expr_cast<IntegerLit>(y);
        if (ix && iy) return ix->value == iy->value;
        const auto* rx = expr_cast<RationalLit>(x);
        const auto* ry = expr_cast<RationalLit>(y);
        if (rx && ry) return rx->numerator == ry->numerator
                            && rx->denominator == ry->denominator;
        return false;
    };
    if (op == BuiltinOp::Hypergeometric1F1 && args.size() == 3U) {
        if (is_zero_expr(args[2])) return ok(make_integer(arena_, BigInt(1)));
        if (symbols_equal(args[0], args[1])) {
            return simplify_expr(arena_.make<FuncCall>(BuiltinOp::Exp,
                std::vector<ExprPtr>{args[2]}));
        }
    }
    if (op == BuiltinOp::Hypergeometric2F1 && args.size() == 4U) {
        if (is_zero_expr(args[3])) return ok(make_integer(arena_, BigInt(1)));
        if (symbols_equal(args[1], args[2])) {
            ExprPtr one_minus_z = arena_.make<Binary>(BinaryOp::Sub,
                make_integer(arena_, BigInt(1)), args[3]);
            ExprPtr neg_a = arena_.make<Unary>(UnaryOp::Neg, args[0]);
            return simplify_expr(arena_.make<Binary>(BinaryOp::Pow,
                one_minus_z, neg_a));
        }
        if (symbols_equal(args[0], args[2])) {
            ExprPtr one_minus_z = arena_.make<Binary>(BinaryOp::Sub,
                make_integer(arena_, BigInt(1)), args[3]);
            ExprPtr neg_b = arena_.make<Unary>(UnaryOp::Neg, args[1]);
            return simplify_expr(arena_.make<Binary>(BinaryOp::Pow,
                one_minus_z, neg_b));
        }
    }

    // Elliptic integrals (F5.9 / Task #20).
    auto half = [&]() {
        return arena_.make<Binary>(BinaryOp::Div,
            arena_.make<Constant>(MathConstant::Pi),
            make_integer(arena_, BigInt(2)));
    };
    if ((op == BuiltinOp::EllipticK || op == BuiltinOp::EllipticE)
        && args.size() == 1U) {
        if (is_zero_expr(args[0])) return ok(half());
        if (op == BuiltinOp::EllipticE) {
            if (const auto* il = expr_cast<IntegerLit>(args[0]);
                il && il->value == BigInt(1)) {
                return ok(make_integer(arena_, BigInt(1)));
            }
        }
    }
    if (op == BuiltinOp::EllipticF && args.size() == 2U) {
        if (is_zero_expr(args[1])) return ok(args[0]);
    }
    if (op == BuiltinOp::EllipticPi && args.size() == 2U) {
        if (is_zero_expr(args[0])) {
            return simplify_expr(arena_.make<FuncCall>(BuiltinOp::EllipticK,
                std::vector<ExprPtr>{args[1]}));
        }
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

}  // namespace cas::symbolic::detail
