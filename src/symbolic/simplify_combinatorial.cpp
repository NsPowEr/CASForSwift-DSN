// F7.5 follow-up — exact-value identities for the combinatorial /
// erfc builtins. Closes the bulk of the special_fn corpus regressions
// where `factorial(n)`, `binomial(n, k)`, and `erfc(0)` left the
// canonical form unevaluated.
//
// Rules implemented (exact; no floating point, no closed pattern
// table beyond well-defined mathematical identities):
//
//   factorial(0)              = 1
//   factorial(n)               for IntegerLit n > 0 -> 1 * 2 * ... * n
//   factorial(n)               for IntegerLit n < 0 -> ComplexInfinity
//                              (pole of Gamma(n+1); matches Γ semantics)
//
//   binomial(n, 0)            = 1
//   binomial(n, n)            = 1                       (n IntegerLit)
//   binomial(n, 1)            = n                       (any n)
//   binomial(n, k)             for IntegerLit n >= k >= 0
//                              -> n! / (k! · (n-k)!)
//   binomial(n, 2)            = n * (n - 1) / 2         (symbolic n)
//
//   erfc(0)                   = 1
//   erfc(+oo)                 = 0
//   erfc(-oo)                 = 2

#include "simplify_impl.hpp"
#include "cas/extended_real.hpp"

namespace cas::symbolic::detail {

namespace {

[[nodiscard]] ExprPtr factorial_bigint(AstArena& arena, const BigInt& n) {
    BigInt acc(1);
    for (BigInt i(2); !(n < i); i = i + BigInt(1)) acc = acc * i;
    return arena.make<IntegerLit>(acc);
}

}  // namespace

Result<ExprPtr> Simplifier::simplify_funcall_combinatorial(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr /*target_before*/)
{
    if (op == BuiltinOp::Factorial && args.size() == 1U) {
        ExprPtr a = args.front();
        if (const auto* lit = expr_cast<IntegerLit>(a)) {
            if (lit->value.is_zero()) return ok(make_integer(arena_, BigInt(1)));
            if (!lit->value.is_negative()) {
                return ok(factorial_bigint(arena_, lit->value));
            }
            // Negative integer: factorial(-k) = Gamma(-k+1) which has
            // a simple pole. Match the Gamma semantics from F7.5.E1.
            return ok(arena_.make<Constant>(MathConstant::ComplexInfinity));
        }
    }

    if (op == BuiltinOp::Binomial && args.size() == 2U) {
        ExprPtr n = args[0];
        ExprPtr k = args[1];
        // binomial(n, 0) = 1 for any n.
        if (is_zero_expr(k)) return ok(make_integer(arena_, BigInt(1)));
        // binomial(n, 1) = n for any n.
        if (const auto* k_lit = expr_cast<IntegerLit>(k);
            k_lit != nullptr && k_lit->value == BigInt(1)) {
            return ok(n);
        }
        // binomial(n, n) = 1 when both are the same IntegerLit.
        if (const auto* n_lit = expr_cast<IntegerLit>(n);
            n_lit != nullptr) {
            if (const auto* k_lit = expr_cast<IntegerLit>(k);
                k_lit != nullptr) {
                if (k_lit->value.is_negative() ||
                    n_lit->value < k_lit->value) {
                    return ok(make_integer(arena_, BigInt(0)));
                }
                if (n_lit->value == k_lit->value) {
                    return ok(make_integer(arena_, BigInt(1)));
                }
                // n! / (k! (n-k)!)
                ExprPtr nf = factorial_bigint(arena_, n_lit->value);
                ExprPtr kf = factorial_bigint(arena_, k_lit->value);
                BigInt nmk = n_lit->value - k_lit->value;
                ExprPtr nmkf = factorial_bigint(arena_, nmk);
                ExprPtr denom_raw = arena_.make<Product>(
                    std::vector<ExprPtr>{kf, nmkf});
                auto denom_s = simplify_expr(denom_raw);
                ExprPtr denom = denom_s.is_ok() ? denom_s.value() : denom_raw;
                ExprPtr quot = arena_.make<Binary>(BinaryOp::Div, nf, denom);
                return simplify_expr(quot);
            }
        }
        // binomial(n, 2) = n*(n-1)/2 with symbolic n.
        if (const auto* k_lit = expr_cast<IntegerLit>(k);
            k_lit != nullptr && k_lit->value == BigInt(2)) {
            ExprPtr n_minus_one = arena_.make<Binary>(
                BinaryOp::Sub, n, make_integer(arena_, BigInt(1)));
            ExprPtr num = arena_.make<Product>(
                std::vector<ExprPtr>{n, n_minus_one});
            ExprPtr quot = arena_.make<Binary>(
                BinaryOp::Div, num, make_integer(arena_, BigInt(2)));
            return simplify_expr(quot);
        }
    }

    if (op == BuiltinOp::Erfc && args.size() == 1U) {
        ExprPtr a = args.front();
        if (is_zero_expr(a)) return ok(make_integer(arena_, BigInt(1)));
        if (is_pos_infinity(a)) return ok(make_integer(arena_, BigInt(0)));
        if (is_neg_infinity(a)) return ok(make_integer(arena_, BigInt(2)));
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

}  // namespace cas::symbolic::detail
