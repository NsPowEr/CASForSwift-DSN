// Risch IBP exp-fold logic for Product simplification.
//
// When fold_exp_products is active (e.g. inside integrate_by_parts), this
// merges exponential products: exp(a) · exp(b) → exp(a+b).
//
// This is non-canonical globally but required in differential fields during
// IBP simplification and differentiation check to prevent infinite loop.

#include "simplify_arithmetic_chain_impl.hpp"

namespace cas::symbolic::detail {

[[nodiscard]] static ExprPtr multiply_by_bigint(ExprPtr expr, const BigInt& val, AstArena& arena) {
    if (val == BigInt(1)) return expr;
    if (val.is_zero()) return make_integer(arena, BigInt(0));
    if (val.is_negative()) {
        BigInt abs_val = BigInt(0) - val;
        ExprPtr abs_coeff = make_integer(arena, abs_val);
        ExprPtr product = (abs_val == BigInt(1)) ? expr : arena.make<Binary>(BinaryOp::Mul, abs_coeff, expr);
        return arena.make<Unary>(UnaryOp::Neg, product);
    } else {
        ExprPtr coeff = make_integer(arena, val);
        return arena.make<Binary>(BinaryOp::Mul, coeff, expr);
    }
}

[[nodiscard]] static std::optional<ExprPtr> try_get_exp_arg(ExprPtr base, const BigInt& exponent, AstArena& arena) {
    if (const auto* fc = expr_cast<FuncCall>(base)) {
        if (fc->func_id == BuiltinOp::Exp && fc->args.size() == 1U) {
            return multiply_by_bigint(fc->args[0], exponent, arena);
        }
    }
    if (const auto* bin = expr_cast<Binary>(base)) {
        if (bin->op == BinaryOp::Pow && is_constant_expr(bin->left, MathConstant::E)) {
            return multiply_by_bigint(bin->right, exponent, arena);
        }
    }
    if (is_constant_expr(base, MathConstant::E)) {
        return make_integer(arena, exponent);
    }
    return std::nullopt;
}

Result<void> Simplifier::fold_exponential_products(
    std::vector<std::pair<ExprPtr, BigInt>>& symbolic,
    ComplexRational& coefficient)
{
    std::vector<ExprPtr> exp_args;
    std::vector<std::size_t> indices_to_remove;
    for (std::size_t i = 0; i < symbolic.size(); ++i) {
        if (auto arg = try_get_exp_arg(symbolic[i].first, symbolic[i].second, arena_)) {
            exp_args.push_back(*arg);
            indices_to_remove.push_back(i);
        }
    }
    if (exp_args.size() > 1) {
        for (auto it = indices_to_remove.rbegin(); it != indices_to_remove.rend(); ++it) {
            symbolic.erase(symbolic.begin() + static_cast<std::ptrdiff_t>(*it));
        }
        ExprPtr sum_exp = arena_.make<Sum>(std::move(exp_args));
        auto sum_exp_s = simplify_expr(sum_exp);
        ExprPtr simplified_sum = sum_exp_s.is_ok() ? sum_exp_s.value() : sum_exp;
        ExprPtr new_exp = arena_.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{simplified_sum});
        auto new_exp_s = simplify_expr(new_exp);
        ExprPtr final_exp = new_exp_s.is_ok() ? new_exp_s.value() : new_exp;

        auto insert_factor = [&](ExprPtr f) {
            if (is_zero_expr(f)) return;
            LiteralComplex comp;
            auto exact = try_get_exact_complex(f, comp);
            if (exact.is_ok() && exact.value()) {
                coefficient = coefficient * comp.value;
                return;
            }
            if (const auto* prod = expr_cast<Product>(f)) {
                for (ExprPtr factor : prod->factors) {
                    symbolic.push_back({factor, BigInt(1)});
                }
            } else {
                symbolic.push_back({f, BigInt(1)});
            }
        };
        insert_factor(final_exp);
        merge_symbolic_factors(symbolic);
    }
    return ok();
}

} // namespace cas::symbolic::detail
