#pragma once
// Internal helper shared (T-048 anti-monolith split) between
// summation_zeilberger_helpers.cpp (gamma expansion / cancellation) and
// summation_zeilberger_shift_ratio.cpp (Pochhammer shift-ratio).  Not part of
// the public CAS API.

#include "cas/ast.hpp"

#include <cstdlib>
#include <vector>

namespace cas::symbolic::zeilberger_detail {

// Flatten a Mul/Div/Pow tree into numerator/denominator factor lists. Integer
// powers (|n| ≤ 32) are expanded into repeated factors; Div swaps the roles of
// its right operand.
inline void collect_product_factors(
    ExprPtr e, std::vector<ExprPtr>& num_out, std::vector<ExprPtr>& den_out) {
    if (!e) return;
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr f : prod->factors) collect_product_factors(f, num_out, den_out);
        return;
    }
    if (const auto* bin = expr_cast<Binary>(e)) {
        if (bin->op == BinaryOp::Mul) {
            collect_product_factors(bin->left,  num_out, den_out);
            collect_product_factors(bin->right, num_out, den_out);
            return;
        }
        if (bin->op == BinaryOp::Div) {
            collect_product_factors(bin->left,  num_out, den_out);
            // Swap roles for the denominator side.
            collect_product_factors(bin->right, den_out, num_out);
            return;
        }
        if (bin->op == BinaryOp::Pow) {
            // Handle Pow(base, n) with integer n: positive → base in num; negative → den.
            long long n = 0;
            bool got_int = false;
            if (const auto* il = expr_cast<IntegerLit>(bin->right);
                il && il->value.bit_length() <= 16U) {
                n = static_cast<long long>(il->value.abs().to_u64()) *
                    (il->value.is_negative() ? -1 : 1);
                got_int = true;
            } else if (const auto* un = expr_cast<Unary>(bin->right);
                       un && un->op == UnaryOp::Neg) {
                if (const auto* il = expr_cast<IntegerLit>(un->operand);
                    il && il->value.bit_length() <= 16U) {
                    n = -static_cast<long long>(il->value.abs().to_u64());
                    got_int = true;
                }
            }
            if (got_int && n != 0 && std::abs(n) <= 32) {
                long long count = std::abs(n);
                for (long long i = 0; i < count; ++i)
                    collect_product_factors(bin->left,
                        (n > 0) ? num_out : den_out,
                        (n > 0) ? den_out : num_out);
                return;
            }
        }
    }
    num_out.push_back(e);
}

}  // namespace cas::symbolic::zeilberger_detail
