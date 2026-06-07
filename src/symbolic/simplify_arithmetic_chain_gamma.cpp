// L3-04 Gamma reflection identity applied during Product simplification.
//
//   Γ(z)·Γ(1−z) = π / sin(πz)        [m = 1 case]
//   Γ(z)·Γ(−z)  = −π / (z · sin(πz)) [m = 0 case]
//
// Fires when two `Γ` factors of `symbolic` (the "base, exponent" pairs that
// `simplify_product_factors` collects in its non-numeric bucket) sum to an
// integer m ∈ {0, 1}.  The pair is rewritten in place as
// π · sin(πz)^{−1}  (and an additional z^{−1} factor with sign flip on the
// numeric coefficient when m = 0).
//
// Reflection is iterated to a fixed point so chains such as
// Γ(z)·Γ(1−z)·Γ(z+1)·Γ(−z) collapse fully.  After the pass, the caller
// re-runs `merge_symbolic_factors` to canonicalise duplicate base entries.

#include "simplify_arithmetic_chain_impl.hpp"

namespace cas::symbolic::detail {

Result<void> Simplifier::apply_gamma_reflection_pairs(
    std::vector<std::pair<ExprPtr, BigInt>>& symbolic,
    ComplexRational& coefficient)
{
    auto gamma_arg = [](const std::pair<ExprPtr, BigInt>& p) -> ExprPtr {
        if (p.second != BigInt(1)) return nullptr;
        const auto* fc = expr_cast<FuncCall>(p.first);
        if (!fc || fc->func_id != BuiltinOp::Gamma || fc->args.size() != 1U)
            return nullptr;
        return fc->args[0];
    };
    bool reflected_any = false;
    bool keep_scanning = true;
    while (keep_scanning) {
        keep_scanning = false;
        std::vector<std::size_t> gamma_idx;
        for (std::size_t i = 0; i < symbolic.size(); ++i)
            if (gamma_arg(symbolic[i]) != nullptr) gamma_idx.push_back(i);
        for (std::size_t a = 0; a < gamma_idx.size() && !keep_scanning; ++a) {
            for (std::size_t b = a + 1; b < gamma_idx.size() && !keep_scanning; ++b) {
                ExprPtr za = gamma_arg(symbolic[gamma_idx[a]]);
                ExprPtr zb = gamma_arg(symbolic[gamma_idx[b]]);
                if (!za || !zb) continue;
                auto sum_simp = simplify_expr(
                    arena_.make<Sum>(std::vector<ExprPtr>{za, zb}));
                if (sum_simp.is_error()) return fail<void>(sum_simp.error());
                const auto* m_lit = expr_cast<IntegerLit>(sum_simp.value());
                if (!m_lit) continue;
                const bool m_is_one  = (m_lit->value == BigInt(1));
                const bool m_is_zero = m_lit->value.is_zero();
                if (!m_is_one && !m_is_zero) continue;
                ExprPtr pi_const = arena_.make<Constant>(MathConstant::Pi);
                auto pi_z_simp = simplify_expr(
                    arena_.make<Product>(std::vector<ExprPtr>{pi_const, za}));
                if (pi_z_simp.is_error()) return fail<void>(pi_z_simp.error());
                auto sin_simp = simplify_expr(
                    arena_.make<FuncCall>(BuiltinOp::Sin,
                        std::vector<ExprPtr>{pi_z_simp.value()}));
                if (sin_simp.is_error()) return fail<void>(sin_simp.error());
                std::size_t ia = gamma_idx[a], ib = gamma_idx[b];
                if (ia > ib) std::swap(ia, ib);
                symbolic.erase(symbolic.begin() + ib);
                symbolic.erase(symbolic.begin() + ia);
                symbolic.push_back({pi_const,          BigInt(1)});
                symbolic.push_back({sin_simp.value(),  BigInt(-1)});
                if (m_is_zero) {
                    coefficient = coefficient * ComplexRational(Rational(BigInt(-1)));
                    symbolic.push_back({za, BigInt(-1)});
                }
                reflected_any = true;
                keep_scanning = true;
            }
        }
    }
    if (reflected_any) merge_symbolic_factors(symbolic);
    return ok();
}

}  // namespace cas::symbolic::detail
