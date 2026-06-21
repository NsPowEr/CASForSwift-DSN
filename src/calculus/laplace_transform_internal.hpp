#ifndef CAS_CALCULUS_LAPLACE_TRANSFORM_INTERNAL_HPP
#define CAS_CALCULUS_LAPLACE_TRANSFORM_INTERNAL_HPP

#include "calculus_internal.hpp"

#include <utility>

// Helpers shared (T-050 anti-monolith split) between the forward transform in
// laplace_transform.cpp and the inverse transform in
// laplace_transform_inverse.cpp.
namespace cas::calculus::laplace_detail {

[[nodiscard]] inline ExprPtr make_int(symbolic::CASContext& ctx, long long v) {
    return ctx.arena().make<IntegerLit>(BigInt(v));
}

// Test e ≡ Symbol(name).
[[nodiscard]] inline bool is_sym(ExprPtr e, const Symbol& sym) {
    const auto* s = expr_cast<Symbol>(e);
    return s && s->name == sym.name;
}

// Factorial as expr (n!).
[[nodiscard]] inline ExprPtr factorial_expr(symbolic::CASContext& ctx, long long n) {
    BigInt acc(1);
    for (long long i = 2; i <= n; ++i) acc = acc * BigInt(i);
    return ctx.arena().make<IntegerLit>(std::move(acc));
}

// Does expr depend on var (symbolic occurrence)?
[[nodiscard]] inline bool depends_on_t(ExprPtr e, const Symbol& v) {
    if (!e) return false;
    if (auto* s = expr_cast<Symbol>(e)) return s->name == v.name;
    if (auto* un = expr_cast<Unary>(e)) return depends_on_t(un->operand, v);
    if (auto* bin = expr_cast<Binary>(e)) return depends_on_t(bin->left, v) || depends_on_t(bin->right, v);
    if (auto* sum = expr_cast<Sum>(e)) { for (auto t : sum->terms) if (depends_on_t(t, v)) return true; return false; }
    if (auto* prod = expr_cast<Product>(e)) { for (auto f : prod->factors) if (depends_on_t(f, v)) return true; return false; }
    if (auto* fc = expr_cast<FuncCall>(e)) { for (auto a : fc->args) if (depends_on_t(a, v)) return true; return false; }
    return false;
}

}  // namespace cas::calculus::laplace_detail

#endif  // CAS_CALCULUS_LAPLACE_TRANSFORM_INTERNAL_HPP
