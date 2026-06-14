#pragma once
// Internal types and functions shared between ode_kovacic_case1.cpp and ode_kovacic.cpp.

#include "ode_kovacic.hpp"
#include "cas/algebra.hpp"
#include "cas/rational.hpp"
#include "calculus_internal.hpp"

namespace cas::calculus {
namespace kovacic_impl {

// Struct for one partial-fraction pole of the invariant r.
struct PFPole {
    ExprPtr  pole;   // location c
    unsigned power;  // order 1 or 2
    ExprPtr  coeff;  // numerator A in  A/(x−c)^power
};

// Pair of candidate rational 1-forms for Case 1.
struct OmegaPair {
    ExprPtr plus;   // ω₊ = Σ k₊ᵢ/(x−cᵢ) + poly
    ExprPtr minus;  // ω₋ = Σ k₋ᵢ/(x−cᵢ) + poly
};

// Compute the Kovacic invariant  r = p²/4 − p'/2 − q
// where  p = a1/a2,  q = a0/a2.
[[nodiscard]] Result<ExprPtr> compute_r(
    ExprPtr a2, ExprPtr a1, ExprPtr a0,
    const Symbol& x, symbolic::CASContext& ctx);

// Attempt Kovacic Case 1: return {ω₊, ω₋} s.t. ωᵢ' + ωᵢ² = r.
// Fails with Unimplemented if Case 2/3 is needed.
[[nodiscard]] Result<OmegaPair> case1_omega(
    ExprPtr r, const Symbol& x, symbolic::CASContext& ctx);

// F4.K2 — Kovacic Case 2 (dihedral D∞ subgroup) — SCAFFOLD ONLY.
// Algorithm body NOT YET implemented; returns Unimplemented with explicit
// diagnostic referencing HC-KV-03 ledger and Kovacic_Case2.md spec.
//
// Entry point reserved so solve_ode_kovacic can route through Case 2
// before final dispatch.  No algorithmic computation performed; no
// silent wrong-answer possible.  See:
//   - Kovacic 1986 §3.2
//   - .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Kovacic_Case2.md
//   - HARDCODE_LEDGER.md HC-KV-03
[[nodiscard]] Result<OmegaPair> case2_omega(
    ExprPtr r, const Symbol& x, symbolic::CASContext& ctx);

// Arithmetic short-hands used by both translation units.
[[nodiscard]] inline ExprPtr kv_int(AstArena& a, long long v) {
    return a.make<IntegerLit>(BigInt(v));
}
[[nodiscard]] inline ExprPtr kv_add(AstArena& a, ExprPtr l, ExprPtr r) {
    return a.make<Binary>(BinaryOp::Add, l, r);
}
[[nodiscard]] inline ExprPtr kv_sub(AstArena& a, ExprPtr l, ExprPtr r) {
    return a.make<Binary>(BinaryOp::Sub, l, r);
}
[[nodiscard]] inline ExprPtr kv_mul(AstArena& a, ExprPtr l, ExprPtr r) {
    return a.make<Binary>(BinaryOp::Mul, l, r);
}
[[nodiscard]] inline ExprPtr kv_div(AstArena& a, ExprPtr l, ExprPtr r) {
    return a.make<Binary>(BinaryOp::Div, l, r);
}
[[nodiscard]] inline ExprPtr kv_neg(AstArena& a, ExprPtr e) {
    return a.make<Unary>(UnaryOp::Neg, e);
}
[[nodiscard]] inline ExprPtr kv_exp(AstArena& a, ExprPtr e) {
    return a.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{e});
}
[[nodiscard]] inline CASError kv_unimpl(std::string msg) {
    return CASError{.kind = CASErrorKind::Unimplemented, .message = std::move(msg), .hint = std::nullopt};
}
[[nodiscard]] inline bool kv_is_zero(ExprPtr e, symbolic::CASContext& ctx) {
    auto r = ctx.simplify(e);
    if (r.is_error()) return false;
    if (auto* il = expr_cast<IntegerLit>(r.value())) return il->value.is_zero();
    if (auto* rl = expr_cast<RationalLit>(r.value())) return rl->numerator.is_zero();
    return false;
}

// Helper functions for Laurent expansion & arithmetic
[[nodiscard]] std::optional<BigInt> bigint_isqrt(const BigInt& n);
[[nodiscard]] std::optional<Rational> as_rational(ExprPtr e);
[[nodiscard]] std::optional<Rational> rational_sqrt(const Rational& r);
[[nodiscard]] Result<ExprPtr> reverse_polynomial(
    ExprPtr poly, const Symbol& x, const Symbol& y, AstArena& a, symbolic::CASContext& ctx);
[[nodiscard]] std::optional<std::vector<Rational>> compute_taylor_rational(
    ExprPtr num, ExprPtr den_other, const Symbol& x, ExprPtr c,
    unsigned terms_needed, symbolic::CASContext& ctx);
[[nodiscard]] std::optional<std::vector<Rational>> compute_laurent_sqrt(
    const std::vector<Rational>& u, unsigned terms_needed);

} // namespace kovacic_impl
} // namespace cas::calculus
