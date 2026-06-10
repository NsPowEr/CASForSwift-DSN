// F7.5.F1 — Extended-real helpers
//
// Promotes ±∞, ComplexInfinity, Indeterminate to first-class constants.
// Provides predicates that accept both the new canonical form
// (Constant(NegInfinity)) and the legacy form (Unary(Neg,
// Constant(Infinity))) during transition, so existing call sites stay
// correct while migrators move them to the canonical form.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Extended_Real_AST.md
// Ledger: HC-F70-A43-EXTENDED-REAL (closing on this migration).

#pragma once

#include "cas/ast.hpp"

namespace cas {

[[nodiscard]] inline bool is_pos_infinity(ExprPtr e) noexcept {
    if (!e) return false;
    if (const auto* c = expr_cast<Constant>(e))
        return c->value == MathConstant::Infinity;
    return false;
}

[[nodiscard]] inline bool is_neg_infinity(ExprPtr e) noexcept {
    if (!e) return false;
    if (const auto* c = expr_cast<Constant>(e))
        return c->value == MathConstant::NegInfinity;
    if (const auto* u = expr_cast<Unary>(e)) {
        if (u->op == UnaryOp::Neg) {
            if (const auto* inner = expr_cast<Constant>(u->operand)) {
                return inner->value == MathConstant::Infinity;
            }
        }
    }
    return false;
}

[[nodiscard]] inline bool is_complex_infinity(ExprPtr e) noexcept {
    if (!e) return false;
    if (const auto* c = expr_cast<Constant>(e))
        return c->value == MathConstant::ComplexInfinity;
    return false;
}

[[nodiscard]] inline bool is_indeterminate(ExprPtr e) noexcept {
    if (!e) return false;
    if (const auto* c = expr_cast<Constant>(e))
        return c->value == MathConstant::Indeterminate;
    return false;
}

[[nodiscard]] inline bool is_signed_infinity(ExprPtr e) noexcept {
    return is_pos_infinity(e) || is_neg_infinity(e);
}

[[nodiscard]] inline bool is_any_infinity(ExprPtr e) noexcept {
    return is_signed_infinity(e) || is_complex_infinity(e);
}

[[nodiscard]] inline ExprPtr make_pos_infinity(AstArena& arena) {
    return arena.make<Constant>(MathConstant::Infinity);
}

[[nodiscard]] inline ExprPtr make_neg_infinity(AstArena& arena) {
    return arena.make<Constant>(MathConstant::NegInfinity);
}

[[nodiscard]] inline ExprPtr make_complex_infinity(AstArena& arena) {
    return arena.make<Constant>(MathConstant::ComplexInfinity);
}

[[nodiscard]] inline ExprPtr make_indeterminate(AstArena& arena) {
    return arena.make<Constant>(MathConstant::Indeterminate);
}

}  // namespace cas
