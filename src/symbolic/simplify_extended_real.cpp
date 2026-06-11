// F7.5.F1 Phase 2 — Extended-real arithmetic propagation.
//
// Implements the semantic rules from
// .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Extended_Real_AST.md
// for Sum, Product, and Pow when one or more operands are extended-real
// (PosInfinity, NegInfinity, ComplexInfinity, Indeterminate).
//
// Closes the residual scope of HC-F70-A43-EXTENDED-REAL by giving the
// simplifier honest, total semantics on the extended real line instead
// of either propagating an `Undefined` error or pattern-matching the
// few cases that previously existed in scattered code.
//
// The helpers are conservative: when an extended-real operand is present
// but the remaining finite operands have unknown sign (free symbols), the
// helpers return std::nullopt so the caller may fall back to its normal
// rewriting path. This avoids fabricating a sign for symbolic factors
// whose sign is not derivable from the AST alone.

#include "simplify_impl.hpp"
#include "cas/extended_real.hpp"
#include <vector>

namespace cas::symbolic::detail {

namespace {

enum class ExtClass {
    Finite,
    PosInf,
    NegInf,
    ComplexInf,
    Indet,
};

[[nodiscard]] ExtClass classify(ExprPtr e) noexcept {
    if (is_pos_infinity(e)) return ExtClass::PosInf;
    if (is_neg_infinity(e)) return ExtClass::NegInf;
    if (is_complex_infinity(e)) return ExtClass::ComplexInf;
    if (is_indeterminate(e)) return ExtClass::Indet;
    return ExtClass::Finite;
}

// True iff e is a literal numeric whose value is strictly positive.
[[nodiscard]] bool is_literal_positive(ExprPtr e) noexcept {
    if (const auto* il = expr_cast<IntegerLit>(e))
        return !il->value.is_negative() && !il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e))
        return !rl->numerator.is_negative() && !rl->numerator.is_zero();
    return false;
}

// True iff e is a literal numeric whose value is strictly negative.
[[nodiscard]] bool is_literal_negative(ExprPtr e) noexcept {
    if (const auto* il = expr_cast<IntegerLit>(e))
        return il->value.is_negative();
    if (const auto* rl = expr_cast<RationalLit>(e))
        return rl->numerator.is_negative();
    return false;
}

[[nodiscard]] bool is_literal_zero(ExprPtr e) noexcept {
    return is_zero_expr(e);
}

[[nodiscard]] bool is_literal_one(ExprPtr e) noexcept {
    return is_one_expr(e);
}

// Sign tracker for product of extended-reals against finite factors with
// known sign. Returns nullopt when any finite factor has unknown sign.
struct SignTracker {
    int flips = 0;            // count of factors contributing a sign flip
    bool any_unknown = false; // true if some factor has indeterminate sign
};

[[nodiscard]] SignTracker tally_finite_signs(const std::vector<ExprPtr>& factors) {
    SignTracker t;
    for (ExprPtr f : factors) {
        if (classify(f) != ExtClass::Finite) continue;
        if (is_literal_positive(f)) continue;
        if (is_literal_negative(f)) { ++t.flips; continue; }
        if (is_literal_zero(f)) continue;  // handled separately
        // Symbolic / not literal — sign unknown.
        t.any_unknown = true;
    }
    return t;
}

}  // namespace

std::optional<ExprPtr> try_simplify_sum_extended_real(
    const std::vector<ExprPtr>& terms,
    AstArena& arena)
{
    bool has_pos = false;
    bool has_neg = false;
    bool has_complex = false;
    bool has_new_enum = false;  // NegInf / ComplexInf / Indet
    bool any_extended = false;

    for (ExprPtr t : terms) {
        switch (classify(t)) {
        case ExtClass::Indet:
            // Indeterminate is absorbing under +.
            return make_indeterminate(arena);
        case ExtClass::PosInf:
            has_pos = true; any_extended = true; break;
        case ExtClass::NegInf:
            has_neg = true; any_extended = true; has_new_enum = true; break;
        case ExtClass::ComplexInf:
            has_complex = true; any_extended = true; has_new_enum = true; break;
        case ExtClass::Finite:
            break;
        }
    }

    if (!any_extended) return std::nullopt;

    // Conservative: when only legacy MathConstant::Infinity appears (no
    // NegInfinity / ComplexInfinity / Indeterminate operand), defer to the
    // existing code paths that depend on the pre-F7.5.F1 semantics
    // (callers such as the limit engine inspect symbolic +∞ sums directly).
    if (!has_new_enum) return std::nullopt;

    // ComplexInfinity has no sign; adding a signed infinity to it leaves
    // the result ambiguous (the limit depends on the direction of approach).
    if (has_complex && (has_pos || has_neg))
        return make_indeterminate(arena);

    if (has_pos && has_neg)
        return make_indeterminate(arena);

    if (has_complex)
        return make_complex_infinity(arena);

    if (has_pos)
        return make_pos_infinity(arena);

    // has_neg
    return make_neg_infinity(arena);
}

std::optional<ExprPtr> try_simplify_product_extended_real(
    const std::vector<ExprPtr>& factors,
    AstArena& arena)
{
    int neg_inf_count = 0;
    bool has_complex_inf = false;
    bool has_literal_zero_factor = false;
    bool has_new_enum = false;  // NegInf / ComplexInf / Indet
    bool any_extended = false;

    for (ExprPtr f : factors) {
        switch (classify(f)) {
        case ExtClass::Indet:
            // Indeterminate is absorbing under ·.
            return make_indeterminate(arena);
        case ExtClass::PosInf:
            any_extended = true; break;
        case ExtClass::NegInf:
            ++neg_inf_count; any_extended = true; has_new_enum = true; break;
        case ExtClass::ComplexInf:
            has_complex_inf = true; any_extended = true; has_new_enum = true; break;
        case ExtClass::Finite:
            if (is_literal_zero(f)) has_literal_zero_factor = true;
            break;
        }
    }

    if (!any_extended) return std::nullopt;
    // Conservative: pure legacy MathConstant::Infinity (no new enum) defers
    // to the existing flow so the limit engine's `Undefined` fallback path
    // for `0 * Infinity` remains observable in callers that depend on it.
    if (!has_new_enum) return std::nullopt;

    // 0 * (anything infinite) = Indeterminate (per IEEE / projective ext-real
    // and the standard limit calculus convention).
    if (has_literal_zero_factor) return make_indeterminate(arena);

    if (has_complex_inf) return make_complex_infinity(arena);

    // Sign tracking for products of signed infinities and finite known-sign
    // literals. Symbolic factors with unknown sign force the result to
    // ComplexInfinity (unsigned), since the directional sign cannot be
    // determined symbolically.
    const SignTracker signs = tally_finite_signs(factors);
    if (signs.any_unknown) return make_complex_infinity(arena);

    const int total_flips = neg_inf_count + signs.flips;
    const bool negative = (total_flips % 2) == 1;
    return negative ? make_neg_infinity(arena) : make_pos_infinity(arena);
}

std::optional<ExprPtr> try_simplify_pow_extended_real(
    ExprPtr base,
    ExprPtr exponent,
    AstArena& arena)
{
    const ExtClass cb = classify(base);
    const ExtClass ce = classify(exponent);

    if (cb == ExtClass::Indet || ce == ExtClass::Indet)
        return make_indeterminate(arena);

    if (cb == ExtClass::Finite && ce == ExtClass::Finite)
        return std::nullopt;

    // Conservative: pure legacy MathConstant::Infinity (no NegInf / ComplexInf
    // operand) defers to the existing simplify_power code path that already
    // handles `+∞^n` with integer exponent (n>0 → +∞, n<0 → 0). This avoids
    // disturbing observed limit-engine semantics on legacy-only inputs.
    const bool has_new_enum =
        (cb == ExtClass::NegInf || cb == ExtClass::ComplexInf ||
         ce == ExtClass::NegInf || ce == ExtClass::ComplexInf);
    if (!has_new_enum) return std::nullopt;

    // 0^0 — Indeterminate. (Captured here for completeness when the caller
    // routes 0^0 through this helper; simplify_power handles the finite case
    // separately.)
    if (cb == ExtClass::Finite && is_literal_zero(base) &&
        ce == ExtClass::Finite && is_literal_zero(exponent))
        return make_indeterminate(arena);

    // 1^(±∞) and 1^(ComplexInfinity) = Indeterminate.
    if (cb == ExtClass::Finite && is_literal_one(base)) {
        if (ce == ExtClass::PosInf || ce == ExtClass::NegInf ||
            ce == ExtClass::ComplexInf)
            return make_indeterminate(arena);
    }

    // ComplexInfinity in the exponent with any base other than 1 is
    // direction-dependent; treat as Indeterminate.
    if (ce == ExtClass::ComplexInf) return make_indeterminate(arena);

    // (±∞)^0 = Indeterminate.
    if ((cb == ExtClass::PosInf || cb == ExtClass::NegInf ||
         cb == ExtClass::ComplexInf) &&
        ce == ExtClass::Finite && is_literal_zero(exponent))
        return make_indeterminate(arena);

    // PosInfinity in the base.
    if (cb == ExtClass::PosInf) {
        if (ce == ExtClass::PosInf) return make_pos_infinity(arena);
        if (ce == ExtClass::NegInf)
            return make_integer(arena, BigInt(0));
        // Finite, nonzero exponent.
        if (is_literal_positive(exponent)) return make_pos_infinity(arena);
        if (is_literal_negative(exponent)) return make_integer(arena, BigInt(0));
        // Symbolic exponent of unknown sign — leave to normal flow.
        return std::nullopt;
    }

    // NegInfinity in the base.
    if (cb == ExtClass::NegInf) {
        if (ce == ExtClass::PosInf) return make_complex_infinity(arena);
        if (ce == ExtClass::NegInf) return make_integer(arena, BigInt(0));
        if (ce == ExtClass::Finite) {
            // (-∞)^n with n a positive integer literal: parity decides sign.
            if (const auto* il = expr_cast<IntegerLit>(exponent)) {
                if (!il->value.is_negative() && !il->value.is_zero()) {
                    const bool odd = (il->value % BigInt(2)) == BigInt(1);
                    return odd ? make_neg_infinity(arena)
                               : make_pos_infinity(arena);
                }
                if (il->value.is_negative()) return make_integer(arena, BigInt(0));
            }
            // (-∞)^positive non-integer rational → ComplexInfinity (branch).
            if (is_literal_positive(exponent)) return make_complex_infinity(arena);
            if (is_literal_negative(exponent)) return make_integer(arena, BigInt(0));
            return std::nullopt;
        }
    }

    // ComplexInfinity in the base.
    if (cb == ExtClass::ComplexInf) {
        // ComplexInfinity^positive = ComplexInfinity.
        // ComplexInfinity^negative = 0.
        if (ce == ExtClass::Finite) {
            if (is_literal_positive(exponent)) return make_complex_infinity(arena);
            if (is_literal_negative(exponent)) return make_integer(arena, BigInt(0));
            return std::nullopt;
        }
        if (ce == ExtClass::PosInf) return make_complex_infinity(arena);
        if (ce == ExtClass::NegInf) return make_integer(arena, BigInt(0));
    }

    return std::nullopt;
}

}  // namespace cas::symbolic::detail
