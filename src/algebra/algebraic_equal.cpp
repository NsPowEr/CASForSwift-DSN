#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/symbolic.hpp"
#include "cas/normal_form.hpp"
#include "algebra_internal.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace cas::symbolic {

// Implemented in algebraic_equal_weierstrass.cpp.
[[nodiscard]] bool weierstrass_zero_diff(ExprPtr diff_expr, CASContext& ctx);
[[nodiscard]] bool trig_exponential_zero_diff(ExprPtr diff_expr, CASContext& ctx);
[[nodiscard]] bool radical_zero_diff(ExprPtr diff_expr, CASContext& ctx);

Result<bool> mathematically_equal(ExprPtr lhs, ExprPtr rhs, CASContext& context) {
    const bool owns_operation = !context.operation_active_;
    if (owns_operation) {
        context.operation_active_ = true;
        // A51: inizializza il budget dell'operazione che stiamo aprendo. Senza,
        // le simplify interne vedono l'operazione gia' attiva e non lo fanno
        // loro, quindi il gate wall-clock misura dall'ultima operazione del
        // contesto — o dall'epoch, se questa e' la prima, e allora scatta
        // subito. Il confronto diventava non deterministico: dipendeva dalla
        // storia del contesto invece che dai suoi operandi.
        context.begin_operation_budget(/*capture_trace=*/false);
    }

    auto finalize = [&]() {
        if (owns_operation) {
            context.operation_active_ = false;
            context.end_operation_budget();
        }
    };

    auto lhs_s = context.simplify(lhs);
    if (lhs_s.is_error()) { finalize(); return fail<bool>(lhs_s.error()); }

    auto rhs_s = context.simplify(rhs);
    if (rhs_s.is_error()) { finalize(); return fail<bool>(rhs_s.error()); }

    // A44: unify the two spellings of the factorial (u! and factorial(u)) to
    // gamma(u+1) before comparing — exact identity, same motivation as the
    // hyperbolic rewrite below (CAS emits Γ(x+1), Maxima emits x!).
    {
        auto lhs_f = algebra::factorial_gamma_normalize(lhs_s.value(), context.arena());
        auto rhs_f = algebra::factorial_gamma_normalize(rhs_s.value(), context.arena());
        if (lhs_f.get() != lhs_s.value().get() || rhs_f.get() != rhs_s.value().get()) {
            auto lhs_f_s = context.simplify(lhs_f);
            auto rhs_f_s = context.simplify(rhs_f);
            if (lhs_f_s.is_ok()) lhs_s = lhs_f_s;
            if (rhs_f_s.is_ok()) rhs_s = rhs_f_s;
        }
    }

    // A43 §4: collapse the interchangeable spellings of the non-elementary
    // family (li→Ei(ln), Shi/Chi→Ei combinations, erfi→erf of an imaginary
    // argument) before comparing. Exact identities, same motivation as the two
    // normalisations around it: our integrator emits `(√π/2)·erfi(x)` where
    // Maxima emits `−i·√π·erf(i·x)/2` — equal, not structurally comparable.
    {
        auto lhs_n = algebra::nonelementary_normalize(lhs_s.value(), context.arena());
        auto rhs_n = algebra::nonelementary_normalize(rhs_s.value(), context.arena());
        if (lhs_n.get() != lhs_s.value().get() || rhs_n.get() != rhs_s.value().get()) {
            auto lhs_n_s = context.simplify(lhs_n);
            auto rhs_n_s = context.simplify(rhs_n);
            if (lhs_n_s.is_ok()) lhs_s = lhs_n_s;
            if (rhs_n_s.is_ok()) rhs_s = rhs_n_s;
        }
    }

    // F7.5.A4: rewrite sech/csch/coth/tanh to canonical cosh/sinh quotients
    // BEFORE structural / algebraic comparison so notational mismatches
    // (CAS cosh(x)^-2 vs Maxima sech(x)^2 etc.) collapse.
    {
        auto lhs_h = algebra::hyperbolic_normalize(lhs_s.value(), context.arena());
        auto rhs_h = algebra::hyperbolic_normalize(rhs_s.value(), context.arena());
        if (lhs_h.get() != lhs_s.value().get() || rhs_h.get() != rhs_s.value().get()) {
            auto lhs_h_s = context.simplify(lhs_h);
            auto rhs_h_s = context.simplify(rhs_h);
            if (lhs_h_s.is_ok()) lhs_s = lhs_h_s;
            if (rhs_h_s.is_ok()) rhs_s = rhs_h_s;
        }
    }

    // A50 (identità 2): espande sin/cos/sinh/cosh di un argomento Sum con
    // fase costante + parte simbolica via la formula di addizione esatta,
    // prima del confronto strutturale — stessa motivazione dei tre
    // normalizzatori sopra (spelling interscambiabili, non canonicalizzazione
    // globale). Chiude D(F)=f sulla famiglia trigonometrica traslata.
    {
        auto lhs_t = algebra::trig_addition_normalize(lhs_s.value(), context.arena());
        auto rhs_t = algebra::trig_addition_normalize(rhs_s.value(), context.arena());
        if (lhs_t.get() != lhs_s.value().get() || rhs_t.get() != rhs_s.value().get()) {
            auto lhs_t_s = context.simplify(lhs_t);
            auto rhs_t_s = context.simplify(rhs_t);
            if (lhs_t_s.is_ok()) lhs_s = lhs_t_s;
            if (rhs_t_s.is_ok()) rhs_s = rhs_t_s;
        }
    }

    if (structural_equal(lhs_s.value(), rhs_s.value())) {
        finalize();
        return ok(true);
    }

    // F7.5 (B.2 algebraic): single square-root extension Q(x)(√p). Proves
    // equality of antiderivatives that differ by a constant but whose
    // derivatives carry √p in a denominator (e.g. ∫√(x²−1): CAS acosh-form vs
    // Maxima log-form). Placed BEFORE the RootOf dispatch: a √(quadratic) is a
    // RootOf of t²−p, so `try_rootof_decision` would otherwise short-circuit
    // these to `false`. Sound: bails on anything outside the single-radical
    // class, the zero test needs both Q[x] coordinates of the numerator to
    // vanish, so it can only ever prove equality.
    {
        auto radical_diff = context.arena().make<Binary>(
            BinaryOp::Sub, lhs_s.value(), rhs_s.value());
        // Run the radical zero test as a fresh top-level operation: its
        // `polynomial_normal_form`/`expand` normalisation is not idempotent under
        // an already active operation (we re-simplified the operands above), which
        // would leave the coordinate polynomials un-reduced and the test would
        // spuriously fail. Restore the flag afterwards.
        const bool saved_active = context.operation_active_;
        context.operation_active_ = false;
        const bool radical_eq = radical_zero_diff(radical_diff, context);
        context.operation_active_ = saved_active;
        if (radical_eq) {
            finalize();
            return ok(true);
        }
    }

    // F7.5.A1 / HC-F75-CYCLOTOMIC-ROOTOF: RootOf-specific decisions.
    // See src/algebra/algebraic_equal_cyclotomic.cpp::try_rootof_decision.
    // RootOf dispatch — try the ORIGINAL operands first. simplify() above may
    // keep a RootOf (e.g. the biquadratic x⁴+1) while rewriting the OTHER side
    // into radical form; the cyclotomic enumerator emits exp-form roots, so a
    // simplified radical `other` would spuriously compare unequal and the call
    // would short-circuit to false. On the originals both sides are still in
    // their given (exp) notation, so Φ_n recognition matches. The simplified
    // forms are tried as a fallback for inputs where simplify EXPOSES a RootOf.
    for (const std::pair<ExprPtr, ExprPtr>& ops :
             {std::pair<ExprPtr, ExprPtr>{lhs, rhs},
              std::pair<ExprPtr, ExprPtr>{lhs_s.value(), rhs_s.value()}}) {
        if (auto rootof_dec = algebra::try_rootof_decision(ops.first, ops.second, context);
            rootof_dec.has_value()) {
            finalize();
            return ok(*rootof_dec);
        }
    }

    // F1.6 bridge: Constant::I ↔ ComplexLit(0,1) (and other purely-imaginary
    // ComplexLit forms) live on different AST kinds but represent the same
    // algebraic value.  A direct simplify of (lhs - rhs) drives both forms
    // through the same Sum coefficient pool (ComplexRational accumulator) and
    // collapses to IntegerLit(0) when equal.
    {
        auto diff_simple = context.simplify(
            context.arena().make<Binary>(BinaryOp::Sub, lhs_s.value(), rhs_s.value()));
        if (diff_simple.is_ok()) {
            if (const auto* il = expr_cast<IntegerLit>(diff_simple.value());
                il != nullptr && il->value.is_zero()) {
                finalize();
                return ok(true);
            }
            if (const auto* rl = expr_cast<RationalLit>(diff_simple.value());
                rl != nullptr && rl->numerator.is_zero()) {
                finalize();
                return ok(true);
            }
            if (const auto* cl = expr_cast<ComplexLit>(diff_simple.value());
                cl != nullptr && cl->re_num.is_zero() && cl->im_num.is_zero()) {
                finalize();
                return ok(true);
            }
        }
    }

    auto diff_expr = context.arena().make<Binary>(BinaryOp::Sub, lhs_s.value(), rhs_s.value());
    auto normal_diff = polynomial_normal_form(diff_expr, context);
    if (normal_diff.is_ok()) {
        if (expr_is<IntegerLit>(normal_diff.value()) && expr_cast<IntegerLit>(normal_diff.value())->value.is_zero()) {
            finalize();
            return ok(true);
        }
    }

    // F7.5 Weierstrass equivalence: when diff depends only on sin/cos of a
    // single var, substitute t = tan(v/2) and check if numerator collapses
    // to 0 as rational(t). Catches `1/(1+sin x)` vs `sin/(cos+1)+1` forms
    // and the wider trig-rational family the linear simplifier can't unify.
    if (weierstrass_zero_diff(diff_expr, context)) {
        finalize();
        return ok(true);
    }

    // F7.5 (B.2): exponential Laurent-in-e^{ix} canonical form for trig
    // POLYNOMIALS — succeeds on linearisable identities (e.g. sin³x·cos²x vs its
    // Fourier expansion) where the t=tan(x/2) rational-of-t numerator does not
    // collapse within budget. Sound: bails (false) on any non-trig-polynomial
    // sub-term, so it can only prove equality, never assert a false one.
    if (trig_exponential_zero_diff(diff_expr, context)) {
        finalize();
        return ok(true);
    }

    auto lhs_parts = algebra::split_num_den(lhs_s.value(), context);
    auto rhs_parts = algebra::split_num_den(rhs_s.value(), context);
    if (lhs_parts.is_ok() && rhs_parts.is_ok()) {
        auto cross_l = algebra::multiply_exprs(lhs_parts.value().numerator, rhs_parts.value().denominator, context);
        auto cross_r = algebra::multiply_exprs(rhs_parts.value().numerator, lhs_parts.value().denominator, context);
        if (cross_l.is_ok() && cross_r.is_ok()) {
            // Sibling of the `structural_equal(lhs_s, rhs_s)` test above, on the
            // cleared-denominator form: if num_l·den_r and num_r·den_l are the
            // SAME tree then the two fractions are equal, under the same
            // non-vanishing-denominator assumption this whole path already
            // makes. Cheap, and it decides before the normalisation below.
            //
            // Not redundant: A48 (simplify does not reach a fixpoint on a Sum
            // that a distribution step nested under a Product — `A − A` with
            // structurally identical `A = (2x−2)·exp(2x)` normalises to
            // `(−2 − 2x + 2x + 2)·exp(2x)`, not to 0) makes the normalisation
            // below fail on exactly the shape this test settles. Measured on
            // the A43 antiderivatives: `2·(2x−2)⁻¹·exp(2x)` vs
            // `(x−1)⁻¹·exp(2x)` — equal, and the cross products are identical
            // trees, yet their difference did not collapse.
            if (structural_equal(cross_l.value(), cross_r.value())) {
                finalize();
                return ok(true);
            }
            auto cross_diff_expr = context.arena().make<Binary>(BinaryOp::Sub, cross_l.value(), cross_r.value());
            auto cross_diff = polynomial_normal_form(cross_diff_expr, context);
            if (cross_diff.is_ok() && expr_is<IntegerLit>(cross_diff.value()) && expr_cast<IntegerLit>(cross_diff.value())->value.is_zero()) {
                finalize();
                return ok(true);
            }
            // B.2: trig RATIONALS — after clearing denominators the cross
            // difference num_l·den_r − num_r·den_l is a trig polynomial, so the
            // exponential Laurent form proves A/B = C/D for the linearisable
            // family (e.g. Weierstrass-form antiderivative derivatives) that the
            // opaque polynomial_normal_form leaves un-reduced.
            if (trig_exponential_zero_diff(cross_diff_expr, context)) {
                finalize();
                return ok(true);
            }
        }
    }

    finalize();
    return ok(false);
}

}  // namespace cas::symbolic
