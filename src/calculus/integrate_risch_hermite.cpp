// integrate_risch_hermite.cpp — Polynomial-quotient integration,
// Hermite reduction, and Rothstein-Trager assembly for the Risch pipeline.
// Bronstein "Symbolic Integration I", Chapters 5-6 (rational part).
//
// Public API (declared in integrate_risch_internal.hpp):
//   integrate_risch_poly_and_rational_part() — steps 3..end of integrate_risch
// integrate_log_polynomial_part() (∫ poly-in-ln(u) dx) moved to
// integrate_risch_logpoly.cpp (anti-monolith split).

#include "integrate_risch_internal.hpp"
#include "integrate_engine.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::calculus {

// Defined in integrate_risch_logpoly.cpp (anti-monolith split); uses
// algebra::PolyExpr, so the declaration needs polynomial_internal.hpp visible
// at the call site (already included above).
[[nodiscard]] Result<ExprPtr> integrate_log_polynomial_part(
    const algebra::PolyExpr& quot,
    ExprPtr u_arg,
    const Symbol& t_top,
    const Symbol& var,
    const DifferentialField& lower_field,
    symbolic::CASContext& context);

// True if `e` contains an elementary transcendental function (exp/ln/trig/
// inverse-trig/…) whose argument depends on `var` — i.e. `e` is not a pure
// rational function of `var`.
[[nodiscard]] static bool has_var_transcendental(ExprPtr e, const Symbol& var) {
    if (!e) return false;
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        bool arg_has_var = false;
        for (ExprPtr a : fc->args) if (depends_on(a, var)) { arg_has_var = true; break; }
        if (arg_has_var) {
            switch (fc->func_id) {
                case BuiltinOp::Exp: case BuiltinOp::Ln:  case BuiltinOp::Log:
                case BuiltinOp::Sin: case BuiltinOp::Cos: case BuiltinOp::Tan:
                case BuiltinOp::Sec: case BuiltinOp::Csc: case BuiltinOp::Cot:
                case BuiltinOp::Asin: case BuiltinOp::Acos: case BuiltinOp::Atan:
                case BuiltinOp::Unknown:  // asinh/acosh/… parsed as Unknown
                    return true;
                default: break;
            }
        }
        for (ExprPtr a : fc->args) if (has_var_transcendental(a, var)) return true;
        return false;
    }
    if (const auto* b = expr_cast<Binary>(e))
        return has_var_transcendental(b->left, var) || has_var_transcendental(b->right, var);
    if (const auto* u = expr_cast<Unary>(e)) return has_var_transcendental(u->operand, var);
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms) if (has_var_transcendental(t, var)) return true;
        return false;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr f : p->factors) if (has_var_transcendental(f, var)) return true;
        return false;
    }
    return false;
}

// HARDCODE-OF-PASSAGE: HC-A38-01 (see HARDCODE_LEDGER.md)
// Decide whether the t_top-polynomial coefficient `coeff` makes the polynomial-
// quotient integration UNSOUND here.
//
// The branch below integrates each coefficient with the lower-field machinery
// (limited_integrate_field first, plain integrate() as fallback).  Both are
// verified/sound in isolation, but A38 found that letting a coefficient carry
// a LOWER-TOWER generator (nested log-in-log, or an unrelated sibling like an
// exp generator) can drive the fallback's recursive integrate() call back into
// a STRUCTURALLY IDENTICAL polynomial-quotient sub-problem — an unbounded
// recursion (reproduced empirically: entry trace cycles kz=1→kz=0→kz=1... with
// no shrinking measure).  Until that recursion is bounded by a real
// termination argument (HC-A38-01), bail here — Unimplemented, never a hang
// or a silent wrong answer (REGOLA ZERO).  Rational integrands and
// single-generator towers never carry a sibling generator, so neither is
// affected; this is exactly why the existing end-to-end regression suite
// (including plain ∫ln(x)dx through the NEW limited_integrate_field-aware
// kz≥1 path) stays green.
[[nodiscard]] static bool coeff_blocks_poly_quotient(
    ExprPtr coeff, const Symbol& t_top, const DifferentialField& field,
    const Symbol& var, symbolic::CASContext& context) {
    if (!coeff) return false;
    bool has_sibling = false;
    for (const auto& ext : field.extensions()) {
        if (ext.t_var.name == t_top.name) continue;
        if (depends_on(coeff, ext.t_var)) { has_sibling = true; break; }
    }
    if (!has_sibling) return false;
    // Map the generator symbols back to their real expressions and simplify; a
    // surviving var-transcendental marks a genuinely non-elementary coefficient.
    auto real = field.from_field_generators(coeff, context);
    if (real.is_error()) return true;  // cannot verify reducibility → safe bail
    ExprPtr c = real.value();
    if (auto s = context.simplify(c); s.is_ok()) c = s.value();
    return has_var_transcendental(c, var);
}

// Steps 3..end of integrate_risch():
//   3.  Decompose gen_expr into P/Q w.r.t. topmost generator t_n.
//   3b. Extract polynomial quotient; handle base and transcendental cases.
//   4.  Hermite reduction: ∫ P/Q dt = A/B + ∫ C/D dt (D square-free).
//   5.  Rothstein-Trager for the square-free remainder.
//   Assembly: combine poly_integral_part + rational + log parts, map back.
Result<ExprPtr> integrate_risch_poly_and_rational_part(
    ExprPtr gen_expr,
    ExprPtr expr_original,
    const Symbol& var,
    const DifferentialField& field,
    ExprPtr poly_integral_part,
    symbolic::CASContext& context) {
    AstArena& arena = context.arena();
    (void)expr_original;  // round-trip verification handled by the per-coefficient guard

    // 3. Decompose into P/Q with respect to the topmost generator t_n
    Symbol t_top = field.extensions().empty() ? var : field.extensions().back().t_var;

    auto rational_res = algebra::apart_num_den(gen_expr, context);
    if (rational_res.is_error()) return fail<ExprPtr>(rational_res.error());
    ExprPtr P = rational_res.value().numerator;
    ExprPtr Q = rational_res.value().denominator;

    // 3b. Polynomial part
    auto P_poly = algebra::parse_polynomial(P, t_top, context);
    auto Q_poly = algebra::parse_polynomial(Q, t_top, context);

    if (P_poly.is_ok() && Q_poly.is_ok() && !algebra::is_zero_poly(Q_poly.value())) {
        if (algebra::poly_degree(P_poly.value()) >= algebra::poly_degree(Q_poly.value())) {
            auto div_res = algebra::divide_poly_with_remainder(P_poly.value(), Q_poly.value(), context);
            if (div_res.is_ok()) {
                const auto& quot = div_res.value().quotient;
                const auto& rem  = div_res.value().remainder;

                std::vector<ExprPtr> int_terms;
                if (field.extensions().empty()) {
                    // Base case: t_top = x, D(t_top) = 1
                    for (std::size_t k = 0; k < quot.size(); ++k) {
                        ExprPtr coeff = quot[k];
                        if (!coeff || algebra::poly_is_zero_expr(coeff)) continue;
                        if (k == 0) {
                            int_terms.push_back(arena.make<Binary>(BinaryOp::Mul,
                                coeff, arena.make<Symbol>(t_top)));
                        } else {
                            ExprPtr kp1 = arena.make<IntegerLit>(BigInt(static_cast<long long>(k + 1)));
                            ExprPtr t_pow = arena.make<Binary>(BinaryOp::Pow,
                                arena.make<Symbol>(t_top), kp1);
                            int_terms.push_back(arena.make<Binary>(BinaryOp::Div,
                                arena.make<Binary>(BinaryOp::Mul, coeff, t_pow), kp1));
                        }
                    }
                } else {
                    const auto& ext = field.extensions().back();
                    const bool handle_log_polynomial = (ext.type == ExtensionType::Logarithmic);

                    if (handle_log_polynomial) {
                        // Soundness precondition: bail when a coefficient carries a
                        // lower-tower generator that doesn't reduce away (nested
                        // log-in-log, or an unrelated sibling like an exp generator).
                        // A38 wired the §5.10 per-degree recursion below to the
                        // tower-recursive parametric solver (limited_integrate_field,
                        // A1+A26) and confirmed it is sound and reachable for
                        // single-generator towers, but for a lower-generator-carrying
                        // coefficient the fallback's recursive integrate() call was
                        // observed to re-enter a structurally identical polynomial-
                        // quotient sub-problem with no shrinking measure (unbounded
                        // recursion, reproduced via trace).  Bounding that recursion
                        // is tracked as a follow-up (see TASKLIST_MASTER.md); until
                        // then, Unimplemented here — never a hang or a silent wrong
                        // answer (REGOLA ZERO).
                        for (std::size_t k = 0; k < quot.size(); ++k) {
                            if (coeff_blocks_poly_quotient(quot[k], t_top, field, var, context)) {
                                return make_unimplemented<ExprPtr>(
                                    "calculus", "integrate_risch",
                                    "log-polynomial coefficient carries an "
                                    "irreducible sibling/lower generator (e.g. exp "
                                    "inside a ln-polynomial coefficient, or a nested "
                                    "log-in-log tower); the field-aware §5.10 recursion "
                                    "cannot yet bound its fallback recursion for this "
                                    "shape",
                                    cas::error::reason_codes::RISCH_LOG_EXTENSION_GENERAL,
                                    "Bound the fallback recursion for lower-generator "
                                    "coefficients in integrate_log_polynomial_part "
                                    "(Bronstein §5)",
                                    "F0.8");
                            }
                        }
                        // The lower field k = Q(x, t_1..t_{n-1}) is where the §5.10
                        // per-degree limited integration problems live (A38).
                        DifferentialField lower_field(field.base_var(),
                            std::vector<DifferentialExtension>(
                                field.extensions().begin(), field.extensions().end() - 1U));
                        auto B_res = integrate_log_polynomial_part(
                            quot, ext.argument, t_top, var, lower_field, context);
                        if (B_res.is_error()) return fail<ExprPtr>(B_res.error());
                        int_terms.push_back(B_res.value());
                    }

                    for (std::size_t k = 0; !handle_log_polynomial && k < quot.size(); ++k) {
                        ExprPtr coeff = quot[k];
                        if (!coeff) continue;
                        // Pre-simplify: divide_poly_with_remainder can leave
                        // Product([0,…]) uncollapsed; normalise here.
                        {
                            auto coeff_tog = algebra::together(coeff, context);
                            if (coeff_tog.is_ok()) {
                                auto coeff_simp = context.simplify(coeff_tog.value());
                                if (coeff_simp.is_ok()) coeff = coeff_simp.value();
                            }
                        }
                        if (algebra::poly_is_zero_expr(coeff)) continue;

                        if (ext.type == ExtensionType::Exponential) {
                            // ∫ a_k * t^k where t = exp(u), Dt = u't
                            if (k == 0) {
                                auto base_int = integrate(coeff, var, context);
                                if (base_int.is_ok()) int_terms.push_back(base_int.value());
                                else return base_int;
                            } else {
                                // Solve Dy + k*u'*y = a_k (Risch DE for exp, §5.8)
                                auto du_res = diff(ext.argument, var, 1U, context);
                                if (du_res.is_error()) return fail<ExprPtr>(du_res.error());
                                ExprPtr du = du_res.value();
                                ExprPtr f_expr = arena.make<Binary>(BinaryOp::Mul,
                                    arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k))), du);
                                
                                DifferentialField sub_field(field.base_var(), std::vector<DifferentialExtension>(
                                    field.extensions().begin(), field.extensions().end() - 1U));
                                auto y_res = solve_risch_de_general(f_expr, coeff, var, sub_field, context);
                                if (y_res.is_ok()) {
                                    ExprPtr t_pow = (k == 1)
                                        ? arena.make<Symbol>(t_top)
                                        : arena.make<Binary>(BinaryOp::Pow,
                                            arena.make<Symbol>(t_top),
                                            arena.make<IntegerLit>(BigInt(k)));
                                    int_terms.push_back(arena.make<Binary>(BinaryOp::Mul,
                                        y_res.value(), t_pow));
                                } else {
                                    // F0.8-MIGRATED
                                    return make_unimplemented<ExprPtr>(
                                        "calculus", "integrate_risch_exponential_extension",
                                        "Risch DE for exponential term has no polynomial solution",
                                        cas::error::reason_codes::RISCH_EXPONENTIAL_DE,
                                        "Implement full Risch DE solver for exponential extensions (Bronstein §5.8)",
                                        "F0.8");
                                }
                            }
                        }
                    }
                }
                if (!int_terms.empty()) {
                    ExprPtr raw = int_terms.size() == 1
                        ? int_terms[0]
                        : arena.make<Sum>(std::move(int_terms));
                    auto simp = context.simplify(raw);
                    if (simp.is_ok()) poly_integral_part = simp.value();
                }
                auto rem_expr = algebra::polynomial_to_expr(rem, t_top, context);
                if (rem_expr.is_ok()) P = rem_expr.value();
            }
        }
    }

    // Re-simplify P after polynomial division (normalise coefficient algebra).
    {
        auto P_tog = algebra::together(P, context);
        if (P_tog.is_ok()) {
            auto P_simp = context.simplify(P_tog.value());
            if (P_simp.is_ok()) P = P_simp.value();
        }
    }

    ExprPtr rational_part;
    ExprPtr rem_P;
    ExprPtr rem_Q;
    Result<ExprPtr> log_part_res = ok(static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(0))));
    bool P_is_zero = expr_is<IntegerLit>(P) && expr_ref<IntegerLit>(P).value.is_zero();
    if (P_is_zero) {
        rational_part = arena.make<IntegerLit>(BigInt(0));
        rem_P = arena.make<IntegerLit>(BigInt(0));
        rem_Q = Q;
    } else {
        // 4. Hermite Reduction: ∫ P/Q dt = A/B + ∫ C/D dt (D square-free)
        auto hermite_res = hermite_reduce(P, Q, t_top, field, context);
        if (hermite_res.is_error()) return fail<ExprPtr>(hermite_res.error());
        rational_part = hermite_res.value().rational_part;
        rem_P = hermite_res.value().remaining_P;
        rem_Q = hermite_res.value().remaining_Q;
        // 5. Rothstein-Trager for the square-free remainder
        log_part_res = integrate_rothstein_trager(rem_P, rem_Q, t_top, field, context);
    }

    if (log_part_res.is_ok()) {
        // Combine rational + log parts; map back from field generators.
        ExprPtr total_gen;
        bool rat_zero = expr_is<IntegerLit>(rational_part)
            && expr_ref<IntegerLit>(rational_part).value.is_zero();
        bool log_zero = expr_is<IntegerLit>(log_part_res.value())
            && expr_ref<IntegerLit>(log_part_res.value()).value.is_zero();
        if (rat_zero && log_zero) {
            total_gen = arena.make<IntegerLit>(BigInt(0));
        } else if (rat_zero) {
            total_gen = log_part_res.value();
        } else if (log_zero) {
            total_gen = rational_part;
        } else {
            total_gen = arena.make<Sum>(std::vector<ExprPtr>{rational_part, log_part_res.value()});
        }

        auto back_res = field.from_field_generators(total_gen, context);
        if (back_res.is_error()) return back_res;

        // poly_integral_part is in field generators too; must map back.
        auto poly_back_res = field.from_field_generators(poly_integral_part, context);
        if (poly_back_res.is_error()) return poly_back_res;
        bool poly_zero = expr_is<IntegerLit>(poly_back_res.value())
            && expr_ref<IntegerLit>(poly_back_res.value()).value.is_zero();
        bool back_zero = expr_is<IntegerLit>(back_res.value())
            && expr_ref<IntegerLit>(back_res.value()).value.is_zero();
        // A27 zero-result guard (REGOLA 0.2): 0 is an antiderivative only of
        // the zero function.  If every component collapsed to 0 while the
        // integrand is not literally 0, some sub-step silently lost the
        // integral (this exact path returned 0 for ∫ 1/(x·ln x·ln ln x) dx) —
        // report Unimplemented instead of a wrong answer.
        if (poly_zero && back_zero) {
            bool integrand_zero = false;
            if (auto is = context.simplify(gen_expr); is.is_ok()) {
                integrand_zero = (expr_is<IntegerLit>(is.value())
                    && expr_ref<IntegerLit>(is.value()).value.is_zero());
            }
            if (!integrand_zero) {
                return make_unimplemented<ExprPtr>(
                    "calculus", "integrate_risch",
                    "rational/log/poly parts all reduced to 0 for a non-zero integrand",
                    cas::error::reason_codes::RISCH_LOG_EXTENSION_GENERAL,
                    "A sub-step (Hermite/Rothstein-Trager) lost the integral; "
                    "the result would be silently wrong",
                    "A27");
            }
        }
        if (poly_zero) return back_res;
        if (back_zero) return context.simplify(poly_back_res.value());
        return context.simplify(arena.make<Sum>(std::vector<ExprPtr>{
            poly_back_res.value(), back_res.value()}));
    }

    // F0.8-MIGRATED
    return make_unimplemented<ExprPtr>(
        "calculus", "integrate_risch",
        "transcendental extension: integrability undecidable with current pipeline",
        cas::error::reason_codes::RISCH_LOG_EXTENSION_GENERAL,
        "Implement full Risch structure theorem for this extension type (Bronstein §5)",
        "F0.8");
}

} // namespace cas::calculus
