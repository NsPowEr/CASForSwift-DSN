// integrate_risch.cpp — Public dispatcher for the Risch integration pipeline.
// Bronstein "Symbolic Integration I", entry point + pre-processing steps.
//
// Steps handled here:
//   0.  Simple pattern matching (exp(x), ln(x)).
//   1.  Build differential extension tower.
//   2.  Map integrand to field generators.
//   2b. Logarithmic-derivative recognition (structure theorem, HPP-007).
//   2b-bis. PolyRischDE wiring for single-generator towers (cap.8 §6.4.1/§6.4.2).
//   2c. Product(f, exp(g)) Risch DE shortcut (cap.6).
//   Steps 3..end delegated to integrate_risch_poly_and_rational_part()
//   (defined in integrate_risch_hermite.cpp).

#include "integrate_risch_internal.hpp"
#include "integrate_engine.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

namespace cas::calculus {

Result<ExprPtr> integrate_risch(ExprPtr expr, const Symbol& var, symbolic::CASContext& context) {
    AstArena& arena = context.arena();

    // 0. Simple pattern matching for fundamental transcendental functions
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->args.size() == 1 && integrate_detail::is_same_symbol(call->args[0], var)) {
            if (call->func_id == BuiltinOp::Exp) {
                return ok(expr);
            }
            if (call->func_id == BuiltinOp::Ln) {
                ExprPtr x = arena.make<Symbol>(var);
                return ok(arena.make<Sum>(std::vector<ExprPtr>{
                    arena.make<Binary>(BinaryOp::Mul, x, expr),
                    arena.make<Unary>(UnaryOp::Neg, x)
                }));
            }
        }
    }

    // 1. Build Differential Extension Tower
    auto field_res = DifferentialField::build(expr, var, context);
    if (field_res.is_error()) return fail<ExprPtr>(field_res.error());
    const auto& field = field_res.value();

    // 2. Map expression to the differential field (generators t_1, ..., t_n)
    auto gen_expr_res = field.to_field_generators(expr, context);
    if (gen_expr_res.is_error()) return fail<ExprPtr>(gen_expr_res.error());
    ExprPtr gen_expr = gen_expr_res.value();

    // 2b. Logarithmic-derivative recognition (Risch structure theorem).
    // If integrand == c · D(g)/g for some generator g in the field tower
    // and c constant, then ∫ = c · ln(g). Handles ∫ 1/(x·ln(x)) dx = ln(ln(x))
    // and similar nested-log cases that Hermite/Rothstein-Trager cannot
    // express (resultant root is non-constant rational function of x).
    //
    // Formal-extraction approach (HPP-007 closure): for each candidate
    // g in the field tower, compute c = integrand / D(ln(g)) symbolically,
    // accept iff c is constant in var AND verifies the round-trip
    // D(c·ln(g)) − integrand = 0.  This eliminates the previous closed
    // set {±1, ±1/2, ±2} which silently missed cases like c=3 from
    // ∫ 3/(x·ln(x)) dx = 3·ln(ln(x)).
    auto verify_and_return = [&](ExprPtr cF) -> Result<ExprPtr> {
        auto DcF = diff(cF, var, 1U, context);
        if (DcF.is_error()) return fail<ExprPtr>(DcF.error());
        ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, DcF.value(), expr);
        auto delta_tog = algebra::together(delta, context);
        if (delta_tog.is_error()) return fail<ExprPtr>(delta_tog.error());
        auto delta_simp = context.simplify(delta_tog.value());
        if (delta_simp.is_error()) return fail<ExprPtr>(delta_simp.error());
        if (expr_is<IntegerLit>(delta_simp.value())
            && expr_ref<IntegerLit>(delta_simp.value()).value.is_zero()) {
            return context.simplify(cF);
        }
        return make_unimplemented<ExprPtr>(
            "calculus", "try_risch_log_candidate",
            "formal constant did not pass round-trip verification",
            cas::error::reason_codes::RISCH_NO_MATCH,
            "Verify the extension tower is well-formed for the integrand",
            "F0.8");
    };
    // Formal-extraction strategy for c = expr/DF.  The simplifier alone cannot
    // cancel mixed factors like x · (x·ln(x))^-1 · ln(x) because it treats the
    // transcendental ln(x) as opaque.  We use apart_num_den to split both expr
    // and DF into (numerator, denominator), cross-multiply to form a single
    // rational, and rely on simplify to cancel the resulting symbolic factors
    // structurally.  When that yields a var-independent expression, it is the
    // candidate constant — round-trip differentiation then verifies it.
    for (const auto& ext : field.extensions()) {
        ExprPtr g;
        if (ext.type == ExtensionType::Logarithmic) {
            g = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{ext.argument});
        } else if (ext.type == ExtensionType::Exponential) {
            g = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{ext.argument});
        } else {
            continue;
        }
        ExprPtr F_unit = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{g});
        auto DF_res = diff(F_unit, var, 1U, context);
        if (DF_res.is_error()) continue;
        ExprPtr DF = DF_res.value();
        // Split both sides into num/den so the cross-multiplied ratio
        //   c = (num(expr) · den(DF)) / (den(expr) · num(DF))
        // exposes the common transcendental factors to simplify().
        auto expr_nd = algebra::apart_num_den(expr,    context);
        auto df_nd   = algebra::apart_num_den(DF,      context);
        if (expr_nd.is_error() || df_nd.is_error()) continue;
        ExprPtr num = arena.make<Binary>(BinaryOp::Mul,
            expr_nd.value().numerator, df_nd.value().denominator);
        ExprPtr den = arena.make<Binary>(BinaryOp::Mul,
            expr_nd.value().denominator, df_nd.value().numerator);
        auto num_simp = context.simplify(num);
        auto den_simp = context.simplify(den);
        if (num_simp.is_error() || den_simp.is_error()) continue;
        ExprPtr c_raw = arena.make<Binary>(BinaryOp::Div,
            num_simp.value(), den_simp.value());
        // Replace each transcendental generator g_j of the differential field
        // with a fresh symbol u_j so the simplifier sees a multivariate
        // rational in {var, u_1, ..., u_k}: polynomial-level GCD cancellation
        // is then reachable (otherwise factors like (x·ln(x))^-1 vs ln(x)
        // remain opaque transcendentals and never cancel structurally).
        ExprPtr c_subst = c_raw;
        std::vector<Symbol> probe_syms;
        probe_syms.reserve(field.extensions().size());
        for (const auto& sub : field.extensions()) {
            ExprPtr gen;
            if (sub.type == ExtensionType::Logarithmic) {
                gen = arena.make<FuncCall>(BuiltinOp::Ln,
                    std::vector<ExprPtr>{sub.argument});
            } else if (sub.type == ExtensionType::Exponential) {
                gen = arena.make<FuncCall>(BuiltinOp::Exp,
                    std::vector<ExprPtr>{sub.argument});
            } else {
                continue;
            }
            Symbol fresh = context.make_fresh_symbol("hpp7");
            probe_syms.push_back(fresh);
            ExprPtr fresh_e = arena.make<Symbol>(fresh);
            c_subst = risch_helpers::deep_replace_expr(c_subst, gen, fresh_e, arena);
        }
        auto c_tog = algebra::together(c_subst, context);
        if (c_tog.is_error()) continue;
        auto c_exp  = algebra::expand(c_tog.value(), context);
        if (c_exp.is_error()) continue;
        auto c_simp = context.simplify(c_exp.value());
        if (c_simp.is_error()) continue;
        if (depends_on(c_simp.value(), var)) continue;  // not a constant in var
        bool depends_on_probe = false;
        for (const Symbol& u : probe_syms) {
            if (depends_on(c_simp.value(), u)) { depends_on_probe = true; break; }
        }
        if (depends_on_probe) continue;
        ExprPtr cF = arena.make<Binary>(BinaryOp::Mul, c_simp.value(), F_unit);
        auto res = verify_and_return(cF);
        if (res.is_ok()) return res;
    }

    // 2b-bis. PolyRischDE wiring (cap.8 §6.4.1/§6.4.2) — antiderivata diretta
    // di una integranda polinomiale in θ = log(u) o θ = exp(u) con u ∈ Q(var).
    //
    // Equivalenza:  ∫ g(x,θ) dx = y(x,θ)  ⟺  y' + 0·y = g.
    // Quindi chiamiamo solve_risch_de_{logarithmic,exponential}_q con f=0
    // e ricaviamo direttamente y come polinomio in θ a coefficienti razionali.
    //
    // Pre-processo cap.9 (structure theorem, sub-case log-factorization):
    // espandi log(u(x)) → Σ m_i · log(p_i) per fondere generatori log
    // algebricamente dipendenti in un singolo log irreducibile prima del scan.
    //
    // Riconoscimento del generatore θ: walk dell'AST per raccogliere tutte le
    // istanze di FuncCall(Ln/Exp).  Wiring fast-path solo per torre con un
    // unico generatore (log(u) o exp(u)) ben definito; torre mista o generatori
    // multipli ricadono su 2c/Hermite/Trager.
    {
        // Pre-processo cap.9 composto: prima log-factorization, poi
        // exp-decomposition.
        ExprPtr preproc = expand_log_args_via_factorization(expr, var, context);
        preproc = expand_exp_args_via_decomposition(preproc, var, context);
        ExprPtr ln_arg = nullptr;
        ExprPtr exp_arg = nullptr;
        bool generator_conflict = false;
        std::function<void(ExprPtr)> scan = [&](ExprPtr e) {
            if (!e || generator_conflict) return;
            if (const auto* fc = expr_cast<FuncCall>(e)) {
                if (fc->func_id == BuiltinOp::Ln && fc->args.size() == 1U) {
                    if (!ln_arg) ln_arg = fc->args[0];
                    else if (!risch_helpers::deep_struct_equal(ln_arg, fc->args[0]))
                        generator_conflict = true;
                } else if (fc->func_id == BuiltinOp::Exp && fc->args.size() == 1U) {
                    if (!exp_arg) exp_arg = fc->args[0];
                    else if (!risch_helpers::deep_struct_equal(exp_arg, fc->args[0]))
                        generator_conflict = true;
                }
                for (ExprPtr a : fc->args) scan(a);
                return;
            }
            if (const auto* bin = expr_cast<Binary>(e)) {
                scan(bin->left); scan(bin->right); return;
            }
            if (const auto* un = expr_cast<Unary>(e)) { scan(un->operand); return; }
            if (const auto* sum = expr_cast<Sum>(e)) {
                for (ExprPtr t : sum->terms) scan(t); return;
            }
            if (const auto* prod = expr_cast<Product>(e)) {
                for (ExprPtr f : prod->factors) scan(f); return;
            }
        };
        scan(preproc);

        if (!generator_conflict) {
            ExprPtr zero = arena.make<IntegerLit>(BigInt(0));
            auto roundtrip_ok = [&](ExprPtr y) -> bool {
                auto Dy = diff(y, var, 1U, context);
                if (Dy.is_error()) return false;
                ExprPtr Dy_exp = expand_log_args_via_factorization(Dy.value(), var, context);
                Dy_exp = expand_exp_args_via_decomposition(Dy_exp, var, context);
                ExprPtr expr_exp = expand_log_args_via_factorization(expr, var, context);
                expr_exp = expand_exp_args_via_decomposition(expr_exp, var, context);
                ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, Dy_exp, expr_exp);
                auto delta_tog = algebra::together(delta, context);
                if (delta_tog.is_error()) return false;
                auto delta_simp = context.simplify(delta_tog.value());
                if (delta_simp.is_error()) return false;
                if (const auto* il = expr_cast<IntegerLit>(delta_simp.value()))
                    return il->value.is_zero();
                if (const auto* rl = expr_cast<RationalLit>(delta_simp.value()))
                    return rl->numerator.is_zero();
                return false;
            };
            if (ln_arg && !exp_arg) {
                auto y = solve_risch_de_logarithmic_q(zero, preproc, ln_arg, var, context);
                if (y.is_ok() && roundtrip_ok(y.value())) {
                    return context.simplify(y.value());
                }
            } else if (exp_arg && !ln_arg) {
                auto y = solve_risch_de_exponential_q(zero, preproc, exp_arg, var, context);
                if (y.is_ok() && roundtrip_ok(y.value())) {
                    return context.simplify(y.value());
                }
            }
        }
    }

    // 2c. Product(f, exp(g)) Risch DE shortcut (Bronstein cap. 6).
    // For integrand = f(x)·exp(g(x)) with f, g polynomial in var, the
    // antiderivative has the form y(x)·exp(g(x)) where y satisfies:
    //   Dy + g'·y = f
    {
        ExprPtr exp_arg;
        std::vector<ExprPtr> poly_factors;
        if (const auto* call = expr_cast<FuncCall>(expr);
            call && call->func_id == BuiltinOp::Exp && call->args.size() == 1U) {
            exp_arg = call->args[0];
        } else if (const auto* prod = expr_cast<Product>(expr)) {
            for (ExprPtr f : prod->factors) {
                if (const auto* c2 = expr_cast<FuncCall>(f);
                    c2 && c2->func_id == BuiltinOp::Exp && c2->args.size() == 1U) {
                    if (exp_arg) { exp_arg = nullptr; break; }
                    exp_arg = c2->args[0];
                } else {
                    poly_factors.push_back(f);
                }
            }
        }
        if (exp_arg) {
            ExprPtr f_poly = poly_factors.empty()
                ? static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(1)))
                : (poly_factors.size() == 1U
                    ? poly_factors[0]
                    : static_cast<ExprPtr>(arena.make<Product>(std::move(poly_factors))));
            auto f_simp = context.simplify(f_poly);
            if (f_simp.is_ok()) f_poly = f_simp.value();
            auto g_simp = context.simplify(exp_arg);
            if (g_simp.is_ok()) exp_arg = g_simp.value();
            auto g_prime = diff(exp_arg, var, 1U, context);
            if (g_prime.is_ok()) {
                DifferentialField sub_field(field.base_var(), std::vector<DifferentialExtension>(
                    field.extensions().begin(), field.extensions().end() - 1U));
                auto y_res = solve_risch_de_general(g_prime.value(), f_poly, var, sub_field, context);
                if (y_res.is_ok()) {
                    ExprPtr exp_g = arena.make<FuncCall>(BuiltinOp::Exp,
                        std::vector<ExprPtr>{exp_arg});
                    ExprPtr antider = arena.make<Binary>(BinaryOp::Mul, y_res.value(), exp_g);
                    // Verify: D(y·exp(g)) = (y'+g'·y)·exp(g) = f·exp(g)
                    // iff y' + g'·y = f (pure rational check, no exp factor).
                    bool de_verified = false;
                    auto dy_res = diff(y_res.value(), var, 1U, context);
                    if (dy_res.is_ok()) {
                        ExprPtr gp_y = arena.make<Binary>(BinaryOp::Mul,
                            g_prime.value(), y_res.value());
                        ExprPtr lhs = arena.make<Sum>(std::vector<ExprPtr>{
                            dy_res.value(), gp_y});
                        ExprPtr residue = arena.make<Binary>(BinaryOp::Sub, lhs, f_poly);
                        auto resid_tog = algebra::together(residue, context);
                        ExprPtr resid_for_simp = resid_tog.is_ok() ? resid_tog.value() : residue;
                        auto resid_simp = context.simplify(resid_for_simp);
                        if (resid_simp.is_ok()
                            && expr_is<IntegerLit>(resid_simp.value())
                            && expr_ref<IntegerLit>(resid_simp.value()).value.is_zero()) {
                            de_verified = true;
                        }
                    }
                    // Fallback: full D(antider)==expr check.
                    if (!de_verified) {
                        auto D_check = diff(antider, var, 1U, context);
                        if (D_check.is_ok()) {
                            ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, D_check.value(), expr);
                            auto delta_tog = algebra::together(delta, context);
                            if (delta_tog.is_ok()) {
                                auto delta_simp = context.simplify(delta_tog.value());
                                if (delta_simp.is_ok()
                                    && expr_is<IntegerLit>(delta_simp.value())
                                    && expr_ref<IntegerLit>(delta_simp.value()).value.is_zero()) {
                                    de_verified = true;
                                }
                            }
                        }
                    }
                    if (de_verified) {
                        return context.simplify(antider);
                    }
                }
            }
        }
    }

    // Steps 3..end: poly quotient extraction + Hermite + Rothstein-Trager.
    ExprPtr poly_integral_part = arena.make<IntegerLit>(BigInt(0));
    return integrate_risch_poly_and_rational_part(
        gen_expr, expr, var, field, poly_integral_part, context);
}

} // namespace cas::calculus
