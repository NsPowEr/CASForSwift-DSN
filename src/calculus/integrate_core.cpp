#include "integrate_engine.hpp"

#include "cas/algebra.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/error_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <utility>
#include <vector>

namespace cas::calculus::integrate_detail {

thread_local std::size_t Integrator::depth_ = 0U;

Integrator::Integrator(symbolic::CASContext& context) noexcept : context_(context), arena_(context.arena()) {}

Result<ExprPtr> Integrator::integrate(ExprPtr expr, const Symbol& var) {
    // HC-F75-A3-HARD-TIMEOUT: poll the interrupt flag at every recursive
    // entry into the integrator. The check is a single atomic load + branch
    // (CASContext::check_interrupt is inline noexcept), so the cost per
    // integrate node is negligible compared to the simplify/algebraic work
    // a single node typically triggers. Polling at the function head — both
    // here and in integrate_once — gives the runner's per-entry SIGALRM
    // handler a deterministic cancellation point even for heavy integrands
    // whose downstream simplify path does not itself poll.
    if (auto irq = context_.check_interrupt(); irq.is_error()) {
        return fail<ExprPtr>(irq.error());
    }
    if (depth_ >= context_.max_integration_depth()) {
        // F0.8-MIGRATED
        return make_unimplemented<ExprPtr>(
            "calculus", "Integrator::integrate",
            "recursion depth exceeded max_integration_depth",
            cas::error::reason_codes::INTEGRATE_DEPTH_EXCEEDED,
            "Increase ctx.max_integration_depth() or simplify the integrand before calling integrate()",
            "F0.8");
    }
    DepthGuard guard(depth_);

    // 0a. Algebraic-extension fast-path: integrate(RootOf(M, ω)) cannot be
    //     reduced to elementary form without the Algebraic Risch machinery
    //     (Bronstein §8).  Until that lands, emit a *symbolic* unevaluated
    //     Integral(RootOf(M, ω), x) so that downstream callers (e.g. Kovacic
    //     Case 3 solve_ode) see a valid AST node instead of an error.
    //     Tracked HC-KV-05 in HARDCODE_LEDGER.md.
    if (expr_is<RootOf>(expr)) {
        return ok(arena_.make<Integral>(expr, var, std::nullopt, std::nullopt));
    }

    // 0. Dirac sifting fast-path: ∫ δ(arg(t)) · f(t) dt = f(a)/|arg'(a)|.
    //    Detect early per evitare che integrate_once tratti δ come integrand
    //    elementare e fallisca.  La routine restituisce Unimplemented se
    //    l'integranda non contiene δ riconoscibile, fall-through naturale.
    {
        auto sifting = try_integrate_dirac_sifting(expr, var, context_);
        if (sifting.is_ok()) {
            return context_.simplify(sifting.value());
        }
    }

    // 1. Try elementary patterns first (substitution, partial fractions, arctan, etc.)
    {
        auto direct = integrate_once(expr, var);
        if (direct.is_ok()) {
            return direct;
        }
    }

    // NEW: Try general u-substitution recognition (CAS-L2-16)
    {
        auto sub_res = integrate_by_substitution(expr, var, context_);
        if (sub_res.is_ok() && sub_res.value().has_value()) {
            return ok(sub_res.value().value());
        }
    }

    // 2. Simplify then retry elementary patterns
    auto simplified = context_.simplify(expr);
    if (simplified.is_ok() && !structural_equal(simplified.value(), expr)) {
        auto direct = integrate_once(simplified.value(), var);
        if (direct.is_ok()) {
            return direct;
        }
    }

    // 3. Risch Algorithm (handles transcendental extensions, rational functions)
    {
        auto risch_result = integrate_risch(expr, var, context_);
        if (risch_result.is_ok()) {
            return context_.simplify(risch_result.value());
        }
    }

    // 4. Weierstrass substitution for rational R(sin x, cos x) integrands
    //    that Risch's trig handling cannot decompose. Generic algorithm —
    //    not a pattern table — driven by structural detection of sin/cos
    //    of var only.
    {
        auto weier_result = integrate_weierstrass(expr, var, context_);
        if (weier_result.is_ok()) {
            return context_.simplify(weier_result.value());
        }
    }

    // F0.8-MIGRATED
    return make_unimplemented<ExprPtr>(
        "calculus", "Integrator::integrate",
        "expression has no applicable integration strategy",
        cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
        "Try supplying assumptions, simplifying, or extending the Risch pipeline",
        "F0.8");
}

Result<bool> Integrator::expressions_match_after_simplify(ExprPtr lhs, ExprPtr rhs) {
    auto simplified_lhs = context_.simplify(lhs);
    if (simplified_lhs.is_error()) {
        return fail<bool>(simplified_lhs.error());
    }
    auto simplified_rhs = context_.simplify(rhs);
    if (simplified_rhs.is_error()) {
        return fail<bool>(simplified_rhs.error());
    }
    return ok(structural_equal(simplified_lhs.value(), simplified_rhs.value()));
}

Result<ExprPtr> Integrator::integrate_rational(ExprPtr expr, const Symbol& var) {
    auto parts = algebra::apart_num_den(expr, context_);
    if (parts.is_error()) return fail<ExprPtr>(parts.error());
    
    ExprPtr P = parts.value().numerator;
    ExprPtr Q = parts.value().denominator;
    
    // Handle polynomial part if deg(P) >= deg(Q)
    auto P_poly = algebra::parse_polynomial(P, var, context_);
    auto Q_poly = algebra::parse_polynomial(Q, var, context_);
    
    ExprPtr poly_integral = make_integer(arena_, 0);
    
    if (P_poly.is_ok() && Q_poly.is_ok() && !algebra::is_zero_poly(Q_poly.value())) {
        if (algebra::poly_degree(P_poly.value()) >= algebra::poly_degree(Q_poly.value())) {
            auto div_res = algebra::divide_poly_with_remainder(P_poly.value(), Q_poly.value(), context_);
            if (div_res.is_ok()) {
                std::vector<ExprPtr> terms;
                for (std::size_t k = 0; k < div_res.value().quotient.size(); ++k) {
                    ExprPtr coeff = div_res.value().quotient[k];
                    if (!coeff) continue;
                    if (const auto* il = expr_cast<IntegerLit>(coeff); il && il->value.is_zero()) continue;
                    
                    ExprPtr kp1 = make_integer(arena_, static_cast<long long>(k + 1));
                    ExprPtr term;
                    if (k == 0) {
                        term = make_product(arena_, {coeff, arena_.make<Symbol>(var)});
                    } else {
                        term = make_product(arena_, {
                            make_binary(arena_, BinaryOp::Div, coeff, kp1),
                            make_binary(arena_, BinaryOp::Pow, arena_.make<Symbol>(var), kp1)
                        });
                    }
                    terms.push_back(term);
                }
                if (!terms.empty()) {
                    auto poly_sum = make_sum(arena_, std::move(terms));
                    auto simplified = context_.simplify(poly_sum);
                    if (simplified.is_ok()) poly_integral = simplified.value();
                }
                
                auto rem_expr = algebra::polynomial_to_expr(div_res.value().remainder, var, context_);
                if (rem_expr.is_ok()) P = rem_expr.value();
            }
        }
    }

    DifferentialField field(var);
    auto hermite = hermite_reduce(P, Q, var, field, context_);
    if (hermite.is_error()) return fail<ExprPtr>(hermite.error());
    
    auto lrt_res = algebra::integrate_rational_lrt(hermite.value().remaining_P, hermite.value().remaining_Q, var, context_);
    if (lrt_res.is_error()) return fail<ExprPtr>(lrt_res.error());
    
    std::vector<ExprPtr> total;
    auto is_expr_zero = [](ExprPtr e) {
        if (!e) return true;
        if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
        return false;
    };

    if (!is_expr_zero(poly_integral)) total.push_back(poly_integral);
    if (!is_expr_zero(hermite.value().rational_part)) total.push_back(hermite.value().rational_part);
    if (!is_expr_zero(lrt_res.value())) total.push_back(lrt_res.value());
    
    if (total.empty()) return ok(make_integer(arena_, 0));
    if (total.size() == 1) return ok(total[0]);
    return ok(make_sum(arena_, std::move(total)));
}

Result<ExprPtr> Integrator::integrate_once(ExprPtr expr, const Symbol& var) {
    // HC-F75-A3-HARD-TIMEOUT: see comment in integrate() above. integrate_once
    // is the per-node recursive entry, so polling here guarantees that any
    // sub-integrand a strategy (by-parts, substitution, Hermite, …) hands
    // back to the dispatcher passes through the same cancellation check.
    if (auto irq = context_.check_interrupt(); irq.is_error()) {
        return fail<ExprPtr>(irq.error());
    }
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot integrate a null expression"));
    }

    if (expr_is<IntegerLit>(expr) || expr_is<RationalLit>(expr) || expr_is<Constant>(expr)) {
        return ok(make_product(arena_, {expr, arena_.make<Symbol>(var)}));
    }
    if (expr_is<DecimalLit>(expr)) {
        // F0.8-MIGRATED
        return make_unimplemented<ExprPtr>(
            "calculus", "Integrator::integrate_once",
            "DecimalLit",
            cas::error::reason_codes::INTEGRATE_DECIMAL_INPUT,
            "Convert DecimalLit to Rational before calling integrate()",
            "F0.8");
    }
    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        if (symbol->name == var.name) {
            return ok(make_product(arena_, {
                make_rational(arena_, 1, 2),
                make_binary(arena_, BinaryOp::Pow, expr, make_integer(arena_, 2)),
            }));
        }
        return ok(make_product(arena_, {expr, arena_.make<Symbol>(var)}));
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            if (matches_reciprocal_sqrt_one_minus_square(unary->operand, var)) {
                return ok(make_function(arena_, "arccos", {arena_.make<Symbol>(var)}));
            }
            auto inner = integrate_once(unary->operand, var);
            if (inner.is_error()) {
                return inner;
            }
            return ok(make_unary(arena_, UnaryOp::Neg, inner.value()));
        }
        // F0.8-MIGRATED
        return make_unimplemented<ExprPtr>(
            "calculus", "Integrator::integrate_once",
            "Unary non-Neg (likely Factorial)",
            cas::error::reason_codes::INTEGRATE_FACTORIAL,
            "Factorial integration is not supported; rewrite the integrand manually",
            "F0.8");
    }
    if (const auto* binary = expr_cast<Binary>(expr)) {
        return integrate_binary(*binary, var);
    }
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        return integrate_function(*call, var);
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(sum->terms.size());
        for (ExprPtr term : sum->terms) {
            auto integrated = integrate_once(term, var);
            if (integrated.is_error()) {
                return integrated;
            }
            terms.push_back(integrated.value());
        }
        return ok(make_sum(arena_, std::move(terms)));
    }
    if (const auto* product = expr_cast<Product>(expr)) {
        return integrate_product(*product, var);
    }

    // F0.8-MIGRATED
    return make_unimplemented<ExprPtr>(
        "calculus", "Integrator::integrate_once",
        "unknown expression kind",
        cas::error::reason_codes::INTEGRATE_UNKNOWN_EXPR,
        "Add a dispatch branch in integrate_once() for this expression node type",
        "F0.8");
}

Result<ExprPtr> Integrator::integrate_binary(const Binary& binary, const Symbol& var) {
    switch (binary.op) {
    case BinaryOp::Add: {
        auto lhs = integrate_once(binary.left, var);
        if (lhs.is_error()) {
            return lhs;
        }
        auto rhs = integrate_once(binary.right, var);
        if (rhs.is_error()) {
            return rhs;
        }
        return ok(make_sum(arena_, {lhs.value(), rhs.value()}));
    }
    case BinaryOp::Sub: {
        auto lhs = integrate_once(binary.left, var);
        if (lhs.is_error()) {
            return lhs;
        }
        auto rhs = integrate_once(binary.right, var);
        if (rhs.is_error()) {
            return rhs;
        }
        return ok(make_sum(arena_, {lhs.value(), make_unary(arena_, UnaryOp::Neg, rhs.value())}));
    }
    case BinaryOp::Mul:
        return integrate_product(Product({binary.left, binary.right}), var);
    case BinaryOp::Div:
        if (is_one(binary.left) && is_same_symbol(binary.right, var)) {
            return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {arena_.make<Symbol>(var)})}));
        }
        if (is_one(binary.left)) {
            if (const auto* denominator_power = expr_cast<Binary>(binary.right);
                denominator_power != nullptr && denominator_power->op == BinaryOp::Pow) {
                if (const auto* exponent = expr_cast<IntegerLit>(denominator_power->right);
                    exponent != nullptr && !exponent->value.is_negative() && !exponent->value.is_zero()) {
                    Binary reciprocal_power{
                        BinaryOp::Pow,
                        denominator_power->left,
                        arena_.make<IntegerLit>(-(exponent->value))};
                    return integrate_power(reciprocal_power, var);
                }
            }
        }
        if (auto numerator = exact_scalar_from_expr(binary.left);
            numerator.has_value() && numerator->numerator() != BigInt(0)) {
            if (auto affine = extract_affine_argument(binary.right, var);
                affine.has_value() && affine->coefficient.numerator() != BigInt(0)) {
                return ok(make_product(arena_, {
                    make_rational(arena_, (*numerator) / affine->coefficient),
                    make_function(arena_, "ln", {make_function(arena_, "abs", {binary.right})}),
                }));
            }
        }
        if (matches_one_plus_square(binary.right, var)) {
            if (auto numerator = extract_affine_argument(binary.left, var); numerator.has_value()) {
                std::vector<ExprPtr> terms;
                if (!numerator->coefficient.numerator().is_zero()) {
                    terms.push_back(make_product(arena_, {
                        make_rational(arena_, numerator->coefficient / Rational(BigInt(2))),
                        make_function(arena_, "ln", {make_function(arena_, "abs", {binary.right})}),
                    }));
                }
                if (!numerator->constant.numerator().is_zero()) {
                    terms.push_back(make_product(arena_, {
                        make_rational(arena_, numerator->constant),
                        make_function(arena_, "arctan", {arena_.make<Symbol>(var)}),
                    }));
                }
                if (!terms.empty()) {
                    return ok(make_sum(arena_, std::move(terms)));
                }
            }
        }
        if (is_one(binary.left) && matches_one_plus_square(binary.right, var)) {
            return ok(make_function(arena_, "arctan", {arena_.make<Symbol>(var)}));
        }
        {
            ExprPtr arctan_base{};
            if (is_one(binary.left) && matches_square_plus_constant_square(binary.right, var, arctan_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_binary(arena_, BinaryOp::Div,
                    make_function(arena_, "arctan", {make_binary(arena_, BinaryOp::Div, x, arctan_base)}),
                    arctan_base));
            }
        }
        if (is_one(binary.left)) {
            const auto* sqrt_call = expr_cast<FuncCall>(binary.right);
            if (sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U &&
                matches_one_minus_square(sqrt_call->args.front(), var)) {
                return ok(make_function(arena_, "arcsin", {arena_.make<Symbol>(var)}));
            }
            ExprPtr arcsin_base{};
            if (sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U &&
                matches_constant_square_minus_variable_square(sqrt_call->args.front(), var, arcsin_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_function(arena_, "arcsin", {make_binary(arena_, BinaryOp::Div, x, arcsin_base)}));
            }
            ExprPtr constant_base{};
            if (sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U &&
                matches_square_plus_constant_square(sqrt_call->args.front(), var, constant_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_function(arena_, "ln", {make_function(arena_, "abs", {make_sum(arena_, {x, binary.right})})}));
            }
            if (matches_square_minus_constant_square(binary.right, var, constant_base)) {
                ExprPtr x = arena_.make<Symbol>(var);
                return ok(make_product(arena_, {
                    make_binary(arena_, BinaryOp::Div, make_integer(arena_, 1), make_product(arena_, {make_integer(arena_, 2), constant_base})),
                    make_function(arena_, "ln", {make_function(arena_, "abs", {make_binary(arena_, BinaryOp::Div, make_sum(arena_, {x, make_unary(arena_, UnaryOp::Neg, constant_base)}), make_sum(arena_, {x, constant_base}))})}),
                }));
            }
        }
        if (is_negative_one(binary.left)) {
            const auto* sqrt_call = expr_cast<FuncCall>(binary.right);
            if (sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U &&
                matches_one_minus_square(sqrt_call->args.front(), var)) {
                return ok(make_function(arena_, "arccos", {arena_.make<Symbol>(var)}));
            }
        }
        // x/sqrt(Q): closed-form via u = Q, du = Q' dx, result = sqrt(Q) or -sqrt(Q)
        if (is_same_symbol(binary.left, var)) {
            const auto* sqrt_call = expr_cast<FuncCall>(binary.right);
            if (sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U) {
                ExprPtr radicand = sqrt_call->args.front();
                ExprPtr cbase{};
                // x/sqrt(1-x^2) = -sqrt(1-x^2)
                if (matches_one_minus_square(radicand, var)) {
                    return ok(make_unary(arena_, UnaryOp::Neg, binary.right));
                }
                // x/sqrt(a^2-x^2) = -sqrt(a^2-x^2)
                if (matches_constant_square_minus_variable_square(radicand, var, cbase)) {
                    return ok(make_unary(arena_, UnaryOp::Neg, binary.right));
                }
                // x/sqrt(x^2+1) = sqrt(x^2+1)
                if (matches_one_plus_square(radicand, var)) {
                    return ok(binary.right);
                }
                // x/sqrt(x^2+a^2) = sqrt(x^2+a^2)
                if (matches_square_plus_constant_square(radicand, var, cbase)) {
                    return ok(binary.right);
                }
                ExprPtr cbase2{};
                // x/sqrt(x^2-a^2) = sqrt(x^2-a^2)
                if (matches_square_minus_constant_square(radicand, var, cbase2)) {
                    return ok(binary.right);
                }
            }
        }
        if (!depends_on(binary.right, var)) {
            auto res = integrate_once(binary.left, var);
            if (res.is_ok()) {
                return ok(make_binary(arena_, BinaryOp::Div, res.value(), binary.right));
            }
        }
        // Pre-pass: cancel poly gcd(num,den) so (x^2+x+1)/(x^3-1) reduces to
        // 1/(x-1) before LRT/PFD emits spurious RootSum forms.
        if (auto dn = algebra::polynomial_degree(binary.left, var, context_),
                dd = algebra::polynomial_degree(binary.right, var, context_);
            dn.is_ok() && dd.is_ok()) {
            auto g = algebra::polynomial_gcd(binary.left, binary.right, var, context_);
            auto dg = g.is_ok() ? algebra::polynomial_degree(g.value(), var, context_)
                                : Result<std::size_t>(fail<std::size_t>(g.error()));
            if (dg.is_ok() && dg.value() > 0U) {
                auto pn = algebra::polynomial_exact_divide(binary.left, g.value(), var, context_);
                auto pd = algebra::polynomial_exact_divide(binary.right, g.value(), var, context_);
                if (pn.is_ok() && pd.is_ok()) {
                    auto simp = context_.simplify(make_binary(arena_, BinaryOp::Div, pn.value(), pd.value()));
                    if (simp.is_ok())
                        if (auto r = integrate_once(simp.value(), var); r.is_ok()) return r;
                }
            }
        }
        if (auto quadratic_integral = integrate_linear_over_quadratic(binary, var);
            quadratic_integral.is_ok()) {
            return quadratic_integral;
        }
        if (auto rational_integral = integrate_via_partial_fractions(make_binary(arena_, BinaryOp::Div, binary.left, binary.right), var);
            rational_integral.is_ok()) {
            return rational_integral;
        }
        // Fallback to integrate_rational if numer/denom are polynomials
        {
            auto P_poly = algebra::parse_polynomial(binary.left, var, context_);
            auto Q_poly = algebra::parse_polynomial(binary.right, var, context_);
            if (P_poly.is_ok() && Q_poly.is_ok()) {
                auto res = integrate_rational(make_binary(arena_, BinaryOp::Div, binary.left, binary.right), var);
                if (res.is_ok()) return res;
            }
        }
        return integrate_once(make_product(arena_, {binary.left, make_binary(arena_, BinaryOp::Pow, binary.right, make_integer(arena_, -1))}), var);
    case BinaryOp::Pow:
        return integrate_power(binary, var);
    case BinaryOp::Mod:
        // F0.8-MIGRATED
        return make_unimplemented<ExprPtr>(
            "calculus", "Integrator::integrate_binary",
            "BinaryOp::Mod expression",
            cas::error::reason_codes::INTEGRATE_MODULO,
            "Modulo integration is undefined in the symbolic domain; reformulate the integrand",
            "F0.8");
    case BinaryOp::Equal:
    case BinaryOp::Less:
    case BinaryOp::Greater:
    case BinaryOp::LessEqual:
    case BinaryOp::GreaterEqual:
        // F0.8-MIGRATED
        return make_unimplemented<ExprPtr>(
            "calculus", "Integrator::integrate_binary",
            "comparison/relational BinaryOp",
            cas::error::reason_codes::INTEGRATE_COMPARISON,
            "Integration of comparison operators is not defined; use Piecewise or conditional integration",
            "F0.8");
    }

    return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unknown binary operator"));
}

Result<ExprPtr> integrate_indefinite_impl(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    return Integrator(ctx).integrate(expr, var);
}

}  // namespace cas::calculus::integrate_detail
