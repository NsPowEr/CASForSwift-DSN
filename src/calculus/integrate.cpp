#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "cas/bigfloat.hpp"
#include "cas/extended_real.hpp"
#include "cas/rational.hpp"
#include "integrate_definite_guards.hpp"
#include "integrate_definite_patterns.hpp"
#include "integrate_engine.hpp"

#include <functional>
#include <optional>

namespace cas::calculus {

namespace {

// Adopt the canonical extended-real predicates; the local copies missed
// Constant(NegInfinity) (only handled Unary(Neg, Constant(Infinity))).
using cas::is_pos_infinity;
using cas::is_neg_infinity;

} // anonymous namespace
Result<ExprPtr> integrate(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    if (ctx.is_caching_enabled()) {
        auto key = symbolic::CASContext::IntegrateKey{expr, var.name};
        if (auto found = ctx.integrate_cache_.get(key)) {
            return ok(*found);
        }
    }

    // A53 — l'integrazione e' UNA operazione, con un budget suo.
    //
    // Prima: le ops si contavano solo in Simplifier/Substituter e ogni
    // `ctx.simplify()` dell'integratore era top-level, quindi azzerava
    // contatore e timer — nessun budget limitava il TOTALE, e il risultato
    // dipendeva dal tempo concesso invece che dall'integranda (misurato: la
    // stessa entry consuma per intero qualunque cap, e cap diversi cambiano
    // entry diverse). L'unico freno era il SIGALRM esterno, cioe' il carico
    // della macchina.
    //
    // Lo scope e' rientrante: le integrazioni annidate (by-parts, sostituzione)
    // NON riaprono l'operazione, quindi il loro costo resta addebitato a quella
    // piu' esterna — che e' il punto: e' il totale a dover essere limitato.
    // Il tetto e' `max_integration_ops` (derivazione dalla misura al getter).
    symbolic::CASContext::OperationScope op_scope(
        ctx, ctx.is_trace_enabled(), ctx.max_integration_ops());
    auto primitive = integrate_detail::integrate_indefinite_impl(expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }
    auto materialized = symbolic::materialize_expr(primitive.value(), ctx.arena());
    if (materialized.is_error()) {
        return materialized;
    }
    if (ctx.is_caching_enabled()) {
        auto key = symbolic::CASContext::IntegrateKey{expr, var.name};
        ctx.integrate_cache_.put(key, materialized.value());
    }
    return materialized;
}

Result<ExprPtr> definite_integral(ExprPtr expr, const Symbol& var, ExprPtr lower, ExprPtr upper, symbolic::CASContext& ctx) {
    ExprPtr normalized_expr = integrate_detail::normalize_definite_integrand(expr, ctx);

    // Extensible pattern table: each matcher returns nullopt to skip, value to commit.
    DefiniteContext dc{
        .integrand = expr,
        .integrand_normalized = normalized_expr,
        .var = var,
        .lower = lower,
        .upper = upper,
        .ctx = ctx,
    };
    for (DefinitePatternFn matcher : definite_patterns()) {
        auto match = matcher(dc);
        if (match.is_error()) return fail<ExprPtr>(match.error());
        if (match.value().has_value()) return ok(match.value().value());
    }

    // Generic infinite-domain fallback: only the Gaussian pattern is currently handled there;
    // anything else over (-inf, +inf) goes to Unimplemented.
    if (is_neg_infinity(lower) && is_pos_infinity(upper)) {
        return fail<ExprPtr>(integrate_detail::make_error(CASErrorKind::Unimplemented,
            "Integrazione su dominio infinito: pattern non riconosciuto."));
    }

    auto pole_check = integrate_detail::reject_poles_in_closed_interval(normalized_expr, var, lower, upper, ctx);
    if (pole_check.is_error()) {
        return fail<ExprPtr>(pole_check.error());
    }

    auto primitive = integrate(normalized_expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }

    auto lower_value = ctx.substitute(primitive.value(), var, lower);
    if (lower_value.is_error()) {
        return lower_value;
    }

    auto upper_value = ctx.substitute(primitive.value(), var, upper);
    if (upper_value.is_error()) {
        return upper_value;
    }

    return ctx.simplify(integrate_detail::make_sum(ctx.arena(), {
        upper_value.value(),
        integrate_detail::make_unary(ctx.arena(), UnaryOp::Neg, lower_value.value()),
    }));
}

}  // namespace cas::calculus
