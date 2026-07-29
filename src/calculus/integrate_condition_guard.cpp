// A53 — confine unico dei tentativi dell'integratore.
//
// Ogni metodo pubblico di `Integrator` e' un TENTATIVO: il chiamante lo esegue
// e, se fallisce, ne prova un altro sullo stesso integrando. Le side-conditions
// (A31) emesse durante un tentativo fallito non appartengono al risultato del
// tentativo che poi riesce; ma dentro un'operazione aperta (A53 apre
// `OperationScope` su calculus::integrate per dare un budget al TOTALE)
// l'accumulatore `side_conditions_` non viene piu' azzerato dalle simplify
// interne, quindi le terrebbe. Misurato: `∫e^{-x²}` usciva con `x>0`, emessa
// dal fallback Meijer/Mellin dentro una sotto-integrazione di by-parts poi
// scartata, mentre il risultato viene dal completamento del quadrato — esatto
// su tutto R.
//
// I wrapper stanno tutti qui invece che accanto alle rispettive `_impl` perche'
// il contratto sia leggibile in un punto solo, e perche' aggiungere un metodo
// pubblico senza il suo wrapper diventi visibile (il file elenca il perimetro).
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Domain_Conditions_Propagation.md

#include "integrate_engine.hpp"

namespace cas::calculus::integrate_detail {

Result<ExprPtr> Integrator::integrate(ExprPtr expr, const Symbol& var) {
    return attempt_with_condition_rollback(context_, [&] { return integrate_impl(expr, var); });
}

Result<bool> Integrator::expressions_match_after_simplify(ExprPtr lhs, ExprPtr rhs) {
    return attempt_with_condition_rollback(
        context_, [&] { return expressions_match_after_simplify_impl(lhs, rhs); }, accept_if_true);
}

Result<ExprPtr> Integrator::try_u_substitution_for_product(const Product& product, const Symbol& var) {
    return attempt_with_condition_rollback(
        context_, [&] { return try_u_substitution_for_product_impl(product, var); });
}

Result<ExprPtr> Integrator::integrate_rational(ExprPtr expr, const Symbol& var) {
    return attempt_with_condition_rollback(context_, [&] { return integrate_rational_impl(expr, var); });
}

Result<ExprPtr> Integrator::integrate_via_partial_fractions(ExprPtr expr, const Symbol& var) {
    return attempt_with_condition_rollback(
        context_, [&] { return integrate_via_partial_fractions_impl(expr, var); });
}

Result<ExprPtr> Integrator::integrate_once(ExprPtr expr, const Symbol& var) {
    return attempt_with_condition_rollback(context_, [&] { return integrate_once_impl(expr, var); });
}

Result<ExprPtr> Integrator::integrate_binary(const Binary& binary, const Symbol& var) {
    return attempt_with_condition_rollback(context_, [&] { return integrate_binary_impl(binary, var); });
}

Result<ExprPtr> Integrator::integrate_linear_over_quadratic(const Binary& quotient, const Symbol& var) {
    return attempt_with_condition_rollback(
        context_, [&] { return integrate_linear_over_quadratic_impl(quotient, var); });
}

Result<ExprPtr> Integrator::integrate_sqrt_quadratic(ExprPtr radicand, const Symbol& var) {
    return attempt_with_condition_rollback(
        context_, [&] { return integrate_sqrt_quadratic_impl(radicand, var); });
}

Result<ExprPtr> Integrator::integrate_xsq_over_sqrt_quadratic(ExprPtr radicand, const Symbol& var) {
    return attempt_with_condition_rollback(
        context_, [&] { return integrate_xsq_over_sqrt_quadratic_impl(radicand, var); });
}

Result<ExprPtr> Integrator::integrate_inverse_sqrt_quadratic(ExprPtr radicand, const Symbol& var) {
    return attempt_with_condition_rollback(
        context_, [&] { return integrate_inverse_sqrt_quadratic_impl(radicand, var); });
}

Result<ExprPtr> Integrator::integrate_monomial_over_sqrt_quadratic(
    long long k, ExprPtr radicand, const Symbol& var) {
    return attempt_with_condition_rollback(
        context_, [&] { return integrate_monomial_over_sqrt_quadratic_impl(k, radicand, var); });
}

Result<ExprPtr> Integrator::integrate_product(const Product& product, const Symbol& var) {
    return attempt_with_condition_rollback(context_, [&] { return integrate_product_impl(product, var); });
}

Result<ExprPtr> Integrator::integrate_power(const Binary& power, const Symbol& var) {
    return attempt_with_condition_rollback(context_, [&] { return integrate_power_impl(power, var); });
}

Result<ExprPtr> Integrator::integrate_inverse_quadratic_power(
    ExprPtr radicand, const BigInt& negative_exponent, const Symbol& var) {
    return attempt_with_condition_rollback(context_, [&] {
        return integrate_inverse_quadratic_power_impl(radicand, negative_exponent, var);
    });
}

Result<ExprPtr> Integrator::integrate_function_direct(const std::string& name, ExprPtr argument) {
    return attempt_with_condition_rollback(
        context_, [&] { return integrate_function_direct_impl(name, argument); });
}

Result<ExprPtr> Integrator::integrate_power_direct(ExprPtr base, ExprPtr exponent, const Symbol& var) {
    return attempt_with_condition_rollback(
        context_, [&] { return integrate_power_direct_impl(base, exponent, var); });
}

Result<ExprPtr> Integrator::integrate_function(const FuncCall& call, const Symbol& var) {
    return attempt_with_condition_rollback(context_, [&] { return integrate_function_impl(call, var); });
}

}  // namespace cas::calculus::integrate_detail
