#pragma once

#include "calculus_internal.hpp"
#include "cas/rational.hpp"

#include <optional>
#include <string>
#include <vector>

namespace cas::calculus::integrate_detail {

struct AffineArgument {
    Rational coefficient;
    Rational constant;
};

struct QuadraticArgument {
    Rational quadratic;
    Rational linear;
    Rational constant;
};

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message);
[[nodiscard]] ExprPtr make_integer(AstArena& arena, long long value);
[[nodiscard]] ExprPtr make_rational(AstArena& arena, long long numerator, long long denominator);
[[nodiscard]] ExprPtr make_rational(AstArena& arena, const Rational& value);
[[nodiscard]] ExprPtr make_unary(AstArena& arena, UnaryOp op, ExprPtr operand);
[[nodiscard]] ExprPtr make_binary(AstArena& arena, BinaryOp op, ExprPtr lhs, ExprPtr rhs);
[[nodiscard]] ExprPtr make_sum(AstArena& arena, std::vector<ExprPtr> terms);
[[nodiscard]] ExprPtr make_product(AstArena& arena, std::vector<ExprPtr> factors);
[[nodiscard]] ExprPtr make_function(AstArena& arena, std::string name, std::vector<ExprPtr> args);
[[nodiscard]] std::string canonical_function_name(const std::string& name);
[[nodiscard]] bool depends_on(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool is_same_symbol(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool is_one(ExprPtr expr);
[[nodiscard]] bool is_negative_one(ExprPtr expr);
[[nodiscard]] std::optional<Rational> exact_scalar_from_expr(ExprPtr expr);
[[nodiscard]] std::optional<AffineArgument> extract_affine_argument(ExprPtr expr, const Symbol& var);
[[nodiscard]] std::optional<QuadraticArgument> extract_quadratic_argument(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool matches_square_of_variable(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool matches_square_plus_constant_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base);
[[nodiscard]] bool matches_square_minus_constant_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base);
[[nodiscard]] bool matches_one_plus_square(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool matches_one_minus_square(ExprPtr expr, const Symbol& var);
[[nodiscard]] bool is_rational_value(ExprPtr expr, long long numerator, long long denominator);
[[nodiscard]] bool matches_constant_square_minus_variable_square(ExprPtr expr, const Symbol& var, ExprPtr& constant_base);
[[nodiscard]] bool matches_reciprocal_sqrt_one_minus_square(ExprPtr expr, const Symbol& var);

// A53 — l'integratore procede per TENTATIVI: prova una strategia, e se fallisce
// ne prova un'altra sullo stesso integrando. Le side-conditions (A31) emesse da
// un tentativo fallito NON appartengono al risultato di quello che riesce, ma
// l'accumulatore del contesto e' monotono dentro un'operazione aperta e le
// terrebbe. Ogni tentativo va quindi eseguito qui dentro: il rollback scatta da
// solo sul ramo d'errore, il commit sul ramo di successo.
//
// Confine noto (dichiarato, non aggirato): copre il tentativo che FALLISCE. Un
// tentativo che riesce e il cui risultato viene poi scartato dal chiamante e'
// coperto solo se il chiamante stesso e' un tentativo — cioe' se anch'esso
// finisce per fallire. Attribuire le condizioni al singolo ExprPtr prodotto
// richiederebbe il canale `Simplified` di A31 §10.4 esteso a integrate().
// `accept` decide se il tentativo ha CONTRIBUITO al risultato: solo allora le
// sue condizioni restano. Serve parametrico perche' l'integratore ha tre modi
// di dire «rifiutato» — `Result` in errore, `ok(false)` di un predicato,
// `ok(nullopt)` di un matcher — e trattarne uno solo lascerebbe passare gli
// altri due (e' il caso misurato: il fallback Meijer/Mellin emette `x>0` e il
// suo risultato viene poi scartato piu' in alto).
template <class Fn, class Accept>
[[nodiscard]] auto attempt_with_condition_rollback(symbolic::CASContext& ctx, Fn&& fn, Accept&& accept)
    -> decltype(fn()) {
    symbolic::CASContext::SideConditionRollback guard(ctx);
    auto result = fn();
    if (accept(result)) {
        guard.commit();
    }
    return result;
}

template <class Fn>
[[nodiscard]] auto attempt_with_condition_rollback(symbolic::CASContext& ctx, Fn&& fn) -> decltype(fn()) {
    return attempt_with_condition_rollback(
        ctx, std::forward<Fn>(fn), [](const auto& r) { return r.is_ok(); });
}

// `ok(false)` e' un rifiuto, non un errore: le condizioni emesse dal confronto
// non appartengono a nulla. `ok(true)` invece le conferma — se il confronto e'
// vero solo dove vale una condizione, quella condizione limita davvero il
// risultato che il chiamante sta per accettare, e cancellarla sarebbe un
// silent-wrong di dominio.
inline constexpr auto accept_if_true = [](const auto& r) { return r.is_ok() && r.value(); };

// Matcher che restituiscono `Result<std::optional<T>>`: `ok(nullopt)` = «questa
// forma non e' mia», stesso trattamento del rifiuto.
inline constexpr auto accept_if_matched = [](const auto& r) { return r.is_ok() && r.value().has_value(); };

class Integrator {
public:
    explicit Integrator(symbolic::CASContext& context) noexcept;
    [[nodiscard]] Result<ExprPtr> integrate(ExprPtr expr, const Symbol& var);

    [[nodiscard]] Result<bool> expressions_match_after_simplify(ExprPtr lhs, ExprPtr rhs);
    [[nodiscard]] Result<ExprPtr> try_u_substitution_for_product(const Product& product, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_rational(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_via_partial_fractions(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_once(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_binary(const Binary& binary, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_linear_over_quadratic(const Binary& quotient, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_sqrt_quadratic(ExprPtr radicand, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_xsq_over_sqrt_quadratic(ExprPtr radicand, const Symbol& var);
    // ∫ dx/√(A x² + B x + C) for any rational A≠0, B, C — general completing-the-square
    // (arcsin / ln family). Generalizes the specific a²−x² / x²±a² matchers.
    [[nodiscard]] Result<ExprPtr> integrate_inverse_sqrt_quadratic(ExprPtr radicand, const Symbol& var);
    // ∫ xᵏ/√(c − d·x²) dx for integer k ≥ 0 via the reduction formula
    //   I_k = [(k−1)·c·I_{k−2} − xᵏ⁻¹·√(c−dx²)] / (k·d)
    // (no closed table). Radicand must be c − d·x² (B=0, A=−d<0, C=c>0) — the form
    // produced by the asin/acos IBP chain; other shapes → Unimplemented (diagnostic).
    [[nodiscard]] Result<ExprPtr> integrate_monomial_over_sqrt_quadratic(
        long long k, ExprPtr radicand, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_product(const Product& product, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_power(const Binary& power, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_inverse_quadratic_power(
        ExprPtr radicand, const BigInt& negative_exponent, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_function_direct(const std::string& name, ExprPtr argument);
    [[nodiscard]] Result<ExprPtr> integrate_power_direct(ExprPtr base, ExprPtr exponent, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_function(const FuncCall& call, const Symbol& var);

private:
    // A53 — ogni metodo pubblico sopra e' un tentativo: il wrapper omonimo (in
    // integrate_condition_guard.cpp) esegue l'`_impl` corrispondente dentro
    // `attempt_with_condition_rollback`, cosi' il contratto «un tentativo che
    // fallisce non lascia side-conditions» vale per costruzione e non per
    // annotazione dei ~60 call-site (che sarebbe violabile per omissione).
    [[nodiscard]] Result<ExprPtr> integrate_impl(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<bool> expressions_match_after_simplify_impl(ExprPtr lhs, ExprPtr rhs);
    [[nodiscard]] Result<ExprPtr> try_u_substitution_for_product_impl(const Product& product, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_rational_impl(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_via_partial_fractions_impl(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_once_impl(ExprPtr expr, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_binary_impl(const Binary& binary, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_linear_over_quadratic_impl(const Binary& quotient, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_sqrt_quadratic_impl(ExprPtr radicand, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_xsq_over_sqrt_quadratic_impl(ExprPtr radicand, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_inverse_sqrt_quadratic_impl(ExprPtr radicand, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_monomial_over_sqrt_quadratic_impl(
        long long k, ExprPtr radicand, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_product_impl(const Product& product, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_power_impl(const Binary& power, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_inverse_quadratic_power_impl(
        ExprPtr radicand, const BigInt& negative_exponent, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_function_direct_impl(const std::string& name, ExprPtr argument);
    [[nodiscard]] Result<ExprPtr> integrate_power_direct_impl(ExprPtr base, ExprPtr exponent, const Symbol& var);
    [[nodiscard]] Result<ExprPtr> integrate_function_impl(const FuncCall& call, const Symbol& var);

    symbolic::CASContext& context_;
    AstArena& arena_;

    struct DepthGuard {
        std::size_t& depth;
        explicit DepthGuard(std::size_t& d) : depth(d) { ++depth; }
        ~DepthGuard() { --depth; }
    };
    static thread_local std::size_t depth_;
};

[[nodiscard]] Result<ExprPtr> integrate_indefinite_impl(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx);

}  // namespace cas::calculus::integrate_detail
