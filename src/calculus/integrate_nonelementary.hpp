#pragma once

// A43 §5 — Fallback per le antiderivate NON elementari (Ei, Si, Ci, Shi, Chi,
// li, Li₂, erfi). Interfaccia interna condivisa fra i due file d'attuazione:
//   integrate_nonelementary.cpp          — dispatcher + famiglia esponenziale
//   integrate_nonelementary_special.cpp  — gaussiana, trigonometrica, li, Li₂
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Nonelementary_Antiderivatives.md

#include "calculus_internal.hpp"
#include "integrate_engine.hpp"
#include "cas/rational.hpp"

#include <optional>
#include <vector>

namespace cas::calculus {

/// @brief A43 §5 — fallback per le antiderivate non elementari.
///
/// Montato in `integrate_core.cpp` rigorosamente DOPO Risch (spec §5: Risch è
/// la procedura di decisione — anticiparlo maschererebbe i casi elementari) e
/// PRIMA del fallback Meijer G, la cui forma ₁F₁ è equivalente ma meno
/// leggibile su questa classe.
[[nodiscard]] Result<ExprPtr> integrate_nonelementary_fallback(
    ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx);

}  // namespace cas::calculus

namespace cas::calculus::nonelementary {

/// @brief Un fattore del prodotto appiattito: `base^exponent`.
struct Factor {
    ExprPtr base;
    Rational exponent;
};

/// @brief Vista moltiplicativa piatta dell'integranda: `coefficient · Π base^exp`.
///
/// Serve perché la STESSA integranda arriva qui in forme strutturali diverse a
/// seconda del percorso che l'ha prodotta — `Pow(Q,-1)`, `Binary(Div,…)`,
/// `Product(N, Pow(D,-1))` — e un riconoscitore legato alla forma mancherebbe i
/// casi equivalenti (lezione misurata in A46: 40 s vs 0.1 s per pura forma).
struct ProductView {
    Rational coefficient{Rational(1)};
    std::vector<Factor> factors;
};

[[nodiscard]] ProductView flatten_product(ExprPtr expr);

/// @brief `true` se `expr` è `exp(u)` o `e^u`; restituisce `u` in `argument`.
[[nodiscard]] bool is_exponential(ExprPtr expr, ExprPtr& argument);

/// @brief Il polo `r` e il fattore di scala di un fattore `(α·x + β)^{-k}`.
///
/// `(α x + β)^{-k} = α^{-k} · (x − r)^{-k}` con `r = −β/α`.
struct LinearPole {
    Rational root;   // r
    Rational scale;  // α^{-k}
    long long order; // k ≥ 1
};

/// @param max_order budget di riduzione: un polo di ordine k costa k−1 passi di
///        integrazione per parti, quindi il tetto naturale è lo stesso della
///        ricorsione dell'integratore (`CASContext::max_integration_depth`),
///        configurabile e non una costante inventata.
[[nodiscard]] std::optional<LinearPole> as_linear_pole(
    const Factor& factor, const Symbol& var, std::size_t max_order);

/// @brief Converte un esponente intero non negativo entro `max_order` in
///        `long long`; `nullopt` se negativo o fuori budget.
[[nodiscard]] std::optional<long long> bounded_degree(const BigInt& value, std::size_t max_order);

// ∫ e^{a·x+b}·(x−r)^{−k} dx per k ≥ 1 (spec §5, riduzione per parti — non tabella).
[[nodiscard]] Result<ExprPtr> exp_over_linear_power(
    const Rational& a, const Rational& b, const Rational& r, long long k,
    const Symbol& var, symbolic::CASContext& ctx);

// ∫ e^{a·x+b}·x^n dx per n ≥ 0 — elementare, stessa ricorsione per parti.
[[nodiscard]] Result<ExprPtr> exp_times_monomial(
    const Rational& a, const Rational& b, long long n,
    const Symbol& var, symbolic::CASContext& ctx);

// Le famiglie non esponenziali (file `_special`).
[[nodiscard]] Result<ExprPtr> try_gaussian(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> try_trig_over_linear(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> try_log_integral(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx);
[[nodiscard]] Result<ExprPtr> try_dilogarithm(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx);

}  // namespace cas::calculus::nonelementary
