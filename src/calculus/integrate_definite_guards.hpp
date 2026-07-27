#pragma once

// A53 (anti-monolito) — guardie di ammissibilita' dell'integrale DEFINITO,
// estratte da integrate.cpp: riconoscimento di poli razionali e di singolarita'
// trigonometriche/algebriche nell'intervallo, piu' la normalizzazione
// dell'integranda. Sono decisioni sul DOMINIO, indipendenti dalla ricerca
// dell'antiderivata, e hanno percio' vita propria.

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

namespace cas::calculus::integrate_detail {

// Errore Undefined se l'intervallo chiuso [lower, upper] contiene un polo della
// funzione razionale, o una singolarita' non razionale (trig/algebrica) che il
// trattamento improprio non copre.
[[nodiscard]] Result<void> reject_poles_in_closed_interval(
    ExprPtr integrand, const Symbol& var, ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx);

// together + simplify sull'integranda; se una delle due fallisce restituisce
// l'espressione di partenza (la normalizzazione e' un'ottimizzazione di forma,
// mai una precondizione di correttezza).
[[nodiscard]] ExprPtr normalize_definite_integrand(ExprPtr expr, symbolic::CASContext& ctx);

}  // namespace cas::calculus::integrate_detail
