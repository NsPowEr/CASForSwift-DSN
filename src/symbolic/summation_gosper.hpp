#pragma once

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include <optional>

namespace cas::symbolic {

/**
 * @brief Algoritmo di Gosper per la summazione indefinita ipergeometrica.
 * 
 * Trova una funzione s(k) tale che s(k+1) - s(k) = term(k).
 * Se il termine non è ipergeometrico o non esiste una somma ipergeometrica,
 * restituisce std::nullopt.
 * 
 * @param term Il termine generale t(k) da sommare
 * @param k La variabile di sommazione
 * @param ctx Il contesto CAS
 * @return La somma indefinita s(k), oppure nullopt se non trovata
 */
[[nodiscard]] Result<std::optional<ExprPtr>> gosper_sum(
    const ExprPtr& term, 
    const Symbol& k, 
    symbolic::CASContext& ctx);

} // namespace cas::symbolic
