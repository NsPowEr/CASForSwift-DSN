#pragma once
// A46 — Rioboo conversion (Bronstein, Symbolic_Integration_I.md §2.8): resa in
// forma REALE della somma formale sui residui prodotta da Lazard-Rioboo-Trager.
//
// Interfaccia condivisa fra il driver LRT (partial_fractions_lrt.cpp) e la
// macchina LogToReal (partial_fractions_logtoreal.cpp).

#include "polynomial_internal.hpp"

#include <utility>
#include <vector>

namespace cas {
namespace symbolic {
class CASContext;
}

namespace algebra {

// A45: un operando *zero* non e' l'identita' moltiplicativa (0·b = 0, non b).
// Un ExprPtr nullo invece significa "fattore assente" e viene saltato.
[[nodiscard]] ExprPtr rioboo_mul(AstArena& arena, ExprPtr a, ExprPtr b);

// Forma reale di Rioboo per un fattore QUADRATICO z² + a·z + b del resultant di
// Rothstein-Trager, i cui due residui sono una coppia coniugata (richiede
// 4b − a² > 0, cioe' discriminante negativo):
//
//   (−a/2)·log(norm) − √(4b−a²)·arctan( √(4b−a²)·Q₁ / (2Q₀ − a·Q₁) )
//
// dove Q₁·z + Q₀ = S(z,x) mod (z² + a·z + b) e
// norm = Q₀² − a·Q₁·Q₀ + b·Q₁².
//
// `a` e `b` sono espressioni qualunque, non necessariamente razionali: e' cosi'
// che il caso quartico (a, b in un'estensione quadratica di Q) riusa la stessa
// forma chiusa del caso razionale.
[[nodiscard]] Result<ExprPtr> rioboo_quadratic_real_form(
    ExprPtr a,
    ExprPtr b,
    const PolyExpr& G_z_x,
    const Symbol& var,
    const Symbol& z_var,
    symbolic::CASContext& ctx);

// Numero di radici reali DISTINTE di un polinomio razionale, per sequenza di
// Sturm valutata a ±∞ (solo segni dei coefficienti di testa): esatto, nessuna
// aritmetica in virgola mobile e nessuna tolleranza.
[[nodiscard]] std::size_t count_real_roots_rational(RatPoly poly);

// Spezza un quartico R(z) senza radici reali nei suoi DUE fattori quadratici
// reali z² + aᵢ·z + bᵢ, quando i residui vivono in un'estensione quadratica di
// Q — cioe' quando il cubico risolvente ha una radice razionale y₀ con
// p²−4q+4y₀ ≥ 0 e y₀²−4s ≥ 0. La fattorizzazione trovata e' SEMPRE verificata
// per espansione prima di essere restituita.
//
// Vettore vuoto = fuori dalla classe supportata (il chiamante torna alla somma
// formale `RootSum`): i residui di un quartico generico non sono esprimibili
// per radicali reali, e restituire una forma chiusa sbagliata sarebbe peggio.
[[nodiscard]] Result<std::vector<std::pair<ExprPtr, ExprPtr>>>
real_quadratic_factors_of_quartic(
    const PolyExpr& R_z,
    const Symbol& z_var,
    symbolic::CASContext& ctx);

}  // namespace algebra
}  // namespace cas
