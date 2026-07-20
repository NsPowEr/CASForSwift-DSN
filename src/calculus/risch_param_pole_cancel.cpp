// risch_param_pole_cancel.cpp — HC-A26-PRIMITIVE-PARAMQ-RATIONAL:
// cancellazione dei poli fra le forzanti di una Risch DE parametrica su un
// monomio PRIMITIVO (Bronstein "Symbolic Integration I", §7.1).
//
// Perché serve.  Nel caso primitivo (Dt ∈ k, cioè t = log(u)) il Teorema 5.1.1
// dà S = k: non esistono polinomi speciali non banali, quindi (formula 5.1)
//
//      k⟨t⟩ = k[t].
//
// Di conseguenza, per QUALSIASI soluzione q ∈ k[t] di
//
//      a·D(q) + b·q = Σ_i c_i·g_i        con a, b ∈ k[t]                (7.3)
//
// il membro sinistro è in k[t].  Quindi anche Σ_i c_i·g_i DEVE essere in k[t]:
// le parti proprie (i poli in t) delle singole forzanti devono CANCELLARSI
// nella combinazione.  Questo è un vincolo Q-lineare sui c_i, non una
// proprietà delle g_i prese una per una — ed è esattamente ciò che il vecchio
// bound di denominatore "D = lcm(denominatori)" sbagliava: assumeva che la
// soluzione avesse i poli delle forzanti, moltiplicava tutto per D e produceva
// f_new = f − D'/D non polinomiale (bail "f_new is not polynomial"), pur
// esistendo una soluzione polinomiale.
//
// Esempio minimo (il caso che ha aperto la voce di ledger): su k = Q(x),
// t = log(x), le forzanti g = (1/x + 1/(x·t), −1/(x·t)) hanno entrambe un polo
// in t = 0, ma la combinazione c₀ = c₁ le cancella e lascia 1/x, la cui
// primitiva q = t è polinomiale.  Il vecchio codice prendeva D = t·x² e
// falliva; questa riduzione trova il vincolo c₀ = c₁, riduce la famiglia alla
// singola forzante 1/x e risolve.
//
// Cosa fa questo file: dato {g_i}, calcola il sottospazio {c : Σ c_i g_i ∈ k[t]}
// e restituisce una base di forzanti ridotte (tutte polinomiali in t) più la
// matrice di cambio base, così che il chiamante possa risolvere il problema
// parametrico ridotto e rimappare le costanti sulla famiglia originale.

#include "calculus_internal.hpp"
#include "risch_parametric_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

using detail::null_space_basis;
using detail::poly_coeffs_q;
using detail::rational_to_expr;
using detail::row_echelon;

// Parte propria di g rispetto a t: g = (parte polinomiale in t) + num/den con
// deg_t(num) < deg_t(den).  Ritorna nullopt se g non è una funzione razionale
// in t trattabile.  `den_is_trivial` segnala che g è già polinomiale in t.
struct ProperPart {
    ExprPtr num;               // numeratore della parte propria (poly in t)
    ExprPtr den;               // denominatore (poly in t, grado ≥ 1)
    bool    den_is_trivial;    // true ⇒ g ∈ k[t], num/den ignorati
};

[[nodiscard]] std::optional<ProperPart> proper_part_in_t(
    ExprPtr g, const Symbol& t, symbolic::CASContext& ctx) {
    auto parts = algebra::apart_num_den(g, ctx);
    if (parts.is_error()) return std::nullopt;
    ExprPtr num = parts.value().numerator;
    ExprPtr den = parts.value().denominator;

    auto den_poly = algebra::parse_polynomial(den, t, ctx);
    if (den_poly.is_error()) return std::nullopt;
    if (algebra::poly_degree(den_poly.value()) == 0U) {
        // Nessun polo in t: g è già in k[t] (i coefficienti restano in k).
        return ProperPart{nullptr, nullptr, true};
    }
    auto num_poly = algebra::parse_polynomial(num, t, ctx);
    if (num_poly.is_error()) return std::nullopt;

    auto div = algebra::divide_poly_with_remainder(num_poly.value(), den_poly.value(), ctx);
    if (div.is_error()) return std::nullopt;
    auto rem_expr = algebra::polynomial_to_expr(div.value().remainder, t, ctx);
    if (rem_expr.is_error()) return std::nullopt;
    if (algebra::poly_is_zero_expr(rem_expr.value())) {
        // La divisione è esatta: g è polinomiale in t a meno di normalizzazione.
        return ProperPart{nullptr, nullptr, true};
    }
    return ProperPart{rem_expr.value(), den, false};
}

// Aggiunge a `rows` le equazioni Q-lineari nei c_i imposte da
// Σ_i c_i · coeff_i ≡ 0, dove coeff_i ∈ k = Q(base_var) sono funzioni
// razionali: si portano a denominatore comune e si annulla ogni potenza di
// base_var del numeratore.
[[nodiscard]] bool append_k_constraint(
    const std::vector<ExprPtr>& coeffs, const Symbol& base_var,
    std::vector<std::vector<Rational>>& rows, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    const std::size_t m = coeffs.size();

    // Denominatore comune (prodotto: non serve l'lcm minimo, solo un comune).
    ExprPtr common_den = arena.make<IntegerLit>(BigInt(1));
    std::vector<ExprPtr> nums(m), dens(m);
    for (std::size_t i = 0; i < m; ++i) {
        auto parts = algebra::apart_num_den(coeffs[i], ctx);
        if (parts.is_error()) return false;
        nums[i] = parts.value().numerator;
        dens[i] = parts.value().denominator;
        common_den = arena.make<Binary>(BinaryOp::Mul, common_den, dens[i]);
    }
    if (auto s = ctx.simplify(common_den); s.is_ok()) common_den = s.value();

    // scaled_i = num_i · (common_den / den_i), polinomio in base_var.
    std::vector<std::vector<Rational>> scaled_coeffs(m);
    std::size_t max_deg = 0;
    for (std::size_t i = 0; i < m; ++i) {
        ExprPtr factor = arena.make<Binary>(BinaryOp::Div, common_den, dens[i]);
        ExprPtr scaled = arena.make<Binary>(BinaryOp::Mul, nums[i], factor);
        if (auto tog = algebra::together(scaled, ctx); tog.is_ok()) scaled = tog.value();
        if (auto ex = algebra::expand(scaled, ctx); ex.is_ok()) scaled = ex.value();
        if (auto s = ctx.simplify(scaled); s.is_ok()) scaled = s.value();
        auto cq = poly_coeffs_q(scaled, base_var, ctx);
        if (!cq) return false;
        scaled_coeffs[i] = std::move(*cq);
        max_deg = std::max(max_deg, scaled_coeffs[i].size());
    }

    for (std::size_t e = 0; e < max_deg; ++e) {
        std::vector<Rational> row(m, Rational(BigInt(0)));
        bool nonzero = false;
        for (std::size_t i = 0; i < m; ++i) {
            if (e < scaled_coeffs[i].size()) {
                row[i] = scaled_coeffs[i][e];
                if (!row[i].numerator().is_zero()) nonzero = true;
            }
        }
        if (nonzero) rows.push_back(std::move(row));
    }
    return true;
}

}  // namespace

Result<std::optional<ParametricForcingReduction>> reduce_parametric_forcing_poles(
    const std::vector<ExprPtr>& g_vec,
    const Symbol& t,
    const Symbol& base_var,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();
    const std::size_t m = g_vec.size();
    if (m == 0U) return ok(std::optional<ParametricForcingReduction>{});

    // 1. Parte propria di ogni forzante.  Se nessuna ha poli in t non c'è nulla
    //    da fare: il chiamante prosegue sul percorso ordinario.
    std::vector<ProperPart> pp;
    pp.reserve(m);
    bool any_pole = false;
    for (ExprPtr g : g_vec) {
        auto p = proper_part_in_t(g, t, ctx);
        if (!p) return ok(std::optional<ParametricForcingReduction>{});
        if (!p->den_is_trivial) any_pole = true;
        pp.push_back(*p);
    }
    if (!any_pole) return ok(std::optional<ParametricForcingReduction>{});

    // 2. Denominatore comune L delle parti proprie e numeratori riscalati
    //    R_i ∈ k[t] con Σ c_i R_i / L = parte propria della combinazione.
    ExprPtr L = arena.make<IntegerLit>(BigInt(1));
    for (const auto& p : pp) {
        if (p.den_is_trivial) continue;
        auto g_common = algebra::polynomial_gcd(L, p.den, t, ctx);
        if (g_common.is_error()) return ok(std::optional<ParametricForcingReduction>{});
        ExprPtr prod = arena.make<Binary>(BinaryOp::Mul, L, p.den);
        L = arena.make<Binary>(BinaryOp::Div, prod, g_common.value());
        if (auto s = ctx.simplify(L); s.is_ok()) L = s.value();
    }

    std::vector<algebra::PolyExpr> R(m);
    std::size_t max_t_deg = 0;
    for (std::size_t i = 0; i < m; ++i) {
        ExprPtr Ri;
        if (pp[i].den_is_trivial) {
            Ri = arena.make<IntegerLit>(BigInt(0));
        } else {
            ExprPtr factor = arena.make<Binary>(BinaryOp::Div, L, pp[i].den);
            Ri = arena.make<Binary>(BinaryOp::Mul, pp[i].num, factor);
            if (auto tog = algebra::together(Ri, ctx); tog.is_ok()) Ri = tog.value();
            if (auto ex = algebra::expand(Ri, ctx); ex.is_ok()) Ri = ex.value();
            if (auto s = ctx.simplify(Ri); s.is_ok()) Ri = s.value();
        }
        auto poly = algebra::parse_polynomial(Ri, t, ctx);
        if (poly.is_error()) return ok(std::optional<ParametricForcingReduction>{});
        R[i] = poly.value();
        max_t_deg = std::max(max_t_deg, R[i].size());
    }

    // 3. Vincolo: Σ_i c_i·R_i ≡ 0 come polinomio in t, coefficiente per
    //    coefficiente (ogni coefficiente vive in k = Q(base_var)).
    std::vector<std::vector<Rational>> rows;
    for (std::size_t e = 0; e < max_t_deg; ++e) {
        std::vector<ExprPtr> coeffs(m);
        for (std::size_t i = 0; i < m; ++i) {
            ExprPtr c = (e < R[i].size()) ? R[i][e] : ExprPtr{};
            coeffs[i] = c ? c : static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(0)));
        }
        if (!append_k_constraint(coeffs, base_var, rows, ctx))
            return ok(std::optional<ParametricForcingReduction>{});
    }

    // 4. Null space del vincolo → base delle combinazioni ammissibili.
    std::vector<std::vector<Rational>> basis;
    if (rows.empty()) {
        for (std::size_t i = 0; i < m; ++i) {
            std::vector<Rational> e_i(m, Rational(BigInt(0)));
            e_i[i] = Rational(BigInt(1));
            basis.push_back(std::move(e_i));
        }
    } else {
        auto pivots = row_echelon(rows, m);
        basis = null_space_basis(rows, pivots, m);
    }

    // 5. Forzanti ridotte g'_j = Σ_i basis[j][i]·g_i — polinomiali in t per
    //    costruzione (le parti proprie si cancellano nel null space).
    ParametricForcingReduction out;
    out.basis = std::move(basis);
    out.g_reduced.reserve(out.basis.size());
    for (const auto& v : out.basis) {
        std::vector<ExprPtr> terms;
        for (std::size_t i = 0; i < m; ++i) {
            if (v[i].numerator().is_zero()) continue;
            terms.push_back(arena.make<Binary>(BinaryOp::Mul,
                rational_to_expr(v[i], arena), g_vec[i]));
        }
        ExprPtr g_red = terms.empty()
            ? static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(0)))
            : (terms.size() == 1U ? terms.front()
                                  : static_cast<ExprPtr>(arena.make<Sum>(std::move(terms))));
        if (auto tog = algebra::together(g_red, ctx); tog.is_ok()) g_red = tog.value();
        if (auto s = ctx.simplify(g_red); s.is_ok()) g_red = s.value();
        out.g_reduced.push_back(g_red);
    }

    return ok(std::optional<ParametricForcingReduction>{std::move(out)});
}

}  // namespace cas::calculus
