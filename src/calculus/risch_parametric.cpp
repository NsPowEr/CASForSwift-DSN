// F5.1 / B9-Task#21 — Risch Bronstein cap.7 Parametric Problems.
//
// Implementa:
//   1. solve_risch_de_parametric_q
//      Risolve  y' + f·y = Σ_i c_i · g_i  per y ∈ Q[var] e c_i ∈ Q.
//      Algoritmo: sistema lineare omogeneo Σ_j eq_j(y_k, c_i) = 0; calcolo
//      base dello spazio nullo via Gauss-Jordan parziale + back-substitution
//      sulle colonne libere.
//
//   2. limited_integration_q
//      Decomposizione  f = g' + Σ_i c_i · h_i  per f, h_i ∈ Q[var]; nel caso
//      polinomiale puro la parte g è l'antiderivata esatta di (f - Σ c_i·h_i)
//      con {c_i} scelti per annullare componenti su {h_i}.
//
// Riferimento: Manuel Bronstein, "Symbolic Integration I — Transcendental
// Functions" (Springer, 2nd ed. 2005), §7.1 "Parametric Risch DE", §7.2
// "Limited Integration".

#include "calculus_internal.hpp"
#include "risch_parametric_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

using detail::as_rational;
using detail::null_space_basis;
using detail::row_echelon;

// Estrae vettore di coefficienti razionali da un PolyExpr ExprPtr-coefficient.
[[nodiscard]] std::optional<std::vector<Rational>>
extract_q_coeffs(const algebra::PolyExpr& poly) {
    std::vector<Rational> out;
    out.reserve(poly.size());
    for (ExprPtr c : poly.coefficients()) {
        if (auto r = as_rational(c)) {
            out.push_back(*r);
        } else {
            return std::nullopt;
        }
    }
    return out;
}

// Ricostruisce Σ_k c_k · var^k come ExprPtr.
[[nodiscard]] ExprPtr build_poly(
    AstArena& arena, const std::vector<Rational>& coeffs, const Symbol& var) {
    std::vector<ExprPtr> terms;
    terms.reserve(coeffs.size());
    ExprPtr var_sym = arena.make<Symbol>(var);
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        const Rational& c = coeffs[k];
        if (c.numerator().is_zero()) continue;
        ExprPtr c_expr = (c.denominator() == BigInt(1))
            ? static_cast<ExprPtr>(arena.make<IntegerLit>(c.numerator()))
            : static_cast<ExprPtr>(arena.make<RationalLit>(c.numerator(), c.denominator()));
        if (k == 0U) { terms.push_back(c_expr); continue; }
        ExprPtr xk = (k == 1U)
            ? var_sym
            : static_cast<ExprPtr>(arena.make<Binary>(BinaryOp::Pow,
                  var_sym,
                  arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k)))));
        if (c.numerator() == BigInt(1) && c.denominator() == BigInt(1)) {
            terms.push_back(xk);
        } else {
            terms.push_back(arena.make<Product>(std::vector<ExprPtr>{c_expr, xk}));
        }
    }
    if (terms.empty()) return arena.make<IntegerLit>(BigInt(0));
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

// Strip trailing zeros (normalizza il grado).
void strip_trailing(std::vector<Rational>& v) {
    while (!v.empty() && v.back().numerator().is_zero()) v.pop_back();
}

// Calcola il bound di grado e_max per y dato f e g_max (degree massimo di
// qualsiasi g_i nella combinazione lineare).  Replica lo schema di
// solve_risch_de_poly_q ma applicato al grado massimo delle forzanti, perché
// la combinazione Σc_i·g_i ha grado al più max(deg g_i).
[[nodiscard]] int compute_degree_bound(int deg_f, int deg_g_max) {
    if (deg_g_max < 0) return 0;        // tutte g_i nulle → y può essere 0.
    if (deg_f < 0) return deg_g_max + 1;  // f ≡ 0  →  y' = g  →  deg(y) = d+1.
    if (deg_f == 0) return deg_g_max;     // f costante non nulla.
    // deg(f·y) = deg_f + e domina deg(y') = e - 1; uguagliando a deg_g_max:
    //   deg_f + e = deg_g_max  →  e = deg_g_max - deg_f.
    return std::max(0, deg_g_max - deg_f);
}

}  // namespace

Result<std::vector<ParametricRischDeQSolution>>
solve_risch_de_parametric_q(
    ExprPtr f_expr,
    const std::vector<ExprPtr>& g_vec,
    const Symbol& var,
    symbolic::CASContext& ctx) {

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<std::vector<ParametricRischDeQSolution>>(
            "calculus", "solve_risch_de_parametric_q",
            msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Parametric Risch DE: extend coefficient extraction or upgrade to "
            "transcendental tower (Bronstein cap.8)",
            "F0.8");
    };

    // 1. Parsing polinomi.
    auto f_poly_res = algebra::parse_polynomial(f_expr, var, ctx);
    if (f_poly_res.is_error())
        return fail_unimpl("f is not a polynomial in the base variable");
    auto f_opt = extract_q_coeffs(f_poly_res.value());
    if (!f_opt.has_value())
        return fail_unimpl("f has non-rational coefficients");
    std::vector<Rational> f = std::move(*f_opt);
    strip_trailing(f);
    const int deg_f = static_cast<int>(f.size()) - 1;

    const std::size_t m = g_vec.size();
    if (m == 0U)
        return fail_unimpl("empty forcing basis g_vec");
    std::vector<std::vector<Rational>> g_list;
    g_list.reserve(m);
    int deg_g_max = -1;
    for (ExprPtr g_expr : g_vec) {
        auto g_poly_res = algebra::parse_polynomial(g_expr, var, ctx);
        if (g_poly_res.is_error())
            return fail_unimpl("g_i is not a polynomial in the base variable");
        auto g_opt = extract_q_coeffs(g_poly_res.value());
        if (!g_opt.has_value())
            return fail_unimpl("g_i has non-rational coefficients");
        std::vector<Rational> g = std::move(*g_opt);
        strip_trailing(g);
        deg_g_max = std::max(deg_g_max,
            static_cast<int>(g.size()) - 1);
        g_list.push_back(std::move(g));
    }

    // 2. Bound di grado e_max per y.
    const int e_max = compute_degree_bound(deg_f, deg_g_max);
    const std::size_t n_y  = static_cast<std::size_t>(e_max) + 1U;
    const std::size_t n_eq = static_cast<std::size_t>(std::max(deg_g_max, e_max + std::max(0, deg_f))) + 1U;
    const std::size_t n_unk = n_y + m;
    const std::size_t n_cols = n_unk;  // sistema omogeneo, no RHS column.

    // 3. Costruzione matrice.  Colonne: y_0,...,y_e, c_1,...,c_m.
    //    Equazione per coefficiente x^j:
    //       (j+1)·y_{j+1}                            (contributo di y')
    //     + Σ_k f_k · y_{j-k}                        (contributo di f·y)
    //     − Σ_i (g_i)_j · c_i                        (rhs spostato a sinistra)
    //     = 0
    std::vector<std::vector<Rational>> M(n_eq,
        std::vector<Rational>(n_cols, Rational(BigInt(0))));
    for (std::size_t j = 0; j < n_eq; ++j) {
        // y' contribution.
        if (j + 1U < n_y) {
            M[j][j + 1U] = M[j][j + 1U] +
                Rational(BigInt(static_cast<std::int64_t>(j + 1U)));
        }
        // f·y contribution.
        for (std::size_t k = 0; k <= static_cast<std::size_t>(std::max(0, deg_f)); ++k) {
            if (deg_f < 0) break;
            if (k > j) break;
            std::size_t idx = j - k;
            if (idx >= n_y) continue;
            M[j][idx] = M[j][idx] + f[k];
        }
        // -Σ c_i · g_i contribution.
        for (std::size_t i = 0; i < m; ++i) {
            const auto& g = g_list[i];
            if (j < g.size()) M[j][n_y + i] = M[j][n_y + i] - g[j];
        }
    }

    // 4. Row echelon + null space basis.
    auto pivots = row_echelon(M, n_cols);
    auto basis = null_space_basis(M, pivots, n_cols);

    // 5. Trasforma ogni vettore base in ParametricRischDeQSolution.
    std::vector<ParametricRischDeQSolution> out;
    out.reserve(basis.size());
    for (const auto& v : basis) {
        ParametricRischDeQSolution sol;
        std::vector<Rational> y_coeffs(n_y);
        for (std::size_t k = 0; k < n_y; ++k) y_coeffs[k] = v[k];
        strip_trailing(y_coeffs);
        sol.y = build_poly(ctx.arena(), y_coeffs, var);
        sol.c.reserve(m);
        for (std::size_t i = 0; i < m; ++i) sol.c.push_back(v[n_y + i]);
        out.push_back(std::move(sol));
    }

    return ok(std::move(out));
}

Result<LimitedIntegrationQSolution>
limited_integration_q(
    ExprPtr f_expr,
    const std::vector<ExprPtr>& residue_basis,
    const Symbol& var,
    symbolic::CASContext& ctx) {

    // Limited integration in Q[x] è il caso speciale parametric con
    //   y' + 0·y = f - Σ c_i · h_i
    // dove 0·y = 0 (f scalare nel termine f del DE).  Equivalente: cerchiamo
    // g, c tale che f - Σ c_i · h_i sia integrabile esattamente in Q[x] e g
    // sia la sua antiderivata polinomiale.
    //
    // Per polinomi puri Q[x], tutto è integrabile esattamente; c_i ottimali
    // sono tutti zero, g = ∫f.  L'utilità della routine emerge quando f e h_i
    // hanno componenti non-integrabili (extension non-Q[x]); qui implementiamo
    // il caso polinomiale come building block consistente con cap.7.

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<LimitedIntegrationQSolution>(
            "calculus", "limited_integration_q",
            msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Limited integration: lift to transcendental field for non-polynomial input",
            "F0.8");
    };

    auto f_poly_res = algebra::parse_polynomial(f_expr, var, ctx);
    if (f_poly_res.is_error())
        return fail_unimpl("f is not a polynomial in the base variable");
    auto f_opt = extract_q_coeffs(f_poly_res.value());
    if (!f_opt.has_value())
        return fail_unimpl("f has non-rational coefficients");
    std::vector<Rational> f = std::move(*f_opt);
    strip_trailing(f);

    // Caso Q[x] puro: antiderivata diretta, residuo nullo.
    std::vector<Rational> g_coeffs(f.size() + 1, Rational(BigInt(0)));
    for (std::size_t k = 0; k < f.size(); ++k) {
        g_coeffs[k + 1] = f[k] /
            Rational(BigInt(static_cast<std::int64_t>(k + 1)));
    }
    LimitedIntegrationQSolution sol;
    sol.g = build_poly(ctx.arena(), g_coeffs, var);
    sol.c.assign(residue_basis.size(), Rational(BigInt(0)));
    return ok(std::move(sol));
}

}  // namespace cas::calculus
