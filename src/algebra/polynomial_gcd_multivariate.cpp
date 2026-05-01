#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <map>
#include <vector>

namespace cas::algebra {

namespace {

// Converte MultivariatePolynomial (variabile singola x) in IntPoly
Result<IntPoly> multivariate_single_var_to_intpoly(const MultivariatePolynomial& poly, const Symbol& var) {
    if (poly.is_zero()) return ok(IntPoly{});

    std::map<std::size_t, BigInt> coeff_map;
    std::size_t max_deg = 0;

    for (const auto& term : poly.terms()) {
        std::size_t deg = 0;
        for (const auto& factor : term.factors) {
            if (factor.first.name == var.name) {
                deg = factor.second;
            } else {
                // Termini con altre variabili non supportati qui
                return fail<IntPoly>(make_error(CASErrorKind::Unimplemented,
                    "gcd_multivariate: atteso polinomio univariato dopo valutazione"));
            }
        }
        coeff_map[deg] += term.coefficient;
        if (deg > max_deg) max_deg = deg;
    }

    IntPoly result;
    result.resize(max_deg + 1, BigInt(0));
    for (const auto& [deg, coeff] : coeff_map) {
        result[deg] = coeff;
    }
    // Normalizza: rimuovi zeri in coda
    result.normalize([](const BigInt& v) { return v.is_zero(); });
    return ok(std::move(result));
}

// Normalizza IntPoly: rende il coefficiente direttivo positivo
void normalize_sign(IntPoly& poly) {
    if (!poly.empty() && poly.leading_coeff().is_negative()) {
        for (auto& c : poly.coefficients()) c = -c;
    }
}

// Ricava i gradi di una variabile in un MultivariatePolynomial
std::size_t degree_in_var(const MultivariatePolynomial& poly, const Symbol& var) {
    std::size_t deg = 0;
    for (const auto& term : poly.terms()) {
        for (const auto& factor : term.factors) {
            if (factor.first.name == var.name) {
                if (factor.second > deg) deg = factor.second;
            }
        }
    }
    return deg;
}

// Interpola un vettore di valori BigInt ai punti 1, 2, ..., n usando Lagrange
// Restituisce coefficienti [a0, a1, ..., a_{n-1}] di a0 + a1*y + ... + a_{n-1}*y^{n-1}
// ma con denominatori (tutto su Rational)
Result<std::vector<Rational>> lagrange_interpolate(const std::vector<BigInt>& values) {
    const std::size_t n = values.size();
    if (n == 0) return ok(std::vector<Rational>{});

    // Costruisce il polinomio di Lagrange come vettore di coefficienti Rational
    // p(y) = Σ_i v[i] * L_i(y),  L_i(y) = Π_{j≠i} (y - (j+1)) / ((i+1) - (j+1))
    std::vector<Rational> result(n, Rational(BigInt(0)));

    for (std::size_t i = 0; i < n; ++i) {
        if (values[i].is_zero()) continue;

        // Calcola L_i(y) come polinomio di grado n-1
        // Inizia con [1] e moltiplica per ogni (y - (j+1)) / ((i+1) - (j+1))
        std::vector<Rational> li(1, Rational(BigInt(1)));

        for (std::size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            // Dividi per scalare ((i+1) - (j+1)) = i - j
            long long denom = static_cast<long long>(i) - static_cast<long long>(j);
            Rational scale(BigInt(1), BigInt(denom > 0 ? denom : -denom));
            if (denom < 0) scale = -scale;

            // Moltiplica li per (y - (j+1)) = y - (j+1)
            // Shift: nuovi coefficienti per li * y
            std::vector<Rational> new_li(li.size() + 1, Rational(BigInt(0)));
            long long eval_pt = static_cast<long long>(j) + 1LL;

            for (std::size_t k = 0; k < li.size(); ++k) {
                new_li[k + 1] = new_li[k + 1] + li[k] * scale;
                new_li[k] = new_li[k] - li[k] * scale * Rational(BigInt(eval_pt));
            }
            li = std::move(new_li);
        }

        // Somma v[i] * L_i(y) al risultato
        Rational vi(values[i]);
        if (li.size() > result.size()) result.resize(li.size(), Rational(BigInt(0)));
        for (std::size_t k = 0; k < li.size(); ++k) {
            result[k] = result[k] + vi * li[k];
        }
    }

    // Rimuovi zeri in coda
    while (!result.empty() && result.back().numerator().is_zero()) {
        result.pop_back();
    }

    return ok(std::move(result));
}

// Verifica se un MultivariatePolynomial G divide P (per valutazione in più punti)
bool multivariate_divides_check(const MultivariatePolynomial& P,
                                const MultivariatePolynomial& G,
                                const std::vector<Symbol>& all_vars,
                                symbolic::CASContext& ctx) {
    // Valuta ad un punto lontano dall'origine per evitare falsi positivi
    std::vector<long long> test_pts = {2LL, 3LL, 5LL};

    for (long long pt : test_pts) {
        // Valuta tutte le variabili a pt
        BigInt p_val(0), g_val(0);
        auto eval_poly = [&](const MultivariatePolynomial& poly) -> BigInt {
            BigInt val(0);
            for (const auto& term : poly.terms()) {
                BigInt term_v = term.coefficient;
                for (const auto& var_sym : all_vars) {
                    // Trova esponente di questa variabile nel termine
                    unsigned int exp = 0;
                    for (const auto& factor : term.factors) {
                        if (factor.first.name == var_sym.name) { exp = factor.second; break; }
                    }
                    if (exp > 0) {
                        for (unsigned int e = 0; e < exp; ++e) term_v = term_v * BigInt(pt);
                    }
                }
                val += term_v;
            }
            return val;
        };

        p_val = eval_poly(P);
        g_val = eval_poly(G);

        if (g_val.is_zero()) continue; // Punto di valutazione degenere, skippa
        if (!(p_val % g_val).is_zero()) return false;
    }
    return true;
}

} // namespace

// Algoritmo: GCD multivariato via evaluation-interpolation
// Funziona per polinomi in 2 variabili con coefficienti interi.
Result<MultivariatePolynomial> gcd_multivariate_eval_interp(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx) {

    if (P.is_zero()) return ok(Q);
    if (Q.is_zero()) return ok(P);

    // Raccogli tutte le variabili
    auto vars_p = P.variables();
    auto vars_q = Q.variables();
    std::vector<Symbol> all_vars = vars_p;
    for (const auto& v : vars_q) {
        bool found = false;
        for (const auto& u : all_vars) if (u.name == v.name) { found = true; break; }
        if (!found) all_vars.push_back(v);
    }
    std::sort(all_vars.begin(), all_vars.end(), [](const Symbol& a, const Symbol& b) {
        return a.name < b.name;
    });

    // Caso univariato: usa gcd_integer_poly_primitive direttamente
    if (all_vars.size() == 1U) {
        auto p_int = multivariate_single_var_to_intpoly(P, all_vars[0]);
        if (p_int.is_error()) return fail<MultivariatePolynomial>(p_int.error());
        auto q_int = multivariate_single_var_to_intpoly(Q, all_vars[0]);
        if (q_int.is_error()) return fail<MultivariatePolynomial>(q_int.error());

        auto g_int = gcd_integer_poly_primitive(std::move(p_int.value()), std::move(q_int.value()));
        normalize_sign(g_int);

        // Converti IntPoly → MultivariatePolynomial
        std::vector<MultivariateTerm> terms;
        for (std::size_t deg = 0; deg < g_int.size(); ++deg) {
            if (!g_int[deg].is_zero()) {
                std::vector<std::pair<Symbol, unsigned int>> factors;
                if (deg > 0) factors.push_back({all_vars[0], static_cast<unsigned int>(deg)});
                terms.push_back({.coefficient = g_int[deg], .factors = std::move(factors)});
            }
        }
        return ok(MultivariatePolynomial(std::move(terms)));
    }

    // Caso bivariato: x = all_vars[0], y = all_vars[1]
    if (all_vars.size() != 2U) {
        return fail<MultivariatePolynomial>(make_error(CASErrorKind::Unimplemented,
            "polynomial_gcd_multivariate: supporta al massimo 2 variabili"));
    }

    const Symbol& x_var = all_vars[0];
    const Symbol& y_var = all_vars[1];

    // Grado massimo del GCD in y: min(deg_y(P), deg_y(Q))
    std::size_t deg_y_p = degree_in_var(P, y_var);
    std::size_t deg_y_q = degree_in_var(Q, y_var);
    std::size_t deg_y_gcd_bound = std::min(deg_y_p, deg_y_q);

    // Grado massimo di GCD in x: min(deg_x(P), deg_x(Q))
    std::size_t deg_x_p = degree_in_var(P, x_var);
    std::size_t deg_x_q = degree_in_var(Q, x_var);
    std::size_t deg_x_gcd_bound = std::min(deg_x_p, deg_x_q);

    // Numero di punti di valutazione necessari
    std::size_t n_points = deg_y_gcd_bound + 2; // +1 per determinare il grado, +1 di sicurezza

    // Raccoglie i GCD univariati ad ogni punto di valutazione
    // Indicizzati: gcd_at[t] = IntPoly del GCD a y = t+1
    std::vector<IntPoly> gcd_at(n_points);
    std::size_t min_gcd_deg = deg_x_gcd_bound + 1; // inizializza alto

    for (std::size_t i = 0; i < n_points; ++i) {
        long long t = static_cast<long long>(i) + 1LL;
        ExprPtr t_expr = ctx.arena().make<IntegerLit>(BigInt(t));

        auto P_t = P.evaluate_at(y_var, t_expr);
        if (P_t.is_error()) return fail<MultivariatePolynomial>(P_t.error());
        auto Q_t = Q.evaluate_at(y_var, t_expr);
        if (Q_t.is_error()) return fail<MultivariatePolynomial>(Q_t.error());

        auto p_int = multivariate_single_var_to_intpoly(P_t.value(), x_var);
        if (p_int.is_error()) return fail<MultivariatePolynomial>(p_int.error());
        auto q_int = multivariate_single_var_to_intpoly(Q_t.value(), x_var);
        if (q_int.is_error()) return fail<MultivariatePolynomial>(q_int.error());

        auto h_t = gcd_integer_poly_primitive(std::move(p_int.value()), std::move(q_int.value()));
        normalize_sign(h_t);
        gcd_at[i] = std::move(h_t);

        if (!gcd_at[i].empty() && gcd_at[i].degree() < min_gcd_deg) {
            min_gcd_deg = gcd_at[i].degree();
        }
    }

    // Filtra: solo punti con GCD di grado == min_gcd_deg
    // (punti "sfortunati" danno GCD di grado minore, che è quello corretto)

    // Interpola: per ogni grado k = 0..min_gcd_deg, il coefficiente di x^k
    // è un polinomio in y di grado ≤ deg_y_gcd_bound.
    std::vector<MultivariateTerm> result_terms;

    for (std::size_t k = 0; k <= min_gcd_deg; ++k) {
        // Raccoglie i valori del coefficiente di x^k ai vari punti di y
        std::vector<BigInt> coeff_values(n_points, BigInt(0));
        for (std::size_t i = 0; i < n_points; ++i) {
            if (!gcd_at[i].empty() && k < gcd_at[i].size()) {
                coeff_values[i] = gcd_at[i][k];
            }
        }

        // Interpola la sequenza coeff_values con y = 1, 2, ..., n_points
        auto interp = lagrange_interpolate(coeff_values);
        if (interp.is_error()) return fail<MultivariatePolynomial>(interp.error());

        const auto& poly_in_y = interp.value();

        // Costruisce i termini del risultato per grado x=k
        for (std::size_t yd = 0; yd < poly_in_y.size(); ++yd) {
            const Rational& c = poly_in_y[yd];
            if (c.numerator().is_zero()) continue;

            // Il coefficiente deve essere intero (denominatore = 1)
            if (c.denominator() != BigInt(1)) {
                // Coefficiente razionale nel GCD → problema con il contenuto
                // Approssima: arrotonda, oppure fallisci
                return fail<MultivariatePolynomial>(make_error(CASErrorKind::InternalError,
                    "gcd_multivariate: coefficiente non intero nell'interpolazione"));
            }

            std::vector<std::pair<Symbol, unsigned int>> factors;
            if (k > 0) factors.push_back({x_var, static_cast<unsigned int>(k)});
            if (yd > 0) factors.push_back({y_var, static_cast<unsigned int>(yd)});

            result_terms.push_back({.coefficient = c.numerator(), .factors = std::move(factors)});
        }
    }

    MultivariatePolynomial G(std::move(result_terms));

    // Verifica: G deve dividere P e Q
    if (!multivariate_divides_check(P, G, all_vars, ctx) ||
        !multivariate_divides_check(Q, G, all_vars, ctx)) {
        return fail<MultivariatePolynomial>(make_error(CASErrorKind::InternalError,
            "gcd_multivariate: candidato GCD non verifica la divisibilità"));
    }

    return ok(std::move(G));
}

// API pubblica: GCD multivariato tra ExprPtr
Result<ExprPtr> polynomial_gcd_multivariate(ExprPtr p, ExprPtr q, symbolic::CASContext& ctx) {
    auto expanded_p = expand_expr_impl(p, ctx);
    if (expanded_p.is_error()) return expanded_p;
    auto expanded_q = expand_expr_impl(q, ctx);
    if (expanded_q.is_error()) return expanded_q;

    auto P = parse_multivariate_polynomial(expanded_p.value(), ctx);
    if (P.is_error()) return fail<ExprPtr>(P.error());
    auto Q = parse_multivariate_polynomial(expanded_q.value(), ctx);
    if (Q.is_error()) return fail<ExprPtr>(Q.error());

    // Prova prima l'algoritmo evaluation-interpolation
    auto G = gcd_multivariate_eval_interp(P.value(), Q.value(), ctx);
    if (G.is_error()) {
        // Fallback: gcd_heuristic (meno affidabile)
        G = gcd_heuristic(P.value(), Q.value());
        if (G.is_error()) return fail<ExprPtr>(G.error());
    }

    auto g_expr = multivariate_to_expr(G.value(), ctx);
    if (g_expr.is_error()) return g_expr;

    return ctx.simplify(g_expr.value());
}

} // namespace cas::algebra
