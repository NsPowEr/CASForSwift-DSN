// A29 anti-monolith split from risch_rde_bronstein.cpp (zero logic changes):
// this TU holds the parametric Risch DE tower descent
// solve_risch_de_parametric_field (Bronstein Symbolic Integration I, §7.1);
// the non-parametric solve_risch_de_field/solve_risch_de_general stay in
// risch_rde_bronstein.cpp.

#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <vector>
#include <algorithm>

namespace cas::calculus {

Result<std::vector<ParametricRischDeQSolution>> solve_risch_de_parametric_field(
    ExprPtr f,
    const std::vector<ExprPtr>& g_vec,
    std::size_t ext_idx,
    const DifferentialField& field,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    // Base case: Q(x)
    if (ext_idx == 0U) {
        auto poly_res = solve_risch_de_parametric_q(f, g_vec, field.base_var(), ctx);
        if (poly_res.is_ok()) return poly_res;
        // The Q[x] solver bailed (typically a rational f or g, which the
        // primitive tower descent produces).  Split by f:
        //   f == 0 → rational limited integration (integrate each g_i, split
        //            rational + log/arctan atoms);
        //   f != 0 → full parametric Risch DE over Q(x) via the P/D ansatz.
        // Both are sound by construction (back-substitution verified).
        // (A26 / HC-A26-PRIMITIVE-PARAMQ-RATIONAL.)
        ExprPtr f_s = f;
        if (auto s = ctx.simplify(f); s.is_ok()) f_s = s.value();
        bool f_is_zero = false;
        if (const auto* il = expr_cast<IntegerLit>(f_s)) f_is_zero = il->value.is_zero();
        if (const auto* rl = expr_cast<RationalLit>(f_s)) f_is_zero = rl->numerator.is_zero();
        if (f_is_zero)
            return solve_param_limited_integration_rational_q(g_vec, field.base_var(), ctx);
        auto rat = solve_param_risch_de_rational_q(f, g_vec, field.base_var(), ctx);
        if (rat.is_ok()) return rat;
        return poly_res;
    }

    // Topmost extension at this level
    const auto& ext = field.extensions()[ext_idx - 1U];
    const Symbol& t = ext.t_var;
    const std::size_t m = g_vec.size();

    // HC-A26-PRIMITIVE-PARAMQ-RATIONAL — cancellazione dei poli fra le forzanti.
    //
    // Caso PRIMITIVO (t = log(u), Dt ∈ k): Teorema 5.1.1 ⇒ S = k, e la formula
    // (5.1) dà k⟨t⟩ = k[t].  Quindi ogni soluzione q sta in k[t] e il membro
    // sinistro a·D(q) + b·q è polinomiale in t: di conseguenza Σ_i c_i·g_i deve
    // essere polinomiale in t ANCHE QUANDO le singole g_i hanno poli — i poli
    // devono cancellarsi nella combinazione.  È un vincolo Q-lineare sui c_i.
    //
    // Il bound di denominatore sottostante (D = lcm dei denominatori) assume
    // invece che la soluzione erediti i poli delle forzanti: moltiplica per D e
    // produce f_new = f − D'/D non polinomiale, facendo fallire il descent con
    // "f_new is not polynomial" anche quando una soluzione polinomiale esiste
    // (repro: g = (1/x + 1/(x·t), −1/(x·t)) su t = log x, soluzione y = t con
    // c₀ = c₁ = 1; misurato D = t·x², f_new = −2/x − 1/(t·x)).
    //
    // Riduciamo perciò la famiglia al sottospazio in cui i poli si cancellano:
    // le forzanti ridotte sono polinomiali in t e il descent ordinario funziona.
    // La ricorsione termina in ≤2 livelli — le forzanti ridotte non hanno poli,
    // quindi la seconda chiamata riceve `nullopt` e prosegue diretta.
    //
    // Vincolo di applicabilità (esplicito, non silenzioso): serve f polinomiale
    // in t.  Con f avente poli occorre il ParamRdeNormalDenominator completo
    // (SplitFactor + WeakNormalizer su torre, Bronstein §6.1/§7.1), fuori scope
    // qui: in quel caso si prosegue sul percorso ordinario invariato.
    if (ext.type == ExtensionType::Logarithmic) {
        bool f_is_poly_in_t = algebra::parse_polynomial(f, t, ctx).is_ok();
        if (f_is_poly_in_t) {
            auto red = reduce_parametric_forcing_poles(g_vec, t, field.base_var(), ctx);
            if (red.is_ok() && red.value().has_value()) {
                const auto& R = red.value().value();
                if (R.g_reduced.empty()) {
                    // Nessuna combinazione cancella i poli ⇒ nessuna soluzione
                    // (esito legittimo, non un errore).
                    return ok(std::vector<ParametricRischDeQSolution>{});
                }
                auto sub = solve_risch_de_parametric_field(f, R.g_reduced, ext_idx, field, ctx);
                if (sub.is_error()) return sub;
                // Rimappa le costanti sulla famiglia originale:
                // c_i = Σ_j d_j · basis[j][i].
                std::vector<ParametricRischDeQSolution> out;
                out.reserve(sub.value().size());
                for (const auto& sol : sub.value()) {
                    std::vector<Rational> c(m, Rational(BigInt(0)));
                    for (std::size_t j = 0; j < R.basis.size() && j < sol.c.size(); ++j) {
                        if (sol.c[j].numerator().is_zero()) continue;
                        for (std::size_t i = 0; i < m; ++i)
                            c[i] = c[i] + sol.c[j] * R.basis[j][i];
                    }
                    out.push_back({sol.y, std::move(c)});
                }
                return ok(std::move(out));
            }
        }
    }

    // Denominator bound: LCM of denominators of f and all g_i
    auto f_parts_res = algebra::apart_num_den(f, ctx);
    if (f_parts_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(f_parts_res.error());
    ExprPtr B = f_parts_res.value().denominator;

    ExprPtr Q = arena.make<IntegerLit>(BigInt(1));
    for (ExprPtr g_expr : g_vec) {
        auto g_parts_res = algebra::apart_num_den(g_expr, ctx);
        if (g_parts_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(g_parts_res.error());
        ExprPtr den = g_parts_res.value().denominator;
        auto gcd_den = algebra::polynomial_gcd(Q, den, t, ctx);
        if (gcd_den.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(gcd_den.error());
        ExprPtr Q_den = arena.make<Binary>(BinaryOp::Mul, Q, den);
        Q = arena.make<Binary>(BinaryOp::Div, Q_den, gcd_den.value());
        if (auto s = ctx.simplify(Q); s.is_ok()) Q = s.value();
    }

    auto gcd_qb = algebra::polynomial_gcd(Q, B, t, ctx);
    if (gcd_qb.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(gcd_qb.error());
    ExprPtr QB = arena.make<Binary>(BinaryOp::Mul, Q, B);
    ExprPtr D = arena.make<Binary>(BinaryOp::Div, QB, gcd_qb.value());
    if (auto s = ctx.simplify(D); s.is_ok()) D = s.value();

    // f_new = f - D'/D,  g_new_i = g_i * D
    auto D_prime_res = field.derive(D, ctx);
    if (D_prime_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(D_prime_res.error());
    ExprPtr D_prime = D_prime_res.value();

    ExprPtr D_prime_over_D = arena.make<Binary>(BinaryOp::Div, D_prime, D);
    ExprPtr f_new = arena.make<Binary>(BinaryOp::Sub, f, D_prime_over_D);
    if (auto s = ctx.simplify(f_new); s.is_ok()) f_new = s.value();

    std::vector<ExprPtr> g_new_vec;
    g_new_vec.reserve(m);
    for (ExprPtr g_expr : g_vec) {
        ExprPtr g_new = arena.make<Binary>(BinaryOp::Mul, g_expr, D);
        if (auto s = ctx.simplify(g_new); s.is_ok()) g_new = s.value();
        g_new_vec.push_back(g_new);
    }

    // Parse f_new and all g_new_i as polynomials in t
    auto f_poly_res = algebra::parse_polynomial(f_new, t, ctx);
    if (f_poly_res.is_error()) {
        return make_unimplemented<std::vector<ParametricRischDeQSolution>>(
            "calculus", "solve_risch_de_parametric_field",
            "f_new is not polynomial",
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Risch DE solver: make sure f_new is polynomial");
    }
    const auto& f_poly = f_poly_res.value();

    std::vector<algebra::PolyExpr> g_polys;
    g_polys.reserve(m);
    int dg_max = -1;
    for (ExprPtr g_new : g_new_vec) {
        auto g_poly_res = algebra::parse_polynomial(g_new, t, ctx);
        if (g_poly_res.is_error()) {
            return make_unimplemented<std::vector<ParametricRischDeQSolution>>(
                "calculus", "solve_risch_de_parametric_field",
                "g_new_i is not polynomial",
                cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
                "Risch DE solver: make sure g_new_i is polynomial");
        }
        g_polys.push_back(g_poly_res.value());
        int dg = static_cast<int>(poly_degree(g_poly_res.value()));
        if (is_zero_poly(g_poly_res.value())) dg = -1;
        dg_max = std::max(dg_max, dg);
    }

    int df = static_cast<int>(poly_degree(f_poly));
    if (is_zero_poly(f_poly)) df = -1;

    // Degree bound N — Bronstein §6.3/§7.1 RdeBoundDegreePrim / ParamRdeBoundDegreePrim.
    // Dopo il clearing per D l'equazione è  D(q) + f_new·q = Σ c_i·g_new_i,
    // cioè "a" = 1 ⇒ d_a = 0 (nessun coefficiente su Dq).  Con d_a = 0 la
    // formula è:
    //   d_b > d_a (df > 0)        → n = max(0, d_c − d_b)          (ramo df>0)
    //   d_b ≤ d_a (df ≤ 0)        → n = max(0, d_c − d_a + 1) = d_c + 1
    // Il ramo df == 0 rientra nel SECONDO caso (d_b = d_a = 0), quindi
    // n = dg_max + 1, NON dg_max: la soluzione può avere grado uno in più
    // della forzante quando il termine di grado massimo si cancella contro
    // f_new·q (es. t = log x, f_new = −1/x, g = 1 → q = x·t di grado 1 con
    // forzante di grado 0).  Il vecchio dg_max troncava a 0 e perdeva q.
    //
    // COMPLETEZZA (d_a = 0, monomio log/primitivo — dimostrazione, non gap).
    // Il branch di cancellazione di testa §6.3.3 scatta solo se
    // α = −f_new = Dz/z è una derivata logaritmica (z ∈ k*).  Sostituendo
    // q = z·h nell'equazione D(q) + f_new·q = c si ha
    //   D(z·h) + f_new·z·h = z·Dh + h·(Dz + f_new·z) = z·Dh   (Dz = −f_new·z),
    // cioè Dh = c/z.  Poiché deg_t(c/z) ≤ dg_max (z ∈ k, grado 0 in t) e
    // l'integrazione primitiva alza il grado in t di ESATTAMENTE uno,
    // deg(q) = deg(h) ≤ dg_max + 1.  Se invece α NON è derivata logaritmica,
    // l'eq. di testa omogenea D(q_n) + f_new·q_n = 0 non ha soluzione k non
    // nulla per n > dg_max, quindi deg(q) ≤ dg_max.  In entrambi i casi
    // N = dg_max + 1 è un bound COMPLETO — l'incremento di cancellazione è già
    // incluso.  Nessuna soluzione di grado dg_max+2 può esistere: richiederebbe
    // η = Dt derivata di un elemento di k, ma ∫η = log(u) ∉ k (è ciò che rende
    // t trascendente).  Il "further increment oltre dg_max+1" temuto dal ledger
    // HC-A26-RDEBOUND-CANCEL-GAP è il ramo generale d_a > 0; con d_a = 0 (a = 1
    // dopo il clearing per D) i due sotto-rami collassano in ≤ dg_max+1.  Non è
    // un gap.  (Sovrastimare N sarebbe comunque SOUND: candidati extra
    // verificati per back-substitution.)
    int N = 0;
    if (ext.type == ExtensionType::Logarithmic) {
        // d_a = 0: d_b>d_a (df>0) → max(0, d_c−d_b); d_b≤d_a (df≤0) → d_c+1.
        N = (df > 0) ? std::max(0, dg_max - df) : dg_max + 1;
    } else { // Exponential
        if (df > 0) {
            N = std::max(0, dg_max - df);
        } else if (df == 0) {
            N = std::max(dg_max, 0);
        } else {
            N = dg_max + 1;
        }
    }
    if (N < 0) N = 0;

    // df > 0 — non-cancellation parametric PolyRischDE (A1).  For the log/exp
    // monomials here, deg_t(f_new) > 0 ⇒ deg(b) > max(0, δ(t)−1), i.e. the
    // "deg(b) is Large Enough" case: Bronstein Symbolic Integration I §7.1,
    // ParamPolyRischDENoCancel1.  Solve  D(q) + f_new·q = Σ c_i·g_new_i  in
    // K[t], then divide the polynomial solutions q by D (as the df≤0 tail does).
    // Sound: solve_param_poly_risch_de_nocancel1 verifies each candidate by
    // field back-substitution.  The residual constant system is solved for any
    // tower K ⊇ Q(x) via ConstantSystem (Bronstein §7.1, Lemma 7.1.2).
    if (df > 0) {
        auto nc = solve_param_poly_risch_de_nocancel1(f_new, g_new_vec, N, t, field, ctx);
        if (nc.is_error()) return nc;
        std::vector<ParametricRischDeQSolution> out_sols;
        out_sols.reserve(nc.value().size());
        for (auto& sol : nc.value()) {
            ExprPtr y = arena.make<Binary>(BinaryOp::Div, sol.y, D);
            auto y_tog = algebra::together(y, ctx);
            if (y_tog.is_ok()) {
                if (auto s = ctx.simplify(y_tog.value()); s.is_ok()) y = s.value();
            }
            out_sols.push_back({y, std::move(sol.c)});
        }
        return ok(std::move(out_sols));
    }

    // df <= 0. f_new is f_0 in lower field.
    ExprPtr f_0 = leading_coefficient(f_poly);
    if (!f_0) f_0 = arena.make<IntegerLit>(BigInt(0));

    // We start at i = N.  Bind a const& so the bounds-checked const operator[]
    // is selected (the non-const overload is unchecked, std-vector style):
    // N = dg_max+1 routinely exceeds deg(g_s), and OOB must read as the zero
    // coefficient (nullptr → 0), not a heap-buffer-overflow.
    std::vector<ExprPtr> H_vec(m);
    for (std::size_t s = 0; s < m; ++s) {
        const algebra::PolyExpr& gp = g_polys[s];
        ExprPtr coeff = gp[static_cast<std::size_t>(N)];
        H_vec[s] = coeff ? coeff : arena.make<IntegerLit>(BigInt(0));
    }

    auto solve_recursive = [&](auto& self, int i, const std::vector<ExprPtr>& H_current) -> Result<std::vector<ParametricRischDeQSolution>> {
        ExprPtr F_eff = f_0;
        if (ext.type == ExtensionType::Exponential && i > 0) {
            auto u_prime_res = field.derive(ext.argument, ctx);
            if (u_prime_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(u_prime_res.error());
            ExprPtr i_coef = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i)));
            ExprPtr term = arena.make<Binary>(BinaryOp::Mul, i_coef, u_prime_res.value());
            F_eff = arena.make<Binary>(BinaryOp::Add, term, f_0);
            if (auto s = ctx.simplify(F_eff); s.is_ok()) F_eff = s.value();
        }

        auto sols_res = solve_risch_de_parametric_field(F_eff, H_current, ext_idx - 1U, field, ctx);
        if (sols_res.is_error()) return sols_res;

        const auto& sols = sols_res.value();
        if (i == 0) {
            return ok(sols);
        }

        const std::size_t num_sols = sols.size();
        std::vector<ExprPtr> H_next(num_sols);

        ExprPtr theta_prime = nullptr;
        if (ext.type == ExtensionType::Logarithmic) {
            auto theta_prime_res = field.derive(arena.make<Symbol>(t.name), ctx);
            if (theta_prime_res.is_error()) return fail<std::vector<ParametricRischDeQSolution>>(theta_prime_res.error());
            theta_prime = theta_prime_res.value();
        }

        for (std::size_t r = 0; r < num_sols; ++r) {
            const auto& sol = sols[r];
            std::vector<ExprPtr> terms;
            for (std::size_t s = 0; s < H_current.size(); ++s) {
                const algebra::PolyExpr& gp = g_polys[s];  // const& → bounds-checked operator[]
                ExprPtr g_val = gp[static_cast<std::size_t>(i - 1)];
                if (!g_val) g_val = arena.make<IntegerLit>(BigInt(0));

                if (sol.c[s].numerator().is_zero()) continue;
                ExprPtr c_val = arena.make<RationalLit>(sol.c[s].numerator(), sol.c[s].denominator());
                terms.push_back(arena.make<Binary>(BinaryOp::Mul, c_val, g_val));
            }
            ExprPtr sum_g = terms.empty() ? arena.make<IntegerLit>(BigInt(0)) :
                            (terms.size() == 1U ? terms[0] : arena.make<Sum>(std::move(terms)));

            if (ext.type == ExtensionType::Logarithmic) {
                // PRIMITIVE-CASE GAP (HC-A26-PRIMITIVE-PARAMQ-RATIONAL): the
                // correction i·y·θ' with θ' = D(t) = u'/u re-introduces a
                // denominator, so H_next can be rational in the lower field.
                // The base case solve_risch_de_parametric_q is polynomial-only
                // (Q[x]) and will return a diagnostic Unimplemented for such a
                // forcing.  Completing this needs ParamRischDE over Q(x)
                // (weak-normalizer + denominator bound, Bronstein §5.12/§6.5).
                ExprPtr i_coef = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i)));
                ExprPtr corr = arena.make<Product>(std::vector<ExprPtr>{i_coef, sol.y, theta_prime});
                ExprPtr val = arena.make<Binary>(BinaryOp::Sub, sum_g, corr);
                if (auto s = ctx.simplify(val); s.is_ok()) val = s.value();
                H_next[r] = val;
            } else {
                if (auto s = ctx.simplify(sum_g); s.is_ok()) sum_g = s.value();
                H_next[r] = sum_g;
            }
        }

        auto next_sols_res = self(self, i - 1, H_next);
        if (next_sols_res.is_error()) return next_sols_res;

        const auto& next_sols = next_sols_res.value();
        std::vector<ParametricRischDeQSolution> combined_sols;
        combined_sols.reserve(next_sols.size());

        for (const auto& next_sol : next_sols) {
            std::vector<ExprPtr> y_i_terms;
            for (std::size_t r = 0; r < num_sols; ++r) {
                const Rational& cr = next_sol.c[r];
                if (cr.numerator().is_zero()) continue;
                ExprPtr cr_e = arena.make<RationalLit>(cr.numerator(), cr.denominator());
                y_i_terms.push_back(arena.make<Binary>(BinaryOp::Mul, cr_e, sols[r].y));
            }
            ExprPtr y_i = y_i_terms.empty() ? arena.make<IntegerLit>(BigInt(0)) :
                          (y_i_terms.size() == 1U ? y_i_terms[0] : arena.make<Sum>(std::move(y_i_terms)));
            if (auto s = ctx.simplify(y_i); s.is_ok()) y_i = s.value();

            ExprPtr ti = (i == 1) ? static_cast<ExprPtr>(arena.make<Symbol>(t.name)) :
                         static_cast<ExprPtr>(arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(t.name),
                             arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(i)))));
            ExprPtr term_ti = arena.make<Binary>(BinaryOp::Mul, y_i, ti);
            ExprPtr y_comb = arena.make<Binary>(BinaryOp::Add, next_sol.y, term_ti);
            if (auto s = ctx.simplify(y_comb); s.is_ok()) y_comb = s.value();

            std::vector<Rational> c_comb(m, Rational(BigInt(0)));
            for (std::size_t r = 0; r < num_sols; ++r) {
                const Rational& cr = next_sol.c[r];
                for (std::size_t s = 0; s < m; ++s) {
                    c_comb[s] = c_comb[s] + cr * sols[r].c[s];
                }
            }

            combined_sols.push_back({y_comb, std::move(c_comb)});
        }

        return ok(std::move(combined_sols));
    };

    auto sols_res = solve_recursive(solve_recursive, N, H_vec);
    if (sols_res.is_error()) return sols_res;

    std::vector<ParametricRischDeQSolution> out_sols;
    out_sols.reserve(sols_res.value().size());
    for (auto& sol : sols_res.value()) {
        ExprPtr y = arena.make<Binary>(BinaryOp::Div, sol.y, D);
        auto y_tog = algebra::together(y, ctx);
        if (y_tog.is_ok()) {
            if (auto s = ctx.simplify(y_tog.value()); s.is_ok()) y = s.value();
        }
        out_sols.push_back({y, std::move(sol.c)});
    }

    return ok(std::move(out_sols));
}

} // namespace cas::calculus
