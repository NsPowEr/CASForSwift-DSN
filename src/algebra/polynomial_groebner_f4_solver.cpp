#include "polynomial_groebner_f4.hpp"

#include "algebra_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"

#include <optional>

namespace cas::algebra {

Result<std::vector<std::vector<ExprPtr>>> solve_nonlinear_system_f4(
    const std::vector<ExprPtr>& equations,
    const std::vector<Symbol>& variables,
    symbolic::CASContext& ctx) {

    const std::size_t n = variables.size();
    const MonomialOrder order = MonomialOrder::Lex;
    std::vector<PolyF4> F;
    for (ExprPtr eq : equations) {
        auto p = expr_to_f4(eq, variables, ctx);
        if (p.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(p.error());
        F.push_back(p.value());
    }

    auto groebner = f4_groebner(std::move(F), order, &ctx);
    if (groebner.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(groebner.error());
    auto G = std::move(groebner.value());
    inter_reduce(G, order);

    if (G.empty()) return ok(std::vector<std::vector<ExprPtr>>{});

    std::optional<PolyF4> pure_last;
    for (const auto& g : G) {
        bool only_last = true;
        for (const auto& [mon, coeff] : g.terms) {
            (void)coeff;
            for (std::size_t i = 0; i < n - 1U; ++i) {
                if (mon[i] > 0U) {
                    only_last = false;
                    break;
                }
            }
            if (!only_last) break;
        }
        if (only_last) {
            pure_last = g;
            break;
        }
    }

    if (!pure_last) {
        return fail<std::vector<std::vector<ExprPtr>>>(make_error(
            CASErrorKind::Unimplemented,
            "Shape lemma solver: pure polynomial not found"));
    }

    auto last_expr_res = f4_to_expr(*pure_last, variables, ctx);
    if (last_expr_res.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(last_expr_res.error());

    auto last_roots_res = solve_polynomial(last_expr_res.value(), variables.back(), ctx);
    if (last_roots_res.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(last_roots_res.error());

    std::vector<std::vector<ExprPtr>> all_solutions;
    for (ExprPtr root : last_roots_res.value()) {
        std::vector<ExprPtr> sol(n);
        sol[n - 1U] = root;
        all_solutions.push_back(std::move(sol));
    }

    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        std::vector<std::vector<ExprPtr>> next_solutions;
        for (const auto& sol : all_solutions) {
            std::vector<ExprPtr> roots_for_i;
            for (const auto& g : G) {
                Monomial lm = g.leading_monomial(order);
                if (lm.size() <= static_cast<std::size_t>(i) || lm[static_cast<std::size_t>(i)] == 0U) {
                    continue;
                }

                bool other_vars_ok = true;
                for (std::size_t k = 0; k < n; ++k) {
                    if (k < static_cast<std::size_t>(i) && lm[k] > 0U) {
                        other_vars_ok = false;
                        break;
                    }
                }
                if (!other_vars_ok) continue;

                auto g_expr_res = f4_to_expr(g, variables, ctx);
                if (g_expr_res.is_error()) continue;
                ExprPtr g_expr = g_expr_res.value();

                for (std::size_t k = static_cast<std::size_t>(i) + 1U; k < n; ++k) {
                    auto sub = symbolic::substitute(g_expr, variables[k], sol[k], ctx);
                    if (sub.is_ok()) g_expr = sub.value();
                }
                auto solved = solve_polynomial(g_expr, variables[static_cast<std::size_t>(i)], ctx);
                if (solved.is_ok() && !solved.value().empty()) {
                    roots_for_i = solved.value();
                    break;
                }
            }
            for (ExprPtr r : roots_for_i) {
                std::vector<ExprPtr> new_sol = sol;
                new_sol[static_cast<std::size_t>(i)] = r;
                next_solutions.push_back(std::move(new_sol));
            }
        }
        all_solutions = std::move(next_solutions);
    }

    return ok(std::move(all_solutions));
}

} // namespace cas::algebra
