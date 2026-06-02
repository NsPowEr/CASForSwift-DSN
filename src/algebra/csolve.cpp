#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"
#include "polynomial_groebner_f4.hpp"
#include "algebra_internal.hpp"

#include <vector>

namespace cas::algebra {

namespace {

// ─── Linear system fast-path (rectangular, parametric) ────────────────────────
//
// Detects systems where every poly is linear in every var. Builds a dense
// Q-coefficient matrix via O(m·n) substitution probes, runs Gauss-Jordan
// with explicit pivot tracking, and produces a particular solution by
// setting free variables to 0.
//
// Rationale: F4 on linear underdetermined systems either rejects or
// degenerates into nonterminating Buchberger iteration (the system has
// infinitely many parametric solutions over Q[x]). Linear systems are
// O(m·n²) by Gauss; routing them through F4 violates Regola Zero
// ("scegliere la soluzione matematicamente corretta e generale").
// This fast-path is the canonical algorithm for the linear subdomain.

[[nodiscard]] bool is_linear_in_vars(const std::vector<ExprPtr>& polys,
                                     const std::vector<Symbol>& vars,
                                     symbolic::CASContext& ctx) {
    for (auto p : polys) {
        for (auto& v : vars) {
            auto deg = polynomial_degree(p, v, ctx);
            if (deg.is_error()) return false;
            if (deg.value() > 1U) return false;
        }
    }
    return true;
}

[[nodiscard]] Result<Rational> expr_to_rational_strict(ExprPtr e, symbolic::CASContext& ctx) {
    if (const auto* i = expr_cast<IntegerLit>(e)) return ok(Rational(i->value));
    if (const auto* r = expr_cast<RationalLit>(e)) {
        return ok(Rational(r->numerator, r->denominator));
    }
    if (const auto* u = expr_cast<Unary>(e); u && u->op == UnaryOp::Neg) {
        auto inner = expr_to_rational_strict(u->operand, ctx);
        if (inner.is_ok()) {
            return ok(-inner.value());
        }
    }
    auto parts = split_num_den(e, ctx);
    if (parts.is_ok()) {
        auto num = expr_to_integer_coefficient(parts.value().numerator);
        auto den = expr_to_integer_coefficient(parts.value().denominator);
        if (num.is_ok() && den.is_ok() && !den.value().is_zero()) {
            return ok(Rational(num.value(), den.value()));
        }
    }
    return fail<Rational>(make_error(
        CASErrorKind::InternalError,
        "csolve linear fast-path: coefficient is not rational"));
}

[[nodiscard]] Result<ExprPtr> solve_linear_rect(
    const std::vector<ExprPtr>& polys,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();
    const std::size_t m = polys.size();
    const std::size_t n = vars.size();

    ExprPtr zero_expr = arena.make<IntegerLit>(0);
    ExprPtr one_expr = arena.make<IntegerLit>(1);

    auto substitute_all = [&](ExprPtr p, const std::vector<ExprPtr>& values) -> Result<ExprPtr> {
        ExprPtr cur = p;
        for (std::size_t l = 0; l < n; ++l) {
            auto r = symbolic::substitute(cur, vars[l], values[l], ctx);
            if (r.is_error()) return r;
            cur = r.value();
        }
        return symbolic::simplify(cur, ctx);
    };

    // Build augmented matrix [A | b] over Q. A is m×n, b is m×1.
    // For each poly_i: poly_i(vars) = sum_j A[i][j]*vars[j] + b_i_const.
    // Equation poly_i = 0 → sum_j A[i][j]*vars[j] = -b_i_const.
    std::vector<std::vector<Rational>> M(m, std::vector<Rational>(n + 1, Rational(BigInt(0))));

    std::vector<ExprPtr> zero_point(n, zero_expr);

    for (std::size_t i = 0; i < m; ++i) {
        // Constant term: poly(vars=0).
        auto const_expr = substitute_all(polys[i], zero_point);
        if (const_expr.is_error()) return fail<ExprPtr>(const_expr.error());
        auto const_q = expr_to_rational_strict(const_expr.value(), ctx);
        if (const_q.is_error()) return fail<ExprPtr>(const_q.error());
        M[i][n] = -const_q.value();  // RHS = -constant

        // Coefficient of vars[j]: poly(vars[j]=1, others=0) - constant.
        for (std::size_t j = 0; j < n; ++j) {
            std::vector<ExprPtr> probe(n, zero_expr);
            probe[j] = one_expr;
            auto probe_expr = substitute_all(polys[i], probe);
            if (probe_expr.is_error()) return fail<ExprPtr>(probe_expr.error());
            auto probe_q = expr_to_rational_strict(probe_expr.value(), ctx);
            if (probe_q.is_error()) return fail<ExprPtr>(probe_q.error());
            M[i][j] = probe_q.value() - const_q.value();
        }
    }

    // Gauss-Jordan with explicit pivot tracking (rectangular).
    std::vector<int> pivot_col(m, -1);
    std::size_t row = 0;
    for (std::size_t col = 0; col < n && row < m; ++col) {
        std::size_t pivot = row;
        while (pivot < m && M[pivot][col].numerator().is_zero()) ++pivot;
        if (pivot == m) continue;  // free column
        std::swap(M[row], M[pivot]);
        pivot_col[row] = static_cast<int>(col);

        Rational lead = M[row][col];
        for (std::size_t j = col; j <= n; ++j) {
            M[row][j] = M[row][j] / lead;
        }
        for (std::size_t k = 0; k < m; ++k) {
            if (k == row) continue;
            Rational f = M[k][col];
            if (f.numerator().is_zero()) continue;
            for (std::size_t j = col; j <= n; ++j) {
                M[k][j] = M[k][j] - f * M[row][j];
            }
        }
        ++row;
    }

    // Consistency check: any row of all-zero coefficients with non-zero RHS → no solution.
    for (std::size_t i = row; i < m; ++i) {
        if (!M[i][n].numerator().is_zero()) {
            return ok(arena.make<Matrix>(0U, n, std::vector<ExprPtr>{}));
        }
    }

    // Particular solution: free vars = 0, pivot vars = augmented column.
    std::vector<ExprPtr> x(n, arena.make<IntegerLit>(0));
    for (std::size_t i = 0; i < row; ++i) {
        if (pivot_col[i] < 0) continue;
        std::size_t col = static_cast<std::size_t>(pivot_col[i]);
        x[col] = make_rational_expr(arena, M[i][n]);
    }

    return ok(arena.make<Matrix>(1U, n, std::move(x)));
}

[[nodiscard]] std::vector<ExprPtr> extract_list(ExprPtr expr) {
    if (!expr) return {};
    if (const auto* m = expr_cast<Matrix>(expr)) {
        return m->elements;
    }
    if (const auto* f = expr_cast<FuncCall>(expr)) {
        if (f->name == "list" || f->name == "vector") {
            return f->args;
        }
    }
    return {expr};
}

[[nodiscard]] std::vector<Symbol> extract_vars(ExprPtr expr) {
    auto list = extract_list(expr);
    std::vector<Symbol> vars;
    for (auto e : list) {
        if (const auto* s = expr_cast<Symbol>(e)) {
            vars.push_back(*s);
        }
    }
    return vars;
}

[[nodiscard]] ExprPtr to_poly_form(ExprPtr eq, symbolic::CASContext& ctx) {
    if (const auto* b = expr_cast<Binary>(eq)) {
        if (b->op == BinaryOp::Equal) {
            return ctx.arena().make<Binary>(BinaryOp::Sub, b->left, b->right);
        }
    }
    return eq;
}

} // namespace

[[nodiscard]] Result<ExprPtr> csolve(
    ExprPtr eqs_expr,
    ExprPtr vars_expr,
    symbolic::CASContext& ctx) {
    
    auto eqs_list = extract_list(eqs_expr);
    auto vars = extract_vars(vars_expr);

    if (eqs_list.empty() || vars.empty()) {
        return ok(ctx.arena().make<Matrix>(0, 0, std::vector<ExprPtr>{}));
    }

    std::vector<ExprPtr> polys;
    for (auto eq : eqs_list) {
        polys.push_back(to_poly_form(eq, ctx));
    }

    // Linear fast-path: rectangular Gauss-Jordan with parametric free-vars=0
    // solution. Avoids F4 nonterminating Buchberger on underdetermined linear
    // systems (e.g. Gosper's polynomial ansatz with one free additive
    // constant). Routes only when every poly is provably linear in every
    // unknown — preserves F4 dispatch for genuinely nonlinear systems.
    if (is_linear_in_vars(polys, vars, ctx)) {
        return solve_linear_rect(polys, vars, ctx);
    }

    auto sol_res = solve_nonlinear_system_f4(polys, vars, ctx);
    // L2-13 fallback: 2-variable resultant elimination when F4 fails.
    // For systems [f(x,y), g(x,y)], compute h(x) = resultant_y(f, g);
    // solve h(x)=0 → x roots; for each x root, back-substitute and
    // solve g(x,y)=0 in y. Robust for low-arity systems where F4
    // struggles with coefficient swell.
    if (sol_res.is_error() && polys.size() == 2U && vars.size() == 2U) {
        ExprPtr f = polys[0];
        ExprPtr g = polys[1];
        // Eliminate vars[1] via resultant wrt vars[1].
        auto res_h = polynomial_resultant(f, g, vars[1], ctx);
        if (res_h.is_ok()) {
            auto x_roots = solve_polynomial(res_h.value(), vars[0], ctx);
            if (x_roots.is_ok()) {
                std::vector<ExprPtr> flat;
                AstArena& arena = ctx.arena();
                for (auto x_val : x_roots.value()) {
                    // Substitute x into g and solve for y.
                    auto g_sub = symbolic::substitute(g, vars[0], x_val, ctx);
                    if (g_sub.is_error()) continue;
                    auto y_roots = solve_polynomial(g_sub.value(), vars[1], ctx);
                    if (y_roots.is_error()) continue;
                    for (auto y_val : y_roots.value()) {
                        flat.push_back(x_val);
                        flat.push_back(y_val);
                    }
                }
                std::size_t n_sols = flat.size() / 2U;
                return ok(arena.make<Matrix>(n_sols, 2U, std::move(flat)));
            }
        }
        return fail<ExprPtr>(sol_res.error());
    }
    if (sol_res.is_error()) return fail<ExprPtr>(sol_res.error());

    const auto& solutions = sol_res.value();
    std::vector<ExprPtr> flat_elements;
    flat_elements.reserve(solutions.size() * vars.size());
    for (const auto& sol : solutions) {
        for (auto val : sol) {
            flat_elements.push_back(val);
        }
    }

    return ok(ctx.arena().make<Matrix>(solutions.size(), vars.size(), std::move(flat_elements)));
}

} // namespace cas::algebra
