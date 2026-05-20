#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"
#include "polynomial_groebner_f4.hpp"
#include "algebra_internal.hpp"

#include <vector>

namespace cas::algebra {

namespace {

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
