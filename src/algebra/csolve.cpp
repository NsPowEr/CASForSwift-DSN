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

    try {
        auto sol_res = solve_nonlinear_system_f4(polys, vars, ctx);
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
    } catch (const std::exception& e) {
        return fail<ExprPtr>(make_error(CASErrorKind::InternalError, std::string("C++ Exception in csolve: ") + e.what()));
    } catch (...) {
        return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unknown Exception in csolve"));
    }
}

} // namespace cas::algebra
