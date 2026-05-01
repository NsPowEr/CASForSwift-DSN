#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"
#include "polynomial_groebner_f4.hpp"

#include <vector>

namespace cas::algebra {

[[nodiscard]] Result<std::vector<ExprPtr>> polynomial_groebner(
    const std::vector<ExprPtr>& equations,
    const std::vector<Symbol>& variables,
    symbolic::CASContext& ctx) {
    
    if (equations.empty()) return ok(std::vector<ExprPtr>{});

    // Converti a PolyF4 usando coefficienti razionali
    std::vector<PolyF4> F;
    F.reserve(equations.size());
    for (ExprPtr eq : equations) {
        auto r = expr_to_f4(eq, variables, ctx);
        if (r.is_error()) return fail<std::vector<ExprPtr>>(r.error());
        if (r.value().terms.empty()) continue; // polinomio zero
        F.push_back(r.value());
    }

    if (F.empty()) return ok(std::vector<ExprPtr>{});

    // Calcola base di Gröbner (GRevLex per efficienza)
    auto G = f4_groebner(F, MonomialOrder::GRevLex);

    // Inter-riduzione -> base ridotta unica
    inter_reduce(G, MonomialOrder::GRevLex);

    // Converti tornando a ExprPtr
    std::vector<ExprPtr> result;
    result.reserve(G.size());
    for (const PolyF4& g : G) {
        auto e = f4_to_expr(g, variables, ctx);
        if (e.is_error()) return fail<std::vector<ExprPtr>>(e.error());
        result.push_back(e.value());
    }
    
    return ok(result);
}

} // namespace cas::algebra
