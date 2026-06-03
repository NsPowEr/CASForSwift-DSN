// F5.7 — Zeilberger driver: wires the symbolic zeilberger_sum into the
// definite-sum pipeline.  Detects the "outer" parameter n from the upper
// bound (e.g., upper = n → n_param = n) and dispatches to zeilberger_sum.
// Handles multi-candidate n detection when upper contains multiple free symbols.

#include "summation_internal.hpp"
#include "../symbolic/summation_zeilberger.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include <optional>
#include <vector>

namespace cas::calculus {

namespace {

// Collect all Symbol nodes that appear in expr.
void collect_symbols(ExprPtr e, std::vector<Symbol>& out) {
    if (!e) return;
    if (const auto* sym = expr_cast<Symbol>(e)) {
        for (const auto& s : out)
            if (s.name == sym->name) return;
        out.push_back(*sym);
        return;
    }
    if (const auto* bin = expr_cast<Binary>(e)) {
        collect_symbols(bin->left, out); collect_symbols(bin->right, out); return;
    }
    if (const auto* un = expr_cast<Unary>(e)) {
        collect_symbols(un->operand, out); return;
    }
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        for (ExprPtr a : fc->args) collect_symbols(a, out); return;
    }
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms) collect_symbols(t, out); return;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr f : p->factors) collect_symbols(f, out); return;
    }
}

}  // namespace

Result<ExprPtr> try_zeilberger_definite(
    ExprPtr term, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {

    // Zeilberger only applies when upper is (or contains) a symbolic parameter
    // distinct from the summation variable var.
    std::vector<Symbol> upper_syms;
    collect_symbols(upper, upper_syms);
    // Remove var itself from candidates.
    upper_syms.erase(
        std::remove_if(upper_syms.begin(), upper_syms.end(),
            [&var](const Symbol& s) { return s.name == var.name; }),
        upper_syms.end());

    if (upper_syms.empty()) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Zeilberger: upper bound contains no free parameter; "
                       "use Gosper or Abramov for purely-numeric-bound sums.",
        });
    }

    // Try each candidate n_param.
    for (const Symbol& n_cand : upper_syms) {
        // Skip if the term does not depend on n_cand — Gosper/Abramov suffice.
        bool term_has_n = false;
        {
            std::vector<Symbol> term_syms;
            collect_symbols(term, term_syms);
            for (const auto& s : term_syms)
                if (s.name == n_cand.name) { term_has_n = true; break; }
        }
        if (!term_has_n) continue;

        auto res = symbolic::zeilberger_sum(term, n_cand, var, lower, upper, ctx);
        if (res.is_error()) continue;
        if (res.value().has_value())
            return ok(*res.value());
    }

    return fail<ExprPtr>(CASError{
        .kind    = CASErrorKind::Unimplemented,
        .message = "Zeilberger: no creative-telescoping recurrence found within "
                   "order J ≤ ctx.max_zeilberger_order() and polynomial degree "
                   "D ≤ ctx.max_zeilberger_poly_degree().  Increase these limits "
                   "or check that the term is proper hypergeometric in k.",
    });
}

}  // namespace cas::calculus
