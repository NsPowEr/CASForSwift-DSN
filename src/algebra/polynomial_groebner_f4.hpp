#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "polynomial_internal.hpp"
#include <vector>
#include <map>

namespace cas {
namespace symbolic {
class CASContext;
}

namespace algebra {

using Monomial = std::vector<unsigned int>;

enum class MonomialOrder {
    Lex,
    GRevLex
};

struct PolyF4 {
    std::map<Monomial, Rational> terms;

    [[nodiscard]] bool is_zero() const { return terms.empty(); }
    
    Monomial leading_monomial(MonomialOrder order = MonomialOrder::Lex) const;
    Rational leading_coefficient(MonomialOrder order = MonomialOrder::Lex) const;
    
    void make_monic(MonomialOrder order = MonomialOrder::Lex);
};

Result<std::vector<PolyF4>> f4_groebner(
    std::vector<PolyF4> G,
    MonomialOrder order = MonomialOrder::GRevLex,
    symbolic::CASContext* ctx = nullptr);

Result<ExprPtr> f4_to_expr(
    const PolyF4& p,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

Result<PolyF4> expr_to_f4(
    ExprPtr expr,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

void inter_reduce(std::vector<PolyF4>& G, MonomialOrder order);

// Returns true iff G is a reduced Groebner basis:
// (1) Groebner: every S-polynomial reduces to 0 mod G.
// (2) Monic: every element has leading coefficient 1.
// (3) Fully inter-reduced: no LM of any g_i divides any term of any other g_j.
// Complexity O(|G|^2 * max_terms); use only for testing/validation.
[[nodiscard]] bool is_reduced_groebner_basis(
    const std::vector<PolyF4>& G,
    MonomialOrder order);

Result<std::vector<std::vector<ExprPtr>>> solve_nonlinear_system_f4(
    const std::vector<ExprPtr>& equations,
    const std::vector<Symbol>& variables,
    symbolic::CASContext& ctx);

} // namespace algebra
} // namespace cas
