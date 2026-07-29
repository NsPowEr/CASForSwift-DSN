#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/algebra.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <string>
#include <vector>

// polynomial_internal_xgcd.cpp — Euclide esteso su polinomi a coefficienti
// simbolici. Estratto da polynomial_internal.cpp, che era a 498 righe (limite
// anti-monolito 500) e non aveva spazio per la guardia di terminazione di A51.

namespace cas::algebra {

[[nodiscard]] Result<PolyXGCDResult> poly_extended_gcd(const PolyExpr& a, const PolyExpr& b, symbolic::CASContext& ctx) {
    if (a.is_zero()) {
        std::vector<ExprPtr> s_coeffs = {poly_make_integer(ctx.arena(), 0)};
        std::vector<ExprPtr> t_coeffs = {poly_make_integer(ctx.arena(), 1)};
        return ok(PolyXGCDResult{b, PolyExpr{s_coeffs}, PolyExpr{t_coeffs}});
    }
    if (b.is_zero()) {
        std::vector<ExprPtr> s_coeffs = {poly_make_integer(ctx.arena(), 1)};
        std::vector<ExprPtr> t_coeffs = {poly_make_integer(ctx.arena(), 0)};
        return ok(PolyXGCDResult{a, PolyExpr{s_coeffs}, PolyExpr{t_coeffs}});
    }

    auto div_res = divide_poly_with_remainder(a, b, ctx);
    if (div_res.is_error()) return fail<PolyXGCDResult>(div_res.error());
    
    auto q = div_res.value().quotient;
    auto r = div_res.value().remainder;

    if (r.is_zero()) {
        std::vector<ExprPtr> s_coeffs = {poly_make_integer(ctx.arena(), 0)};
        std::vector<ExprPtr> t_coeffs = {poly_make_integer(ctx.arena(), 1)};
        return ok(PolyXGCDResult{b, PolyExpr{s_coeffs}, PolyExpr{t_coeffs}});
    }

    // A51 — terminazione STRUTTURALE, non per fiducia nei coefficienti: Euclide
    // termina perche' deg(r) < deg(b), ma qui la riduzione a zero del leader
    // simbolico passa da `simplify` e, se quella non riesce, il grado apparente
    // non cala e la ricorsione non finisce (misurato: stack-overflow ASan).
    if (r.degree() >= b.degree()) {
        return fail<PolyXGCDResult>(make_error(CASErrorKind::Unimplemented,
            "poly_extended_gcd: il resto non decresce di grado ("
                + std::to_string(r.degree()) + " >= " + std::to_string(b.degree())
                + "): coefficienti simbolici non ridotti, Euclide non terminerebbe"));
    }

    auto ext_res = poly_extended_gcd(b, r, ctx);
    if (ext_res.is_error()) return fail<PolyXGCDResult>(ext_res.error());
    
    auto g = ext_res.value().gcd;
    auto s1 = ext_res.value().s;
    auto t1 = ext_res.value().t;

    auto q_t1_res = poly_multiply(q, t1, ctx);
    if (q_t1_res.is_error()) return fail<PolyXGCDResult>(q_t1_res.error());
    
    auto s_res = poly_subtract(s1, q_t1_res.value(), ctx);
    if (s_res.is_error()) return fail<PolyXGCDResult>(s_res.error());

    return ok(PolyXGCDResult{g, t1, s_res.value()});
}

}  // namespace cas::algebra
