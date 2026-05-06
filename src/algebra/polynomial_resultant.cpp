#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <utility>

namespace cas {
namespace algebra {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}

[[nodiscard]] ExprPtr mul_expr(AstArena& arena, ExprPtr a, ExprPtr b) {
    return arena.make<Binary>(BinaryOp::Mul, a, b);
}
[[nodiscard]] ExprPtr pow_expr(AstArena& arena, ExprPtr base, ExprPtr exp) {
    return arena.make<Binary>(BinaryOp::Pow, base, exp);
}
[[nodiscard]] ExprPtr div_expr(AstArena& arena, ExprPtr a, ExprPtr b) {
    return arena.make<Binary>(BinaryOp::Div, a, b);
}

[[nodiscard]] Result<PolyExpr> pseudo_remainder_poly_expr(PolyExpr A, const PolyExpr& B, symbolic::CASContext& ctx) {
    if (is_zero_poly(B))
        return fail<PolyExpr>(make_error(CASErrorKind::InvalidArgument, "Divisor cannot be zero"));

    std::size_t m = poly_degree(A);
    std::size_t n = poly_degree(B);
    if (m < n) return ok(std::move(A));

    ExprPtr b_n = leading_coefficient(B);
    PolyExpr R = A;

    for (std::size_t step = 0; step <= m - n; ++step) {
        std::size_t deg_r = poly_degree(R);
        if (deg_r == 0 && is_zero_poly(R)) break;

        if (deg_r == m - step) {
            ExprPtr lc_r = leading_coefficient(R);
            for (auto& coeff : R.coefficients()) {
                if (coeff) {
                    auto s = poly_simplify_expr(mul_expr(ctx.arena(), coeff, b_n), ctx);
                    if (s.is_error()) return fail<PolyExpr>(s.error());
                    coeff = s.value();
                }
            }
            PolyExpr term;
            term.resize(deg_r - n + 1, nullptr);
            term[deg_r - n] = lc_r;

            auto sub_term_res = poly_multiply(B, term, ctx);
            if (sub_term_res.is_error()) return fail<PolyExpr>(sub_term_res.error());

            auto sub_res = poly_subtract(R, sub_term_res.value(), ctx);
            if (sub_res.is_error()) return fail<PolyExpr>(sub_res.error());
            R = sub_res.value();
        } else {
            for (auto& coeff : R.coefficients()) {
                if (coeff) {
                    auto s = poly_simplify_expr(mul_expr(ctx.arena(), coeff, b_n), ctx);
                    if (s.is_error()) return fail<PolyExpr>(s.error());
                    coeff = s.value();
                }
            }
        }
    }
    normalize_poly(R);
    return ok(std::move(R));
}

} // namespace

[[nodiscard]] Result<ExprPtr> polynomial_resultant(ExprPtr p, ExprPtr q, const Symbol& var, symbolic::CASContext& ctx) {
    auto res_f = parse_polynomial(p, var, ctx);
    if (res_f.is_error()) return fail<ExprPtr>(res_f.error());
    auto res_g = parse_polynomial(q, var, ctx);
    if (res_g.is_error()) return fail<ExprPtr>(res_g.error());

    PolyExpr f = res_f.value();
    PolyExpr g = res_g.value();

    if (is_zero_poly(f) || is_zero_poly(g))
        return ok(poly_make_integer(ctx.arena(), 0));

    std::size_t n = poly_degree(f);
    std::size_t m = poly_degree(g);

    int sign_correction = 1;
    if (n < m) {
        std::swap(f, g);
        std::swap(n, m);
        if ((n * m) % 2 != 0) sign_correction = -1;
    }

    std::size_t d = n - m;
    ExprPtr b = poly_make_integer(ctx.arena(), ((d + 1) % 2 != 0) ? -1 : 1);

    auto prem_res = pseudo_remainder_poly_expr(f, g, ctx);
    if (prem_res.is_error()) return fail<ExprPtr>(prem_res.error());
    PolyExpr h = prem_res.value();

    for (auto& coeff : h.coefficients()) {
        if (coeff) {
            auto s = poly_simplify_expr(mul_expr(ctx.arena(), coeff, b), ctx);
            if (s.is_error()) return fail<ExprPtr>(s.error());
            coeff = s.value();
        }
    }
    normalize_poly(h);

    ExprPtr lc = leading_coefficient(g);
    ExprPtr S_last;
    {
        auto s = poly_simplify_expr(pow_expr(ctx.arena(), lc, poly_make_integer(ctx.arena(), d)), ctx);
        if (s.is_error()) return fail<ExprPtr>(s.error());
        ExprPtr lc_pow = s.value();
        auto c_res = poly_simplify_expr(mul_expr(ctx.arena(), poly_make_integer(ctx.arena(), -1), lc_pow), ctx);
        if (c_res.is_error()) return fail<ExprPtr>(c_res.error());
        S_last = c_res.value();
    }
    ExprPtr c = S_last;

    while (!is_zero_poly(h)) {
        std::size_t k = poly_degree(h);
        f = g;
        g = h;
        m = k;
        std::size_t d_next = poly_degree(f) - k;

        {
            auto cpow = poly_simplify_expr(pow_expr(ctx.arena(), c, poly_make_integer(ctx.arena(), d_next)), ctx);
            if (cpow.is_error()) return fail<ExprPtr>(cpow.error());
            auto neg_lc = poly_simplify_expr(mul_expr(ctx.arena(), poly_make_integer(ctx.arena(), -1), lc), ctx);
            if (neg_lc.is_error()) return fail<ExprPtr>(neg_lc.error());
            auto b_res = poly_simplify_expr(mul_expr(ctx.arena(), neg_lc.value(), cpow.value()), ctx);
            if (b_res.is_error()) return fail<ExprPtr>(b_res.error());
            b = b_res.value();
        }

        auto next_prem = pseudo_remainder_poly_expr(f, g, ctx);
        if (next_prem.is_error()) return fail<ExprPtr>(next_prem.error());
        h = next_prem.value();

        for (auto& coeff : h.coefficients()) {
            if (coeff) {
                auto s = poly_simplify_expr(div_expr(ctx.arena(), coeff, b), ctx);
                if (s.is_error()) return fail<ExprPtr>(s.error());
                coeff = s.value();
            }
        }
        normalize_poly(h);

        lc = leading_coefficient(g);

        if (d_next > 1) {
            auto q_res = poly_simplify_expr(pow_expr(ctx.arena(), c, poly_make_integer(ctx.arena(), d_next - 1)), ctx);
            if (q_res.is_error()) return fail<ExprPtr>(q_res.error());
            auto neg_lc_res = poly_simplify_expr(mul_expr(ctx.arena(), poly_make_integer(ctx.arena(), -1), lc), ctx);
            if (neg_lc_res.is_error()) return fail<ExprPtr>(neg_lc_res.error());
            auto neg_lc_pow = poly_simplify_expr(pow_expr(ctx.arena(), neg_lc_res.value(), poly_make_integer(ctx.arena(), d_next)), ctx);
            if (neg_lc_pow.is_error()) return fail<ExprPtr>(neg_lc_pow.error());
            auto c_res = poly_simplify_expr(div_expr(ctx.arena(), neg_lc_pow.value(), q_res.value()), ctx);
            if (c_res.is_error()) return fail<ExprPtr>(c_res.error());
            c = c_res.value();
        } else {
            auto c_res = poly_simplify_expr(mul_expr(ctx.arena(), poly_make_integer(ctx.arena(), -1), lc), ctx);
            if (c_res.is_error()) return fail<ExprPtr>(c_res.error());
            c = c_res.value();
        }

        auto s_last_res = poly_simplify_expr(mul_expr(ctx.arena(), poly_make_integer(ctx.arena(), -1), c), ctx);
        if (s_last_res.is_error()) return fail<ExprPtr>(s_last_res.error());
        S_last = s_last_res.value();
    }

    if (poly_degree(g) > 0)
        return ok(poly_make_integer(ctx.arena(), 0));

    ExprPtr raw = S_last;
    if (sign_correction == -1) {
        auto neg = poly_simplify_expr(mul_expr(ctx.arena(), poly_make_integer(ctx.arena(), -1), S_last), ctx);
        if (neg.is_error()) return fail<ExprPtr>(neg.error());
        raw = neg.value();
    }
    // Normalize: simplify any residual rational factors
    auto normalized = ctx.simplify(raw);
    return normalized.is_ok() ? normalized : ok(raw);
}

[[nodiscard]] Result<ExprPtr> polynomial_discriminant(ExprPtr p, const Symbol& var, symbolic::CASContext& ctx) {
    if (!p) return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Discriminant requires non-null polynomial"));

    auto res_f = parse_polynomial(p, var, ctx);
    if (res_f.is_error()) return fail<ExprPtr>(res_f.error());
    PolyExpr f = std::move(res_f.value());

    if (is_zero_poly(f)) return ok(poly_make_integer(ctx.arena(), 0));
    std::size_t n = poly_degree(f);
    if (n < 2) return ok(poly_make_integer(ctx.arena(), 0));

    // Disc(f) = (-1)^(n(n-1)/2) * (1/an) * Res(f, f')
    ExprPtr an = leading_coefficient(f);
    
    // Compute derivative f'
    auto deriv_res = cas::calculus::diff(p, var, 1U, ctx);
    if (deriv_res.is_error()) return fail<ExprPtr>(deriv_res.error());
    
    auto res_val = polynomial_resultant(p, deriv_res.value(), var, ctx);
    if (res_val.is_error()) return res_val;

    // (-1)^(n(n-1)/2): parity from (n mod 4): sign flips when n≡2 or n≡3 (mod 4)
    const int sign = ((n % 4U == 2U) || (n % 4U == 3U)) ? -1 : 1;

    auto signed_res = poly_simplify_expr(mul_expr(ctx.arena(), poly_make_integer(ctx.arena(), sign), res_val.value()), ctx);
    if (signed_res.is_error()) return fail<ExprPtr>(signed_res.error());

    auto disc = poly_simplify_expr(div_expr(ctx.arena(), signed_res.value(), an), ctx);
    if (disc.is_error()) return fail<ExprPtr>(disc.error());

    return disc;
}

}  // namespace algebra
}  // namespace cas
