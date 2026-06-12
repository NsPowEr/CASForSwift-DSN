// Kovacic Case 1 Laurent series calculations and helper arithmetic functions.
// Ref: Kovacic J.J. (1986) §3, "An Algorithm for Solving Second Order Linear
// Homogeneous Differential Equations" J. Symbolic Computation 2, 3–43.

#include "ode_kovacic_internal.hpp"
#include <optional>
#include <vector>
#include <algorithm>

namespace cas::calculus {
namespace kovacic_impl {

std::optional<BigInt> bigint_isqrt(const BigInt& n) {
    if (n.is_negative()) return std::nullopt;
    if (n.is_zero()) return BigInt(0);

    const std::size_t bits = n.bit_length();
    const std::size_t init_bits = (bits + 1U) / 2U;
    BigInt x = BigInt(1).shift_left_bits(init_bits);

    const std::size_t max_iter = 2U * bits + 64U;
    for (std::size_t iter = 0U; iter < max_iter; ++iter) {
        BigInt q = n / x;
        BigInt x_next = (x + q).shift_right_bits(1U);
        if (x_next >= x) {
            BigInt c = (x_next < x) ? x_next : x;
            while (!c.is_zero() && c * c > n) c = c - BigInt(1);
            if (c * c == n) return c;
            return std::nullopt;
        }
        x = x_next;
    }
    return std::nullopt;
}

std::optional<Rational> as_rational(ExprPtr e) {
    if (auto* il = expr_cast<IntegerLit>(e)) return Rational(il->value, BigInt(1));
    if (auto* rl = expr_cast<RationalLit>(e)) return Rational(rl->numerator, rl->denominator);
    return std::nullopt;
}

std::optional<Rational> rational_sqrt(const Rational& r) {
    if (r < Rational(BigInt(0))) return std::nullopt;
    auto sp = bigint_isqrt(r.numerator());
    auto sq = bigint_isqrt(r.denominator());
    if (!sp || !sq) return std::nullopt;
    return Rational(*sp, *sq);
}

Result<ExprPtr> reverse_polynomial(
    ExprPtr poly, const Symbol& x, const Symbol& y, AstArena& a, symbolic::CASContext& ctx) {
    
    auto coeffs_res = algebra::univariate_coefficients(poly, x, ctx);
    if (coeffs_res.is_error()) return fail<ExprPtr>(coeffs_res.error());
    
    auto coeffs = coeffs_res.value();
    std::reverse(coeffs.begin(), coeffs.end());
    
    ExprPtr result = kv_int(a, 0);
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        if (kv_is_zero(coeffs[i], ctx)) continue;
        ExprPtr term = coeffs[i];
        if (i > 0) {
            ExprPtr ypow = (i == 1)
                ? static_cast<ExprPtr>(a.make<Symbol>(y.name))
                : a.make<Binary>(BinaryOp::Pow, a.make<Symbol>(y.name), kv_int(a, i));
            term = kv_mul(a, term, ypow);
        }
        result = kv_is_zero(result, ctx) ? term : kv_add(a, result, term);
    }
    return ctx.simplify(result);
}

std::optional<std::vector<Rational>> compute_taylor_rational(
    ExprPtr num, ExprPtr den_other, const Symbol& x, ExprPtr c,
    unsigned terms_needed, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();
    Symbol y(ctx.make_fresh_symbol("y"));
    
    auto shift_val = kv_add(a, a.make<Symbol>(y.name), c);
    auto num_shifted_res = ctx.substitute(num, x, shift_val);
    auto den_shifted_res = ctx.substitute(den_other, x, shift_val);
    if (num_shifted_res.is_error() || den_shifted_res.is_error())
        return std::nullopt;
        
    auto num_simp = ctx.simplify(num_shifted_res.value());
    auto den_simp = ctx.simplify(den_shifted_res.value());
    if (num_simp.is_error() || den_simp.is_error())
        return std::nullopt;
        
    auto num_coeffs_res = algebra::univariate_coefficients(num_simp.value(), y, ctx);
    auto den_coeffs_res = algebra::univariate_coefficients(den_simp.value(), y, ctx);
    if (num_coeffs_res.is_error() || den_coeffs_res.is_error())
        return std::nullopt;
        
    const auto& num_coeffs = num_coeffs_res.value();
    const auto& den_coeffs = den_coeffs_res.value();
    
    std::vector<Rational> num_r;
    for (auto coeff : num_coeffs) {
        auto cs = ctx.simplify(coeff);
        if (cs.is_error()) return std::nullopt;
        auto cr = as_rational(cs.value());
        if (!cr) return std::nullopt;
        num_r.push_back(*cr);
    }
    
    std::vector<Rational> den_r;
    for (auto coeff : den_coeffs) {
        auto cs = ctx.simplify(coeff);
        if (cs.is_error()) return std::nullopt;
        auto cr = as_rational(cs.value());
        if (!cr) return std::nullopt;
        den_r.push_back(*cr);
    }
    
    if (den_r.empty()) return std::nullopt;
    
    Rational d_0 = den_r[0];
    if (d_0.numerator().is_zero()) return std::nullopt;
    
    std::vector<Rational> u(terms_needed);
    for (unsigned k = 0; k < terms_needed; ++k) {
        Rational a_k = (k < num_r.size()) ? num_r[k] : Rational(BigInt(0));
        Rational sum(BigInt(0));
        for (unsigned j = 0; j < k; ++j) {
            Rational d_kj = (k - j < den_r.size()) ? den_r[k - j] : Rational(BigInt(0));
            sum = sum + u[j] * d_kj;
        }
        u[k] = (a_k - sum) / d_0;
    }
    
    return u;
}

std::optional<std::vector<Rational>> compute_laurent_sqrt(
    const std::vector<Rational>& u, unsigned terms_needed) {
    
    if (u.empty()) return std::nullopt;
    
    auto s0 = rational_sqrt(u[0]);
    if (!s0) return std::nullopt;
    
    std::vector<Rational> s(terms_needed);
    s[0] = *s0;
    
    Rational two_s0 = Rational(BigInt(2)) * s[0];
    
    for (unsigned k = 1; k < terms_needed; ++k) {
        Rational sum(BigInt(0));
        for (unsigned i = 1; i < k; ++i) {
            sum = sum + s[i] * s[k - i];
        }
        s[k] = (u[k] - sum) / two_s0;
    }
    
    return s;
}

} // namespace kovacic_impl
} // namespace cas::calculus
