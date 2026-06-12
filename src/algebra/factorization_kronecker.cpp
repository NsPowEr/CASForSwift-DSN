#include "polynomial_internal.hpp"
#include "algebra_internal.hpp"
#include "cas/numtheory.hpp"
#include <vector>
#include <algorithm>

namespace cas::algebra {

namespace {

// Interpolate a polynomial of degree at most x.size() - 1 in Q[x]
// ascending coefficients
std::vector<Rational> interpolate_lagrange(const std::vector<BigInt>& x, const std::vector<BigInt>& y) {
    std::size_t n = x.size();
    std::vector<Rational> result(n, Rational(0));
    
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<Rational> term{Rational(1)};
        Rational denom(1);
        
        for (std::size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            std::vector<Rational> factor{Rational(-x[j]), Rational(1)};
            std::vector<Rational> next_term(term.size() + 1, Rational(0));
            for (std::size_t k = 0; k < term.size(); ++k) {
                next_term[k] += term[k] * factor[0];
                next_term[k + 1] += term[k] * factor[1];
            }
            while (next_term.size() > 1 && next_term.back().numerator().is_zero()) {
                next_term.pop_back();
            }
            term = next_term;
            denom *= Rational(x[i] - x[j]);
        }
        
        Rational coeff = Rational(y[i]) / denom;
        for (std::size_t k = 0; k < term.size(); ++k) {
            result[k] += term[k] * coeff;
        }
    }
    
    while (result.size() > 1 && result.back().numerator().is_zero()) {
        result.pop_back();
    }
    return result;
}

bool has_integer_coefficients(const std::vector<Rational>& p) {
    for (const auto& coeff : p) {
        if (coeff.denominator() != BigInt(1)) {
            return false;
        }
    }
    return true;
}

IntPoly to_integer_poly(const std::vector<Rational>& p) {
    std::vector<BigInt> coeffs;
    for (const auto& coeff : p) {
        coeffs.push_back(coeff.numerator());
    }
    return IntPoly(coeffs);
}

BigInt eval_int_poly_at_int(const IntPoly& p, const BigInt& x) {
    if (p.empty()) return BigInt(0);
    BigInt result(0);
    BigInt x_pow(1);
    for (const BigInt& coeff : p.coefficients()) {
        result += coeff * x_pow;
        x_pow *= x;
    }
    return result;
}

bool search_tuples(
    std::size_t index,
    std::vector<BigInt>& current_y,
    const std::vector<BigInt>& x,
    const std::vector<std::vector<BigInt>>& D,
    const IntPoly& f,
    symbolic::CASContext& ctx,
    std::vector<IntPoly>& out_factors) {
    
    if (index == x.size()) {
        auto g_rational = interpolate_lagrange(x, current_y);
        if (has_integer_coefficients(g_rational)) {
            IntPoly g = to_integer_poly(g_rational);
            if (g.degree() >= 1 && g.degree() < f.degree()) {
                normalize_integer_poly(g);
                BigInt cont = integer_content(g);
                divide_integer_coefficients_by_scalar(g, cont);
                if (g.leading_coeff().is_negative()) {
                    for (auto& c : g.coefficients()) c = -c;
                }
                
                auto q_res = exact_divide_integer_poly(f, g, ctx);
                if (q_res.is_ok()) {
                    auto f1_res = factorize_kronecker(g, ctx);
                    auto f2_res = factorize_kronecker(q_res.value(), ctx);
                    if (f1_res.is_ok() && f2_res.is_ok()) {
                        out_factors = f1_res.value();
                        out_factors.insert(out_factors.end(), f2_res.value().begin(), f2_res.value().end());
                        return true;
                    }
                }
            }
        }
        return false;
    }
    
    for (const auto& val : D[index]) {
        current_y[index] = val;
        if (search_tuples(index + 1, current_y, x, D, f, ctx, out_factors)) {
            return true;
        }
    }
    return false;
}

} // namespace

Result<std::vector<IntPoly>> factorize_kronecker(
    const IntPoly& f,
    symbolic::CASContext& ctx) {
    
    if (f.empty() || f.degree() == 0U) {
        return ok(std::vector<IntPoly>{});
    }
    if (f.degree() == 1U) {
        return ok(std::vector<IntPoly>{f});
    }

    std::size_t r = f.degree() / 2;
    
    // Choose distinct points: 0, 1, -1, 2, -2, 3, -3, 4, -4...
    std::vector<BigInt> x_points;
    x_points.push_back(BigInt(0));
    for (std::size_t i = 1; x_points.size() <= r + 1; ++i) {
        x_points.push_back(BigInt(static_cast<long long>(i)));
        x_points.push_back(BigInt(-static_cast<long long>(i)));
    }
    x_points.resize(r + 1);

    for (std::size_t s = 1; s <= r; ++s) {
        std::vector<BigInt> x(s + 1);
        std::copy(x_points.begin(), x_points.begin() + s + 1, x.begin());
        
        // Compute f(x_i) and their divisors
        std::vector<std::vector<BigInt>> D(s + 1);
        bool found_zero_root = false;
        BigInt zero_root_val;
        
        for (std::size_t i = 0; i <= s; ++i) {
            BigInt val = eval_int_poly_at_int(f, x[i]);
            if (val.is_zero()) {
                found_zero_root = true;
                zero_root_val = x[i];
                break;
            }
            
            auto div_res = numtheory::divisors(val);
            if (div_res.is_error()) {
                return fail<std::vector<IntPoly>>(div_res.error());
            }
            
            std::vector<BigInt> divs;
            for (const auto& d : div_res.value()) {
                divs.push_back(d);
                divs.push_back(-d);
            }
            D[i] = std::move(divs);
        }
        
        if (found_zero_root) {
            // (x - zero_root_val) is a factor
            IntPoly linear_fac(std::vector<BigInt>{-zero_root_val, BigInt(1)});
            auto q_res = exact_divide_integer_poly(f, linear_fac, ctx);
            if (q_res.is_ok()) {
                auto f1_res = factorize_kronecker(linear_fac, ctx);
                auto f2_res = factorize_kronecker(q_res.value(), ctx);
                if (f1_res.is_ok() && f2_res.is_ok()) {
                    std::vector<IntPoly> result = f1_res.value();
                    result.insert(result.end(), f2_res.value().begin(), f2_res.value().end());
                    return ok(result);
                }
            }
        }
        
        std::vector<BigInt> current_y(s + 1, BigInt(0));
        std::vector<IntPoly> factors;
        if (search_tuples(0, current_y, x, D, f, ctx, factors)) {
            return ok(factors);
        }
    }

    // If no factor was found, f is irreducible
    return ok(std::vector<IntPoly>{f});
}

} // namespace cas::algebra
