#include "polynomial_internal.hpp"
#include "cas/numtheory.hpp"
#include <vector>
#include <algorithm>
#include <random>

namespace cas::algebra {

static BigInt mod_bigint(const BigInt& a, const BigInt& m) {
    BigInt r = a % m;
    if (r.is_negative()) r += m.abs();
    return r;
}

static void poly_mod(IntPoly& p, const BigInt& m) {
    for (auto& c : p.coefficients()) c = mod_bigint(c, m);
    p.normalize([](const BigInt& v) { return v.is_zero(); });
}

static IntPoly poly_sub_mod(const IntPoly& a, const IntPoly& b, const BigInt& m) {
    std::size_t sz = std::max(a.size(), b.size());
    IntPoly r;
    r.resize(sz, BigInt(0));
    for (std::size_t i = 0; i < sz; ++i) {
        BigInt ai = (i < a.size()) ? a[i] : BigInt(0);
        BigInt bi = (i < b.size()) ? b[i] : BigInt(0);
        r[i] = mod_bigint(ai - bi, m);
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

static IntPoly poly_add_mod(const IntPoly& a, const IntPoly& b, const BigInt& m) {
    std::size_t sz = std::max(a.size(), b.size());
    IntPoly r;
    r.resize(sz, BigInt(0));
    for (std::size_t i = 0; i < sz; ++i) {
        BigInt ai = (i < a.size()) ? a[i] : BigInt(0);
        BigInt bi = (i < b.size()) ? b[i] : BigInt(0);
        r[i] = mod_bigint(ai + bi, m);
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

static IntPoly poly_mul_mod(const IntPoly& a, const IntPoly& b, const BigInt& m) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly r;
    r.resize(a.size() + b.size() - 1, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].is_zero()) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] = mod_bigint(r[i + j] + a[i] * b[j], m);
        }
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

static std::pair<IntPoly, IntPoly> poly_div_rem_mod(const IntPoly& a, const IntPoly& b, const BigInt& p) {
    if (b.is_zero()) return {IntPoly{}, a};
    IntPoly q, r = a;
    poly_mod(r, p);
    auto inv_lb_res = numtheory::modular_inverse(mod_bigint(b.leading_coeff(), p), p);
    if (inv_lb_res.is_error()) return {IntPoly{}, r};
    BigInt inv_lb = inv_lb_res.value();
    
    while (!r.is_zero() && r.degree() >= b.degree()) {
        std::size_t deg_diff = r.degree() - b.degree();
        BigInt factor = mod_bigint(r.leading_coeff() * inv_lb, p);
        if (q.size() <= deg_diff) q.resize(deg_diff + 1, BigInt(0));
        q[deg_diff] = factor;
        for (std::size_t i = 0; i < b.size(); ++i) {
            std::size_t target_idx = i + deg_diff;
            r[target_idx] = mod_bigint(r[target_idx] - factor * b[i], p);
        }
        r.normalize([](const BigInt& v) { return v.is_zero(); });
    }
    q.normalize([](const BigInt& v) { return v.is_zero(); });
    return {q, r};
}

static IntPoly poly_gcd_mod(IntPoly a, IntPoly b, const BigInt& p) {
    poly_mod(a, p);
    poly_mod(b, p);
    while (!b.is_zero()) {
        a = poly_div_rem_mod(a, b, p).second;
        std::swap(a, b);
    }
    if (!a.is_zero()) {
        BigInt inv = numtheory::modular_inverse(a.leading_coeff(), p).value();
        for (auto& c : a.coefficients()) c = mod_bigint(c * inv, p);
    }
    return a;
}

static IntPoly poly_pow_mod_poly(IntPoly base, BigInt exponent, const IntPoly& f, const BigInt& p) {
    IntPoly res{std::vector<BigInt>{BigInt(1)}};
    while (!exponent.is_zero()) {
        if (!(exponent % BigInt(2)).is_zero()) {
            res = poly_mul_mod(res, base, p);
            res = poly_div_rem_mod(res, f, p).second;
        }
        base = poly_mul_mod(base, base, p);
        base = poly_div_rem_mod(base, f, p).second;
        exponent /= BigInt(2);
    }
    return res;
}

std::vector<std::pair<IntPoly, std::size_t>> distinct_degree_factorization(IntPoly f, const BigInt& p) {
    std::vector<std::pair<IntPoly, std::size_t>> result;
    IntPoly x{std::vector<BigInt>{BigInt(0), BigInt(1)}};
    IntPoly w = x;
    std::size_t d = 1;
    while (f.degree() >= 2 * d) {
        // w = w^p mod f
        w = poly_pow_mod_poly(w, p, f, p);
        
        // g = gcd(f, w - x)
        IntPoly w_minus_x = poly_sub_mod(w, x, p);
        IntPoly g = poly_gcd_mod(f, w_minus_x, p);
        if (g.degree() > 0) {
            result.push_back({g, d});
            f = poly_div_rem_mod(f, g, p).first;
            w = poly_div_rem_mod(w, f, p).second;
        }
        d++;
    }
    if (f.degree() > 0) {
        result.push_back({f, f.degree()});
    }
    return result;
}

static std::vector<IntPoly> equal_degree_factorization(IntPoly f, std::size_t d, const BigInt& p) {
    if (f.degree() == d) return {f};
    
    std::vector<IntPoly> factors;
    std::size_t poly_seed = f.size();
    for (const auto& coeff : f.coefficients()) {
        std::uint64_t cv = coeff.to_u64();
        poly_seed ^= cv + 0x9e3779b9ULL + (poly_seed << 6U) + (poly_seed >> 2U);
    }
    std::mt19937 rng(static_cast<std::uint32_t>(poly_seed));
    
    while (f.degree() > d) {
        // Choose random a(x) with deg(a) < deg(f)
        IntPoly a;
        a.resize(f.degree());
        for (std::size_t i = 0; i < f.degree(); ++i) {
            a[i] = BigInt(static_cast<long long>(rng() % static_cast<unsigned int>(p.to_u64())));
        }
        a.normalize([](const BigInt& v) { return v.is_zero(); });
        if (a.is_zero()) continue;
        
        IntPoly g = poly_gcd_mod(f, a, p);
        if (g.degree() > 0) {
            if (g.degree() % d == 0) {
                for (auto&& fact : equal_degree_factorization(g, d, p)) factors.push_back(fact);
                f = poly_div_rem_mod(f, g, p).first;
                continue;
            }
        }
        
        // b = a^((p^d-1)/2) - 1 mod f
        if (p > BigInt(2)) {
            BigInt exponent(1);
            for(std::size_t i=0; i<d; ++i) exponent *= p;
            exponent = (exponent - BigInt(1)) / BigInt(2);
            
            IntPoly b = poly_pow_mod_poly(a, exponent, f, p);
            b = poly_sub_mod(b, IntPoly{std::vector<BigInt>{BigInt(1)}}, p);
            
            g = poly_gcd_mod(f, b, p);
            if (g.degree() > 0 && g.degree() < f.degree()) {
                for (auto&& fact : equal_degree_factorization(g, d, p)) factors.push_back(fact);
                f = poly_div_rem_mod(f, g, p).first;
            }
        } else {
            // p = 2 case: trace polynomial T(a) = a + a^2 + a^4 + ... + a^(2^(d-1)) mod f
            IntPoly tr = a;
            IntPoly current_a = a;
            for (std::size_t i = 1; i < d; ++i) {
                // current_a = current_a^2 mod f
                current_a = poly_mul_mod(current_a, current_a, p);
                current_a = poly_div_rem_mod(current_a, f, p).second;
                tr = poly_add_mod(tr, current_a, p);
            }
            
            g = poly_gcd_mod(f, tr, p);
            if (g.degree() > 0 && g.degree() < f.degree()) {
                for (auto&& fact : equal_degree_factorization(g, d, p)) factors.push_back(fact);
                f = poly_div_rem_mod(f, g, p).first;
            }
        }
    }
    factors.push_back(f);
    return factors;
}

Result<std::vector<IntPoly>> factor_polynomial_mod_p(IntPoly f, const BigInt& p) {
    if (f.is_zero()) return ok(std::vector<IntPoly>{});
    
    BigInt lc = f.leading_coeff();
    auto inv_lc_res = numtheory::modular_inverse(mod_bigint(lc, p), p);
    if (inv_lc_res.is_error()) return fail<std::vector<IntPoly>>(inv_lc_res.error());
    
    BigInt inv_lc = inv_lc_res.value();
    for (auto& c : f.coefficients()) c = mod_bigint(c * inv_lc, p);
    f.normalize([](const BigInt& v) { return v.is_zero(); });

    auto ddf = distinct_degree_factorization(f, p);
    std::vector<IntPoly> all_factors;
    for (auto& [g, d] : ddf) {
        if (g.degree() == d) {
            all_factors.push_back(g);
        } else {
            auto edf = equal_degree_factorization(g, d, p);
            all_factors.insert(all_factors.end(), edf.begin(), edf.end());
        }
    }
    return ok(all_factors);
}

} // namespace cas::algebra
