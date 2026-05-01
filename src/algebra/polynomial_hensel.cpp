#include "polynomial_internal.hpp"
#include "cas/numtheory.hpp"
#include <vector>
#include <algorithm>
#include <tuple>

namespace cas::algebra {

// Helper for modular arithmetic on BigInt
static BigInt mod_bigint(const BigInt& a, const BigInt& m) {
    if (m.is_zero()) return a;
    BigInt r = a % m;
    if (r.is_negative()) {
        r += m.abs();
    }
    return r;
}

// Helper for polynomial modular reduction
static void poly_mod(IntPoly& p, const BigInt& m) {
    for (auto& c : p.coefficients()) {
        c = mod_bigint(c, m);
    }
    p.normalize([](const BigInt& v) { return v.is_zero(); });
}

// Polynomial addition modulo m
static IntPoly poly_add_mod(const IntPoly& a, const IntPoly& b, const BigInt& m) {
    std::size_t n = std::max(a.size(), b.size());
    IntPoly r;
    r.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        BigInt sum(0);
        if (i < a.size()) sum += a[i];
        if (i < b.size()) sum += b[i];
        r[i] = mod_bigint(sum, m);
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Polynomial subtraction modulo m
static IntPoly poly_sub_mod(const IntPoly& a, const IntPoly& b, const BigInt& m) {
    std::size_t n = std::max(a.size(), b.size());
    IntPoly r;
    r.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        BigInt diff(0);
        if (i < a.size()) diff += a[i];
        if (i < b.size()) diff -= b[i];
        r[i] = mod_bigint(diff, m);
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Polynomial multiplication modulo m
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

// Polynomial division with remainder modulo p (p must be prime)
static std::pair<IntPoly, IntPoly> poly_div_rem_mod(const IntPoly& a, const IntPoly& b, const BigInt& p) {
    if (b.is_zero()) return {IntPoly{}, a};
    
    IntPoly q, r = a;
    poly_mod(r, p);
    if (r.is_zero()) return {IntPoly{}, IntPoly{}};

    auto inv_lb_res = numtheory::modular_inverse(mod_bigint(b.leading_coeff(), p), p);
    if (inv_lb_res.is_error()) return {IntPoly{}, r}; // Should not happen if p is prime and b normalized
    const BigInt inv_lb = inv_lb_res.value();
    
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

// Extended GCD for polynomials modulo p
static std::tuple<IntPoly, IntPoly, IntPoly> poly_gcd_ext_mod(const IntPoly& a, const IntPoly& b, const BigInt& p) {
    if (a.is_zero()) {
        IntPoly g = b;
        poly_mod(g, p);
        if (!g.is_zero()) {
            BigInt inv = numtheory::modular_inverse(g.leading_coeff(), p).value();
            for (auto& c : g.coefficients()) c = mod_bigint(c * inv, p);
            return {g, IntPoly{std::vector<BigInt>{BigInt(0)}}, IntPoly{std::vector<BigInt>{inv}}};
        }
        return {g, IntPoly{}, IntPoly{}};
    }
    
    IntPoly r0 = a, r1 = b;
    poly_mod(r0, p);
    poly_mod(r1, p);
    
    IntPoly s0{std::vector<BigInt>{BigInt(1)}}, s1{std::vector<BigInt>{BigInt(0)}};
    IntPoly t0{std::vector<BigInt>{BigInt(0)}}, t1{std::vector<BigInt>{BigInt(1)}};
    
    while (!r1.is_zero()) {
        auto [q, rem] = poly_div_rem_mod(r0, r1, p);
        r0 = r1;
        r1 = rem;
        
        IntPoly next_s = poly_sub_mod(s0, poly_mul_mod(q, s1, p), p);
        s0 = s1;
        s1 = next_s;
        
        IntPoly next_t = poly_sub_mod(t0, poly_mul_mod(q, t1, p), p);
        t0 = t1;
        t1 = next_t;
    }
    
    if (!r0.is_zero()) {
        BigInt inv = numtheory::modular_inverse(r0.leading_coeff(), p).value();
        for (auto& c : r0.coefficients()) c = mod_bigint(c * inv, p);
        for (auto& c : s0.coefficients()) c = mod_bigint(c * inv, p);
        for (auto& c : t0.coefficients()) c = mod_bigint(c * inv, p);
    }
    
    return {r0, s0, t0};
}

// Univariate Hensel Lift: lifts f = g*h mod p to f = G*H mod p^k
Result<std::pair<IntPoly, IntPoly>> hensel_lift(
    const IntPoly& f, 
    const IntPoly& g, 
    const IntPoly& h, 
    const BigInt& p, 
    std::size_t k) {
    
    BigInt modulus = p;
    IntPoly G = g, H = h;
    poly_mod(G, p);
    poly_mod(H, p);
    
    // Find a, b such that aG + bH = 1 mod p
    auto [gcd_p, a, b] = poly_gcd_ext_mod(G, H, p);
    
    BigInt target_modulus(1);
    for(std::size_t i=0; i<k; ++i) target_modulus *= p;
    
    while (modulus < target_modulus) {
        BigInt next_modulus = modulus * p; // Linear lifting for simplicity and robustness
        if (next_modulus > target_modulus) next_modulus = target_modulus;

        // E = (f - GH) / modulus mod p
        IntPoly GH = poly_mul_mod(G, H, next_modulus);
        IntPoly diff = poly_sub_mod(f, GH, next_modulus);
        
        IntPoly e;
        e.resize(diff.size());
        for (std::size_t i = 0; i < diff.size(); ++i) {
            e[i] = mod_bigint(diff[i] / modulus, p);
        }
        e.normalize([](const BigInt& v) { return v.is_zero(); });
        
        if (!e.is_zero()) {
            // Solve G*delta_h + H*delta_g = e mod p
            // delta_g = (b*e) mod G
            auto [qg, dg] = poly_div_rem_mod(poly_mul_mod(b, e, p), G, p);
            // delta_h = (a*e + qg*H) mod p
            IntPoly dh = poly_add_mod(poly_mul_mod(a, e, p), poly_mul_mod(qg, H, p), p);
            
            // G = G + dg * modulus
            if (G.size() < dg.size()) G.resize(dg.size(), BigInt(0));
            for (std::size_t i = 0; i < dg.size(); ++i) {
                G[i] = mod_bigint(G[i] + dg[i] * modulus, next_modulus);
            }
            
            // H = H + dh * modulus
            if (H.size() < dh.size()) H.resize(dh.size(), BigInt(0));
            for (std::size_t i = 0; i < dh.size(); ++i) {
                H[i] = mod_bigint(H[i] + dh[i] * modulus, next_modulus);
            }
        }
        
        G.normalize([](const BigInt& v) { return v.is_zero(); });
        H.normalize([](const BigInt& v) { return v.is_zero(); });
        modulus = next_modulus;
    }
    
    return ok(std::make_pair(G, H));
}

// Multi-factor Hensel lifting
Result<std::vector<IntPoly>> hensel_lift_multi(
    const IntPoly& f,
    const std::vector<IntPoly>& factors,
    const BigInt& p,
    std::size_t k) {
    
    if (factors.size() <= 1) return ok(factors);
    
    // Split factors into two groups
    std::size_t mid = factors.size() / 2;
    std::vector<IntPoly> left_factors(factors.begin(), factors.begin() + mid);
    std::vector<IntPoly> right_factors(factors.begin() + mid, factors.end());
    
    // Compute product of each group mod p
    IntPoly g{std::vector<BigInt>{BigInt(1)}}, h{std::vector<BigInt>{BigInt(1)}};
    for (const auto& poly : left_factors) g = poly_mul_mod(g, poly, p);
    for (const auto& poly : right_factors) h = poly_mul_mod(h, poly, p);
    
    // Lift f = g*h mod p to f = G*H mod p^k
    auto lift_res = hensel_lift(f, g, h, p, k);
    if (lift_res.is_error()) return fail<std::vector<IntPoly>>(lift_res.error());
    
    auto [G, H] = lift_res.value();
    
    // Recursively lift groups
    auto left_lifted = hensel_lift_multi(G, left_factors, p, k);
    if (left_lifted.is_error()) return fail<std::vector<IntPoly>>(left_lifted.error());
    
    auto right_lifted = hensel_lift_multi(H, right_factors, p, k);
    if (right_lifted.is_error()) return fail<std::vector<IntPoly>>(right_lifted.error());
    
    std::vector<IntPoly> result = left_lifted.value();
    result.insert(result.end(), right_lifted.value().begin(), right_lifted.value().end());
    return ok(result);
}

} // namespace cas::algebra
