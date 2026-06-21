// polynomial_gcd_zippel_crt.cpp — multi-prime CRT + Farey reconstruction for the
// Zippel-Prony sparse GCD (F3.1 / T-006). Spec: Zippel_Sparse_Interpolation.md §146.
//
// Invoked only when the single-prime center-lift fails its exact-division
// certificate (the true gcd's x_1-leading coefficient is not a unit, so the
// monic-over-F_p coefficients are rationals a_i/lc). We CRT-combine the monic
// gcd mod p across a sequence of primes, Farey-reconstruct each coefficient to an
// exact rational, clear denominators to the primitive integer gcd, and CERTIFY by
// exact division. The certificate guarantees no wrong result ever escapes; if the
// derived prime budget is exhausted without certifying, we return a diagnostic
// Unimplemented and the dispatcher falls back to Brown's modular GCD.

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_gcd_zippel_internal.hpp"

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra::zippel_detail {
namespace {

BigInt pos_mod(const BigInt& a, const BigInt& m) {
    BigInt r = a % m;
    if (r.is_negative()) r += m.abs();
    return r;
}

// floor(sqrt(n)) for n ≥ 0 (Newton).
BigInt isqrt(const BigInt& n) {
    if (n.is_negative() || n.is_zero()) return BigInt(0);
    if (n < BigInt(4)) return BigInt(1);
    BigInt x(1);
    for (std::size_t i = 0; i < (n.bit_length() + 2U) / 2U; ++i) x = x + x;  // 2^⌈b/2⌉
    while (true) {
        BigInt y = (x + n / x) / BigInt(2);
        if (!(y < x)) break;
        x = y;
    }
    while (x * x > n) x = x - BigInt(1);
    return x;
}

// Incremental Garner CRT: given x ≡ a (mod M) and x ≡ b (mod p) with gcd(M,p)=1,
// return x in [0, M·p).
Result<BigInt> crt_step(const BigInt& a, const BigInt& M, const BigInt& b, const BigInt& p) {
    auto inv = numtheory::modular_inverse(pos_mod(M, p), p);
    if (inv.is_error()) return fail<BigInt>(inv.error());
    BigInt t = pos_mod((pos_mod(b, p) - pos_mod(a, p)) * inv.value(), p);
    return ok(a + M * t);
}

// Wang rational reconstruction: find num/den ≡ u (mod m) in lowest terms with
// |num| ≤ N and 0 < den ≤ N. Returns false when no such rational exists. The
// congruence is re-verified exactly before returning.
bool rational_reconstruct(const BigInt& u, const BigInt& m, const BigInt& N,
                          BigInt& num, BigInt& den) {
    BigInt uu = pos_mod(u, m);
    if (uu.is_zero()) { num = BigInt(0); den = BigInt(1); return true; }
    BigInt r0 = m, r1 = uu, t0(0), t1(1);
    while (r1 > N) {
        BigInt q = r0 / r1;
        BigInt r2 = r0 - q * r1; r0 = r1; r1 = r2;
        BigInt t2 = t0 - q * t1; t0 = t1; t1 = t2;
    }
    if (t1.is_zero()) return false;
    BigInt nn = r1, dd = t1;
    if (dd.is_negative()) { nn = -nn; dd = -dd; }
    if (dd > N) return false;
    if (cas::gcd(nn.abs(), dd) != BigInt(1)) return false;
    if (!pos_mod(nn - uu * dd, m).is_zero()) return false;
    num = nn; den = dd;
    return true;
}

}  // namespace

Result<MultivariatePolynomial> gcd_zippel_prony_crt(
    const ZSparsePoly& spP, const ZSparsePoly& spQ,
    const std::vector<Symbol>& vars, const BigInt& first_prime,
    symbolic::CASContext& ctx, std::size_t* out_samples) {

    const MultivariatePolynomial P_poly = from_sparse_z(spP, vars);
    const MultivariatePolynomial Q_poly = from_sparse_z(spQ, vars);

    // Prime budget from a Mignotte-style bound: |gcd coeff| ≤ 2^dtot · cmax, and
    // Farey needs M > 2·B²; primes are ~30-bit. No magic constant — derived from
    // the actual input coefficient sizes and total degree.
    BigInt cmax(1);
    std::size_t dtot = 0;
    auto scan = [&](const ZSparsePoly& sp) {
        for (const auto& [m, c] : sp) {
            if (c.abs() > cmax) cmax = c.abs();
            std::size_t td = 0; for (auto e : m) td += e;
            if (td > dtot) dtot = td;
        }
    };
    scan(spP); scan(spQ);
    const std::size_t bits_needed = 2U * (dtot + cmax.bit_length()) + 64U;
    const std::size_t max_primes = bits_needed / 29U + 6U;
    const std::size_t prime_cap = max_primes * 4U + 8U;  // headroom for unlucky skips

    std::map<ZMonomial, BigInt> acc;   // residue mod M, per monomial
    BigInt M(1);
    std::set<ZMonomial> support;
    bool have = false;
    BigInt p = first_prime;
    std::size_t total_samples = 0, good = 0, tried = 0;

    while (good < max_primes && tried < prime_cap) {
        auto np = numtheory::next_prime(p);
        if (np.is_error()) break;
        p = np.value();
        ++tried;

        std::size_t s = 0;
        auto gp = zippel_gcd_modp(spP, spQ, vars, p, ctx, &s);
        total_samples += s;
        if (gp.is_error()) continue;             // unlucky prime → try next
        const ZSparsePoly& gm = gp.value();

        std::set<ZMonomial> sup;
        for (const auto& [m, c] : gm) sup.insert(m);
        if (!have) {
            support = sup;
            acc.clear();
            for (const auto& [m, c] : gm) acc[m] = c;
            M = p; have = true; ++good;
        } else {
            if (sup != support) continue;        // skeleton differs → unlucky, skip
            std::map<ZMonomial, BigInt> next;
            bool ok_step = true;
            for (const auto& [m, a] : acc) {
                auto cs = crt_step(a, M, gm.at(m), p);
                if (cs.is_error()) { ok_step = false; break; }
                next[m] = cs.value();
            }
            if (!ok_step) continue;
            acc = std::move(next);
            M = M * p; ++good;
        }
        if (good < 2U) continue;                 // need ≥2 primes before reconstructing

        // Farey-reconstruct each coefficient to a rational, then clear denominators.
        const BigInt N = isqrt(M / BigInt(2));
        std::map<ZMonomial, std::pair<BigInt, BigInt>> rats;
        BigInt denom_lcm(1);
        bool all_ok = true;
        for (const auto& [m, u] : acc) {
            BigInt num, den;
            if (!rational_reconstruct(u, M, N, num, den)) { all_ok = false; break; }
            rats[m] = {num, den};
            denom_lcm = denom_lcm / cas::gcd(denom_lcm, den) * den;
        }
        if (!all_ok) continue;

        ZSparsePoly g_int;
        for (const auto& [m, nd] : rats) {
            BigInt scaled = nd.first * (denom_lcm / nd.second);
            if (!scaled.is_zero()) g_int[m] = scaled;
        }
        if (g_int.empty()) continue;

        MultivariatePolynomial g_cand = from_sparse_z(g_int, vars);
        const BigInt content = g_cand.integer_content();   // primitive part
        if (!content.is_zero() && content.abs() != BigInt(1)) {
            ZSparsePoly prim;
            for (const auto& [m, c] : g_int) prim[m] = c / content;
            g_cand = from_sparse_z(prim, vars);
        }
        {   // sign-normalize: lex-leading positive
            ZSparsePoly sp = to_sparse_z(g_cand, vars);
            if (!sp.empty()) {
                if (std::prev(sp.end())->second.is_negative())
                    for (auto& [m, c] : sp) c = -c;
                g_cand = from_sparse_z(sp, vars);
            }
        }
        if (certify_divides(P_poly, g_cand, vars) &&
            certify_divides(Q_poly, g_cand, vars)) {
            if (out_samples) *out_samples = total_samples;
            return ok(std::move(g_cand));
        }
        // Certificate not yet satisfied — accumulate more primes.
    }

    if (out_samples) *out_samples = total_samples;
    return make_unimplemented<MultivariatePolynomial>(
        "algebra", "gcd_zippel_prony_crt",
        "primes=" + std::to_string(good) + ",tried=" + std::to_string(tried),
        "ZIPPEL_CRT_EXHAUSTED",
        "Multi-prime CRT did not certify within the derived prime budget", "F3.1");
}

}  // namespace cas::algebra::zippel_detail
