#include "polynomial_internal.hpp"
#include "algebra_internal.hpp"
#include "cas/numtheory.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

namespace cas::algebra {

namespace {

[[nodiscard]] std::size_t poly_hash_fnv1a(const IntPoly& f) {
    constexpr std::size_t kFnvPrime = 0x00000100000001B3ULL;
    constexpr std::size_t kFnvBasis = 0xcbf29ce484222325ULL;
    std::size_t h = kFnvBasis;
    for (std::size_t i = 0; i < f.size(); ++i) {
        const std::string& dec = f[i].decimal();
        for (unsigned char c : dec) {
            h ^= static_cast<std::size_t>(c);
            h *= kFnvPrime;
        }
        h ^= 0xABU;
        h *= kFnvPrime;
    }
    return h;
}

[[nodiscard]] bool is_small_prime(std::uint64_t n) {
    if (n < 2U) return false;
    if (n == 2U) return true;
    if ((n % 2U) == 0U) return false;
    for (std::uint64_t d = 3U; d <= n / d; d += 2U) {
        if ((n % d) == 0U) return false;
    }
    return true;
}

[[nodiscard]] std::vector<BigInt> get_prime_candidates(const IntPoly& f, std::size_t count) {
    static constexpr int kPool[] = {
        13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
        53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
        101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
        151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199
    };
    static constexpr std::size_t kPoolSize = sizeof(kPool) / sizeof(kPool[0]);
    
    std::vector<BigInt> result;
    BigInt lc = f.leading_coeff();
    
    std::size_t start = poly_hash_fnv1a(f) % kPoolSize;
    for (std::size_t i = 0; i < kPoolSize && result.size() < count; ++i) {
        BigInt bp(kPool[(start + i) % kPoolSize]);
        if (!(lc % bp).is_zero()) {
            result.push_back(bp);
        }
    }
    
    std::uint64_t candidate = 201U;
    while (result.size() < count) {
        if (is_small_prime(candidate)) {
            BigInt bp(static_cast<long long>(candidate));
            if (!(lc % bp).is_zero()) {
                result.push_back(bp);
            }
        }
        candidate += 2U;
    }
    return result;
}

// Is f squarefree mod p ("lucky" prime)?  True iff gcd(f mod p, f' mod p) is a
// constant, i.e. p ∤ disc(f).  On a lucky prime the mod-p factorisation lifts
// faithfully to the Z-factorisation, so the exhaustive Hensel recombination is a
// COMPLETE factor search: a "no factor" result is then a proof of irreducibility.
[[nodiscard]] bool integer_poly_squarefree_mod_p(const IntPoly& f, long long p) {
    std::vector<long long> a;
    for (const auto& c : f.coefficients()) {
        long long m = static_cast<long long>((c % BigInt(p)).to_u64());
        if ((c % BigInt(p)).is_negative()) m = (m + p) % p;
        a.push_back(((m % p) + p) % p);
    }
    while (!a.empty() && a.back() == 0) a.pop_back();
    if (a.size() <= 1) return true;
    std::vector<long long> b(a.size() - 1, 0);
    for (std::size_t i = 1; i < a.size(); ++i) b[i-1] = (a[i] * (long long)i) % p;
    while (!b.empty() && b.back() == 0) b.pop_back();
    if (b.empty()) return false;
    auto inv = [p](long long x){ long long r=1,bs=((x%p)+p)%p,e=p-2; while(e>0){ if(e&1) r=r*bs%p; bs=bs*bs%p; e>>=1;} return r; };
    while (!b.empty()) {
        long long il = inv(b.back());
        while (a.size() >= b.size()) {
            long long f2 = a.back()*il%p;
            std::size_t sh = a.size()-b.size();
            for (std::size_t i=0;i<b.size();++i) a[i+sh]=((a[i+sh]-f2*b[i])%p+p)%p;
            while(!a.empty() && a.back()==0) a.pop_back();
        }
        std::swap(a,b);
    }
    return a.size() <= 1;
}

} // namespace

Result<std::vector<IntPoly>> factorize_univariate_hensel_or_kronecker(
    const IntPoly& f,
    symbolic::CASContext& ctx) {

    if (f.empty() || f.degree() == 0U) {
        return ok(std::vector<IntPoly>{});
    }
    if (f.degree() == 1U) {
        return ok(std::vector<IntPoly>{f});
    }

    const std::size_t max_attempts = ctx.max_hensel_lift_attempts();
    std::vector<BigInt> primes = get_prime_candidates(f, max_attempts);

    std::size_t attempts = 0;
    std::size_t bad_primes = 0;

    for (const auto& p : primes) {
        attempts++;
        
        // 1. Factor mod p
        const std::size_t deg_p_product = f.degree() * static_cast<std::size_t>(p.to_u64());
        auto mod_factors_res = (deg_p_product <= ctx.max_berlekamp_matrix_size())
            ? berlekamp_factor_mod_p(f, p, ctx.max_berlekamp_matrix_size())
            : factor_polynomial_mod_p(f, p, &ctx);
        if (mod_factors_res.is_error()) {
            mod_factors_res = factor_polynomial_mod_p(f, p, &ctx);
        }

        if (mod_factors_res.is_error()) {
            bad_primes++;
            continue;
        }

        const auto& mod_factors = mod_factors_res.value();
        
        // If irreducible modulo p, it's irreducible over Z
        if (mod_factors.size() == 1U) {
            return ok(std::vector<IntPoly>{f});
        }

        // 2. Determine target modulus p^k via Mignotte bound
        std::size_t n = f.degree();
        BigInt pk = p;
        std::size_t k = 1;
        BigInt two_pow_n = BigInt(1).shift_left_bits(n);
        BigInt norm2_sq(0);
        for (const auto& c : f.coefficients()) norm2_sq += c * c;
        while (pk < two_pow_n * norm2_sq * BigInt(2)) {
            pk *= p;
            k++;
        }

        // 3. Hensel lift
        auto lifted_res = hensel_lift_multi(f, mod_factors, p, k);
        if (lifted_res.is_error()) {
            bad_primes++;
            continue;
        }

        const auto& lifted = lifted_res.value();

        // 4. Try recombination.  van Hoeij (fast, but in this implementation not
        // provably complete — it can miss a real factor) is always backed by the
        // EXHAUSTIVE Hensel subset search, which IS complete up to degree n/2.
        // So after this block, an empty found_factor means the exhaustive search
        // found nothing.
        std::optional<IntPoly> found_factor;
        if (lifted.size() >= ctx.van_hoeij_threshold()) {
            found_factor = van_hoeij_knapsack_factor(
                f, lifted, pk, ctx.lll_delta(), ctx.van_hoeij_lll_threshold(), &ctx);
            if (!found_factor.has_value()) {
                found_factor = find_factor_by_hensel_recombination(
                    f, mod_factors, p, k, f.degree() / 2U);
            }
        } else {
            found_factor = find_factor_by_hensel_recombination(
                f, mod_factors, p, k, f.degree() / 2U);
        }

        if (found_factor.has_value()) {
            // Validate and divide
            auto q_res = exact_divide_integer_poly(f, found_factor.value(), ctx);
            if (q_res.is_ok()) {
                // Recursively factor both found_factor and the quotient
                auto f1_res = factorize_univariate_hensel_or_kronecker(found_factor.value(), ctx);
                auto f2_res = factorize_univariate_hensel_or_kronecker(q_res.value(), ctx);
                if (f1_res.is_ok() && f2_res.is_ok()) {
                    std::vector<IntPoly> result = f1_res.value();
                    result.insert(result.end(), f2_res.value().begin(), f2_res.value().end());
                    return ok(result);
                }
            }
        }

        // No factor found.  On a LUCKY prime (f squarefree mod p) the exhaustive
        // Hensel search above is a complete factor search, so "no factor" is a
        // PROOF that f is irreducible over Z — return it directly instead of
        // treating the good prime as bad and falling through to the exponential
        // Kronecker fallback (the SD3 Swinnerton-Dyer hang, HC-F8-SD3-VANHOEIJ-SLOW).
        // We trust ONLY the exhaustive Hensel verdict on a lucky prime, never van
        // Hoeij's (which can miss factors), and require squarefree-mod-p so the
        // mod-p factorisation lifts faithfully.
        if (integer_poly_squarefree_mod_p(f, static_cast<long long>(p.to_u64()))) {
            return ok(std::vector<IntPoly>{f});
        }

        // Unlucky prime (p | disc): the verdict is unreliable, try another prime.
        bad_primes++;
    }

    // If we failed after all attempts, compute bad prime rate
    double bad_prime_rate = (attempts > 0) ? (static_cast<double>(bad_primes) / attempts) : 0.0;

    if (bad_prime_rate > 0.5 || attempts == 0) {
        if (f.degree() <= ctx.kronecker_max_degree()) {
            return factorize_kronecker(f, ctx);
        } else {
            return fail<std::vector<IntPoly>>(make_error(
                CASErrorKind::Unimplemented,
                "Hensel lifting failed after max attempts and degree exceeds Kronecker limit"));
        }
    }

    return fail<std::vector<IntPoly>>(make_error(
        CASErrorKind::Unimplemented,
        "Hensel lifting failed to find factors"));
}

} // namespace cas::algebra
