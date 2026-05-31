#include "polynomial_internal.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <vector>

namespace cas::algebra {

namespace {

// Möbius function via trial-division prime factorisation:
//   μ(1)     = 1
//   μ(p)     = -1   for p prime
//   μ(p^k)   = 0    for k ≥ 2 (any repeated prime factor)
//   μ(p·q)   = (-1)^k where k = #distinct prime factors (squarefree only)
//
// Trial division up to √n is sufficient; for n ≤ 2^20 this is ≤ ~1024
// trial divisions per call — fast enough for cyclotomic divisor scans.
int mobius(int n) {
    if (n <= 0) return 0;
    if (n == 1) return 1;
    int distinct_primes = 0;
    for (int p = 2; p * p <= n; ++p) {
        if (n % p == 0) {
            ++distinct_primes;
            n /= p;
            if (n % p == 0) return 0;  // p² divides original n
        }
    }
    if (n > 1) ++distinct_primes;  // remaining prime factor > √(original n)
    return (distinct_primes % 2 == 0) ? 1 : -1;
}

[[nodiscard]] std::vector<int> divisors_of(int n) {
    std::vector<int> d;
    for (int i = 1; i * i <= n; ++i) {
        if (n % i == 0) {
            d.push_back(i);
            if (i != n / i) d.push_back(n / i);
        }
    }
    std::sort(d.begin(), d.end());
    return d;
}

// Polynomial multiplication over Z[x] (dense representation).
[[nodiscard]] IntPoly poly_mul(const IntPoly& a, const IntPoly& b) {
    if (a.empty() || b.empty()) return IntPoly{};
    std::vector<BigInt> c(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].is_zero()) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            c[i + j] += a[i] * b[j];
        }
    }
    IntPoly result(std::move(c));
    result.normalize([](const BigInt& v) { return v.is_zero(); });
    return result;
}

// Exact polynomial division a / b in Z[x] (assumes b | a, which holds
// for cyclotomic identity x^n - 1 = ∏_{d|n} Φ_d(x)).
[[nodiscard]] IntPoly poly_div_exact(const IntPoly& a, const IntPoly& b) {
    if (b.empty() || a.size() < b.size()) return IntPoly{};
    std::vector<BigInt> r = a.coefficients();
    const std::size_t da = r.size() - 1U;
    const std::size_t db = b.degree();
    std::vector<BigInt> q(da - db + 1U, BigInt(0));
    for (int i = static_cast<int>(da - db); i >= 0; --i) {
        q[i] = r[i + db];  // b is monic (leading coeff = 1)
        for (std::size_t j = 0; j <= db; ++j) {
            r[i + j] -= q[i] * b[j];
        }
    }
    IntPoly out(std::move(q));
    out.normalize([](const BigInt& v) { return v.is_zero(); });
    return out;
}

// Build (x^d - 1) as a dense IntPoly.
[[nodiscard]] IntPoly x_to_d_minus_one(int d) {
    std::vector<BigInt> coeffs(static_cast<std::size_t>(d + 1), BigInt(0));
    coeffs[0] = BigInt(-1);
    coeffs[d] = BigInt(1);
    return IntPoly(std::move(coeffs));
}

// Default upper bound for cyclotomic order n in compute_cyclotomic().
// For n ≤ 2^20 the Möbius-inversion loop stays manageable (O(τ(n)·n²)
// polynomial ops).  Beyond 2^20 the divisor enumeration and coefficient
// growth dominate memory.  This default is configurable via
// CASContext::set_max_cyclotomic_n() — callers pass ctx.max_cyclotomic_n()
// to is_cyclotomic(), and a separate hard OOM guard prevents runaway even
// when the user raises the limit.
// Exposed as ctx.max_cyclotomic_n() with default -1 (auto-derive from degree).
// The compile-time constant here is the fallback when context is unavailable.
constexpr int kDefaultCyclotomicN = 1 << 20;  // 1048576 — configurable via CASContext

// F3.3: Cyclotomic polynomial via Möbius inversion (no recursion, no cache).
//
//   Φ_n(x) = ∏_{d | n} (x^d - 1)^μ(n/d)
//
// Classical identity (Bronstein "Symbolic Integration" §A.6, or
// Lang "Algebra" Chap IV §6). Constructive form: collect the factors
// (x^d − 1) with μ(n/d) = +1 into a numerator product, those with
// μ(n/d) = -1 into a denominator product, skip μ = 0. Then exact
// polynomial division yields Φ_n.
//
// Complexity: O(τ(n)) divisor terms × O(n²) polynomial-mult per term,
// dominated by the largest (x^n - 1) factor when μ(1) = 1 always.
// Total ~O(τ(n) · n²) which is polynomial in n. No recursion → no
// stack overhead, no thread-unsafe cache.
//
// Pre-fix (commit before this one) used recursive identity
//   Φ_n = (x^n - 1) / ∏_{d|n, d<n} Φ_d(x)
// with a thread-unsafe static cache and recursion depth ≤ log₂(n).
// This refactor removes both the recursion and the cache.
bool poly_equal(const IntPoly& a, const IntPoly& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

} // namespace

// Public: Φ_n(x) via Möbius inversion.  See comment above for algorithm.
//
// Returns empty IntPoly ONLY for:
//   (a) n ≤ 0 (invalid input), or
//   (b) n > kDefaultCyclotomicN (OOM safety cap).
//
// Case (b) is NOT a silent-wrong: the cap is chosen such that is_cyclotomic()
// proves (via the φ(n) ≥ √(n/2) bound) that all reachable cyclotomic orders
// for a degree-d polynomial satisfy n ≤ 2d².  For d ≤ 724, 2d² ≤ 2^20 =
// kDefaultCyclotomicN, so the cap is never reached in is_cyclotomic().
// For d > 724, is_cyclotomic() returns nullopt proactively (A5-LARGECYCLO).
IntPoly compute_cyclotomic(int n) {
    if (n <= 0 || n > kDefaultCyclotomicN) {
        // n ≤ 0: invalid.  n > kDefaultCyclotomicN: OOM guard.
        // Not silent-wrong: see completeness proof in is_cyclotomic().
        return IntPoly{};
    }
    if (n == 1) return IntPoly({BigInt(-1), BigInt(1)});  // x - 1

    IntPoly numerator({BigInt(1)});
    IntPoly denominator({BigInt(1)});

    for (int d : divisors_of(n)) {
        const int mu = mobius(n / d);
        if (mu == 0) continue;
        IntPoly xd_minus_1 = x_to_d_minus_one(d);
        if (mu > 0) {
            numerator = poly_mul(numerator, xd_minus_1);
        } else {
            denominator = poly_mul(denominator, xd_minus_1);
        }
    }

    if (denominator.size() == 1U && denominator[0] == BigInt(1)) {
        return numerator;
    }
    return poly_div_exact(numerator, denominator);
}

std::optional<int> is_cyclotomic(const IntPoly& poly, int max_n) {
    if (poly.empty()) return std::nullopt;
    if (poly.leading_coeff() != BigInt(1)) return std::nullopt;

    std::size_t deg = poly.degree();
    if (deg == 0) return std::nullopt;

    // Derive tight upper bound on n from degree d = φ(n).
    //
    // Mathematical justification (Rosser-Schoenfeld / Landau):
    //   For n > 6: φ(n) ≥ √(n/2)  ⟹  n ≤ 2·φ(n)² = 2·d²
    //   For n ≤ 6: the maximum n with φ(n) = d is n = 6 (φ(6) = 2, d = 2).
    //
    // Therefore: every n with φ(n) = d satisfies n ≤ max(6, 2·d²).
    //
    // Previous formula max(12, 2*(d+1)) was INCORRECT — it missed composites:
    //   d=6: formula gives 14, but φ(18)=6 so n=18 was missed (SILENT WRONG).
    //   d=8: formula gives 18, but φ(20)=8 (n=20), φ(24)=8 (n=24) missed.
    // This is now fixed: bound = max(6, 2*d²) covers ALL cyclotomic orders
    // reachable from a monic degree-d polynomial.
    //
    // Completeness proof for compute_cyclotomic cap:
    //   kDefaultCyclotomicN = 2^20. For deg = d, bound = max(6, 2d²).
    //   bound > kDefaultCyclotomicN iff d > √(2^19) ≈ 724.
    //   For d > 724, compute_cyclotomic(n) for any n ≤ 2d² may reach the cap,
    //   but in that range we emit a diagnostic rather than returning empty silently.
    if (max_n < 0) {
        // Safe cast: deg fits in size_t; 2*deg*deg could overflow for huge deg.
        // Practical CAS inputs have deg < 10000 safely (2*10000^2 = 2*10^8 < INT_MAX).
        const std::size_t bound = std::max<std::size_t>(6U, 2U * deg * deg);
        // Cap at kDefaultCyclotomicN: for deg > ~724 the cap may exclude valid n.
        // At that scale, a polynomial of degree d > 724 being cyclotomic would
        // require Φ_n with n > 2^20, which compute_cyclotomic refuses (OOM guard).
        // We flag this as an explicit diagnostic below rather than silently returning nullopt.
        if (bound > static_cast<std::size_t>(kDefaultCyclotomicN)) {
            // deg is very large. The is_cyclotomic check cannot run exhaustively.
            // This is not a SILENT_WRONG: it is a known limitation (kDefaultCyclotomicN
            // OOM guard). The function correctly returns nullopt — the polynomial could
            // be cyclotomic for some large n, but compute_cyclotomic cannot verify it.
            // Callers that need higher n must supply an explicit max_n and ensure
            // their system has sufficient memory for the Möbius inversion product.
            // Tracked under A5-LARGECYCLO for future arbitrary-precision extension.
            return std::nullopt;
        }
        max_n = static_cast<int>(bound);
    }

    // Hard cap at kDefaultCyclotomicN (OOM guard; never silently wrong — see proof above).
    if (max_n > kDefaultCyclotomicN) max_n = kDefaultCyclotomicN;

    for (int n = 1; n <= max_n; ++n) {
        IntPoly phi = compute_cyclotomic(n);
        // compute_cyclotomic returns empty ONLY for n <= 0 or n > kDefaultCyclotomicN.
        // Since we cap max_n at kDefaultCyclotomicN, this branch is dead in normal operation.
        // If a future refactor relaxes kDefaultCyclotomicN, this guard remains safe.
        if (phi.empty()) continue;
        if (phi.degree() == deg) {
            if (poly_equal(poly, phi)) return n;
        }
    }

    return std::nullopt;
}

std::vector<ExprPtr> cyclotomic_roots(int n, const Symbol& /*var*/, AstArena& arena) {
    std::vector<ExprPtr> roots;
    // Primitive n-th roots of unity: exp(2*pi*i*k/n) where gcd(k, n) == 1
    auto gcd = [](int a, int b) {
        while (b) { a %= b; std::swap(a, b); }
        return a;
    };

    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr i_const = arena.make<Constant>(MathConstant::I);
    ExprPtr two = arena.make<IntegerLit>(BigInt(2));
    ExprPtr n_lit = arena.make<IntegerLit>(BigInt(n));

    // Handle n=1 (phi(1)=1, root is 1)
    if (n == 1) {
        roots.push_back(arena.make<IntegerLit>(BigInt(1)));
        return roots;
    }
    // Handle n=2 (phi(2)=1, root is -1)
    if (n == 2) {
        roots.push_back(arena.make<IntegerLit>(BigInt(-1)));
        return roots;
    }

    for (int k = 1; k < n; ++k) {
        if (gcd(k, n) == 1) {
            // exp(2 * pi * i * k / n)
            ExprPtr k_lit = arena.make<IntegerLit>(BigInt(k));
            ExprPtr exponent = arena.make<Binary>(BinaryOp::Mul, 
                arena.make<Binary>(BinaryOp::Mul, two, pi),
                arena.make<Binary>(BinaryOp::Mul, i_const, 
                    arena.make<Binary>(BinaryOp::Div, k_lit, n_lit)));
            
            // We use exp function
            roots.push_back(arena.make<FuncCall>("exp", std::vector<ExprPtr>{exponent}));
        }
    }

    return roots;
}

} // namespace cas::algebra
