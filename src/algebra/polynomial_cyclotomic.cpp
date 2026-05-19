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

// F2.1.c: hard cap on recursion depth and search range. The recursion
// descends along divisors of n, so depth ≤ log2(n). For n ≤ 2^20 the
// recursion is safe; beyond that the inner divisor-enumeration loop
// `for (int d = 1; d < n; ++d)` and the cache map become the dominant
// cost (O(n) time, O(n) memory). We refuse arguments above kMaxN to
// prevent silent OOM / runaway.
constexpr int kCyclotomicMaxN = 1 << 20;  // 1048576

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
IntPoly compute_cyclotomic(int n) {
    if (n <= 0 || n > kCyclotomicMaxN) {
        // Caller (`is_cyclotomic`) checks empty and skips.
        return IntPoly{};
    }
    if (n == 1) return IntPoly({BigInt(-1), BigInt(1)});  // x - 1

    IntPoly numerator({BigInt(1)});    // (x^d - 1)^(+1) accumulator
    IntPoly denominator({BigInt(1)});  // (x^d - 1)^(-1) accumulator

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

bool poly_equal(const IntPoly& a, const IntPoly& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

} // namespace

std::optional<int> is_cyclotomic(const IntPoly& poly, int max_n) {
    if (poly.empty()) return std::nullopt;
    if (poly.leading_coeff() != BigInt(1)) return std::nullopt;

    std::size_t deg = poly.degree();
    if (deg == 0) return std::nullopt;

    // If max_n not specified, derive from degree.
    // For φ(n) = d: n = p (prime, p-1 = d), n = 2p, or composites up to ~2*(d+1).
    // We use 2*(deg+1) as the base and add a small constant for composites like n=4 (d=2).
    if (max_n < 0) {
        max_n = static_cast<int>(std::max<std::size_t>(12U, 2U * (deg + 1U)));
    }
    // Cap at kCyclotomicMaxN to prevent runaway scan; for deg > ~500
    // a Möbius-formula direct construction is faster than enumeration.
    // Tracked as FE-004 for follow-up.
    if (max_n > kCyclotomicMaxN) max_n = kCyclotomicMaxN;

    for (int n = 1; n <= max_n; ++n) {
        IntPoly phi = compute_cyclotomic(n);
        if (phi.empty()) continue;  // skip cap'd-out entries
        if (phi.degree() == deg) {
            if (poly_equal(poly, phi)) return n;
        }
    }

    return std::nullopt;
}

std::vector<ExprPtr> cyclotomic_roots(int n, const Symbol& var, AstArena& arena) {
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
