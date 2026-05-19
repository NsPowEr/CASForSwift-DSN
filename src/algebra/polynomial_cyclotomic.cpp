#include "polynomial_internal.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include <map>
#include <cmath>

namespace cas::algebra {

namespace {

int mobius(int n) {
    if (n == 1) return 1;
    int p = 0;
    for (int i = 2; i <= n; i++) {
        if (n % i == 0) {
            p++;
            if ((n / i) % i == 0) return 0;
            n /= i;
        }
    }
    return (p % 2 == 0) ? 1 : -1;
}

// F2.1.c: hard cap on recursion depth and search range. The recursion
// descends along divisors of n, so depth ≤ log2(n). For n ≤ 2^20 the
// recursion is safe; beyond that the inner divisor-enumeration loop
// `for (int d = 1; d < n; ++d)` and the cache map become the dominant
// cost (O(n) time, O(n) memory). We refuse arguments above kMaxN to
// prevent silent OOM / runaway.
constexpr int kCyclotomicMaxN = 1 << 20;  // 1048576

IntPoly compute_cyclotomic(int n) {
    // Phi_n(x) = product_{d|n} (x^d - 1)^mu(n/d)
    // We can also use recursive division: x^n - 1 = product_{d|n} Phi_d(x)
    // So Phi_n(x) = (x^n - 1) / product_{d|n, d<n} Phi_d(x)

    static std::map<int, IntPoly> cache;
    if (cache.count(n)) return cache[n];
    if (n <= 0 || n > kCyclotomicMaxN) {
        // Caller must guard against this; we return empty to signal
        // failure rather than crash. is_cyclotomic checks this below.
        return IntPoly{};
    }

    // Simple cases
    if (n == 1) return cache[1] = IntPoly({BigInt(-1), BigInt(1)}); // x - 1
    
    // x^n - 1
    std::vector<BigInt> coeffs(n + 1, BigInt(0));
    coeffs[0] = BigInt(-1);
    coeffs[n] = BigInt(1);
    IntPoly res(std::move(coeffs));

    for (int d = 1; d < n; ++d) {
        if (n % d == 0) {
            // res = res / compute_cyclotomic(d)
            IntPoly phid = compute_cyclotomic(d);
            
            // Polynomial division (exact)
            std::vector<BigInt> q;
            std::vector<BigInt> r = res.coefficients();
            std::size_t deg_r = r.size() - 1;
            std::size_t deg_d = phid.degree();
            
            std::vector<BigInt> quot(deg_r - deg_d + 1, BigInt(0));
            for (int i = static_cast<int>(deg_r - deg_d); i >= 0; --i) {
                quot[i] = r[i + deg_d] / phid.leading_coeff();
                for (std::size_t j = 0; j <= deg_d; ++j) {
                    r[i + j] -= quot[i] * phid[j];
                }
            }
            res = IntPoly(std::move(quot));
            res.normalize([](const BigInt& b) { return b.is_zero(); });
        }
    }

    return cache[n] = res;
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
