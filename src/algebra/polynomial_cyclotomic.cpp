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

IntPoly compute_cyclotomic(int n) {
    // Phi_n(x) = product_{d|n} (x^d - 1)^mu(n/d)
    // We can also use recursive division: x^n - 1 = product_{d|n} Phi_d(x)
    // So Phi_n(x) = (x^n - 1) / product_{d|n, d<n} Phi_d(x)
    
    static std::map<int, IntPoly> cache;
    if (cache.count(n)) return cache[n];

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

std::optional<int> is_cyclotomic(const IntPoly& poly) {
    if (poly.empty()) return std::nullopt;
    if (poly.leading_coeff() != BigInt(1)) return std::nullopt;

    std::size_t deg = poly.degree();
    if (deg == 0) return std::nullopt;

    // Candidates n for deg = phi(n)
    // We only check up to n=100 or something reasonable.
    // Small phi(n) lookup:
    // n=1:1, 2:1, 3:2, 4:2, 5:4, 6:2, 7:6, 8:4, 9:6, 10:4, 12:4, ...
    for (int n = 1; n <= 100; ++n) {
        IntPoly phi = compute_cyclotomic(n);
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
