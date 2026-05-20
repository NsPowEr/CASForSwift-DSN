#include "cas/bigint.hpp"
#include "cas/rational.hpp"
#include "polynomial_internal.hpp"
#include <vector>
#include <algorithm>

namespace cas::algebra {

static Rational dot_product(const LatticeVector& a, const LatticeVector& b) {
    Rational res(BigInt(0));
    std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) {
        res += a[i] * b[i];
    }
    return res;
}

// LLL reduction for a basis b
void lll_reduction(LatticeMatrix& b, double delta_val) {
    std::size_t n = b.size();
    if (n == 0) return;
    std::size_t m = b[0].size();
    
    const long long kDeltaDen = 1000000LL;
    const long long kDeltaNum = static_cast<long long>(delta_val * static_cast<double>(kDeltaDen) + 0.5);
    Rational delta = Rational(BigInt(kDeltaNum), BigInt(kDeltaDen));

    LatticeMatrix b_star(n, LatticeVector(m));
    std::vector<std::vector<Rational>> mu(n, std::vector<Rational>(n));
    std::vector<Rational> B(n);
    
    auto recompute_gs = [&]() {
        for (std::size_t i = 0U; i < n; ++i) {
            b_star[i] = b[i];
            for (std::size_t j = 0U; j < i; ++j) {
                if (B[j].numerator().is_zero()) {
                    mu[i][j] = Rational(BigInt(0));
                    continue;
                }
                mu[i][j] = dot_product(b[i], b_star[j]) / B[j];
                for (std::size_t k_idx = 0U; k_idx < m; ++k_idx) {
                    b_star[i][k_idx] -= mu[i][j] * b_star[j][k_idx];
                }
            }
            B[i] = dot_product(b_star[i], b_star[i]);
        }
    };
    
    std::size_t k = 1;
    while (k < n) {
        recompute_gs();
        
        // Size reduction
        for (int j = (int)k - 1; j >= 0; --j) {
            if (mu[k][j].numerator().abs() * BigInt(2) > mu[k][j].denominator().abs()) {
                BigInt q = mu[k][j].round();
                if (!q.is_zero()) {
                    Rational r_q(q);
                    for (std::size_t l = 0; l < m; ++l) {
                        b[k][l] -= r_q * b[j][l];
                    }
                    recompute_gs();
                }
            }
        }
        
        // Lovász condition
        if (B[k] >= (delta - mu[k][k - 1] * mu[k][k - 1]) * B[k - 1]) {
            k++;
        } else {
            std::swap(b[k], b[k - 1]);
            k = std::max((std::size_t)1, k - 1);
        }
    }
}

// Identify potential factor of f in Z[x] given modular factor g mod p^k
std::optional<IntPoly> find_factor_lll(
    const IntPoly& f,
    const IntPoly& g,
    const BigInt& pk,
    std::size_t max_deg,
    double delta) {
    
    std::size_t d = g.degree();
    if (max_deg < d) max_deg = d;
    
    // Lattice basis for degree m
    for (std::size_t m = d; m <= max_deg; ++m) {
        LatticeMatrix basis(m + 1, LatticeVector(m + 1, Rational(BigInt(0))));
        
        // p^k * x^i for 0 <= i < d
        for (std::size_t i = 0; i < d; ++i) {
            basis[i][i] = Rational(pk);
        }
        
        // g(x) * x^{i-d} for d <= i <= m
        for (std::size_t i = d; i <= m; ++i) {
            for (std::size_t j = 0; j < g.size(); ++j) {
                if (i - d + j <= m) {
                    basis[i][i - d + j] = Rational(g[j]);
                }
            }
        }
        
        lll_reduction(basis, delta);
        
        // Shortest vector is in basis[0]
        IntPoly h;
        h.resize(m + 1);
        bool all_integers = true;
        for (std::size_t i = 0; i <= m; ++i) {
            if (!basis[0][i].is_integer()) {
                all_integers = false;
                break;
            }
            h[i] = basis[0][i].numerator();
        }
        
        if (all_integers) {
            h = primitive_integer_poly(std::move(h));
            if (!h.is_zero() && h.degree() > 0 && h.degree() < f.degree()) {
                auto remainder = pseudo_remainder_integer_poly(f, h);
                normalize_integer_poly(remainder);
                if (remainder.is_zero()) {
                    return h;
                }
            }
        }
    }
    
    return std::nullopt;
}

} // namespace cas::algebra
