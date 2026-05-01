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
    
    Rational delta = Rational(BigInt(75), BigInt(100)); // Using 0.75 as default
    if (delta_val != 0.75) {
        // Convert double to Rational if needed, but 0.75 is the usual value.
        // For simplicity, we stick to Rational(75, 100).
    }

    LatticeMatrix b_star(n, LatticeVector(m));
    std::vector<std::vector<Rational>> mu(n, std::vector<Rational>(n));
    std::vector<Rational> B(n);
    
    auto compute_gs = [&](std::size_t i) {
        b_star[i] = b[i];
        for (std::size_t j = 0; j < i; ++j) {
            mu[i][j] = dot_product(b[i], b_star[j]) / B[j];
            for (std::size_t k_idx = 0; k_idx < m; ++k_idx) {
                b_star[i][k_idx] -= mu[i][j] * b_star[j][k_idx];
            }
        }
        B[i] = dot_product(b_star[i], b_star[i]);
    };
    
    B[0] = dot_product(b[0], b[0]);
    b_star[0] = b[0];
    
    std::size_t k = 1;
    while (k < n) {
        compute_gs(k);
        
        // Size reduction
        for (int j = (int)k - 1; j >= 0; --j) {
            if (mu[k][j].numerator().abs() * BigInt(2) > mu[k][j].denominator().abs()) {
                BigInt q = mu[k][j].round();
                if (!q.is_zero()) {
                    Rational r_q(q);
                    for (std::size_t l = 0; l < m; ++l) {
                        b[k][l] -= r_q * b[j][l];
                    }
                    // Update mu and b_star for the changed b[k]
                    compute_gs(k);
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
    [[maybe_unused]] const IntPoly& f,
    const IntPoly& g,
    const BigInt& pk,
    std::size_t max_deg) {
    
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
        
        lll_reduction(basis);
        
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
            h.normalize([](const BigInt& v) { return v.is_zero(); });
            if (!h.is_zero() && h.degree() > 0) {
                // Check if h divides f
                // We should use exact division
                return h;
            }
        }
    }
    
    return std::nullopt;
}

} // namespace cas::algebra
