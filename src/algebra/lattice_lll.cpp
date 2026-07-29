#include "cas/bigint.hpp"
#include "cas/rational.hpp"
#include "polynomial_internal.hpp"
#include <vector>
#include <algorithm>

namespace cas::algebra {

namespace {

// Rational Gram-Schmidt LLL.  O(n^2 m) recompute per swap with Rational
// coefficients — correct but slow, and used ONLY as a fallback for
// rank-deficient (linearly dependent) input bases, which the fraction-free
// integer path below cannot represent (zero Gram determinants).  All lattices
// produced by the factorisation pipeline (p^k·I block + modular-factor rows)
// are full rank, so this branch never runs on the performance-critical path.
Rational lll_dot(const LatticeVector& a, const LatticeVector& b) {
    Rational res(BigInt(0));
    std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) res += a[i] * b[i];
    return res;
}

void lll_reduction_rational(LatticeMatrix& b, double delta_val) {
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
                if (B[j].numerator().is_zero()) { mu[i][j] = Rational(BigInt(0)); continue; }
                mu[i][j] = lll_dot(b[i], b_star[j]) / B[j];
                for (std::size_t k_idx = 0U; k_idx < m; ++k_idx) {
                    b_star[i][k_idx] -= mu[i][j] * b_star[j][k_idx];
                }
            }
            B[i] = lll_dot(b_star[i], b_star[i]);
        }
    };

    recompute_gs();
    std::size_t k = 1;
    while (k < n) {
        for (int j = (int)k - 1; j >= 0; --j) {
            if (mu[k][j].numerator().abs() * BigInt(2) > mu[k][j].denominator().abs()) {
                BigInt q = mu[k][j].round();
                if (!q.is_zero()) {
                    Rational r_q(q);
                    for (std::size_t l = 0; l < m; ++l) b[k][l] -= r_q * b[j][l];
                    for (int l = 0; l < j; ++l) mu[k][l] -= r_q * mu[j][l];
                    mu[k][j] -= r_q;
                }
            }
        }
        if (B[k] >= (delta - mu[k][k - 1] * mu[k][k - 1]) * B[k - 1]) {
            k++;
        } else {
            std::swap(b[k], b[k - 1]);
            recompute_gs();
            k = std::max((std::size_t)1, k - 1);
        }
    }
}

}  // namespace

// LLL reduction for a basis b.
//
// Fraction-free integer LLL (Cohen, "A Course in Computational Algebraic Number
// Theory", Algorithm 2.6.7; de Weger 1987).  All Gram-Schmidt data is kept as
// integers: the Gram determinants d[i] = prod_{j<i} ||b*_j||^2 and the scaled
// coefficients lam[i][j] = d[j+1]*mu[i][j].  Every quantity is Hadamard-bounded,
// so there is NO Rational denominator growth.  The previous Rational
// Gram-Schmidt recomputed the full orthogonalisation per swap and its b_star
// fractions exploded on van Hoeij knapsack lattices, hanging deg-16 tower
// factorisations and the SD3 Swinnerton-Dyer gate
// (HC-F8-FACTORIZATIONTOWER-PERF / HC-F8-SD3-VANHOEIJ-SLOW).
//
// The input basis is integer-valued by construction (p^k * x^i and modular
// factor rows); the basis stays integral throughout (size reduction subtracts
// integer multiples) and is written back unchanged in type.
void lll_reduction(LatticeMatrix& b_io, double delta_val) {
    const std::size_t n = b_io.size();
    if (n == 0) return;
    const std::size_t m = b_io[0].size();

    // delta = dnum / dden  (Lovász parameter, typically 3/4 .. 0.99).
    const BigInt dden(1000000LL);
    const BigInt dnum(static_cast<long long>(delta_val * 1000000.0 + 0.5));

    // Integer working copy of the basis.
    std::vector<std::vector<BigInt>> b(n, std::vector<BigInt>(m, BigInt(0)));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t c = 0; c < m; ++c) {
            b[i][c] = b_io[i][c].numerator();
        }
    }

    auto dot = [&](const std::vector<BigInt>& u, const std::vector<BigInt>& v) {
        BigInt s(0);
        for (std::size_t c = 0; c < m; ++c) s = s + u[c] * v[c];
        return s;
    };

    std::vector<BigInt> d(n + 1, BigInt(0));
    d[0] = BigInt(1);
    std::vector<std::vector<BigInt>> lam(n, std::vector<BigInt>(n, BigInt(0)));

    // Initial fraction-free Gram-Schmidt (Cohen 2.6.3 integer variant).
    // A zero Gram determinant d[i] signals a linearly dependent (rank-deficient)
    // basis, which the integer recurrence cannot divide through — delegate the
    // whole reduction to the Rational fallback, which handles zero b* vectors.
    bool rank_deficient = false;
    for (std::size_t i = 0; i < n && !rank_deficient; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            BigInt u = dot(b[i], b[j]);
            for (std::size_t t = 0; t < j; ++t) {
                if (d[t].is_zero()) { rank_deficient = true; break; }
                u = (d[t + 1] * u - lam[i][t] * lam[j][t]) / d[t];
            }
            if (rank_deficient) break;
            if (j < i) lam[i][j] = u;
            else {
                d[i + 1] = u;
                if (u.is_zero()) { rank_deficient = true; break; }
            }
        }
    }
    if (rank_deficient) {
        lll_reduction_rational(b_io, delta_val);
        return;
    }

    // RED(kx, l): size-reduce b[kx] against b[l]  (l < kx).
    auto reduce = [&](std::size_t kx, std::size_t l) {
        if (BigInt(2) * lam[kx][l].abs() <= d[l + 1]) return;
        const BigInt q = Rational(lam[kx][l], d[l + 1]).round();
        if (q.is_zero()) return;
        for (std::size_t c = 0; c < m; ++c) b[kx][c] = b[kx][c] - q * b[l][c];
        lam[kx][l] = lam[kx][l] - q * d[l + 1];
        for (std::size_t i = 0; i < l; ++i) lam[kx][i] = lam[kx][i] - q * lam[l][i];
    };

    std::size_t k = 1;
    while (k < n) {
        reduce(k, k - 1);

        // Lovász: no swap iff  dden*(d[k+1]*d[k-1] + lam[k][k-1]^2) >= dnum*d[k]^2.
        const BigInt lhs = dden * (d[k + 1] * d[k - 1] + lam[k][k - 1] * lam[k][k - 1]);
        const BigInt rhs = dnum * (d[k] * d[k]);
        if (!(lhs < rhs)) {
            for (int l = static_cast<int>(k) - 2; l >= 0; --l) {
                reduce(k, static_cast<std::size_t>(l));
            }
            ++k;
        } else {
            // SWAP(k): exchange b[k], b[k-1] and update d/lam in integers.
            std::swap(b[k], b[k - 1]);
            for (std::size_t j = 0; j + 1 < k; ++j) std::swap(lam[k][j], lam[k - 1][j]);
            const BigInt lambda = lam[k][k - 1];
            const BigInt B_new = (d[k - 1] * d[k + 1] + lambda * lambda) / d[k];
            for (std::size_t i = k + 1; i < n; ++i) {
                const BigInt t = lam[i][k];
                lam[i][k] = (d[k + 1] * lam[i][k - 1] - lambda * t) / d[k];
                lam[i][k - 1] = (B_new * t + lambda * lam[i][k]) / d[k + 1];
            }
            d[k] = B_new;
            k = (k > 1) ? (k - 1) : 1;
        }
    }

    // Write the reduced integer basis back into the Rational matrix.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t c = 0; c < m; ++c) {
            b_io[i][c] = Rational(b[i][c]);
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
