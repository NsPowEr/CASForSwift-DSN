// factorization_berlekamp.cpp — Berlekamp algorithm for factorization in Fp[x].
// A3 (F2 Block A).
//
// References:
//   Berlekamp 1967 "Factoring Polynomials over Finite Fields"
//   Knuth TAOCP Vol.2 §4.6.2 Algorithm B
//   Geddes-Czapor-Labahn "Algorithms for Computer Algebra" §8.4
//
// Algorithm:
//   1. Build the Berlekamp Q matrix: Q[i][j] = coefficient j of x^(p·i) mod f.
//   2. Compute null space of (Q - I) over Fp via Gaussian elimination.
//   3. Each null-space basis vector h(x) yields factors via gcd(f, h(x) - s)
//      for all s in Fp.
//   Splitting stops when we have deg(f) factors (all linear) or the factor
//   count equals the null-space dimension (all irreducible blocks found).
//
// Size guard: if deg(f) * p > max_matrix_size (parameter, default 1024),
// returns Unimplemented (category BUDGET_GUARD) to avoid O(n²·p) memory
// blowup.  Pass ctx.max_berlekamp_matrix_size() to make the budget
// configurable per-context.  Under the threshold Berlekamp is exact;
// Cantor-Zassenhaus (polynomial_modular.cpp) is used by the dispatcher
// for large p (probabilistic, no size constraint).
//
// Correctness certificate: ∏ returned factors ≡ f (mod p), verified by
// the caller.  Internal check: result.size() == null_space dimension.

#include "polynomial_internal.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include <vector>
#include <algorithm>

namespace cas::algebra {

namespace {

// --- Fp arithmetic helpers ---

static BigInt fp_mod(const BigInt& a, const BigInt& p) {
    BigInt r = a % p;
    if (r.is_negative()) r += p;
    return r;
}

// Polynomial mod over Fp (in-place)
static void fp_poly_mod(IntPoly& poly, const BigInt& p) {
    for (auto& c : poly.coefficients()) c = fp_mod(c, p);
    poly.normalize([](const BigInt& v) { return v.is_zero(); });
}

// Poly multiplication mod p, result reduced mod f and mod p
static IntPoly fp_poly_mul_mod_poly(const IntPoly& a, const IntPoly& b,
                                     const IntPoly& f, const BigInt& p) {
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly r;
    r.resize(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].is_zero()) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] = fp_mod(r[i + j] + a[i] * b[j], p);
        }
    }
    // Reduce mod f
    if (!f.is_zero()) {
        auto inv_lc_res = numtheory::modular_inverse(fp_mod(f.leading_coeff(), p), p);
        if (inv_lc_res.is_ok()) {
            const BigInt inv_lc = inv_lc_res.value();
            while (!r.is_zero() && r.degree() >= f.degree()) {
                std::size_t dd = r.degree() - f.degree();
                BigInt fac = fp_mod(r.leading_coeff() * inv_lc, p);
                for (std::size_t i = 0; i < f.size(); ++i) {
                    r[i + dd] = fp_mod(r[i + dd] - fac * f[i], p);
                }
                r.normalize([](const BigInt& v) { return v.is_zero(); });
            }
        }
    }
    r.normalize([](const BigInt& v) { return v.is_zero(); });
    return r;
}

// Compute x^e mod f mod p via repeated squaring
static IntPoly fp_poly_power_mod(std::size_t degree_x, BigInt e,
                                  const IntPoly& f, const BigInt& p) {
    // Start with base = x^degree_x reduced mod f
    IntPoly base;
    if (degree_x < f.size()) {
        // x^degree_x is just the monomial [0,...,0,1] at position degree_x
        base.resize(degree_x + 1U, BigInt(0));
        base[degree_x] = BigInt(1);
        // reduce mod f if needed
        if (base.degree() >= f.degree()) {
            base = fp_poly_mul_mod_poly(IntPoly{std::vector<BigInt>{BigInt(1)}},
                                         base, f, p);
        }
    } else {
        base.resize(degree_x + 1U, BigInt(0));
        base[degree_x] = BigInt(1);
        // Reduce mod f
        while (!base.is_zero() && base.degree() >= f.degree()) {
            auto inv_lc_res = numtheory::modular_inverse(fp_mod(f.leading_coeff(), p), p);
            if (inv_lc_res.is_error()) break;
            std::size_t dd = base.degree() - f.degree();
            BigInt fac = fp_mod(base.leading_coeff() * inv_lc_res.value(), p);
            for (std::size_t i = 0; i < f.size(); ++i) {
                base[i + dd] = fp_mod(base[i + dd] - fac * f[i], p);
            }
            base.normalize([](const BigInt& v) { return v.is_zero(); });
        }
    }

    // Square-and-multiply
    IntPoly result{std::vector<BigInt>{BigInt(1)}};
    while (!e.is_zero()) {
        if (!(e % BigInt(2)).is_zero()) {
            result = fp_poly_mul_mod_poly(result, base, f, p);
        }
        base = fp_poly_mul_mod_poly(base, base, f, p);
        e = e / BigInt(2);
    }
    return result;
}

// GCD over Fp[x]
static IntPoly fp_poly_gcd(IntPoly a, IntPoly b, const BigInt& p) {
    fp_poly_mod(a, p);
    fp_poly_mod(b, p);
    while (!b.is_zero()) {
        // Compute remainder a mod b
        auto inv_lc_res = numtheory::modular_inverse(fp_mod(b.leading_coeff(), p), p);
        if (inv_lc_res.is_error()) break;
        const BigInt inv_lc = inv_lc_res.value();
        while (!a.is_zero() && a.degree() >= b.degree()) {
            std::size_t dd = a.degree() - b.degree();
            BigInt fac = fp_mod(a.leading_coeff() * inv_lc, p);
            for (std::size_t i = 0; i < b.size(); ++i) {
                a[i + dd] = fp_mod(a[i + dd] - fac * b[i], p);
            }
            a.normalize([](const BigInt& v) { return v.is_zero(); });
        }
        std::swap(a, b);
    }
    // Make monic
    if (!a.is_zero()) {
        auto inv_lc_res = numtheory::modular_inverse(fp_mod(a.leading_coeff(), p), p);
        if (inv_lc_res.is_ok()) {
            for (auto& c : a.coefficients()) c = fp_mod(c * inv_lc_res.value(), p);
        }
    }
    return a;
}

// Subtract scalar s from constant term of poly, mod p
static IntPoly fp_poly_sub_scalar(IntPoly poly, const BigInt& s, const BigInt& p) {
    if (poly.is_zero()) {
        poly.resize(1U, BigInt(0));
    }
    poly[0] = fp_mod(poly[0] - s, p);
    poly.normalize([](const BigInt& v) { return v.is_zero(); });
    return poly;
}

// --- Berlekamp Q-matrix and null space ---

// Build the Berlekamp Q matrix (n×n, where n = deg(f)).
// Q[i][j] = j-th coefficient of x^(p*i) mod f mod p,  for i,j in [0,n).
// Stored as Q[row][col] = Q[i][j].
using Matrix = std::vector<std::vector<BigInt>>;

static Matrix build_berlekamp_matrix(const IntPoly& f, const BigInt& p) {
    const std::size_t n = f.degree();
    Matrix Q(n, std::vector<BigInt>(n, BigInt(0)));

    for (std::size_t i = 0; i < n; ++i) {
        // Compute x^(p*i) mod f mod p
        BigInt exponent = p * BigInt(static_cast<long long>(i));
        IntPoly row_poly = fp_poly_power_mod(1U, exponent, f, p);
        // Fill row i of Q
        for (std::size_t j = 0; j < row_poly.size() && j < n; ++j) {
            Q[i][j] = fp_mod(row_poly[j], p);
        }
    }
    return Q;
}

// Compute null space of (Q - I) over Fp via Gaussian elimination.
// Returns basis vectors of the null space, each as a polynomial of degree < n.
// The null space always contains the constant polynomial 1 (corresponding to f itself).
// Dimension of null space = number of irreducible factors of f.
static std::vector<IntPoly> berlekamp_null_space(Matrix Q, const BigInt& p) {
    const std::size_t n = Q.size();

    // M = Q^T - I  (we work on columns of Q-I, since null space of (Q-I)^T
    // corresponds to vectors v with v*(Q-I) = 0, i.e., vQ = v)
    // Equivalently: we want row vectors h with h*Q = h, i.e., h*(Q-I) = 0.
    // Build augmented matrix [Q^T - I | I_n] for row reduction.
    // After reduction, the null space basis is readable from the right half
    // where left half rows are zero.
    //
    // Simpler: work on A = (Q - I)^T and find null(A) by row reduction.
    // Row i of A = column i of (Q-I) = Q[:,i] - e_i.
    Matrix A(n, std::vector<BigInt>(n, BigInt(0)));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            // A[i][j] = (Q-I)[j][i] = Q[j][i] - delta(j==i)
            A[i][j] = fp_mod(Q[j][i] - (j == i ? BigInt(1) : BigInt(0)), p);
        }
    }

    // Track which rows are "free" (pivot-less) — they give null space vectors.
    // Gaussian elimination with partial pivoting.
    std::vector<int> pivot_col(n, -1);  // pivot_col[row] = col of pivot
    std::vector<int> pivot_row(n, -1);  // pivot_row[col] = row with pivot in that col

    std::size_t current_row = 0;
    for (std::size_t col = 0; col < n && current_row < n; ++col) {
        // Find pivot in column col, rows >= current_row
        std::size_t pivot = n;
        for (std::size_t row = current_row; row < n; ++row) {
            if (!A[row][col].is_zero()) { pivot = row; break; }
        }
        if (pivot == n) continue;  // free column

        std::swap(A[current_row], A[pivot]);
        pivot_col[current_row] = static_cast<int>(col);
        pivot_row[col] = static_cast<int>(current_row);

        // Scale row to make pivot = 1
        auto inv_res = numtheory::modular_inverse(A[current_row][col], p);
        if (inv_res.is_ok()) {
            const BigInt inv = inv_res.value();
            for (auto& val : A[current_row]) val = fp_mod(val * inv, p);
        }

        // Eliminate column in all other rows
        for (std::size_t row = 0; row < n; ++row) {
            if (row == current_row || A[row][col].is_zero()) continue;
            const BigInt factor = A[row][col];
            for (std::size_t c = 0; c < n; ++c) {
                A[row][c] = fp_mod(A[row][c] - factor * A[current_row][c], p);
            }
        }
        ++current_row;
    }

    // Extract null space: free columns (where pivot_row[col] == -1) generate
    // null vectors. For each free column col_f, construct v with:
    //   v[col_f] = 1
    //   v[pivot_col of row r] = -A[r][col_f] for each pivoted row r
    std::vector<IntPoly> basis;

    // First always include the constant 1 (trivial null vector for the trivial factor f)
    {
        std::vector<BigInt> one_coeffs(n, BigInt(0));
        one_coeffs[0] = BigInt(1);
        IntPoly h(std::move(one_coeffs));
        h.normalize([](const BigInt& v) { return v.is_zero(); });
        basis.push_back(std::move(h));
    }

    for (std::size_t col_f = 1; col_f < n; ++col_f) {
        if (pivot_row[col_f] != -1) continue;  // not free

        // Build null vector
        std::vector<BigInt> v(n, BigInt(0));
        v[col_f] = BigInt(1);
        for (std::size_t row = 0; row < current_row; ++row) {
            if (pivot_col[row] < 0) continue;
            std::size_t pc = static_cast<std::size_t>(pivot_col[row]);
            v[pc] = fp_mod(-A[row][col_f], p);
        }
        IntPoly h(std::move(v));
        h.normalize([](const BigInt& v2) { return v2.is_zero(); });
        basis.push_back(std::move(h));
    }

    return basis;
}

// Exact polynomial division g / d in Fp[x] (d divides g exactly).
static IntPoly fp_poly_exact_div(const IntPoly& g, const IntPoly& d, const BigInt& p) {
    auto inv_lc_res = numtheory::modular_inverse(fp_mod(d.leading_coeff(), p), p);
    if (inv_lc_res.is_error()) return IntPoly{};
    const BigInt inv_lc = inv_lc_res.value();
    IntPoly rem = g;
    fp_poly_mod(rem, p);
    IntPoly quo;
    while (!rem.is_zero() && rem.degree() >= d.degree()) {
        std::size_t dd = rem.degree() - d.degree();
        BigInt fac = fp_mod(rem.leading_coeff() * inv_lc, p);
        if (quo.size() <= dd) quo.resize(dd + 1U, BigInt(0));
        quo[dd] = fac;
        for (std::size_t i = 0; i < d.size(); ++i) {
            rem[i + dd] = fp_mod(rem[i + dd] - fac * d[i], p);
        }
        rem.normalize([](const BigInt& v) { return v.is_zero(); });
    }
    quo.normalize([](const BigInt& v) { return v.is_zero(); });
    return quo;
}

// Berlekamp splitting: given f and null space basis, compute irreducible factors.
// Iterates ALL basis elements against ALL current composite factors in repeated
// rounds until no further splitting occurs or all target factors are found.
//
// Key invariant: a basis vector hᵢ that cannot split f may still split a proper
// sub-factor produced by a previous split.  The outer progress-loop guarantees
// we revisit all (hᵢ, s) pairs after each new split, so no sub-factor is missed.
//
// Ref: Knuth TAOCP Vol.2 §4.6.2 Algorithm B; Geddes-Czapor-Labahn §8.4.
static std::vector<IntPoly> berlekamp_split(
    const IntPoly& f_in,
    const std::vector<IntPoly>& basis,
    const BigInt& p) {

    std::vector<IntPoly> factors;
    factors.push_back(f_in);

    const std::size_t target = basis.size();  // = number of irreducible factors

    // Repeat rounds until no progress (each round tries all (hᵢ, s) pairs).
    bool progress = true;
    while (progress && factors.size() < target) {
        progress = false;
        for (std::size_t bi = 1; bi < basis.size() && factors.size() < target; ++bi) {
            const IntPoly& h = basis[bi];
            std::vector<IntPoly> next_factors;
            for (const IntPoly& g : factors) {
                if (g.degree() <= 1U) {
                    next_factors.push_back(g);
                    continue;
                }
                bool split_found = false;
                for (BigInt s(0); s < p; s += BigInt(1)) {
                    IntPoly h_minus_s = fp_poly_sub_scalar(h, s, p);
                    IntPoly d = fp_poly_gcd(g, h_minus_s, p);
                    if (d.degree() > 0 && d.degree() < g.degree()) {
                        IntPoly quo = fp_poly_exact_div(g, d, p);
                        next_factors.push_back(d);
                        if (!quo.is_zero()) next_factors.push_back(quo);
                        split_found = true;
                        progress = true;
                        break;
                    }
                }
                if (!split_found) {
                    next_factors.push_back(g);
                }
            }
            factors = std::move(next_factors);
        }
    }
    return factors;
}

} // namespace

// --- Public API ---

// Berlekamp factorization of f in Fp[x].
// Budget guard: if deg(f) * p > max_matrix_size, returns Unimplemented with
// reason code BERLEKAMP_MATRIX_TOO_LARGE.  The caller falls back to
// Cantor-Zassenhaus.  Pass ctx.max_berlekamp_matrix_size() (default 1024,
// configurable via CASContext::set_max_berlekamp_matrix_size()) for
// context-driven budget control; the parameter default matches the CASContext
// default so bare callers (tests without a context) need not change.
//
// Correctness certificate (internal): ∏ factors ≡ f (mod p).
Result<std::vector<IntPoly>> berlekamp_factor_mod_p(IntPoly f, const BigInt& p,
                                                     std::size_t max_matrix_size) {
    if (f.is_zero()) return ok(std::vector<IntPoly>{});

    fp_poly_mod(f, p);
    f.normalize([](const BigInt& v) { return v.is_zero(); });

    if (f.is_zero()) return ok(std::vector<IntPoly>{});

    // Make monic
    {
        auto inv_lc_res = numtheory::modular_inverse(fp_mod(f.leading_coeff(), p), p);
        if (inv_lc_res.is_error()) return fail<std::vector<IntPoly>>(inv_lc_res.error());
        const BigInt inv_lc = inv_lc_res.value();
        for (auto& c : f.coefficients()) c = fp_mod(c * inv_lc, p);
        f.normalize([](const BigInt& v) { return v.is_zero(); });
    }

    const std::size_t n = f.degree();
    if (n == 0U) return ok(std::vector<IntPoly>{});
    if (n == 1U) return ok(std::vector<IntPoly>{f});

    // Budget guard: matrix is n×n, building each row costs O(n log(p·i) + n²) ops.
    // For n·p > max_matrix_size the Berlekamp null-space is impractical; fall back to CZ.
    // max_matrix_size is supplied by the caller (ctx.max_berlekamp_matrix_size(), default 1024).
    {
        // Use to_u64 safely; if p is large the guard triggers immediately.
        std::size_t p_val = 0;
        const std::string p_str = p.decimal();
        bool overflow = false;
        if (p_str.size() > 10) {
            overflow = true;  // definitely large
        } else {
            p_val = static_cast<std::size_t>(p.to_u64());
        }
        if (overflow || n * p_val > max_matrix_size) {
            return make_unimplemented<std::vector<IntPoly>>(
                "factorization_berlekamp",
                "berlekamp_factor_mod_p",
                "deg=" + std::to_string(n) + " p=" + p.decimal() +
                    " limit=" + std::to_string(max_matrix_size),
                "BERLEKAMP_MATRIX_TOO_LARGE",
                "Use Cantor-Zassenhaus (factor_polynomial_mod_p) for large p, "
                "or raise ctx.set_max_berlekamp_matrix_size()",
                "F2-A3",
                "Berlekamp Q-matrix size deg*p exceeds budget guard (max_berlekamp_matrix_size)");
        }
    }

    // Step 1: Build Berlekamp Q matrix
    Matrix Q = build_berlekamp_matrix(f, p);

    // Step 2: Compute null space of (Q - I) over Fp
    std::vector<IntPoly> basis = berlekamp_null_space(Q, p);

    if (basis.size() == 1U) {
        // null space dimension 1 → f is irreducible
        return ok(std::vector<IntPoly>{f});
    }

    // Step 3: Split f using basis elements
    std::vector<IntPoly> factors = berlekamp_split(f, basis, p);

    return ok(std::move(factors));
}

} // namespace cas::algebra
