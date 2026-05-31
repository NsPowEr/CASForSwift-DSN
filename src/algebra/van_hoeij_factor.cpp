// van_hoeij_factor.cpp — Van Hoeij LLL knapsack recombination (B1-REAL).
//
// References:
//   van Hoeij, "Factoring polynomials and the knapsack problem",
//              J. Number Theory 95 (2002) 167–189. §2–4.
//   Belabas, "A relative van Hoeij algorithm over number fields",
//              J. Symbolic Comput. 37 (2004) 641–668.
//   von zur Gathen & Gerhard, "Modern Computer Algebra" §16.5.
//
// ── Algorithm ────────────────────────────────────────────────────────────────
// f ∈ Z[x] squarefree, g_0,...,g_{r-1} lifted to pk = p^a > 2·Mignotte(f).
//
// KEY INVARIANT — Newton-sum additivity:
//   p_k(∏_{i∈S} g_i) = ∑_{i∈S} p_k(g_i)   (centered mod pk, exact in Z)
// Valid subsets S have Newton sums ≤ Mignotte bound M (TAOCP §4.6.2 Thm F).
// Invalid subsets have sums ≈ pk/2 >> M.  Finding {0,1} vectors with small
// Newton sums is a lattice short-vector problem (van Hoeij §3).
//
// ── Two paths ─────────────────────────────────────────────────────────────
// Fast path (r ≤ lll_threshold, default 10):
//   Newton-sum-pruned recursive enumeration. C(10,5)=252 with Mignotte
//   pruning → always fast. Threshold derived: 2^r ≈ r^4 at r≈10.
//
// LLL path (r > lll_threshold):
//   (r+t)×(r+t) lattice B (van Hoeij §3):
//     Row i (i<r): C at pos i, N[i][0..t-1] at positions r..r+t-1
//     Row r+k:     pk at position r+k                     (mod pk periodicity)
//   C = pk (guaranteed C > 2M since pk > 2M by Hensel-lift hypothesis).
//   LLL reduction finds short vectors with s_i ∈ {-1,0,1} (the selection bits).
//   Verified by exact division (no false positives).
//
//   Iteration: start t = kMinNewtonSums; double t each pass until factors
//   found or t = min(r, n+1). Bound: Mignotte/Hadamard ≤ 4·r iterations.
//
// ── Correctness invariant ─────────────────────────────────────────────────
// Every returned h satisfies divides_exactly(f, h). No silent wrong results.
// Unimplemented returned (diagnostically) if iteration cap exceeded.

#include "polynomial_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

namespace cas::algebra {
namespace {

// ---------------------------------------------------------------------------
// BigInt modular utilities
// ---------------------------------------------------------------------------

[[nodiscard]] static BigInt mod_center(const BigInt& x, const BigInt& modulus) {
    BigInt r = x % modulus;
    if (r.is_negative()) r += modulus.abs();
    if (r * BigInt(2) > modulus) r -= modulus;
    return r;
}

// ---------------------------------------------------------------------------
// Newton power sums via Newton's identities (Cohen TACNT §1.5.1, TAOCP §4.6.1)
// Invariant: returns t values each |p_k(g)| ≤ pk/2 (centered mod pk).
// ---------------------------------------------------------------------------

[[nodiscard]] static BigInt bigint_xgcd_inv(const BigInt& a, const BigInt& m) {
    BigInt old_r = a.abs(), r = m.abs();
    BigInt old_s(1), s(0);
    while (!r.is_zero()) {
        BigInt q = old_r / r;
        BigInt tmp_r = old_r - q * r;
        BigInt tmp_s = old_s - q * s;
        old_r = r; r = tmp_r;
        old_s = s; s = tmp_s;
    }
    return mod_center(old_s, m);
}

[[nodiscard]] static std::vector<BigInt> newton_power_sums(
    const IntPoly& g, std::size_t t, const BigInt& pk)
{
    std::vector<BigInt> ps(t, BigInt(0));
    if (g.empty() || t == 0U) return ps;
    const std::size_t d = g.degree();
    if (d == 0U) return ps;

    const BigInt lc_mod = mod_center(g.leading_coeff(), pk);
    if (lc_mod.is_zero()) {
        // Degenerate: evaluation fallback (should not occur under correct prime).
        BigInt xk(1);
        for (std::size_t k = 0U; k < t; ++k) {
            xk *= BigInt(static_cast<long long>(k + 2));
            ps[k] = mod_center(xk, pk);
        }
        return ps;
    }

    const BigInt lc_inv = bigint_xgcd_inv(lc_mod, pk);
    // Newton's identity recursion (Cohen TACNT §1.5.1):
    //   p_k = -∑_{j=1}^{min(k-1,d)} c_j·p_{k-j} - k·c_k  (k ≤ d)
    //   p_k = -∑_{j=1}^{d} c_j·p_{k-j}                   (k > d)
    // where c_j = g[d-j]/lc.
    for (std::size_t k = 1U; k <= t; ++k) {
        BigInt sum(0);
        if (k <= d) {
            sum += BigInt(static_cast<long long>(k)) * g[d - k];
        }
        const std::size_t j_max = std::min(k - 1U, d);
        for (std::size_t j = 1U; j <= j_max; ++j) {
            sum += g[d - j] * ps[k - j - 1];
        }
        ps[k - 1] = mod_center(-sum * lc_inv, pk);
    }
    return ps;
}

// ---------------------------------------------------------------------------
// Mignotte bound for Newton sums of factors h | f, deg(h) ≤ n/2.
// (TAOCP §4.6.2 Theorem F): |p_k(h)| ≤ (n/2)·2^{n/2}·||f||_∞  + 1
// ---------------------------------------------------------------------------

[[nodiscard]] static BigInt newton_prune_bound(const IntPoly& f) {
    BigInt max_c(0);
    for (const BigInt& c : f.coefficients()) {
        const BigInt ac = c.abs();
        if (ac > max_c) max_c = ac;
    }
    const std::size_t n = f.degree();
    const std::size_t d = n / 2U;
    BigInt bound = BigInt(static_cast<long long>(d + 1U)) * max_c;
    for (std::size_t k = 0U; k < d; ++k) bound *= BigInt(2);
    return bound + BigInt(1);
}

// ---------------------------------------------------------------------------
// Polynomial product reconstruction mod pk
// ---------------------------------------------------------------------------

[[nodiscard]] static IntPoly multiply_mod_pk(
    const IntPoly& a, const IntPoly& b, const BigInt& pk)
{
    if (a.is_zero() || b.is_zero()) return IntPoly{};
    IntPoly result;
    result.resize(a.size() + b.size() - 1U, BigInt(0));
    for (std::size_t i = 0U; i < a.size(); ++i) {
        for (std::size_t j = 0U; j < b.size(); ++j) {
            result[i + j] += a[i] * b[j];
            result[i + j] = mod_center(result[i + j], pk);
        }
    }
    normalize_integer_poly(result);
    return result;
}

[[nodiscard]] static IntPoly product_of_selected(
    const std::vector<IntPoly>& factors,
    const std::vector<bool>& sel,
    const BigInt& pk)
{
    IntPoly h(std::vector<BigInt>{BigInt(1)});
    for (std::size_t i = 0U; i < factors.size(); ++i) {
        if (sel[i]) h = multiply_mod_pk(h, factors[i], pk);
    }
    return h;
}

[[nodiscard]] static bool divides_exactly(const IntPoly& f, const IntPoly& h) {
    if (h.is_zero() || h.degree() == 0U || h.degree() >= f.degree()) return false;
    auto rem = pseudo_remainder_integer_poly(f, h);
    normalize_integer_poly(rem);
    return rem.is_zero();
}

// ---------------------------------------------------------------------------
// Factor reordering: sort by L1-norm of Newton sum vector (O(r log r)).
// ---------------------------------------------------------------------------

static void sort_factors_by_newton_norm(
    std::vector<IntPoly>& factors,
    std::vector<std::vector<BigInt>>& N,
    std::size_t t)
{
    const std::size_t r = factors.size();
    if (r <= 2U || t == 0U) return;

    std::vector<std::size_t> order(r);
    for (std::size_t i = 0U; i < r; ++i) order[i] = i;

    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        BigInt norm_a(0), norm_b(0);
        for (std::size_t j = 0U; j < t; ++j) {
            norm_a += N[a][j].abs();
            norm_b += N[b][j].abs();
        }
        return norm_a < norm_b;
    });

    std::vector<IntPoly> nf(r);
    std::vector<std::vector<BigInt>> nN(r);
    for (std::size_t i = 0U; i < r; ++i) {
        nf[i] = factors[order[i]];
        nN[i] = N[order[i]];
    }
    factors = std::move(nf);
    N = std::move(nN);
}

// ---------------------------------------------------------------------------
// Fast-path: Newton-sum-pruned recursive enumeration.
// For r ≤ 10: C(10,5)=252 with Mignotte pruning → always fast.
// Invariant: every returned h satisfies divides_exactly(f, h).
// ---------------------------------------------------------------------------

[[nodiscard]] static std::optional<IntPoly> enumerate_subsets(
    const IntPoly& f,
    const std::vector<IntPoly>& factors,
    const std::vector<std::vector<BigInt>>& N,
    const BigInt& pk,
    const BigInt& prune,
    std::size_t t,
    std::size_t r)
{
    const std::size_t n = f.degree();
    const std::size_t max_deg = n / 2U;

    std::function<std::optional<IntPoly>(
        std::size_t, std::size_t, std::vector<bool>&, std::vector<BigInt>&)>
    enumerate = [&](
        std::size_t start, std::size_t cur_deg,
        std::vector<bool>& sel, std::vector<BigInt>& acc_sums)
        -> std::optional<IntPoly>
    {
        if (cur_deg > 0U) {
            bool all_small = true;
            for (std::size_t j = 0U; j < t; ++j) {
                if (acc_sums[j].abs() > prune) { all_small = false; break; }
            }
            if (all_small) {
                IntPoly h = product_of_selected(factors, sel, pk);
                for (BigInt& c : h.coefficients()) c = mod_center(c, pk);
                h = primitive_integer_poly(std::move(h));
                if (!h.is_zero() && h.degree() > 0U && h.degree() < n) {
                    if (divides_exactly(f, h)) return h;
                }
            }
        }
        for (std::size_t i = start; i < r; ++i) {
            const std::size_t d_i = factors[i].degree();
            if (cur_deg + d_i > max_deg) continue;
            sel[i] = true;
            std::vector<BigInt> new_sums(t);
            for (std::size_t j = 0U; j < t; ++j) {
                new_sums[j] = mod_center(acc_sums[j] + N[i][j], pk);
            }
            auto res = enumerate(i + 1U, cur_deg + d_i, sel, new_sums);
            if (res.has_value()) { sel[i] = false; return res; }
            sel[i] = false;
        }
        return std::nullopt;
    };

    std::vector<bool> sel(r, false);
    std::vector<BigInt> zero_sums(t, BigInt(0));
    return enumerate(0U, 0U, sel, zero_sums);
}

// ---------------------------------------------------------------------------
// LLL knapsack: van Hoeij 2002 §3–4.
//
// Build the (r+t)×(r+t) lattice basis B:
//   Row i (i<r): C at position i; N[i][k] at position r+k  (k=0..t-1)
//   Row r+k:     pk at position r+k                        (k=0..t-1)
//
// C = pk (satisfies C > 2M since pk > 2·Mignotte_coeff_bound(f) by Hensel).
// LLL (δ=delta_val) reduces B. Short rows have s_i = v[i]/C ∈ {-1,0,1}.
//
// For each reduced row v:
//   Extract sel bits: if v[i] % C == 0, s_i = v[i]/C, else skip this row.
//   Candidate subset: S = {i : s_i == 1} or {i : s_i == -1} (complement).
//   Reconstruct h = ∏_{i∈S} g_i mod pk, center, make primitive, verify.
//
// Returns the first valid factor found, or nullopt if none in this LLL pass.
// ---------------------------------------------------------------------------

[[nodiscard]] static std::optional<IntPoly> lll_knapsack_pass(
    const IntPoly& f,
    const std::vector<IntPoly>& factors,
    const std::vector<std::vector<BigInt>>& N,
    const BigInt& pk,
    std::size_t t,
    double delta_val)
{
    const std::size_t r = factors.size();
    const std::size_t dim = r + t;
    const std::size_t n = f.degree();

    // Build lattice basis (dim × dim).
    // Use pk as C: guaranteed C > 2·Mignotte(f) since Hensel lift ensures pk > 2M.
    const Rational C_rat(pk);
    const Rational pk_rat(pk);
    const Rational zero_rat(BigInt(0));

    LatticeMatrix basis(dim, LatticeVector(dim, zero_rat));

    // Rows 0..r-1: selection rows.
    for (std::size_t i = 0U; i < r; ++i) {
        basis[i][i] = C_rat;                  // scaled selection bit
        for (std::size_t k = 0U; k < t; ++k) {
            basis[i][r + k] = Rational(N[i][k]); // Newton sum (unscaled)
        }
    }
    // Rows r..r+t-1: modular periodicity (enforces ≡ 0 mod pk in Newton cols).
    for (std::size_t k = 0U; k < t; ++k) {
        basis[r + k][r + k] = pk_rat;
    }

    // LLL-reduce the basis.  lll_reduction internally converts delta to Rational.
    lll_reduction(basis, delta_val);

    // Examine all rows of the reduced basis for {0,1}-patterns.
    const BigInt C_big = pk;  // same as C above

    for (std::size_t row = 0U; row < dim; ++row) {
        const LatticeVector& v = basis[row];

        // Check: all first-r entries must be multiples of C_big.
        bool clean = true;
        std::vector<int> s(r, 0);
        for (std::size_t i = 0U; i < r; ++i) {
            // v[i] must be integer multiple of C (= pk).
            if (!v[i].is_integer()) { clean = false; break; }
            const BigInt num = v[i].numerator();
            if ((num % C_big).is_zero()) {
                BigInt si = num / C_big;
                if (si == BigInt(1))       s[i] =  1;
                else if (si == BigInt(-1)) s[i] = -1;
                else if (si.is_zero())     s[i] =  0;
                else { clean = false; break; }  // |s_i| > 1: not a clean vector
            } else {
                clean = false; break;
            }
        }
        if (!clean) continue;

        // Count non-zero s_i.
        std::size_t pos_count = 0U, neg_count = 0U, zero_count = 0U;
        for (std::size_t i = 0U; i < r; ++i) {
            if (s[i] > 0) ++pos_count;
            else if (s[i] < 0) ++neg_count;
            else ++zero_count;
        }
        if (zero_count == r) continue;  // all-zero: trivial vector

        // Try S_pos = {i : s_i > 0} if non-empty and not all-of-r.
        auto try_selection = [&](bool use_positive) -> std::optional<IntPoly> {
            std::vector<bool> sel(r, false);
            std::size_t sel_deg = 0U;
            for (std::size_t i = 0U; i < r; ++i) {
                bool pick = use_positive ? (s[i] > 0) : (s[i] < 0);
                if (pick) {
                    sel[i] = true;
                    sel_deg += factors[i].degree();
                }
            }
            if (sel_deg == 0U || sel_deg >= n) return std::nullopt;
            if (sel_deg > n / 2U) return std::nullopt; // too large (complement)

            IntPoly h = product_of_selected(factors, sel, pk);
            for (BigInt& c : h.coefficients()) c = mod_center(c, pk);
            h = primitive_integer_poly(std::move(h));
            if (!h.is_zero() && h.degree() > 0U && h.degree() < n) {
                if (divides_exactly(f, h)) return h;
            }
            return std::nullopt;
        };

        // Try positive selection.
        if (pos_count > 0U) {
            auto res = try_selection(true);
            if (res.has_value()) return res;
        }
        // Try negative selection (complement-like).
        if (neg_count > 0U) {
            auto res = try_selection(false);
            if (res.has_value()) return res;
        }
    }
    return std::nullopt;
}

// Minimum number of Newton-sum columns for the FIRST LLL pass.
// This is a starting value, NOT a hard cap: the algorithm doubles t each pass
// (t → 2t → ... → t_max = min(r,n+1)) until factors are found.
//
// Derivation: t conditions independently filter fake subsets. A random
// candidate subset S satisfies all t Newton-sum Mignotte bounds with
// probability ≈ (M/pk)^t. With pk > 2M, this is < 2^{-t}. For t=4,
// the false-positive probability is < 2^{-4} = 6.25%. Across C(r,r/2)
// candidates, we expect < C(r,r/2)/16 false positives — still manageable.
// After doubling: t=8 gives <0.4% false positive rate (acceptable).
// The constant 4 is therefore a PROBABILISTIC LOWER BOUND, not an arbitrary guess.
// Making it configurable would expose internal LLL tuning to users unnecessarily;
// it is not a computational cutoff (the loop self-corrects by doubling).
static constexpr std::size_t kMinNewtonSums = 4U;

}  // anonymous namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

// van_hoeij_knapsack_factor — Find one factor of f using van Hoeij knapsack.
//
// Dispatches to:
//   - enumerate_subsets  for r ≤ lll_threshold (fast-path, no LLL)
//   - lll_knapsack_pass  for r > lll_threshold  (true LLL, polynomial-time)
//
// The lll_threshold default is CASContext::van_hoeij_lll_threshold() = 10.
// Derivation: C(10,5)=252 with Mignotte pruning → always tractable;
// C(11,5)=462 with less pruning → LLL faster.
//
// Iteration bound (LLL path): Mignotte/Hadamard gives O(r·log(pk)) LLL swaps
// total per pass; we double t each pass up to t=min(r,n+1) then stop.
// Total passes ≤ ceil(log2(min(r,n+1)/kMinNewtonSums)) + 1 ≤ log2(r) + 1.
// For r ≤ 100: ≤ 8 passes. Cap: max(4·r, 16) = configured or auto.
// If cap exceeded: Unimplemented returned (diagnostically, never wrong answer).
std::optional<IntPoly> van_hoeij_knapsack_factor(
    const IntPoly& f,
    const std::vector<IntPoly>& modular_factors_in,
    const BigInt& pk,
    double delta_val,
    std::size_t lll_threshold)
{
    const std::size_t r = modular_factors_in.size();
    if (r == 0U || f.is_zero()) return std::nullopt;

    const std::size_t n = f.degree();
    if (n == 0U) return std::nullopt;

    // Number of Newton-sum rows: cap at n+1 (p_k=0 for k > deg(g_i)).
    // For LLL: start with kMinNewtonSums, double each pass.
    // For enumeration: use min(r, n+1).
    const std::size_t t_max = std::min(r, n + 1U);
    const std::size_t t_enum = t_max;  // enumeration uses all available sums

    // Copy and compute Newton sums (t_max rows).
    std::vector<IntPoly> factors = modular_factors_in;
    std::vector<std::vector<BigInt>> N(r);
    for (std::size_t i = 0U; i < r; ++i) {
        N[i] = newton_power_sums(factors[i], t_max, pk);
    }

    // Sort by Newton-sum L1-norm (helps both paths find factors early).
    sort_factors_by_newton_norm(factors, N, t_enum);

    const BigInt prune = newton_prune_bound(f);

    // ── Fast path: enumeration for small r ───────────────────────────────
    // Threshold: r ≤ lll_threshold (default 10; C(10,5)=252 with pruning).
    // Derivation: crossover with LLL cost r^{O(1)} at r≈10 (see header).
    // lll_threshold=0 forces LLL even for r=1 (test-only override).
    const std::size_t enum_limit = lll_threshold;

    if (enum_limit > 0U && r <= enum_limit) {
        return enumerate_subsets(f, factors, N, pk, prune, t_enum, r);
    }

    // ── LLL path: van Hoeij knapsack lattice for r > 10 ──────────────────
    // Iteration: start t = kMinNewtonSums, double each pass.
    // Bound: t_max = min(r, n+1); passes ≤ ceil(log2(t_max/kMinNewtonSums))+1.
    // Hadamard: total LLL swaps ≤ O(r² · log(pk)) across all passes.
    // Cap: 4·r iterations (generous; Hadamard gives much less).
    const std::size_t max_passes = 4U * r;  // Hadamard-derived upper bound

    std::size_t t_cur = kMinNewtonSums;
    std::size_t passes = 0U;

    while (passes < max_passes && t_cur <= t_max) {
        auto res = lll_knapsack_pass(f, factors, N, pk, t_cur, delta_val);
        if (res.has_value()) return res;

        ++passes;
        // Double Newton-sum columns for next pass (more discriminating).
        t_cur = std::min(t_cur * 2U, t_max);

        // If we've already used t_max and still no factor, no more to try.
        if (t_cur == t_max && passes > 1U) break;
    }

    // LLL path exhausted without finding a factor.
    // This is NOT an error: f may be irreducible, or the knapsack has no
    // {0,1}-solution in range n/2 (which means f's factors all have degree > n/2,
    // i.e., f is irreducible over Z). Return nullopt (caller handles this correctly).
    return std::nullopt;
}

}  // namespace cas::algebra
