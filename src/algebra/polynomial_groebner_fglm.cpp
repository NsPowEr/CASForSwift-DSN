// polynomial_groebner_fglm.cpp
//
// FGLM order-conversion algorithm (Faugère-Gianni-Lazard-Mora 1993).
// Reference: Faugère, Gianni, Lazard, Mora, "Efficient computation of
// zero-dimensional Gröbner bases by change of ordering",
// J. Symbolic Computation 16 (1993) 329-344.
//
// Algorithm summary (see fglm_convert below for pseudocode mapping):
//   Input : reduced GRevLex Groebner basis G of a zero-dimensional ideal I ⊆ Q[x_1..x_n].
//   Step 1: Enumerate standard monomials (monomial staircase) under GRevLex.
//           A monomial m is standard iff no lm(g) divides m for any g ∈ G.
//           Since I is zero-dimensional, the staircase is finite and its
//           cardinality D = dim_Q(Q[x]/I) (Hilbert series constant term).
//   Step 2: Build multiplication matrices M_xi (D×D over Q):
//           M_xi[j][k] = coord_k of NF(xi * sigma_j, G, GRevLex),
//           where sigma_0..sigma_{D-1} are the standard monomials in some order.
//   Step 3: Enumerate monomials in Lex order. For each new monomial t,
//           compute its coordinate vector v(t) = M_{x1}^{e1} ... M_{xn}^{en} * e_0
//           (or update by one step: v(xi*t') = M_xi * v(t') for t = xi * t').
//   Step 4: Attempt to express v(t) as Q-linear combination of already-known
//           vectors. If dependent: emit Lex basis polynomial
//               h = t - sum c_k * sigma_k,
//           where the c_k come from the dependency. The sigma_k are from the
//           already-found "Lex staircase".
//           If independent: add sigma_k := t to the Lex staircase.
//   Step 5: Terminate when Lex staircase size equals D.
//
// Complexity: O(D^3 + D^2 * |G|) field operations. Dramatically faster than
// computing the Lex basis directly for most systems (Bayer-Stillman 1987).
//
// REGOLA ZERO compliance:
// - All arithmetic over Rational (exact, no double/float).
// - D is bounded by the ideal's true vector space dimension; no artificial cap.
//   A configurable guard `ctx.fglm_max_dimension()` (default 512) causes an
//   explicit Unimplemented diagnostic if exceeded — hardware safety, never silent.
// - Standard monomial enumeration terminates in finite time because I is
//   zero-dimensional (staircase is finite by Hilbert basis theorem + 0-dim cond).

#include "polynomial_groebner_fglm.hpp"
#include "polynomial_groebner_fglm_internal.hpp"
#include "polynomial_groebner_f4.hpp"
#include "algebra_internal.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <optional>
#include <vector>
#include <map>
#include <set>

namespace cas::algebra {

// ── Internal helpers ────────────────────────────────────────────────────────

namespace {

// Check whether monomial m is a standard monomial w.r.t. basis G under order.
// m is standard iff no leading monomial of any g ∈ G divides m.
[[nodiscard]] bool is_standard_monomial(
    const Monomial& m,
    const std::vector<PolyF4>& G,
    MonomialOrder order)
{
    for (const PolyF4& g : G) {
        const Monomial lmg = g.leading_monomial(order);
        if (lmg.empty()) continue;
        // Check divisibility: lmg | m iff lmg[i] <= m[i] for all i
        bool divides = true;
        for (std::size_t i = 0; i < lmg.size(); ++i) {
            if (lmg[i] > m[i]) { divides = false; break; }
        }
        if (divides) return false;
    }
    return true;
}

// Comparators for monomial ordering traversal
struct LexLt {
    bool operator()(const Monomial& a, const Monomial& b) const {
        // Lex order: a < b iff a precedes b in lex ordering
        for (std::size_t i = 0; i < a.size(); ++i) {
            if (a[i] < b[i]) return true;
            if (a[i] > b[i]) return false;
        }
        return false;
    }
};

// Total-degree ascending, then lex tie-break — used for BFS enumeration
// over monomials (both GRevLex staircase and Lex staircase use this to
// enumerate monomials in a controlled finite-first order).
struct TotalDegreeLexLt {
    bool operator()(const Monomial& a, const Monomial& b) const {
        unsigned int da = 0; for (unsigned int e : a) da += e;
        unsigned int db = 0; for (unsigned int e : b) db += e;
        if (da != db) return da < db;
        return a < b; // lex tie-break
    }
};

// normal_form_vector: given monomial m, reduce t*G_grevlex and return
// the coordinate vector in terms of the standard_monomials basis (Q^D).
// The result is: NF(m, G, GRevLex) expressed as linear combination over standard_mons.
// If NF has a term not in standard_mons (shouldn't happen if G is complete), error.
[[nodiscard]] Result<std::vector<Rational>> normal_form_vector(
    const Monomial& m,
    const std::vector<PolyF4>& G,
    const std::vector<Monomial>& standard_monomials,
    MonomialOrder grevlex_order)
{
    // Build poly = {m -> 1}
    PolyF4 poly;
    poly.terms[m] = Rational(1);

    // Reduce by G (fully)
    bool changed = true;
    while (changed && !poly.is_zero()) {
        changed = false;
        for (auto it = poly.terms.begin(); it != poly.terms.end(); ) {
            const Monomial& curr = it->first;
            bool reduced = false;
            for (const PolyF4& g : G) {
                const Monomial lmg = g.leading_monomial(grevlex_order);
                if (lmg.empty()) continue;
                bool divides = true;
                for (std::size_t i = 0; i < lmg.size(); ++i) {
                    if (lmg[i] > curr[i]) { divides = false; break; }
                }
                if (!divides) continue;
                // Eliminate: poly -= factor * shift(g)
                Rational factor = it->second / g.leading_coefficient(grevlex_order);
                Monomial shift(curr.size());
                for (std::size_t i = 0; i < curr.size(); ++i) shift[i] = curr[i] - lmg[i];
                for (const auto& [gm, gc] : g.terms) {
                    Monomial nm(gm.size());
                    for (std::size_t i = 0; i < gm.size(); ++i) nm[i] = gm[i] + shift[i];
                    poly.terms[nm] = poly.terms[nm] - factor * gc;
                    if (poly.terms[nm].numerator().is_zero()) poly.terms.erase(nm);
                }
                changed = true;
                reduced = true;
                it = poly.terms.begin();
                break;
            }
            if (!reduced) ++it;
        }
    }

    // Build coordinate vector
    const std::size_t D = standard_monomials.size();
    std::vector<Rational> vec(D, Rational(0));

    // Build a map: standard monomial -> index (built once ideally, passed in)
    // Here we build it inline for simplicity (O(D) lookup per coefficient)
    for (const auto& [term, coeff] : poly.terms) {
        bool found = false;
        for (std::size_t k = 0; k < D; ++k) {
            if (standard_monomials[k] == term) {
                vec[k] = coeff;
                found = true;
                break;
            }
        }
        if (!found) {
            return fail<std::vector<Rational>>(make_error(
                CASErrorKind::Unimplemented,
                "FGLM: normal form has term outside standard monomial basis — "
                "input GRevLex basis may not be zero-dimensional"));
        }
    }
    return ok(std::move(vec));
}

} // namespace

// ── Public: enumerate_standard_monomials ────────────────────────────────────

std::vector<Monomial> enumerate_standard_monomials(
    const std::vector<PolyF4>& G_grevlex,
    std::size_t n_vars,
    std::size_t max_dim)
{
    // BFS over monomials in total-degree order; collect standard ones.
    // Terminate when we have collected max_dim standard monomials or when
    // the frontier is exhausted (zero-dimensional: finite staircase).
    // The zero monomial (constant 1) is always standard for non-trivial ideals.

    std::vector<Monomial> result;
    std::set<Monomial> visited;
    std::set<Monomial, TotalDegreeLexLt> frontier;

    Monomial one(n_vars, 0);
    frontier.insert(one);
    visited.insert(one);

    // Maximum degree to consider: bounded by D (unknown a priori) but we use
    // max_dim as a hard cap on staircase size, not on degree directly.
    // Degree cap: for a zero-dimensional ideal with D standard monomials,
    // the max degree is at most D (since 1, x1, x1^2,... alone gives D monomials).
    // We bound degree at max_dim to prevent infinite expansion in non-0-dim case.
    const unsigned int degree_cap = static_cast<unsigned int>(max_dim);

    while (!frontier.empty() && result.size() < max_dim) {
        auto it = frontier.begin();
        Monomial m = *it;
        frontier.erase(it);

        if (is_standard_monomial(m, G_grevlex, MonomialOrder::GRevLex)) {
            result.push_back(m);
            // Generate children: multiply by each variable
            for (std::size_t i = 0; i < n_vars; ++i) {
                Monomial child = m;
                child[i]++;
                // Degree check: don't expand beyond cap
                unsigned int deg = 0;
                for (unsigned int e : child) deg += e;
                if (deg > degree_cap) continue;
                if (visited.find(child) == visited.end()) {
                    visited.insert(child);
                    frontier.insert(child);
                }
            }
        }
        // If non-standard, do NOT expand its children (they are also non-standard
        // or will be caught by their own reduction path).
        // Actually we need to still add children even of non-standard monomials
        // to avoid missing standard ones. But wait — if lm(g)|m for some g, it
        // doesn't imply lm(g)|m*xi. So we must also expand non-standard monomials.
        // CORRECTION: We must explore all neighbors. The staircase is not an
        // order ideal in general for arbitrary orders. Re-add expansion below.
        else {
            // Expand non-standard monomial too (its children might be standard
            // via different reductions), but do NOT add it to result.
            for (std::size_t i = 0; i < n_vars; ++i) {
                Monomial child = m;
                child[i]++;
                unsigned int deg = 0;
                for (unsigned int e : child) deg += e;
                if (deg > degree_cap) continue;
                if (visited.find(child) == visited.end()) {
                    visited.insert(child);
                    frontier.insert(child);
                }
            }
        }
    }
    return result;
}

// ── Public: build_multiplication_matrix ─────────────────────────────────────

Result<std::vector<std::vector<Rational>>> build_multiplication_matrix(
    const std::vector<PolyF4>& G_grevlex,
    const std::vector<Monomial>& standard_monomials,
    std::size_t var_index,
    MonomialOrder grevlex_order)
{
    // M_xi[j][k] = coefficient of sigma_k in NF(xi * sigma_j, G, GRevLex).
    const std::size_t D = standard_monomials.size();
    std::vector<std::vector<Rational>> M(D, std::vector<Rational>(D, Rational(0)));

    for (std::size_t j = 0; j < D; ++j) {
        // Compute xi * sigma_j
        Monomial xi_sigma = standard_monomials[j];
        xi_sigma[var_index]++;
        // Compute normal form vector
        auto nf = normal_form_vector(xi_sigma, G_grevlex, standard_monomials, grevlex_order);
        if (nf.is_error()) return fail<std::vector<std::vector<Rational>>>(nf.error());
        M[j] = std::move(nf.value());
    }
    // Transpose: M_xi[i][j] = M[j][i] since above computed row = column j of result
    // Wait: we defined M_xi[j][k] = coord of sigma_k in NF(xi*sigma_j).
    // Above: for row j, M[j][k] = coord of sigma_k in NF(xi*sigma_j). Correct.
    return ok(std::move(M));
}

// ── Public: fglm_convert ────────────────────────────────────────────────────
//
// FGLM pseudocode (matching implementation):
//
//  1. Enumerate standard_monomials = {sigma_0..sigma_{D-1}} under GRevLex.
//  2. For each variable xi, build multiplication matrix M_xi (D×D over Q):
//       M_xi[j][k] = coord of sigma_k in NF(xi*sigma_j, G_grevlex).
//  3. Initialise:
//       lex_staircase = [] (monomials independent in Lex ordering)
//       known_vecs = []    (corresponding coordinate vectors in Q^D)
//       lex_basis = []     (output Lex Groebner basis elements)
//       v_one = e_0        (vector for monomial 1 = sigma_0 = identity in Q^D)
//  4. For each monomial t in Lex order (BFS by Lex degree, then lex):
//       a. Compute v(t): if t = 1, use v_one;
//          else t = xi * t', where t' = t/xi (previous Lex staircase monomial),
//          v(t) = M_xi * v(t').
//       b. Try linear_dependency(known_vecs, v(t)).
//          If DEPENDENT with coefficients c_0..c_{|lex_staircase|-1}:
//            Emit h = t - sum_{k} c_k * lex_staircase[k]  (Lex Groebner element)
//            Add h to lex_basis.
//          If INDEPENDENT:
//            Add t to lex_staircase; add v(t) to known_vecs.
//  5. Terminate when |lex_staircase| == D (staircase complete).
//     Return lex_basis after making each element monic w.r.t. Lex.

Result<std::vector<PolyF4>> fglm_convert(
    const std::vector<PolyF4>& G_grevlex,
    std::size_t n_vars,
    symbolic::CASContext* ctx)
{
    const std::size_t max_dim = ctx ? ctx->fglm_max_dimension() : std::size_t{512U};

    if (G_grevlex.empty()) return ok(std::vector<PolyF4>{});

    // Step 1: Enumerate standard monomials under GRevLex
    std::vector<Monomial> standard_monomials =
        enumerate_standard_monomials(G_grevlex, n_vars, max_dim);
    const std::size_t D = standard_monomials.size();

    if (D == 0) {
        // Trivial ideal (basis contains constant) → Lex basis is also {1}
        return ok(G_grevlex);
    }

    if (D >= max_dim) {
        return fail<std::vector<PolyF4>>(make_error(
            CASErrorKind::Unimplemented,
            "FGLM: standard monomial dimension " + std::to_string(D) +
            " reached fglm_max_dimension cap (" + std::to_string(max_dim) +
            "). Ideal may not be zero-dimensional, or increase ctx.fglm_max_dimension()."));
    }

    // Build index map: standard monomial -> position in Q^D
    std::map<Monomial, std::size_t> std_mon_index;
    for (std::size_t k = 0; k < D; ++k) std_mon_index[standard_monomials[k]] = k;

    // Step 2: Build multiplication matrices M_xi for each variable
    std::vector<std::vector<std::vector<Rational>>> mult_matrices(n_vars);
    for (std::size_t xi = 0; xi < n_vars; ++xi) {
        auto M = build_multiplication_matrix(
            G_grevlex, standard_monomials, xi, MonomialOrder::GRevLex);
        if (M.is_error()) return fail<std::vector<PolyF4>>(M.error());
        mult_matrices[xi] = std::move(M.value());
    }

    // Step 3: Initialise FGLM traversal
    // v_one: coordinate vector of the constant 1 in Q^D
    // sigma_0 should be the zero monomial (constant 1), find it
    std::vector<Rational> v_one(D, Rational(0));
    {
        Monomial zero_mon(n_vars, 0);
        auto it = std_mon_index.find(zero_mon);
        if (it == std_mon_index.end()) {
            return fail<std::vector<PolyF4>>(make_error(
                CASErrorKind::Unimplemented,
                "FGLM: constant 1 not in standard monomials — "
                "ideal appears to have no solutions"));
        }
        v_one[it->second] = Rational(1);
    }

    // Lex staircase: monomials known to be independent
    std::vector<Monomial> lex_staircase;
    std::vector<std::vector<Rational>> known_vecs;
    // Map: Lex staircase monomial -> its coordinate vector (for later lookup)
    std::map<Monomial, std::vector<Rational>> lex_staircase_vecs;

    std::vector<PolyF4> lex_basis;

    // Step 4: BFS over monomials in Lex order
    // We traverse monomials in increasing Lex order. For each t:
    // - If t = 1: use v_one
    // - If t = xi * t' for some t' in lex_staircase:
    //     v(t) = M_xi * v(t')
    // The key insight: to apply M_xi, t' must be in the staircase (known vector).
    // We track: for each staircase monomial, its vector.

    // BFS: start from {1}, expand by multiplying by each variable.
    // We need Lex ordering; use a priority queue with Lex comparator.
    struct LexQueueEntry {
        Monomial mon;
        std::size_t parent_xi; // which variable was multiplied
        Monomial parent_mon;   // the parent monomial = mon/x_{parent_xi}
    };

    // Use a set ordered by Lex for deterministic traversal
    std::set<Monomial, LexLt> lex_frontier;
    std::set<Monomial> lex_visited;
    Monomial one_mon(n_vars, 0);
    lex_frontier.insert(one_mon);
    lex_visited.insert(one_mon);

    // Loop until frontier exhausted. After staircase reaches D, all subsequent
    // monomials are necessarily dependent (D vectors span the quotient space),
    // and each emits a Lex Groebner basis element. We MUST process them all to
    // obtain the complete reduced Lex basis — terminating at D yields an
    // incomplete basis missing generators for border monomials (e.g. x² when
    // staircase = {1, y, x, xy} on input `<x²+y²-1, x²-y>` — solver
    // shape-lemma recovery then has no equation linking x to roots of the
    // pure-last polynomial).
    while (!lex_frontier.empty()) {
        // HC-F70-A33: poll interrupt at FGLM BFS outer iteration.
        if (ctx) { if (auto chk = ctx->check_interrupt(); chk.is_error()) return fail<std::vector<PolyF4>>(chk.error()); }
        Monomial t = *lex_frontier.begin();
        lex_frontier.erase(lex_frontier.begin());

        // Compute v(t)
        std::vector<Rational> vt;
        bool vt_valid = false;

        if (t == one_mon) {
            vt = v_one;
            vt_valid = true;
        } else {
            // Find some i such that t[i] > 0 and t/xi is in lex_staircase
            for (std::size_t xi = 0; xi < n_vars; ++xi) {
                if (t[xi] == 0) continue;
                Monomial parent = t;
                parent[xi]--;
                auto pit = lex_staircase_vecs.find(parent);
                if (pit == lex_staircase_vecs.end()) continue;
                // v(t) = M_xi * v(parent)
                vt = mat_vec_mul(mult_matrices[xi], pit->second);
                vt_valid = true;
                break;
            }
            // If no parent in staircase, this monomial's vector can't be computed
            // (its predecessor was declared dependent/emitted as basis element).
            // We skip it — it cannot be a staircase element either.
            if (!vt_valid) continue;
        }

        // Try linear dependency
        auto dep = linear_dependency(known_vecs, vt);
        if (dep.has_value()) {
            // DEPENDENT: emit Lex Groebner basis element
            // h = t - sum_{k=0}^{|lex_staircase|-1} dep[k] * lex_staircase[k]
            PolyF4 h;
            h.terms[t] = Rational(1);
            const auto& coeffs = dep.value();
            for (std::size_t k = 0; k < coeffs.size(); ++k) {
                if (coeffs[k].numerator().is_zero()) continue;
                Monomial sk = lex_staircase[k];
                h.terms[sk] = h.terms[sk] - coeffs[k];
                if (h.terms[sk].numerator().is_zero()) h.terms.erase(sk);
            }
            if (!h.is_zero()) {
                h.make_monic(MonomialOrder::Lex);
                lex_basis.push_back(std::move(h));
            }
            // Do NOT add t to staircase. Do NOT expand children of t.
        } else {
            // INDEPENDENT: add to staircase
            lex_staircase.push_back(t);
            known_vecs.push_back(vt);
            lex_staircase_vecs[t] = vt;

            // Expand children in Lex order
            for (std::size_t xi = 0; xi < n_vars; ++xi) {
                Monomial child = t;
                child[xi]++;
                // Lex degree guard: don't go beyond max_dim total degree
                unsigned int deg = 0;
                for (unsigned int e : child) deg += e;
                if (deg > static_cast<unsigned int>(max_dim)) continue;
                if (lex_visited.find(child) == lex_visited.end()) {
                    lex_visited.insert(child);
                    lex_frontier.insert(child);
                }
            }
        }
    }

    if (lex_staircase.size() < D) {
        // Not all D dimensions covered; staircase incomplete
        // This can happen if the BFS degree cap was hit.
        return fail<std::vector<PolyF4>>(make_error(
            CASErrorKind::Unimplemented,
            "FGLM: Lex staircase incomplete (" +
            std::to_string(lex_staircase.size()) + "/" + std::to_string(D) +
            "). Increase ctx.fglm_max_dimension() or check ideal is zero-dimensional."));
    }

    // Inter-reduce the Lex basis (final reduced Groebner basis)
    inter_reduce(lex_basis, MonomialOrder::Lex);
    for (auto& g : lex_basis) g.make_monic(MonomialOrder::Lex);

    return ok(std::move(lex_basis));
}

} // namespace cas::algebra
