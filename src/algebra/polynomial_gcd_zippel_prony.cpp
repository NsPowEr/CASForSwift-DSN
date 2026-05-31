// polynomial_gcd_zippel_prony.cpp — REAL Zippel sparse interpolation GCD (F3.1).
// Ref: Zippel "Probabilistic algorithms for sparse polynomials", EUROSAM 1979.
// Survey: Geddes-Czapor-Labahn §7.6.
//
// Pipeline:
//   Phase 1 (skeleton): one anchor evaluation (a_2..a_n) ∈ (Fp*)^{n-1}, run
//     univariate gcd_fp in x_1, record deg_x1(g_anchor).  Extract skeleton:
//     for each x_1-degree k, the monomial set {m_j(x_2..x_n)} appearing in
//     P-coeff-of-x_1^k ∪ Q-coeff-of-x_1^k (strict superset of true gcd skeleton).
//   Phase 2 (Prony Vandermonde): for each k, T = |skeleton_k| extra random
//     evaluations; assemble Vandermonde system m_j(point_t), solve via Gauss
//     elimination over Fp.  Verify recovered α_{k,j} at one fresh point.
//   Lift center-mod-p to Z; certify divisibility g|P and g|Q in Z[x_1..x_n];
//   on cert failure emit Unimplemented (multi-prime CRT extension required).
//
// Sample-count budget: T = max_k |skeleton_k| + extra where
//   extra = ceil(log2(1/δ)) with δ = ctx.gcd_error_probability() (default 1e-3
//   ⇒ extra ≥ 10) — Schwartz-Zippel.

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cas::algebra {

// Defined in polynomial_gcd_brown_modular.cpp.
using ZMonomial   = std::vector<unsigned int>;
using ZSparsePoly = std::map<ZMonomial, BigInt>;

namespace {

static BigInt pos_mod_z(const BigInt& a, const BigInt& m) {
    BigInt r = a % m;
    if (r.is_negative()) r += m.abs();
    return r;
}

// bigint_gcd_z removed (unused); arithmetic helpers live in polynomial_gcd_brown_modular.cpp.

static ZSparsePoly to_sparse_z(const MultivariatePolynomial& p,
                                const std::vector<Symbol>& vars) {
    ZSparsePoly sp;
    for (const auto& term : p.terms()) {
        ZMonomial mono(vars.size(), 0U);
        for (const auto& [sym, exp] : term.factors)
            for (std::size_t i = 0; i < vars.size(); ++i)
                if (vars[i].name == sym.name) mono[i] = exp;
        sp[mono] += term.coefficient;
        if (sp[mono].is_zero()) sp.erase(mono);
    }
    return sp;
}

static MultivariatePolynomial from_sparse_z(const ZSparsePoly& sp,
                                             const std::vector<Symbol>& vars) {
    std::vector<MultivariateTerm> terms;
    for (const auto& [m, c] : sp) {
        if (c.is_zero()) continue;
        std::vector<std::pair<Symbol, unsigned int>> f;
        for (std::size_t i = 0; i < vars.size(); ++i)
            if (m[i] > 0U) f.emplace_back(vars[i], m[i]);
        terms.push_back(MultivariateTerm{ .coefficient = c, .factors = std::move(f) });
    }
    return MultivariatePolynomial(std::move(terms));
}

static std::vector<Symbol> collect_vars_z(const MultivariatePolynomial& p,
                                           const MultivariatePolynomial& q) {
    auto vp = p.variables();
    auto vq = q.variables();
    std::vector<Symbol> all = vp;
    for (const auto& s : vq)
        if (std::none_of(all.begin(), all.end(),
                         [&](const Symbol& a){ return a.name == s.name; }))
            all.push_back(s);
    std::sort(all.begin(), all.end(),
              [](const Symbol& a, const Symbol& b){ return a.name < b.name; });
    return all;
}

static std::size_t deg_in_var_z(const ZSparsePoly& sp, std::size_t vi) {
    std::size_t d = 0U;
    for (const auto& [m, _] : sp)
        if (vi < m.size()) d = std::max<std::size_t>(d, m[vi]);
    return d;
}

static std::size_t total_degree_z(const ZSparsePoly& sp) {
    std::size_t d = 0;
    for (const auto& [m, _] : sp) {
        std::size_t td = 0; for (auto e : m) td += e;
        d = std::max(d, td);
    }
    return d;
}

// Evaluate full multivariate sparse poly at (a_2..a_n) leaving x_1 free.
// Returns dense vector of coefficients (BigInt mod p), indexed by degree in x_1.
static std::vector<BigInt> eval_skeleton_to_univariate_mod_p(
    const ZSparsePoly& sp, const std::vector<BigInt>& a,
    const BigInt& p) {
    // a[i] = value for variable index (i+1).  Variable 0 = x_1 stays free.
    std::size_t d1 = deg_in_var_z(sp, 0U);
    std::vector<BigInt> coeffs(d1 + 1U, BigInt(0));
    for (const auto& [mono, c] : sp) {
        BigInt val = pos_mod_z(c, p);
        std::size_t k1 = (0U < mono.size()) ? mono[0] : 0U;
        for (std::size_t i = 1; i < mono.size(); ++i) {
            unsigned int e = mono[i];
            if (e == 0U) continue;
            BigInt ai = (i - 1U < a.size()) ? a[i - 1U] : BigInt(0);
            BigInt acc(1);
            for (unsigned int t = 0; t < e; ++t) acc = pos_mod_z(acc * ai, p);
            val = pos_mod_z(val * acc, p);
        }
        if (k1 >= coeffs.size()) coeffs.resize(k1 + 1U, BigInt(0));
        coeffs[k1] = pos_mod_z(coeffs[k1] + val, p);
    }
    return coeffs;
}

// Dense univariate gcd in Fp[x] given dense BigInt coefficient vectors mod p.
static std::vector<BigInt> univariate_gcd_fp_dense(
    std::vector<BigInt> a, std::vector<BigInt> b, const BigInt& p) {
    while (!a.empty() && a.back().is_zero()) a.pop_back();
    while (!b.empty() && b.back().is_zero()) b.pop_back();
    while (!b.empty()) {
        auto inv = numtheory::modular_inverse(pos_mod_z(b.back(), p), p);
        if (inv.is_error()) return {};
        BigInt ilc = inv.value();
        std::vector<BigInt> rem = a;
        while (rem.size() >= b.size() && !rem.empty()) {
            std::size_t dd = rem.size() - b.size();
            BigInt factor = pos_mod_z(rem.back() * ilc, p);
            for (std::size_t i = 0; i < b.size(); ++i)
                rem[i + dd] = pos_mod_z(rem[i + dd] - factor * b[i], p);
            while (!rem.empty() && rem.back().is_zero()) rem.pop_back();
        }
        a = std::move(b);
        b = std::move(rem);
    }
    if (a.empty()) return {};
    // Monic.
    auto inv = numtheory::modular_inverse(pos_mod_z(a.back(), p), p);
    if (inv.is_error()) return {};
    BigInt ilc = inv.value();
    for (auto& c : a) c = pos_mod_z(c * ilc, p);
    return a;
}


}  // namespace

// ── REAL Zippel Prony GCD ────────────────────────────────────────────────────
[[nodiscard]] Result<MultivariatePolynomial> gcd_zippel_prony(
    const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx, std::size_t* out_samples_used) {

    if (P.is_zero()) return ok(Q);
    if (Q.is_zero()) return ok(P);

    const std::vector<Symbol> vars = collect_vars_z(P, Q);
    const std::size_t n = vars.size();
    if (n < 2U) {
        // Zippel sparse is for multivariate inputs; univariate → not applicable.
        return make_unimplemented<MultivariatePolynomial>(
            "algebra", "gcd_zippel_prony",
            "n_vars=" + std::to_string(n),
            "ZIPPEL_PRONY_NOT_MULTIVARIATE",
            "Use univariate path for n<2 (gcd_integer_poly_dispatch)",
            "F3.1");
    }

    ZSparsePoly spP = to_sparse_z(P, vars);
    ZSparsePoly spQ = to_sparse_z(Q, vars);

    // Pick anchor prime near 2^30 with input-dependent hash.
    std::size_t h = (spP.size() * 2654435761ULL) ^ (spQ.size() * 40503ULL);
    long long start = 1073741827LL + static_cast<long long>(h % 65537ULL);
    auto np0 = numtheory::next_prime(BigInt(start - 1));
    if (np0.is_error())
        return make_unimplemented<MultivariatePolynomial>(
            "algebra", "gcd_zippel_prony", "next_prime_error",
            "ZIPPEL_PRONY_PRIME_INIT", "next_prime failed", "F3.1");
    BigInt p = np0.value();

    // Sample count derived from Schwartz-Zippel.
    // T_needed = ceil(log(1/δ) / log(1/(1 - d/q)))  with δ = ctx.gcd_error_probability(),
    // d = total degree bound, q = p.
    const std::size_t dtot = std::max(total_degree_z(spP), total_degree_z(spQ));
    const double delta = ctx.gcd_error_probability();
    const double q_d   = static_cast<double>(dtot) /
                         std::max(1.0, static_cast<double>(p.bit_length()) * 0.30103);
    // (d / q) ≈ d / 10^{bits·log10(2)} — but we want d/q.  Use a robust lower
    // bound: assume q ≥ 2^30 ⇒ d/q ≤ d / 10^9 ≈ 0.  We thus need very few extra
    // samples beyond the skeleton size.  Cap at 4 redundant samples.
    (void)q_d;
    const std::size_t extra_samples =
        std::max<std::size_t>(2U,
            (delta > 0.0 && delta < 1.0)
                ? static_cast<std::size_t>(std::ceil(std::log2(1.0 / delta)))
                : 10U);

    // ── Phase 1: Skeleton anchor.
    // Pick random a = (a_2..a_n) in (Fp*)^{n-1}.  Use deterministic input-hash
    // seed for reproducibility.
    std::size_t samples_used = 0;
    auto eval_full_at = [&](const std::vector<BigInt>& a, std::vector<BigInt>& Pa_out,
                             std::vector<BigInt>& Qa_out) {
        Pa_out = eval_skeleton_to_univariate_mod_p(spP, a, p);
        Qa_out = eval_skeleton_to_univariate_mod_p(spQ, a, p);
        ++samples_used;
    };

    // Deterministic node generator over Fp* (skip 0).
    auto gen_node = [&](std::size_t k) {
        long long base = static_cast<long long>(7U + ((k * 2654435761ULL) % 1000003ULL));
        return BigInt(base);
    };

    std::vector<BigInt> a0(n - 1U);
    for (std::size_t i = 0; i < n - 1U; ++i) a0[i] = gen_node(i + 1U);

    std::vector<BigInt> Pa, Qa;
    eval_full_at(a0, Pa, Qa);
    auto g_anchor = univariate_gcd_fp_dense(Pa, Qa, p);
    if (g_anchor.empty())
        return make_unimplemented<MultivariatePolynomial>(
            "algebra", "gcd_zippel_prony", "empty_anchor",
            "ZIPPEL_PRONY_ANCHOR_VOID",
            "Anchor evaluation produced an empty gcd — unlucky point",
            "F3.1");

    const std::size_t d1 = g_anchor.size() - 1U;  // deg in x_1

    // Skeleton support: for each x_1-degree k, the SET of monomials in x_2..x_n
    // that appear in (P's coeff of x_1^k) ∪ (Q's coeff of x_1^k).  Strict
    // superset of the true gcd's x_1^k-coefficient support (Zippel stability).
    std::vector<std::set<ZMonomial>> skeleton(d1 + 1U);
    auto extract_x1_coeffs = [&](const ZSparsePoly& sp,
                                   std::vector<std::set<ZMonomial>>& sk) {
        for (const auto& [m, c] : sp) {
            std::size_t k = (m.size() > 0U) ? m[0] : 0U;
            if (k > d1) continue;
            ZMonomial mr = m;
            if (!mr.empty()) mr[0] = 0U;
            sk[k].insert(mr);
        }
    };
    extract_x1_coeffs(spP, skeleton);
    extract_x1_coeffs(spQ, skeleton);

    // ── Phase 2: Prony Vandermonde recovery per x_1-degree.
    // For each k, T_k = |skeleton[k]| evaluations resolve α_{k,j}.
    // We perform max_T = max_k T_k + extra_samples total full samples, reusing
    // them across all k.
    std::size_t max_T = 0;
    for (auto& s : skeleton) max_T = std::max(max_T, s.size());
    if (max_T == 0U)
        return make_unimplemented<MultivariatePolynomial>(
            "algebra", "gcd_zippel_prony", "empty_skeleton",
            "ZIPPEL_PRONY_EMPTY_SKELETON",
            "Skeleton extraction yielded no monomials — input too small",
            "F3.1");
    const std::size_t T = max_T + extra_samples;

    // Generate T distinct random points b^(t) ∈ (Fp*)^{n-1}.
    std::vector<std::vector<BigInt>> points(T);
    std::vector<std::vector<BigInt>> Pa_vals(T), Qa_vals(T);
    for (std::size_t t = 0; t < T; ++t) {
        points[t].resize(n - 1U);
        for (std::size_t i = 0; i < n - 1U; ++i)
            points[t][i] = gen_node((t + 1U) * (n - 1U) + i + 1U);
        eval_full_at(points[t], Pa_vals[t], Qa_vals[t]);
    }
    // Per-point gcd in Fp[x_1].
    std::vector<std::vector<BigInt>> gp_vals(T);
    for (std::size_t t = 0; t < T; ++t) {
        gp_vals[t] = univariate_gcd_fp_dense(Pa_vals[t], Qa_vals[t], p);
        // Degree-stability invariant: reject if degree differs from anchor.
        if (gp_vals[t].size() != g_anchor.size()) {
            return make_unimplemented<MultivariatePolynomial>(
                "algebra", "gcd_zippel_prony",
                "anchor_deg=" + std::to_string(d1) + ",this_deg=" +
                    std::to_string(gp_vals[t].empty() ? 0 : gp_vals[t].size() - 1),
                "ZIPPEL_PRONY_UNSTABLE_DEGREE",
                "Random eval gave a different gcd degree — bad anchor; retry with new prime",
                "F3.1");
        }
    }

    // Recover g_p ∈ Fp[x_1..x_n] mod p, monomial by monomial.
    ZSparsePoly g_mod_p;
    for (std::size_t k = 0; k <= d1; ++k) {
        const auto& sk = skeleton[k];
        if (sk.empty()) continue;
        std::vector<ZMonomial> monos(sk.begin(), sk.end());
        const std::size_t Tk = monos.size();
        // Build Vandermonde-like nodes: at sample t, evaluate each candidate
        // monomial m_j ∈ x_2..x_n at points[t].  Row t, column j = m_j(points[t]) mod p.
        // RHS row t = gp_vals[t][k]  (k-th coefficient in x_1 of univariate gcd).
        std::vector<std::vector<BigInt>> A(Tk, std::vector<BigInt>(Tk + 1U, BigInt(0)));
        for (std::size_t t = 0; t < Tk; ++t) {
            for (std::size_t j = 0; j < Tk; ++j) {
                const ZMonomial& m = monos[j];
                BigInt v(1);
                for (std::size_t i = 1; i < m.size(); ++i) {
                    unsigned int e = m[i];
                    if (e == 0U) continue;
                    BigInt ai = points[t][i - 1U];
                    BigInt acc(1);
                    for (unsigned int s = 0; s < e; ++s) acc = pos_mod_z(acc * ai, p);
                    v = pos_mod_z(v * acc, p);
                }
                A[t][j] = v;
            }
            A[t][Tk] = (k < gp_vals[t].size()) ? pos_mod_z(gp_vals[t][k], p) : BigInt(0);
        }
        // Gaussian-eliminate.
        bool singular = false;
        for (std::size_t col = 0; col < Tk && !singular; ++col) {
            std::size_t piv = col;
            while (piv < Tk && A[piv][col].is_zero()) ++piv;
            if (piv == Tk) { singular = true; break; }
            if (piv != col) std::swap(A[piv], A[col]);
            auto inv = numtheory::modular_inverse(A[col][col], p);
            if (inv.is_error()) { singular = true; break; }
            BigInt ilc = inv.value();
            for (std::size_t j = col; j <= Tk; ++j)
                A[col][j] = pos_mod_z(A[col][j] * ilc, p);
            for (std::size_t r = 0; r < Tk; ++r) {
                if (r == col || A[r][col].is_zero()) continue;
                BigInt f = A[r][col];
                for (std::size_t j = col; j <= Tk; ++j)
                    A[r][j] = pos_mod_z(A[r][j] - f * A[col][j], p);
            }
        }
        if (singular) {
            return make_unimplemented<MultivariatePolynomial>(
                "algebra", "gcd_zippel_prony",
                "k=" + std::to_string(k) + ",Tk=" + std::to_string(Tk),
                "ZIPPEL_PRONY_SINGULAR_VANDERMONDE",
                "Random nodes produced singular Vandermonde — retry with fresh nodes/prime",
                "F3.1");
        }
        // Verification by re-evaluation at one fresh point (catches Zippel
        // "bad start" with prob ≤ δ).
        std::vector<BigInt> verify_node(n - 1U);
        for (std::size_t i = 0; i < n - 1U; ++i)
            verify_node[i] = gen_node((T + 17U) * (n - 1U) + i + 1U);
        BigInt predicted(0);
        for (std::size_t j = 0; j < Tk; ++j) {
            const ZMonomial& m = monos[j];
            BigInt mval(1);
            for (std::size_t i = 1; i < m.size(); ++i) {
                unsigned int e = m[i];
                if (e == 0U) continue;
                BigInt ai = verify_node[i - 1U];
                BigInt acc(1);
                for (unsigned int s = 0; s < e; ++s) acc = pos_mod_z(acc * ai, p);
                mval = pos_mod_z(mval * acc, p);
            }
            predicted = pos_mod_z(predicted + pos_mod_z(A[j][Tk], p) * mval, p);
        }
        std::vector<BigInt> Pv, Qv;
        Pv = eval_skeleton_to_univariate_mod_p(spP, verify_node, p);
        Qv = eval_skeleton_to_univariate_mod_p(spQ, verify_node, p);
        ++samples_used;
        auto gv = univariate_gcd_fp_dense(Pv, Qv, p);
        BigInt observed = (k < gv.size()) ? pos_mod_z(gv[k], p) : BigInt(0);
        if (observed != predicted) {
            return make_unimplemented<MultivariatePolynomial>(
                "algebra", "gcd_zippel_prony",
                "k=" + std::to_string(k) + ",observed!=predicted",
                "ZIPPEL_PRONY_VERIFICATION_FAILED",
                "Skeleton hypothesis rejected by fresh eval — Zippel bad start",
                "F3.1");
        }
        // Emit recovered monomials.
        for (std::size_t j = 0; j < Tk; ++j) {
            BigInt alpha = A[j][Tk];
            if (alpha.is_zero()) continue;
            ZMonomial m = monos[j];
            if (m.empty()) m.assign(n, 0U);
            if (m.size() < n) m.resize(n, 0U);
            m[0] = static_cast<unsigned int>(k);
            g_mod_p[m] = alpha;
        }
    }

    if (g_mod_p.empty())
        return make_unimplemented<MultivariatePolynomial>(
            "algebra", "gcd_zippel_prony", "empty_lift",
            "ZIPPEL_PRONY_EMPTY_LIFT",
            "Prony recovery yielded empty polynomial — fall back to Brown's modular",
            "F3.1");

    // Center-lift each coefficient from Fp to Z.
    ZSparsePoly g_lifted;
    for (const auto& [m, c] : g_mod_p) {
        BigInt cv = c;
        BigInt two_c = cv + cv;
        if (two_c > p) cv = cv - p;
        if (!cv.is_zero()) g_lifted[m] = cv;
    }

    MultivariatePolynomial g_cand = from_sparse_z(g_lifted, vars);

    // Certificate: direct sparse exact division in Z[x_1..x_n].
    auto certify = [&](const MultivariatePolynomial& a, const MultivariatePolynomial& b)
        -> bool {
        ZSparsePoly sa = to_sparse_z(a, vars);
        ZSparsePoly sb = to_sparse_z(b, vars);
        if (sb.empty()) return false;
        if (sa.empty()) return true;
        const std::size_t nv = vars.size();
        ZSparsePoly rem = sa;
        auto [dlm, dlc] = *std::prev(sb.end());
        const std::size_t budget = (rem.size() + 1U) * (sb.size() + 1U) + 16U;
        std::size_t steps = 0;
        while (!rem.empty()) {
            if (++steps > budget) return false;
            auto [rlm, rlc] = *std::prev(rem.end());
            for (std::size_t i = 0; i < nv; ++i) {
                std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
                std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
                if (ri < di) return false;
            }
            if ((rlc % dlc) != BigInt(0)) return false;
            BigInt qc = rlc / dlc;
            ZMonomial qm(nv, 0U);
            for (std::size_t i = 0; i < nv; ++i) {
                std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
                std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
                qm[i] = static_cast<unsigned int>(ri - di);
            }
            for (const auto& [dm, dc] : sb) {
                ZMonomial nm(nv, 0U);
                for (std::size_t i = 0; i < nv; ++i)
                    nm[i] = qm[i] + ((i < dm.size()) ? dm[i] : 0U);
                rem[nm] -= qc * dc;
                if (rem[nm].is_zero()) rem.erase(nm);
            }
        }
        return true;
    };

    if (!certify(P, g_cand) || !certify(Q, g_cand)) {
        if (out_samples_used) *out_samples_used = samples_used;
        // Single-prime Prony cert insufficient; dispatcher falls back to Brown's.
        return make_unimplemented<MultivariatePolynomial>(
            "algebra", "gcd_zippel_prony",
            "single-prime certificate failed; samples=" + std::to_string(samples_used),
            "ZIPPEL_PRONY_SINGLE_PRIME_CERT_FAILED",
            "Single-prime Prony lift insufficient — multi-prime CRT extension required",
            "F3.1");
    }

    // Normalize sign (lex-leading positive).
    {
        auto sp = to_sparse_z(g_cand, vars);
        if (!sp.empty()) {
            auto last = std::prev(sp.end());
            if (last->second.is_negative())
                for (auto& [_, c] : sp) c = -c;
            g_cand = from_sparse_z(sp, vars);
        }
    }
    if (out_samples_used) *out_samples_used = samples_used;
    return ok(std::move(g_cand));
}

// Overload without probe pointer.
[[nodiscard]] Result<MultivariatePolynomial> gcd_zippel_prony(
    const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx) {
    return gcd_zippel_prony(P, Q, ctx, nullptr);
}

}  // namespace cas::algebra
