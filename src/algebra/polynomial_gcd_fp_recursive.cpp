// polynomial_gcd_fp_recursive.cpp — Helpers + Brown-in-Fp recursive gcd.
// Implements the public API in polynomial_gcd_fp_internal.hpp.

#include "polynomial_gcd_fp_internal.hpp"
#include "cas/numtheory.hpp"

#include <algorithm>
#include <set>

namespace cas::algebra::fp_helpers {

BigInt pos_mod(const BigInt& a, const BigInt& m) {
    BigInt r = a % m;
    if (r.is_negative()) r += m.abs();
    return r;
}

BigInt bigint_gcd(BigInt a, BigInt b) {
    a = a.abs(); b = b.abs();
    while (!b.is_zero()) { BigInt t = a % b; a = std::move(b); b = std::move(t); }
    return a;
}

BigInt centered_repr(const BigInt& r, const BigInt& M) {
    BigInt two_r = r + r;
    return (two_r > M) ? (r - M) : r;
}

BSparsePoly reduce_sparse_mod_p(const BSparsePoly& sp, const BigInt& p) {
    BSparsePoly r;
    for (const auto& [mono, c] : sp) {
        BigInt cr = pos_mod(c, p);
        if (!cr.is_zero()) r[mono] = cr;
    }
    return r;
}

BigInt sparse_inf_norm(const BSparsePoly& sp) {
    BigInt best(0);
    for (const auto& [_, c] : sp) {
        BigInt ac = c.abs();
        if (ac > best) best = ac;
    }
    return best;
}

std::size_t deg_in_var(const BSparsePoly& sp, std::size_t vi) {
    std::size_t d = 0U;
    for (const auto& [m, _] : sp)
        if (vi < m.size()) d = std::max<std::size_t>(d, m[vi]);
    return d;
}

BSparsePoly eval_var_mod_p(const BSparsePoly& sp, std::size_t var_idx,
                            const BigInt& val, const BigInt& p) {
    BSparsePoly out;
    for (const auto& [mono, c] : sp) {
        BigInt new_c = c;
        if (var_idx < mono.size() && mono[var_idx] > 0U) {
            BigInt acc(1);
            for (unsigned int i = 0; i < mono[var_idx]; ++i)
                acc = pos_mod(acc * val, p);
            new_c = pos_mod(new_c * acc, p);
        }
        if (new_c.is_zero()) continue;
        BMonomial m2 = mono;
        if (var_idx < m2.size()) m2[var_idx] = 0U;
        out[m2] += new_c;
        out[m2] = pos_mod(out[m2], p);
        if (out[m2].is_zero()) out.erase(m2);
    }
    return out;
}

namespace {

// Univariate Fp gcd between two sparse polys with only one active variable.
static std::optional<BSparsePoly> univariate_sparse_gcd_fp(
    const BSparsePoly& A, const BSparsePoly& B,
    std::size_t var_idx, const BigInt& p) {
    auto to_dense = [&](const BSparsePoly& sp) -> std::vector<BigInt> {
        std::size_t d = deg_in_var(sp, var_idx);
        std::vector<BigInt> v(d + 1U, BigInt(0));
        for (const auto& [mono, c] : sp) {
            std::size_t k = (var_idx < mono.size()) ? mono[var_idx] : 0U;
            v[k] = pos_mod(v[k] + c, p);
        }
        while (!v.empty() && v.back().is_zero()) v.pop_back();
        return v;
    };
    auto a = to_dense(A);
    auto b = to_dense(B);
    while (!b.empty()) {
        BigInt inv_lc = pos_mod(b.back(), p);
        auto inv = numtheory::modular_inverse(inv_lc, p);
        if (inv.is_error()) return std::nullopt;
        BigInt ilc = inv.value();
        std::vector<BigInt> rem = a;
        while (!rem.empty() && rem.size() >= b.size()) {
            std::size_t dd = rem.size() - b.size();
            BigInt factor = pos_mod(rem.back() * ilc, p);
            for (std::size_t i = 0; i < b.size(); ++i)
                rem[i + dd] = pos_mod(rem[i + dd] - factor * b[i], p);
            while (!rem.empty() && rem.back().is_zero()) rem.pop_back();
        }
        a = std::move(b);
        b = std::move(rem);
    }
    if (a.empty()) return BSparsePoly{};
    auto inv = numtheory::modular_inverse(pos_mod(a.back(), p), p);
    if (inv.is_error()) return std::nullopt;
    BigInt ilc = inv.value();
    BSparsePoly out;
    for (std::size_t k = 0; k < a.size(); ++k) {
        BigInt c = pos_mod(a[k] * ilc, p);
        if (c.is_zero()) continue;
        BMonomial m;
        m.assign(var_idx + 1U, 0U);
        m[var_idx] = static_cast<unsigned int>(k);
        out[m] = c;
    }
    return out;
}

// Lagrange interpolation in Fp[x_var] from (b_i, val_i) pairs.
static std::optional<BSparsePoly> lagrange_interp_fp(
    const std::vector<BigInt>& bs, const std::vector<BSparsePoly>& vals,
    std::size_t var_idx, std::size_t target_deg, const BigInt& p) {
    const std::size_t n = bs.size();
    if (n != target_deg + 1U || n == 0U) return std::nullopt;
    std::set<BMonomial> all_monos;
    for (const auto& v : vals)
        for (const auto& [m, _] : v) all_monos.insert(m);
    BSparsePoly result;
    std::vector<std::vector<BigInt>> Li_polys(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<BigInt> num = { BigInt(1) };
        BigInt den(1);
        for (std::size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            std::vector<BigInt> nn(num.size() + 1U, BigInt(0));
            for (std::size_t k = 0; k < num.size(); ++k) {
                nn[k + 1U] = pos_mod(nn[k + 1U] + num[k], p);
                nn[k]      = pos_mod(nn[k]     - num[k] * bs[j], p);
            }
            num = std::move(nn);
            den = pos_mod(den * pos_mod(bs[i] - bs[j], p), p);
        }
        auto inv_d = numtheory::modular_inverse(den, p);
        if (inv_d.is_error()) return std::nullopt;
        for (auto& c : num) c = pos_mod(c * inv_d.value(), p);
        Li_polys[i] = std::move(num);
    }
    for (const auto& m : all_monos) {
        for (std::size_t i = 0; i < n; ++i) {
            auto it = vals[i].find(m);
            if (it == vals[i].end() || it->second.is_zero()) continue;
            const BigInt& vi = it->second;
            for (std::size_t k = 0; k < Li_polys[i].size(); ++k) {
                if (Li_polys[i][k].is_zero()) continue;
                BMonomial m2 = m;
                if (m2.size() <= var_idx) m2.resize(var_idx + 1U, 0U);
                m2[var_idx] = static_cast<unsigned int>(k);
                result[m2] = pos_mod(result[m2] + vi * Li_polys[i][k], p);
                if (result[m2].is_zero()) result.erase(m2);
            }
        }
    }
    return result;
}

// Decompose `sp` as a polynomial in `var_idx` over sub-vars: returns a map
// from y-degree → coefficient (a BSparsePoly in n-var space but with
// var_idx-exp = 0 throughout).  Useful for content/pp computation.
static std::map<std::size_t, BSparsePoly> layers_by_var(
        const BSparsePoly& sp, std::size_t var_idx) {
    std::map<std::size_t, BSparsePoly> out;
    for (const auto& [m, c] : sp) {
        std::size_t k = (var_idx < m.size()) ? m[var_idx] : 0U;
        BMonomial m2 = m;
        if (var_idx < m2.size()) m2[var_idx] = 0U;
        out[k][m2] = c;
    }
    return out;
}

// Reassemble layer-map (k → sub-coefficient poly) into a BSparsePoly with the
// k stored back into the var_idx slot.
static BSparsePoly reassemble_layers(
        const std::map<std::size_t, BSparsePoly>& layers,
        std::size_t var_idx, const BigInt& p) {
    BSparsePoly out;
    for (const auto& [k, layer] : layers) {
        for (const auto& [m, c] : layer) {
            if (c.is_zero()) continue;
            BMonomial nm = m;
            if (nm.size() <= var_idx) nm.resize(var_idx + 1U, 0U);
            nm[var_idx] = static_cast<unsigned int>(k);
            out[nm] = pos_mod(c, p);
            if (out[nm].is_zero()) out.erase(nm);
        }
    }
    return out;
}

// Multiply two sparse polys (in arbitrary n-var space) mod p.
static BSparsePoly mul_mod_p(const BSparsePoly& A, const BSparsePoly& B,
                               const BigInt& p) {
    BSparsePoly out;
    if (A.empty() || B.empty()) return out;
    std::size_t n = 0;
    for (const auto& [m, _] : A) n = std::max(n, m.size());
    for (const auto& [m, _] : B) n = std::max(n, m.size());
    for (const auto& [ma, ca] : A) {
        for (const auto& [mb, cb] : B) {
            BMonomial nm(n, 0U);
            for (std::size_t i = 0; i < n; ++i) {
                unsigned int ea = (i < ma.size()) ? ma[i] : 0U;
                unsigned int eb = (i < mb.size()) ? mb[i] : 0U;
                nm[i] = ea + eb;
            }
            BigInt v = pos_mod(out[nm] + ca * cb, p);
            if (v.is_zero()) out.erase(nm);
            else out[nm] = v;
        }
    }
    return out;
}

// Exact division of A by B in Fp[x_1..x_n] (sparse, lex-leading).  Returns
// nullopt if not exactly divisible.  Used for content extraction in
// sparse_gcd_fp recursive Brown algorithm.
static std::optional<BSparsePoly> exact_div_fp(
        const BSparsePoly& A, const BSparsePoly& B, const BigInt& p) {
    if (B.empty()) return std::nullopt;
    if (A.empty()) return BSparsePoly{};
    std::size_t n = 0;
    for (const auto& [m, _] : A) n = std::max(n, m.size());
    for (const auto& [m, _] : B) n = std::max(n, m.size());
    BSparsePoly rem = A;
    BSparsePoly quo;
    auto [dlm, dlc] = *std::prev(B.end());
    auto inv = numtheory::modular_inverse(pos_mod(dlc, p), p);
    if (inv.is_error()) return std::nullopt;
    const BigInt ilc = inv.value();
    const std::size_t budget = (rem.size() + 1U) * (B.size() + 1U) + 16U;
    std::size_t steps = 0;
    while (!rem.empty()) {
        if (++steps > budget) return std::nullopt;
        auto [rlm, rlc] = *std::prev(rem.end());
        BMonomial qm(n, 0U);
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            if (ri < di) return std::nullopt;
            qm[i] = static_cast<unsigned int>(ri - di);
        }
        BigInt qc = pos_mod(rlc * ilc, p);
        quo[qm] = qc;
        for (const auto& [bm, bc] : B) {
            BMonomial nm(n, 0U);
            for (std::size_t i = 0; i < n; ++i)
                nm[i] = qm[i] + ((i < bm.size()) ? bm[i] : 0U);
            BigInt v = pos_mod(rem[nm] - pos_mod(qc * bc, p), p);
            if (v.is_zero()) rem.erase(nm);
            else rem[nm] = v;
        }
    }
    return quo;
}

}  // namespace

std::optional<BSparsePoly> sparse_gcd_fp(
    const BSparsePoly& A, const BSparsePoly& B,
    const std::vector<std::size_t>& active_vars,
    const BigInt& p, symbolic::CASContext& ctx, std::size_t depth) {
    if (A.empty()) return B;
    if (B.empty()) return A;
    if (depth > ctx.max_gcd_recursion_depth() + active_vars.size())
        return std::nullopt;
    if (active_vars.size() <= 1U) {
        if (active_vars.empty()) {
            BSparsePoly r; r[BMonomial{}] = BigInt(1); return r;
        }
        return univariate_sparse_gcd_fp(A, B, active_vars.front(), p);
    }
    const std::size_t eval_idx = active_vars.back();
    const std::vector<std::size_t> sub_active(active_vars.begin(), active_vars.end() - 1U);

    // ── Trailing eval_idx-power extraction (common monomial factor) ─────────
    // When both A and B are divisible by eval_idx^k (in particular when the
    // true gcd contains a factor eval_idx^k itself), the Brown eval-interp
    // loop CANNOT recover that factor: evaluating at eval_idx=0 yields 0
    // (skipped), and evaluating at eval_idx=b≠0 reduces eval_idx^k to b^k,
    // a unit in Fp; after monic normalization the eval_idx^k factor is lost.
    //
    // Fix: extract v_A = min eval_idx-exp across A's monomials (resp. v_B),
    // shift both polys down by m = min(v_A, v_B), and at the end multiply
    // pp_gcd by eval_idx^m to restore the common monomial factor.
    auto min_var_exp = [&](const BSparsePoly& X) -> std::size_t {
        std::size_t m = std::numeric_limits<std::size_t>::max();
        for (const auto& [mono, _] : X) {
            std::size_t e = (eval_idx < mono.size()) ? mono[eval_idx] : 0U;
            if (e < m) m = e;
            if (m == 0U) break;
        }
        return (m == std::numeric_limits<std::size_t>::max()) ? 0U : m;
    };
    const std::size_t vA = min_var_exp(A);
    const std::size_t vB = min_var_exp(B);
    const std::size_t v_common = std::min(vA, vB);
    BSparsePoly A_local;
    BSparsePoly B_local;
    auto shift_down = [&](const BSparsePoly& X, std::size_t k) {
        BSparsePoly out;
        for (const auto& [mono, c] : X) {
            BMonomial nm = mono;
            if (eval_idx < nm.size()) {
                nm[eval_idx] = static_cast<unsigned int>(nm[eval_idx] - k);
            }
            out[nm] = c;
        }
        return out;
    };
    const BSparsePoly& A_eff = (v_common > 0U)
        ? (A_local = shift_down(A, v_common), A_local) : A;
    const BSparsePoly& B_eff = (v_common > 0U)
        ? (B_local = shift_down(B, v_common), B_local) : B;
    auto multiply_by_eval_power = [&](BSparsePoly X) -> BSparsePoly {
        if (v_common == 0U) return X;
        BSparsePoly out;
        for (const auto& [m, c] : X) {
            BMonomial nm = m;
            if (eval_idx >= nm.size()) nm.resize(eval_idx + 1U, 0U);
            nm[eval_idx] = static_cast<unsigned int>(
                nm[eval_idx] + v_common);
            out[nm] = c;
        }
        return out;
    };

    // ── Content / primitive-part extraction (GCL §7.2 Algorithm 7.1) ────────
    // Decompose A = cont_A(sub_vars) · pp_A(eval_idx, sub_vars), similarly B.
    // Then gcd(A,B) = gcd(cont_A, cont_B) · gcd(pp_A, pp_B).
    // Without this step, when the true gcd has a polynomial leading coefficient
    // in a sub-variable (e.g. lc_z(gcd) = xy as a function of y when we recurse
    // on z), each per-sample monic normalization during Lagrange interpolation
    // collapses the y-dependent factor to a scalar and the polynomial structure
    // is lost.  This is the resolution of F3.1-BROWN-FP-RECURSIVE-NONMONIC.
    //
    // To preserve correctness when content extraction unexpectedly cancels real
    // structure (e.g. for fully-primitive inputs), we use content extraction
    // only when at least one of contA, contB is NON-trivial (i.e. the input is
    // imprimitive in eval_idx).  For primitive inputs we fall through to the
    // classical evaluation/Lagrange path that the original implementation used.
    //
    // Reference: Geddes-Czapor-Labahn "Algorithms for Computer Algebra"
    // §7.2 (Primitive PRS) + §7.4 (modular univariate over Fp); SymPy
    // `sympy.polys.modulargcd._modgcd_p` content step
    // (https://github.com/sympy/sympy/blob/master/sympy/polys/modulargcd.py).
    auto is_trivial_const1 = [](const BSparsePoly& sp) {
        if (sp.size() != 1U) return false;
        const auto& [m, c] = *sp.begin();
        return std::all_of(m.begin(), m.end(), [](unsigned int e){ return e == 0U; })
               && c == BigInt(1);
    };
    auto compute_content = [&](const BSparsePoly& X) -> std::optional<BSparsePoly> {
        auto layers = layers_by_var(X, eval_idx);
        if (layers.size() <= 1U) {
            // Empty polynomial: return 1.
            if (layers.empty()) { BSparsePoly one; one[BMonomial{}] = BigInt(1); return one; }
            // Single-layer (eval_idx-free) polynomial: the content IS the polynomial
            // itself (gcd of a single-element set = that element).  Return the actual
            // layer so that the cont_g step can detect a shared factor between an
            // eval_idx-free A and the eval_idx-content of B.
            // NOTE: returning trivial 1 here was an old short-cut that broke the case
            //   gcd(x(w-x), x(wy-xy+1)) w.r.t. y → should find factor x but got 1.
            return layers.begin()->second;
        }
        BSparsePoly cont;
        bool first = true;
        for (const auto& [k, layer] : layers) {
            if (layer.empty()) continue;
            if (first) { cont = layer; first = false; }
            else {
                auto g = sparse_gcd_fp(cont, layer, sub_active, p, ctx, depth + 1U);
                if (!g.has_value()) return std::nullopt;
                cont = std::move(*g);
                if (is_trivial_const1(cont)) return cont;  // early exit.
            }
        }
        return cont;
    };
    auto contA_opt = compute_content(A_eff);
    auto contB_opt = compute_content(B_eff);
    // Content extraction is OPTIONAL — if it fails (sub-recursion couldn't
    // resolve a content gcd) we fall back to the classical primitive-input
    // evaluation path.  This preserves original behavior for primitive inputs
    // while gaining the F3.1-BROWN-FP-RECURSIVE-NONMONIC fix for imprimitive
    // inputs.
    BSparsePoly contA = contA_opt.value_or(BSparsePoly{});
    if (contA.empty()) { contA[BMonomial{}] = BigInt(1); }
    BSparsePoly contB = contB_opt.value_or(BSparsePoly{});
    if (contB.empty()) { contB[BMonomial{}] = BigInt(1); }
    const bool both_primitive = is_trivial_const1(contA) && is_trivial_const1(contB);

    BSparsePoly cont_g;
    BSparsePoly ppA = A_eff, ppB = B_eff;
    if (!both_primitive) {
        auto cont_g_opt = sparse_gcd_fp(contA, contB, sub_active, p, ctx, depth + 1U);
        if (!cont_g_opt.has_value()) {
            // Cont-gcd recursion failed: fall back to primitive treatment.
            cont_g.clear(); cont_g[BMonomial{}] = BigInt(1);
        } else {
            cont_g = std::move(*cont_g_opt);
            auto strip_content = [&](const BSparsePoly& X, const BSparsePoly& cont)
                    -> std::optional<BSparsePoly> {
                if (cont.empty() || is_trivial_const1(cont)) return X;
                auto layers = layers_by_var(X, eval_idx);
                std::map<std::size_t, BSparsePoly> pp_layers;
                for (const auto& [k, layer] : layers) {
                    if (layer.empty()) { pp_layers[k] = {}; continue; }
                    auto q = exact_div_fp(layer, cont, p);
                    if (!q.has_value()) return std::nullopt;
                    pp_layers[k] = std::move(*q);
                }
                return reassemble_layers(pp_layers, eval_idx, p);
            };
            auto ppA_opt = strip_content(A_eff, contA);
            auto ppB_opt = strip_content(B_eff, contB);
            if (!ppA_opt.has_value() || !ppB_opt.has_value()) {
                // Strip failed: fall back to primitive treatment.
                ppA = A_eff; ppB = B_eff;
                cont_g.clear(); cont_g[BMonomial{}] = BigInt(1);
            } else {
                ppA = std::move(*ppA_opt);
                ppB = std::move(*ppB_opt);
            }
        }
    } else {
        cont_g.clear();
        cont_g[BMonomial{}] = BigInt(1);
    }

    // Divisibility fast-path: if ppB | ppA → gcd = ppB; if ppA | ppB → gcd = ppA.
    // Necessary for two failure modes of the Brown eval-interp loop:
    //   (a) ppA == ppB: eval gives equal scalars at every point; gcd normalises to 1.
    //   (b) ppB = eval_var (e.g. y): eval at b=0 is empty (skipped); eval at b≠0
    //       gives a unit; gcd = 1 at every usable point; interpolation gives 1.
    // Both are covered by the divisibility check (equality ⊂ divisibility).
    if (auto q = exact_div_fp(ppB, ppA, p); q.has_value()) {
        BSparsePoly fp = is_trivial_const1(cont_g) ? ppA : mul_mod_p(ppA, cont_g, p);
        return multiply_by_eval_power(std::move(fp));
    }
    if (auto q = exact_div_fp(ppA, ppB, p); q.has_value()) {
        BSparsePoly fp = is_trivial_const1(cont_g) ? ppB : mul_mod_p(ppB, cont_g, p);
        return multiply_by_eval_power(std::move(fp));
    }

    std::size_t dA = deg_in_var(ppA, eval_idx);
    std::size_t dB = deg_in_var(ppB, eval_idx);
    std::size_t deg_bound = std::min(dA, dB);

    // Three subcases when deg_bound == 0:
    //   (a) both dA == 0 AND dB == 0 → neither depends on eval_idx; pp_gcd =
    //       gcd over sub_active.
    //   (b) dA == 0 and dB > 0 (or symmetric) → ppA is eval_idx-free; gcd(P,Q)
    //       divides P so it's also eval_idx-free; reduce to gcd(P, cont_eval(Q))
    //       where cont_eval(Q) = gcd of Q's eval_idx-layer coefficients.
    // Without this, the Brown eval-interp loop fails (samples of the lower-deg
    // poly hit zero coefficients or destroy a leading factor) and Lagrange
    // recovers only the unit gcd, losing real common factors.
    if (deg_bound == 0U) {
        // Compute the gcd-of-layers of the non-trivial side, then gcd with
        // the eval_idx-free side over sub_active.  If both are eval_idx-free,
        // this reduces to gcd(ppA, ppB) over sub_active directly.
        auto layers_of = [&](const BSparsePoly& X) {
            std::map<std::size_t, BSparsePoly> out;
            for (const auto& [m, c] : X) {
                std::size_t k = (eval_idx < m.size()) ? m[eval_idx] : 0U;
                BMonomial m2 = m;
                if (eval_idx < m2.size()) m2[eval_idx] = 0U;
                out[k][m2] = c;
            }
            return out;
        };
        auto layer_content = [&](const BSparsePoly& X) -> std::optional<BSparsePoly> {
            auto layers = layers_of(X);
            if (layers.empty()) {
                BSparsePoly r; r[BMonomial{}] = BigInt(1); return r;
            }
            if (layers.size() == 1U) return layers.begin()->second;
            BSparsePoly acc;
            bool first = true;
            for (auto& [k, layer] : layers) {
                if (layer.empty()) continue;
                if (first) { acc = layer; first = false; continue; }
                auto g = sparse_gcd_fp(acc, layer, sub_active, p, ctx, depth + 1U);
                if (!g.has_value()) return std::nullopt;
                acc = std::move(*g);
                if (is_trivial_const1(acc)) return acc;
            }
            return acc;
        };
        BSparsePoly Aside = (dA == 0U) ? ppA : layer_content(ppA).value_or(BSparsePoly{});
        BSparsePoly Bside = (dB == 0U) ? ppB : layer_content(ppB).value_or(BSparsePoly{});
        if (Aside.empty()) { Aside[BMonomial{}] = BigInt(1); }
        if (Bside.empty()) { Bside[BMonomial{}] = BigInt(1); }
        auto sub_gcd_opt = sparse_gcd_fp(Aside, Bside, sub_active, p, ctx, depth + 1U);
        BSparsePoly sub_gcd;
        if (sub_gcd_opt.has_value()) {
            sub_gcd = std::move(*sub_gcd_opt);
        } else {
            sub_gcd[BMonomial{}] = BigInt(1);
        }
        BSparsePoly combined;
        if (is_trivial_const1(cont_g)) combined = std::move(sub_gcd);
        else if (is_trivial_const1(sub_gcd)) combined = cont_g;
        else combined = mul_mod_p(sub_gcd, cont_g, p);
        return multiply_by_eval_power(std::move(combined));
    }

    // ── Original Brown evaluation / Lagrange interp loop ────────────────────
    std::vector<BigInt> bs;
    std::vector<BSparsePoly> gvals;
    std::optional<std::vector<std::size_t>> stable_subdegs;
    const std::size_t max_tries = (deg_bound + 1U) * 4U + 8U;
    for (std::size_t tries = 0; bs.size() < deg_bound + 1U && tries < max_tries; ++tries) {
        BigInt b(static_cast<long long>(tries));
        if (b >= p) break;
        BSparsePoly Ae = eval_var_mod_p(ppA, eval_idx, b, p);
        BSparsePoly Be = eval_var_mod_p(ppB, eval_idx, b, p);
        if (Ae.empty() || Be.empty()) continue;
        auto g_opt = sparse_gcd_fp(Ae, Be, sub_active, p, ctx, depth + 1U);
        if (!g_opt.has_value()) continue;
        BSparsePoly ge = std::move(*g_opt);
        std::vector<std::size_t> sd;
        sd.reserve(sub_active.size());
        for (auto vi : sub_active) sd.push_back(deg_in_var(ge, vi));
        if (!stable_subdegs.has_value()) {
            stable_subdegs = sd;
        } else {
            bool lower = false, higher = false;
            for (std::size_t k = 0; k < sd.size(); ++k) {
                if (sd[k] < (*stable_subdegs)[k]) { lower = true; }
                if (sd[k] > (*stable_subdegs)[k]) { higher = true; }
            }
            if (higher && !lower) { stable_subdegs = sd; bs.clear(); gvals.clear(); }
            else if (lower) continue;
            if (lower && higher) continue;
        }
        bs.push_back(b);
        gvals.push_back(std::move(ge));
    }
    if (bs.size() < deg_bound + 1U) return std::nullopt;
    auto pp_gcd_opt = lagrange_interp_fp(bs, gvals, eval_idx, deg_bound, p);
    if (!pp_gcd_opt.has_value()) return std::nullopt;
    BSparsePoly pp_gcd = std::move(*pp_gcd_opt);
    BSparsePoly final_result = is_trivial_const1(cont_g)
        ? std::move(pp_gcd)
        : mul_mod_p(pp_gcd, cont_g, p);
    return multiply_by_eval_power(std::move(final_result));
}

}  // namespace cas::algebra::fp_helpers
