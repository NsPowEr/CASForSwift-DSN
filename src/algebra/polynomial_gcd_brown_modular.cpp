// polynomial_gcd_brown_modular.cpp — REAL Brown's modular multivariate GCD (F3.1).
//
// Algorithm (Geddes-Czapor-Labahn §7.4–7.5):
//   for each prime p (skipped if p | lc(P) or p | lc(Q) in main variable):
//       g_p = gcd in Fp[x1,...,xn] via recursive evaluation/interpolation IN Fp
//             (delegated to polynomial_gcd_fp_recursive.cpp).
//       degree-stability: discard old gcd data on a smaller degree (unlucky prime)
//       CRT-combine coefficients of g_p across primes
//   stop when ∏p > 2·Mignotte_bound; reconstruct via centered representation;
//   certify divisibility g | P and g | Q in Z[x1,...,xn]; if certificate fails,
//   add more primes; if budget exhausted → explicit Unimplemented diagnostic.
//
// Reference: Geddes-Czapor-Labahn "Algorithms for Computer Algebra" §7.4 (modular GCD),
// §7.5 (recursion in variables); Brown (1971) "On Euclid's Algorithm and the
// Computation of Polynomial Greatest Common Divisors".

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include "polynomial_gcd_fp_internal.hpp"

#include <algorithm>
#include <cstddef>
#include <set>
#include <string>
#include <vector>

namespace cas::algebra {

using fp_helpers::BMonomial;
using fp_helpers::BSparsePoly;
using fp_helpers::pos_mod;
using fp_helpers::bigint_gcd;
using fp_helpers::centered_repr;
using fp_helpers::reduce_sparse_mod_p;
using fp_helpers::sparse_inf_norm;
using fp_helpers::deg_in_var;
using fp_helpers::sparse_gcd_fp;

namespace {

static BSparsePoly to_sparse(const MultivariatePolynomial& p,
                              const std::vector<Symbol>& vars) {
    BSparsePoly sp;
    for (const auto& term : p.terms()) {
        BMonomial mono(vars.size(), 0U);
        for (const auto& [sym, exp] : term.factors)
            for (std::size_t i = 0; i < vars.size(); ++i)
                if (vars[i].name == sym.name) mono[i] = exp;
        sp[mono] += term.coefficient;
        if (sp[mono].is_zero()) sp.erase(mono);
    }
    return sp;
}

static MultivariatePolynomial from_sparse(const BSparsePoly& sp,
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

static std::vector<Symbol> collect_vars(const MultivariatePolynomial& p,
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

// Multivariate Mignotte bound: 2^{min total deg} * min(||P||_∞, ||Q||_∞).
static BigInt multivar_mignotte_bound(const BSparsePoly& P, const BSparsePoly& Q) {
    std::size_t dP = 0, dQ = 0;
    for (const auto& [m, _] : P) {
        std::size_t td = 0; for (auto e : m) td += e;
        dP = std::max(dP, td);
    }
    for (const auto& [m, _] : Q) {
        std::size_t td = 0; for (auto e : m) td += e;
        dQ = std::max(dQ, td);
    }
    BigInt nP = sparse_inf_norm(P);
    BigInt nQ = sparse_inf_norm(Q);
    BigInt mn = (nP < nQ) ? nP : nQ;
    if (mn.is_zero()) mn = BigInt(1);
    return BigInt(1).shift_left_bits(std::min(dP, dQ)) * mn;
}

// Exact divisibility test in Z[x1,...,xn] via sparse lex pseudo-division.
[[nodiscard]] static bool divides_sparse_z(
    const BSparsePoly& dividend, const BSparsePoly& divisor,
    std::size_t n_vars) {
    if (divisor.empty()) return false;
    if (dividend.empty()) return true;
    BSparsePoly rem = dividend;
    auto [dlm, dlc] = *std::prev(divisor.end());
    const std::size_t budget = (rem.size() + 1U) * (divisor.size() + 1U) + 16U;
    std::size_t steps = 0;
    while (!rem.empty()) {
        if (++steps > budget) return false;
        auto [rlm, rlc] = *std::prev(rem.end());
        for (std::size_t i = 0; i < n_vars; ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            if (ri < di) return false;
        }
        if ((rlc % dlc) != BigInt(0)) return false;
        BigInt qc = rlc / dlc;
        BMonomial qm(n_vars, 0U);
        for (std::size_t i = 0; i < n_vars; ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            qm[i] = static_cast<unsigned int>(ri - di);
        }
        for (const auto& [dm, dc] : divisor) {
            BMonomial nm(n_vars, 0U);
            for (std::size_t i = 0; i < n_vars; ++i)
                nm[i] = qm[i] + ((i < dm.size()) ? dm[i] : 0U);
            rem[nm] -= qc * dc;
            if (rem[nm].is_zero()) rem.erase(nm);
        }
    }
    return true;
}

}  // namespace

[[nodiscard]] Result<MultivariatePolynomial> gcd_brown_modular(
    const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx, std::vector<BigInt>* out_primes_used);

[[nodiscard]] Result<MultivariatePolynomial> gcd_brown_modular(
    const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx, std::vector<BigInt>* out_primes_used) {
    if (P.is_zero()) return ok(Q);
    if (Q.is_zero()) return ok(P);
    const std::vector<Symbol> vars = collect_vars(P, Q);
    if (vars.empty()) {
        BigInt g(0);
        for (const auto& t : P.terms()) g = bigint_gcd(g, t.coefficient.abs());
        for (const auto& t : Q.terms()) g = bigint_gcd(g, t.coefficient.abs());
        return ok(MultivariatePolynomial({{ .coefficient = g, .factors = {} }}));
    }

    // Early-out: if either input is an integer constant (no variables in any
    // term), the gcd reduces to gcd(P, Q) of the content scalars times Q
    // (resp. P).  This avoids feeding a constant into sparse_gcd_fp where the
    // multivariate-gcd-mod-p machinery can stall.
    auto poly_content_if_constant = [](const MultivariatePolynomial& R)
            -> std::optional<BigInt> {
        BigInt c(0);
        for (const auto& t : R.terms()) {
            if (!t.factors.empty()) return std::nullopt;
            c = bigint_gcd(c, t.coefficient.abs());
        }
        return c;
    };
    if (auto cP = poly_content_if_constant(P); cP.has_value()) {
        BigInt cQ(0);
        for (const auto& t : Q.terms()) cQ = bigint_gcd(cQ, t.coefficient.abs());
        BigInt g = bigint_gcd(cP.value(), cQ);
        if (g.is_zero()) g = BigInt(1);
        return ok(MultivariatePolynomial({{ .coefficient = g, .factors = {} }}));
    }
    if (auto cQ = poly_content_if_constant(Q); cQ.has_value()) {
        BigInt cP(0);
        for (const auto& t : P.terms()) cP = bigint_gcd(cP, t.coefficient.abs());
        BigInt g = bigint_gcd(cQ.value(), cP);
        if (g.is_zero()) g = BigInt(1);
        return ok(MultivariatePolynomial({{ .coefficient = g, .factors = {} }}));
    }

    const std::size_t n = vars.size();
    BSparsePoly spP = to_sparse(P, vars);
    BSparsePoly spQ = to_sparse(Q, vars);

    BigInt cont(0);
    for (const auto& [_, c] : spP) cont = bigint_gcd(cont, c.abs());
    BigInt contQ(0);
    for (const auto& [_, c] : spQ) contQ = bigint_gcd(contQ, c.abs());
    BigInt cont_g = bigint_gcd(cont, contQ);
    if (cont_g.is_zero()) cont_g = BigInt(1);

    auto make_prim = [](const BSparsePoly& sp, const BigInt& c) {
        BSparsePoly r;
        for (const auto& [m, v] : sp) r[m] = v / c;
        return r;
    };
    BSparsePoly ppP = make_prim(spP, cont);
    BSparsePoly ppQ = make_prim(spQ, contQ);

    const std::size_t main_var = n - 1U;

    // ─ Polynomial-content pre-extraction w.r.t. main_var (Geddes §7.4 setup) ─
    // Decompose ppP = cont_main(ppP) · pp_main(ppP), similarly ppQ, where
    // cont_main is the (n-1)-variable gcd of the main_var-layer coefficients.
    // Then gcd(ppP, ppQ) = gcd(cont_main_P, cont_main_Q) · gcd(pp_main_P, pp_main_Q).
    // The modular machinery below handles only the primitive part w.r.t. main_var;
    // without this pre-extraction, content factors common to every main_var
    // layer (e.g. P = x·(yz+1)(w+...), gcd_content_z = x) get stripped by
    // remove_spurious_main_var_factor leading to under-divided results.
    auto layers_main = [&](const BSparsePoly& sp) {
        std::map<std::size_t, BSparsePoly> out;
        for (const auto& [m, c] : sp) {
            std::size_t k = (main_var < m.size()) ? m[main_var] : 0U;
            BMonomial m2 = m;
            if (main_var < m2.size()) m2[main_var] = 0U;
            out[k][m2] = c;
        }
        return out;
    };
    auto content_main_z = [&](const BSparsePoly& sp) -> MultivariatePolynomial {
        if (n < 2U) return MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}});
        auto layers = layers_main(sp);
        if (layers.size() <= 1U) {
            return MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}});
        }
        MultivariatePolynomial acc;
        bool first = true;
        for (auto& [k, layer] : layers) {
            if (layer.empty()) continue;
            MultivariatePolynomial cur = sub_sparse_to_mv(layer, vars, main_var);
            if (first) { acc = std::move(cur); first = false; continue; }
            auto rec = gcd_brown_modular(acc, cur, ctx);
            if (rec.is_error()) {
                return MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}});
            }
            acc = std::move(rec.value());
            if (acc.terms().size() == 1U && acc.terms()[0].factors.empty()
                && acc.terms()[0].coefficient == BigInt(1)) {
                return acc;
            }
        }
        return acc;
    };

    MultivariatePolynomial cont_main_P_mv, cont_main_Q_mv;
    MultivariatePolynomial cont_main_gcd_mv =
        MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}});
    bool have_main_content_gcd = false;
    if (n >= 2U) {
        cont_main_P_mv = content_main_z(ppP);
        cont_main_Q_mv = content_main_z(ppQ);
        auto is_one = [](const MultivariatePolynomial& m){
            return m.terms().size() == 1U && m.terms()[0].factors.empty()
                && m.terms()[0].coefficient == BigInt(1);
        };
        if (!is_one(cont_main_P_mv) || !is_one(cont_main_Q_mv)) {
            auto cg = gcd_brown_modular(cont_main_P_mv, cont_main_Q_mv, ctx);
            if (cg.is_ok() && !is_one(cg.value())) {
                cont_main_gcd_mv = std::move(cg.value());
                BSparsePoly cont_sp = mv_to_sub_sparse(cont_main_gcd_mv, vars, main_var);
                if (!cont_sp.empty()) {
                    auto strip = [&](const BSparsePoly& X) -> std::optional<BSparsePoly> {
                        auto layers = layers_main(X);
                        std::map<std::size_t, BSparsePoly> out_layers;
                        for (const auto& [k, layer] : layers) {
                            if (layer.empty()) { out_layers[k] = {}; continue; }
                            BSparsePoly quo;
                            if (!exact_divide_sparse_z(layer, cont_sp, n, quo))
                                return std::nullopt;
                            out_layers[k] = std::move(quo);
                        }
                        BSparsePoly res;
                        for (const auto& [k, layer] : out_layers) {
                            for (const auto& [m, c] : layer) {
                                BMonomial nm = m;
                                if (nm.size() <= main_var) nm.resize(main_var + 1U, 0U);
                                nm[main_var] = static_cast<unsigned int>(k);
                                res[nm] = c;
                            }
                        }
                        return res;
                    };
                    auto sp_pp_P = strip(ppP);
                    auto sp_pp_Q = strip(ppQ);
                    if (sp_pp_P.has_value() && sp_pp_Q.has_value()) {
                        ppP = std::move(*sp_pp_P);
                        ppQ = std::move(*sp_pp_Q);
                        have_main_content_gcd = true;
                    }
                }
            }
        }
    }
    const std::size_t dP_main = deg_in_var(ppP, main_var);
    const std::size_t dQ_main = deg_in_var(ppQ, main_var);
    BigInt lcP(0), lcQ(0);
    for (const auto& [m, c] : ppP) {
        std::size_t k = (main_var < m.size()) ? m[main_var] : 0U;
        if (k == dP_main) lcP = bigint_gcd(lcP, c.abs());
    }
    for (const auto& [m, c] : ppQ) {
        std::size_t k = (main_var < m.size()) ? m[main_var] : 0U;
        if (k == dQ_main) lcQ = bigint_gcd(lcQ, c.abs());
    }

    // ─ Lc-poly-scaling (GCL §7.4.2 Algorithm 7.2) ───────────────────────────
    // Compute lc_bound_poly L = gcd(lc_main(ppP), lc_main(ppQ)) as polynomial
    // in n-1 sub-variables (Z-coefficient).  When n=1 (univariate input) L
    // collapses to a scalar already captured by lcP/lcQ above — skip recursion.
    // When ALL coefficients of L are constant w.r.t. sub-vars (i.e. L has at
    // most one term with empty monomial), L is a scalar and lc-scaling reduces
    // to the legacy scalar path.
    BSparsePoly lc_bound_poly;       // poly in n-1 sub-vars (n-var indexed, main_var=0)
    bool use_lc_scaling = false;
    if (n >= 2U) {
        BSparsePoly lcP_poly = extract_lc_in_var(ppP, main_var);
        BSparsePoly lcQ_poly = extract_lc_in_var(ppQ, main_var);
        // Heuristic: skip recursive call when both lcs are pure scalars
        // (single empty monomial) → existing scalar lcP/lcQ already handle it.
        auto is_scalar = [](const BSparsePoly& sp){
            if (sp.size() != 1U) return sp.empty();
            return std::all_of(sp.begin()->first.begin(), sp.begin()->first.end(),
                               [](unsigned int e){ return e == 0U; });
        };
        if (!is_scalar(lcP_poly) || !is_scalar(lcQ_poly)) {
            MultivariatePolynomial mP = sub_sparse_to_mv(lcP_poly, vars, main_var);
            MultivariatePolynomial mQ = sub_sparse_to_mv(lcQ_poly, vars, main_var);
            // Recursive gcd — dimension n-1 < n, well-founded.
            auto rec = gcd_brown_modular(mP, mQ, ctx);
            if (rec.is_ok()) {
                lc_bound_poly = mv_to_sub_sparse(rec.value(), vars, main_var);
                if (!lc_bound_poly.empty()) use_lc_scaling = true;
            }
            // If recursion fails, fall through to scalar path (no lc-scaling).
        }
    }

    const BigInt M_target = multivar_mignotte_bound(ppP, ppQ);
    BigInt M_need = (M_target + M_target) + BigInt(1);
    // When lc-scaling active, candidate lc-in-main grows by ||L||_∞ factor.
    if (use_lc_scaling) {
        BigInt L_norm = fp_helpers::sparse_inf_norm(lc_bound_poly);
        if (!L_norm.is_zero()) M_need = M_need * L_norm + L_norm;
    }

    BSparsePoly crt_solution;
    BigInt M_acc(1);
    std::optional<std::vector<std::size_t>> stable_subdegs;
    std::set<BMonomial> stable_monos;
    bool have_lucky = false;

    std::size_t h = lcP.bit_length() * 2654435761ULL;
    h ^= ppP.size() * 40503ULL ^ ppQ.size() * 12347ULL;
    long long start = 1073741827LL + static_cast<long long>(h % 65537ULL);
    auto np0 = numtheory::next_prime(BigInt(start - 1));
    BigInt cur_p = np0.is_ok() ? np0.value() : BigInt(1073741827LL);

    const std::size_t max_primes = std::min<std::size_t>(ctx.max_gcd_total_calls(), 128U);
    std::size_t primes_used = 0;
    std::size_t cert_fail_count = 0;
    std::size_t consecutive_useless = 0;
    const std::size_t useless_streak_cap = 16U;
    std::vector<BigInt> primes_log;

    std::vector<std::size_t> active;
    for (std::size_t i = 0; i < n; ++i) active.push_back(i);

    std::size_t prime_attempts = 0;
    const std::size_t max_attempts = max_primes * 4U;
    while (primes_used < max_primes && prime_attempts < max_attempts) {
        ++prime_attempts;
        if (consecutive_useless >= useless_streak_cap) {
            return make_unimplemented<MultivariatePolynomial>(
                "algebra", "gcd_brown_modular",
                "useless_primes_in_a_row=" + std::to_string(consecutive_useless),
                "GCD_BROWN_MODULAR_UNLUCKY_STREAK",
                "Repeated primes produce empty/unusable modular gcds — fallback path",
                "F3.1");
        }
        const BigInt p = cur_p;
        auto np_next = numtheory::next_prime(cur_p);
        if (np_next.is_ok()) cur_p = np_next.value(); else break;

        if (!lcP.is_zero() && (lcP % p).is_zero()) { ++consecutive_useless; continue; }
        if (!lcQ.is_zero() && (lcQ % p).is_zero()) { ++consecutive_useless; continue; }

        BSparsePoly Pp = reduce_sparse_mod_p(ppP, p);
        BSparsePoly Qp = reduce_sparse_mod_p(ppQ, p);
        if (Pp.empty() || Qp.empty()) { ++consecutive_useless; continue; }

        auto g_opt = sparse_gcd_fp(Pp, Qp, active, p, ctx, 0U);
        if (!g_opt.has_value()) {
            ++consecutive_useless; continue;
        }
        BSparsePoly gp = std::move(*g_opt);
        gp = reduce_sparse_mod_p(gp, p);
        if (gp.empty()) { ++consecutive_useless; continue; }

        // ─ Lc-poly-scaling per-prime step (GCL §7.4.2 Alg 7.2 line 6) ────────
        // sparse_gcd_fp returns gp monic in main_var → lc_main(gp) is a scalar
        // (typically 1, possibly a non-unit if the implementation differs).
        // To force lc_main(scaled_gp) ≡ Lp (so all primes share the SAME lc),
        // we multiply gp polynomially by (Lp / lcg) where lcg = lc_main(gp).
        //
        // Case A (common): lcg is a constant scalar c.  Then scale gp by Lp/c:
        //                  gp_new = (Lp as n-var poly) * (c^{-1} mod p) * gp.
        //                  Result: lc_main(gp_new) = Lp.
        // Case B: lcg is a true polynomial in n-1 sub-vars (sparse_gcd_fp did not
        //         normalize to monic).  We require Lp ≡ q * lcg (mod p) for some
        //         scalar q ∈ Fp*; otherwise this prime is unlucky.
        if (use_lc_scaling) {
            BSparsePoly Lp = reduce_lc_bound_mod_p(lc_bound_poly, p);
            if (Lp.empty()) { ++consecutive_useless; continue; }
            BSparsePoly lcg = extract_lc_in_var(gp, main_var);
            if (lcg.empty()) { ++consecutive_useless; continue; }

            // Detect whether lcg is a pure scalar (single empty-monomial term).
            auto is_pure_scalar = [](const BSparsePoly& sp) -> std::optional<BigInt> {
                if (sp.size() != 1U) return std::nullopt;
                auto it = sp.begin();
                bool all_zero = std::all_of(it->first.begin(), it->first.end(),
                                            [](unsigned int e){ return e == 0U; });
                if (!all_zero) return std::nullopt;
                return it->second;
            };

            BigInt q_inv;  // multiplicative scaling factor (∈ Fp*).
            BSparsePoly Lp_eff;  // polynomial we must multiply gp by.
            if (auto cs = is_pure_scalar(lcg); cs.has_value()) {
                // Case A: gp monic-up-to-scalar.  Multiply by c^{-1} * Lp.
                auto inv_c = numtheory::modular_inverse(pos_mod(cs.value(), p), p);
                if (inv_c.is_error()) { ++consecutive_useless; continue; }
                q_inv = inv_c.value();
                Lp_eff = Lp;
            } else {
                // Case B: lcg is a true polynomial.  Require Lp = q · lcg (mod p).
                // (Otherwise the prime is unlucky for lc-scaling — content
                // structure of Lp and lcg disagrees.)
                auto u_opt = compute_lc_scalar_ratio(Lp, lcg, p);
                if (!u_opt.has_value()) { ++consecutive_useless; continue; }
                // We have u s.t. u·lcg = Lp.  So scaling by u (scalar) already
                // makes lc_main(u·gp) = Lp — pure scalar multiplication suffices.
                scale_by_lc(gp, u_opt.value(), p);
                consecutive_useless = 0;
                goto lc_scaling_done;
            }
            // Case A path: scale gp by q_inv (scalar), then polynomial-mult by Lp.
            scale_by_lc(gp, q_inv, p);
            gp = multiply_sparse_mod_p(Lp_eff, gp, p, n);
            if (gp.empty()) { fprintf(stderr, "[DBG-USELESS] gp*Lp empty p=%s\n", p.decimal().c_str()); ++consecutive_useless; continue; }
        }
        lc_scaling_done:;
        consecutive_useless = 0;

        std::vector<std::size_t> sd(n, 0U);
        for (std::size_t i = 0; i < n; ++i) sd[i] = deg_in_var(gp, i);
        if (have_lucky) {
            bool higher = false, lower = false;
            for (std::size_t i = 0; i < n; ++i) {
                if (sd[i] < (*stable_subdegs)[i]) { lower = true; break; }
                if (sd[i] > (*stable_subdegs)[i]) { higher = true; break; }
            }
            if (higher) continue;
            if (lower) {
                stable_subdegs = sd;
                stable_monos.clear();
                crt_solution.clear();
                M_acc = BigInt(1);
                primes_used = 0;
            }
        } else {
            stable_subdegs = sd;
            have_lucky = true;
        }
        for (const auto& [m, _] : gp) stable_monos.insert(m);

        // CRT-merge per monomial.
        const BigInt M_before = M_acc;
        bool merge_failed = false;
        for (const auto& mono : stable_monos) {
            BigInt ri(0);
            auto it = gp.find(mono);
            if (it != gp.end()) ri = it->second;
            BigInt cur = crt_solution.count(mono) ? crt_solution[mono] : BigInt(0);
            BigInt cur_mod_p = pos_mod(cur, p);
            BigInt delta = pos_mod(ri - cur_mod_p, p);
            auto inv_m = numtheory::modular_inverse(pos_mod(M_before, p), p);
            if (inv_m.is_error()) { merge_failed = true; break; }
            BigInt t = pos_mod(delta * inv_m.value(), p);
            BigInt new_val = cur + M_before * t;
            if (new_val.is_zero()) crt_solution.erase(mono);
            else crt_solution[mono] = new_val;
        }
        if (merge_failed) continue;
        M_acc = M_before * p;
        ++primes_used;
        primes_log.push_back(p);

        if (M_acc > M_need) {
            BSparsePoly cand;
            for (const auto& mono : stable_monos) {
                BigInt v = crt_solution.count(mono) ? crt_solution[mono] : BigInt(0);
                BigInt cv = centered_repr(v, M_acc);
                if (!cv.is_zero()) cand[mono] = cv;
            }
            // Helper: finalize and return cand (apply cont_g + sign normalization).
            // When main-var polynomial content was pre-extracted, the candidate is
            // pp_main(gcd); we multiply by cont_main_gcd_mv to recover full gcd.
            auto try_certify_and_return = [&](BSparsePoly trial)
                    -> std::optional<MultivariatePolynomial> {
                if (trial.empty()) return std::nullopt;
                // Strip integer content (a scalar BigInt) — does NOT touch poly lc.
                BigInt cc(0);
                for (const auto& [_, c] : trial) cc = bigint_gcd(cc, c.abs());
                if (cc.is_zero()) cc = BigInt(1);
                for (auto& [_, c] : trial) c = c / cc;
                if (!divides_sparse_z(ppP, trial, n)) return std::nullopt;
                if (!divides_sparse_z(ppQ, trial, n)) return std::nullopt;
                for (auto& [_, c] : trial) c *= cont_g;
                if (!trial.empty()) {
                    auto last = std::prev(trial.end());
                    if (last->second.is_negative())
                        for (auto& [_, c] : trial) c = -c;
                }
                MultivariatePolynomial result = from_sparse(trial, vars);
                if (have_main_content_gcd) result = result * cont_main_gcd_mv;
                return result;
            };

            // Attempt 1: certify cand directly (works when L = lc_main(true_gcd)).
            if (auto r = try_certify_and_return(cand); r.has_value()) {
                if (out_primes_used) *out_primes_used = primes_log;
                return ok(*r);
            }

            // Attempt 2 (GCL §7.4.2 spurious-factor removal): when lc-scaling is
            // active, cand has lc_main = L but the true gcd's lc_main divides L.
            // The quotient h = L / lc_main(true_gcd) is a polynomial in the n-1
            // sub-vars that appears as a *uniform* factor in every main-var
            // layer of cand.  Computing the main-var content (= gcd of layers)
            // and dividing it out recovers the true gcd — UNLESS the true gcd
            // itself has non-trivial main-var content; in that case the divided
            // candidate will fail certification and we fall through to the
            // M_need-doubling path.  Both outcomes are safe (the cert is the
            // final arbiter).
            // Skip attempt 2 when main-var content was already pre-extracted —
            // ppP/ppQ are pp_main, so cand should NOT carry a spurious main-var
            // content factor; stripping would over-divide (e.g. cand = yz+1
            // could be reduced to scalar 1 which would falsely certify).
            if (use_lc_scaling && !have_main_content_gcd) {
                BSparsePoly cand_clean;
                if (remove_spurious_main_var_factor(
                        cand, vars, main_var, ctx, cand_clean)) {
                    if (auto r2 = try_certify_and_return(cand_clean);
                        r2.has_value()) {
                        if (out_primes_used) *out_primes_used = primes_log;
                        return ok(*r2);
                    }
                }
            }

            M_need = M_need + M_need;
            if (++cert_fail_count >= 4U) {
                return make_unimplemented<MultivariatePolynomial>(
                    "algebra", "gcd_brown_modular",
                    "cert_failures=" + std::to_string(cert_fail_count) +
                        ",lc_scaling=" + std::string(use_lc_scaling ? "on" : "off"),
                    "GCD_BROWN_MODULAR_CERT_REPEATEDLY_FAILED",
                    "Modular candidate not divisible — even with lc-poly-scaling and content removal",
                    "F3.1");
            }
        }
    }

    return make_unimplemented<MultivariatePolynomial>(
        "algebra", "gcd_brown_modular",
        "primes_used=" + std::to_string(primes_used) + ",budget=" + std::to_string(max_primes),
        "GCD_BROWN_MODULAR_PRIME_BUDGET_EXHAUSTED",
        "Increase ctx.max_gcd_total_calls() or input is adversarial for current bound",
        "F3.1");
}

[[nodiscard]] Result<MultivariatePolynomial> gcd_brown_modular(
    const MultivariatePolynomial& P, const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx) {
    return gcd_brown_modular(P, Q, ctx, nullptr);
}

}  // namespace cas::algebra
