// polynomial_gcd_fp_recursive.cpp — Helpers + Brown-in-Fp recursive gcd.
// Implements the public API in polynomial_gcd_fp_internal.hpp.

#include "polynomial_gcd_fp_internal.hpp"
#include "cas/numtheory.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <optional>
#include <vector>

namespace cas::algebra::fp_helpers {

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
    auto is_trivial_const1 = [](const BSparsePoly& sp) {
        if (sp.size() != 1U) return false;
        const auto& [m, c] = *sp.begin();
        return std::all_of(m.begin(), m.end(), [](unsigned int e){ return e == 0U; })
               && c == BigInt(1);
    };
    auto compute_content = [&](const BSparsePoly& X) -> std::optional<BSparsePoly> {
        auto layers = layers_by_var(X, eval_idx);
        if (layers.size() <= 1U) {
            if (layers.empty()) { BSparsePoly one; one[BMonomial{}] = BigInt(1); return one; }
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
                if (is_trivial_const1(cont)) return cont;
            }
        }
        return cont;
    };
    auto contA_opt = compute_content(A_eff);
    auto contB_opt = compute_content(B_eff);
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

    // Divisibility fast-path
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

    if (deg_bound == 0U) {
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
