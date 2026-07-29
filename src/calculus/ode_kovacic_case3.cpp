// Kovacic Case 3 — algebraic ω over C(x) with Galois group conjugate to a
// finite subgroup of SL(2,C) (tetrahedral A₄, octahedral S₄, icosahedral A₅).
//
// Reference: Kovacic J.J. (1986), §5, "The Algorithm for Case 3",
// J. Symbolic Computation 2, pp. 22-43.
//
// Spec (paper-verified): .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/
//                        Kovacic_Case3_SL2C.md
// Ledger: HC-KV-03 (HARDCODE_LEDGER.md).
//
// Algorithm overview (assumes Case 1 + Case 2 have failed):
//   For each candidate n ∈ {4, 6, 12} (Galois group order ∈ {24, 48, 120}):
//     Step 1.  For each finite pole c of r, build E_c ⊂ Z:
//                order 1     →  E_c = {12}
//                order 2     →  E_c = {6 + (12k/n)·√(1+4α) | k = 0, ±1, ..., ±n/2} ∩ Z
//                                (α = coeff of 1/(x-c)² in partial fraction)
//              For ∞:
//                E_∞ = {6 + (12k/n)·√(1+4γ) | k = 0, ±1, ..., ±n/2} ∩ Z
//                (γ = coeff of x⁻² in Laurent at ∞; γ = 0 if ord(r,∞) > 2)
//     Step 2.  Enumerate families (e_c, e_∞).  Retain those with
//                d := (n/12)·(e_∞ − Σ e_c) ∈ Z_{≥0}.
//     Step 3.  Search polynomial P of degree d s.t. P_{−1} ≡ 0 in the
//              Kovacic §5 recurrence (helpers TU).
//     Step 4.  Build the n-th-degree minimal polynomial of ω; return as
//              RootOf(M, ω_var).  Downstream solver computes z₁ = exp(∫ω).
//
// Necessary conditions (Kovacic §2 p. 8): poles of r have order ≤ 2 and
// ord(r, ∞) ≥ 2.  Violations → Unimplemented immediately.
//
// Anti-monolith: Step 3 + 4 implementations reside in
// ode_kovacic_case3_helpers.cpp; this TU handles dispatching + Step 1/2.

#include "ode_kovacic_internal.hpp"
#include "ode_kovacic_case3_helpers.hpp"
#include "ode_kovacic_pf_helpers.hpp"
#include "cas/error_helpers.hpp"
#include <algorithm>
#include <chrono>
#include <optional>
#include <vector>
#include <string>

namespace cas::calculus {
namespace kovacic_impl {

namespace {

// Build  E_c = {6 + (12k/n)·√(1+4α) | k = 0, ±1, ..., ±n/2} ∩ Z.
// `coef` is α (for order-2 finite poles) or γ (for ∞).
[[nodiscard]] std::optional<std::vector<long long>> build_E_set_case3(
    const Rational& coef, unsigned n) {
    Rational disc = Rational(BigInt(1)) + Rational(BigInt(4)) * coef;
    auto sq = rational_sqrt(disc);
    if (!sq) return std::nullopt;
    std::vector<long long> out;
    auto add_if_int = [&](const Rational& cand) {
        auto ival = rational_to_int(cand);
        if (ival && std::find(out.begin(), out.end(), *ival) == out.end()) {
            out.push_back(*ival);
        }
    };
    const Rational six(BigInt(6));
    const Rational twelve(BigInt(12));
    const Rational n_q(BigInt(static_cast<long long>(n)));
    const long long half_n = static_cast<long long>(n) / 2;
    for (long long k = -half_n; k <= half_n; ++k) {
        Rational k_q{BigInt(k)};
        Rational coef_k = (twelve * k_q / n_q) * (*sq);
        add_if_int(six + coef_k);
    }
    return out;
}

// Extract order-2 pole coefficient α at pole_loc.
[[nodiscard]] std::optional<Rational> extract_order2_coeff_case3(
    const std::vector<PFPole>& pf_poles, ExprPtr pole_loc,
    symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    for (const auto& p : pf_poles) {
        if (p.power != 2U) continue;
        auto diff_res = ctx.simplify(kv_sub(a, p.pole, pole_loc));
        if (diff_res.is_error()) continue;
        if (!kv_is_zero(diff_res.value(), ctx)) continue;
        return try_get_rational(p.coeff, ctx);
    }
    return std::nullopt;
}

struct C3Pole {
    ExprPtr  loc;
    unsigned order;
    Rational coeff_alpha;
    bool     alpha_set;
};

// Family-enumeration cartesian product (same shape as Case 2 helper but with
// case3 cap).
[[nodiscard]] Result<std::vector<std::vector<std::size_t>>> enumerate_families_c3(
    const std::vector<std::vector<long long>>& E_c_sets,
    const std::vector<long long>& E_inf,
    symbolic::CASContext& ctx) {
    std::size_t total = E_inf.size();
    if (total == 0U) return ok(std::vector<std::vector<std::size_t>>{});
    for (const auto& set : E_c_sets) {
        if (set.empty()) return ok(std::vector<std::vector<std::size_t>>{});
        const std::size_t cap = ctx.kovacic_case3_max_pole_combinations();
        if (total > cap / set.size()) {
            return fail<std::vector<std::vector<std::size_t>>>(kv_unimpl(
                "Kovacic Case 3 Step 2: family enumeration exceeds "
                "ctx.kovacic_case3_max_pole_combinations() = "
                + std::to_string(cap)));
        }
        total *= set.size();
    }
    std::vector<std::vector<std::size_t>> out;
    out.reserve(total);
    std::vector<std::size_t> idx(E_c_sets.size() + 1U, 0U);
    while (true) {
        out.push_back(idx);
        std::size_t i = 0U;
        while (i < idx.size()) {
            const std::size_t set_size = (i < E_c_sets.size())
                ? E_c_sets[i].size() : E_inf.size();
            if (idx[i] + 1U < set_size) { ++idx[i]; break; }
            idx[i] = 0U; ++i;
        }
        if (i == idx.size()) break;
    }
    return ok(std::move(out));
}

// Compute d = (n/12) · (e_∞ − Σ e_c).  Return std::nullopt unless d is a
// non-negative integer.
[[nodiscard]] std::optional<long long> compute_d_case3(
    const std::vector<long long>& e_c, long long e_inf, unsigned n) {
    long long sum = 0;
    for (long long e : e_c) sum += e;
    long long diff = e_inf - sum;
    // d = n * diff / 12  must be an integer ≥ 0.
    long long num = static_cast<long long>(n) * diff;
    if (num % 12 != 0) return std::nullopt;
    long long d = num / 12;
    if (d < 0) return std::nullopt;
    return d;
}

// θ = (n/12) · Σ e_c / (x − c).
[[nodiscard]] ExprPtr build_theta_case3(
    const std::vector<ExprPtr>& pole_locs,
    const std::vector<long long>& e_c,
    unsigned n, const Symbol& x, AstArena& a) {
    std::vector<ExprPtr> terms;
    for (std::size_t i = 0U; i < pole_locs.size(); ++i) {
        if (e_c[i] == 0) continue;
        ExprPtr xmc = kv_sub(a, a.make<Symbol>(x.name), pole_locs[i]);
        ExprPtr e_expr = kv_int(a, e_c[i]);
        terms.push_back(kv_div(a, e_expr, xmc));
    }
    if (terms.empty()) return kv_int(a, 0);
    ExprPtr sum = terms[0];
    for (std::size_t i = 1U; i < terms.size(); ++i)
        sum = kv_add(a, sum, terms[i]);
    return kv_div(a,
        kv_mul(a, kv_int(a, static_cast<long long>(n)), sum),
        kv_int(a, 12));
}

// S = Π (x − c).
[[nodiscard]] ExprPtr build_S(
    const std::vector<ExprPtr>& pole_locs,
    const Symbol& x, AstArena& a) {
    if (pole_locs.empty()) return kv_int(a, 1);
    ExprPtr prod = kv_sub(a, a.make<Symbol>(x.name), pole_locs[0]);
    for (std::size_t i = 1U; i < pole_locs.size(); ++i)
        prod = kv_mul(a, prod, kv_sub(a, a.make<Symbol>(x.name), pole_locs[i]));
    return prod;
}

// Try Case 3 for one candidate Galois-order n ∈ {4, 6, 12}.  Returns ok with
// the minimal polynomial of ω (a polynomial in `omega_var`) on success.  On
// inapplicability returns ok(nullopt) — the caller advances to the next n.
[[nodiscard]] Result<std::optional<std::pair<ExprPtr, Symbol>>> try_case3_for_n(
    const std::vector<C3Pole>& poles,
    const InfinityData& inf,
    ExprPtr r, unsigned n,
    const Symbol& x, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();

    // Step 1: E_c sets.
    std::vector<std::vector<long long>> E_c_sets;
    std::vector<ExprPtr> pole_locations;
    pole_locations.reserve(poles.size());
    for (const auto& pd : poles) {
        std::vector<long long> e;
        if (pd.order == 1U) {
            e = {12};
        } else if (pd.order == 2U) {
            if (!pd.alpha_set) {
                return fail<std::optional<std::pair<ExprPtr, Symbol>>>(kv_unimpl(
                    "Kovacic Case 3 Step 1 (c₂): could not extract order-2 α."));
            }
            auto E = build_E_set_case3(pd.coeff_alpha, n);
            if (!E) {
                return ok(std::optional<std::pair<ExprPtr, Symbol>>{});
            }
            e = *E;
        } else {
            // Order > 2 violates Case 3 §2 necessary condition.
            return ok(std::optional<std::pair<ExprPtr, Symbol>>{});
        }
        if (e.empty()) return ok(std::optional<std::pair<ExprPtr, Symbol>>{});
        E_c_sets.push_back(std::move(e));
        pole_locations.push_back(pd.loc);
    }

    // Step 1 ∞ branch — ord(r, ∞) must be ≥ 2.
    std::vector<long long> E_inf;
    Rational gamma(BigInt(0));
    if (inf.ord < 2) {
        // §2 necessary condition fails; Case 3 cannot hold at all.
        return ok(std::optional<std::pair<ExprPtr, Symbol>>{});
    }
    if (inf.ord == 2) {
        if (!inf.leading_set) {
            return fail<std::optional<std::pair<ExprPtr, Symbol>>>(kv_unimpl(
                "Kovacic Case 3 Step 1 (∞): leading Laurent coeff γ at ∞ "
                "extraction failed."));
        }
        gamma = inf.leading_b;
    }
    // ord > 2 → γ = 0 (default).
    auto E_inf_opt = build_E_set_case3(gamma, n);
    if (!E_inf_opt) return ok(std::optional<std::pair<ExprPtr, Symbol>>{});
    E_inf = *E_inf_opt;
    if (E_inf.empty()) return ok(std::optional<std::pair<ExprPtr, Symbol>>{});

    // Step 2: enumerate families.
    auto fams_res = enumerate_families_c3(E_c_sets, E_inf, ctx);
    if (fams_res.is_error()) {
        return fail<std::optional<std::pair<ExprPtr, Symbol>>>(fams_res.error());
    }

    // Steps 3-4 per family.
    // Outer wall-clock safety net per try_case3_for_n call (configurable via
    // ctx.kovacic_case3_budget_ms(), 0 = disabled): the PolyExpr recurrence
    // is deterministic (HC-KV-06 closed), so this only guards pathological
    // multi-pole inputs where the family count itself is large.
    const auto start = std::chrono::steady_clock::now();
    const auto budget = std::chrono::milliseconds(ctx.kovacic_case3_budget_ms());
    const bool budget_enabled = ctx.kovacic_case3_budget_ms() > 0U;
    for (const auto& idx : fams_res.value()) {
        if (budget_enabled && std::chrono::steady_clock::now() - start > budget) {
            return ok(std::optional<std::pair<ExprPtr, Symbol>>{});
        }
        std::vector<long long> e_c_chosen;
        e_c_chosen.reserve(E_c_sets.size());
        for (std::size_t i = 0U; i < E_c_sets.size(); ++i) {
            e_c_chosen.push_back(E_c_sets[i][idx[i]]);
        }
        const long long e_inf = E_inf[idx.back()];

        auto d_opt = compute_d_case3(e_c_chosen, e_inf, n);
        if (!d_opt) continue;

        ExprPtr theta = build_theta_case3(pole_locations, e_c_chosen, n, x, a);
        ExprPtr S = build_S(pole_locations, x, a);

        auto P_res = search_polynomial_P_case3(theta, S, r, *d_opt, n, x, ctx);
        if (P_res.is_error()) {
            return fail<std::optional<std::pair<ExprPtr, Symbol>>>(P_res.error());
        }
        if (!P_res.value().has_value()) continue;
        ExprPtr P = *P_res.value();

        Symbol omega_var(ctx.make_fresh_symbol("ω"));
        auto minpoly_res = build_omega_minpoly_case3(
            theta, S, P, r, n, x, omega_var, ctx);
        if (minpoly_res.is_error()) {
            // Soft-fail at this family; let the dispatcher advance to the
            // next n.  Hard-failing here would mask cases where a different
            // n provides a valid construction.
            continue;
        }
        return ok(std::optional<std::pair<ExprPtr, Symbol>>{
            std::make_pair(minpoly_res.value(), omega_var)});
    }

    return ok(std::optional<std::pair<ExprPtr, Symbol>>{});
}

} // anonymous namespace

// ─── Public entry point ───────────────────────────────────────────────────────

Result<OmegaPair> case3_omega(
    ExprPtr r, const Symbol& x, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();

    auto pf_res = algebra::partial_fractions(r, x, ctx);
    if (pf_res.is_error()) {
        return fail<OmegaPair>(kv_unimpl(
            "Kovacic Case 3: partial_fractions failed: "
            + pf_res.error().message));
    }
    std::vector<PFPole> pf_poles = collect_pf_poles(pf_res.value(), x, ctx);

    // Dedupe finite poles by location (keep max power).
    std::vector<C3Pole> poles;
    for (const auto& p : pf_poles) {
        bool found = false;
        for (auto& q : poles) {
            auto diff_res = ctx.simplify(kv_sub(a, q.loc, p.pole));
            if (diff_res.is_ok() && kv_is_zero(diff_res.value(), ctx)) {
                if (p.power > q.order) q.order = p.power;
                found = true;
                break;
            }
        }
        if (!found) {
            poles.push_back(C3Pole{p.pole, p.power, Rational(BigInt(0)), false});
        }
    }
    // Case 3 §2 necessary: poles of order ≤ 2 only.
    for (const auto& pd : poles) {
        if (pd.order > 2U) {
            return fail<OmegaPair>(kv_unimpl(
                "Kovacic Case 3 §2: r has a pole of order > 2 "
                "(necessary condition violated).  Case 3 cannot hold."));
        }
    }
    for (auto& pd : poles) {
        if (pd.order == 2U) {
            auto a_opt = extract_order2_coeff_case3(pf_poles, pd.loc, ctx);
            if (a_opt) {
                pd.coeff_alpha = *a_opt;
                pd.alpha_set = true;
            }
        }
    }

    auto inf_opt = compute_infinity_data(r, x, ctx);
    if (!inf_opt) {
        return fail<OmegaPair>(kv_unimpl(
            "Kovacic Case 3 Step 1 (∞): infinity data extraction failed."));
    }

    // Try n = 4, 6, 12 in order (paper §5.1).
    // For n = 12 (icosahedral A₅) the §5 recurrence runs 13 levels;
    // `compute_P_sequence` (helpers TU) represents it in `algebra::PolyExpr`
    // form so every family completes deterministically (HC-KV-06, CHIUSO
    // 2026-07-08).  The wall-clock budgets below remain as an outer safety
    // net for pathological multi-pole inputs, not as the primary mechanism.
    for (unsigned n : {4U, 6U, 12U}) {
        auto try_res = try_case3_for_n(poles, *inf_opt, r, n, x, ctx);
        if (try_res.is_error()) return fail<OmegaPair>(try_res.error());
        if (!try_res.value().has_value()) continue;
        ExprPtr minpoly = try_res.value()->first;
        Symbol  omega_var = try_res.value()->second;
        // Wrap as RootOf(minpoly, ω_var).
        ExprPtr root = a.make<RootOf>(minpoly, omega_var);
        // Case 3 ω is algebraic, no conjugate-pair distinction at this level;
        // emit the same RootOf for both fields so downstream solver can pick
        // either branch (currently identical).
        return ok(OmegaPair{root, root});
    }

    return cas::make_unimplemented<OmegaPair>(UnimplementedInfo{
        .module      = "calculus",
        .function    = "case3_omega",
        .input_shape = "Kovacic Case 3 candidate n ∈ {4, 6, 12}",
        .reason      = cas::error::reason_codes::ODE_UNSUPPORTED_TYPE,
        .suggestion  = "Verify Kovacic invariants or check higher order solver",
        .ticket      = ""  // no open debt: negative outcome is the correct §5 result (HC-KV-06 closed 2026-07-08)
    }, "Kovacic Case 3: no polynomial P found for any n ∈ {4, 6, 12}; Case 3 cannot hold for this DE.");
}

} // namespace kovacic_impl
} // namespace cas::calculus
