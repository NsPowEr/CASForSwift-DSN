// Kovacic Case 2 — algebraic ω of degree 2 over C(x), Galois group ⊆ D⁺.
//
// Reference: Kovacic J.J. (1986), §4, "An Algorithm for Solving Second
// Order Linear Homogeneous Differential Equations", J. Symbolic
// Computation 2, 3-43.  Algorithm pp. 18-19, proof pp. 20-22.
//
// Spec (paper-verified): .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/
//                        Kovacic_Case2.md
// Ledger: HC-KV-03 (HARDCODE_LEDGER.md).
//
// Algorithm overview (assumes Case 1 has failed):
//   Step 1.  For each finite pole c of r, build E_c ⊂ Z:
//              order 1     →  E_c = {4}
//              order 2     →  E_c = {2 + k·√(1+4b) | k = 0, ±2} ∩ Z
//                              (b = coeff of 1/(x-c)² in partial fraction)
//              order v > 2 →  E_c = {v}
//            For ∞:
//              ord_∞ > 2   →  E_∞ = {0, 2, 4}
//              ord_∞ = 2   →  E_∞ = {2 + k·√(1+4b) | k = 0, ±2} ∩ Z
//                              (b = leading Laurent coeff at ∞)
//              ord_∞ < 2   →  E_∞ = {ord_∞}
//   Step 2.  Enumerate families (e_c, e_∞).  Discard families with all-even
//            coordinates or where d := (e_∞ - Σe_c)/2 ∉ Z_{≥0}.
//   Step 3.  Search polynomial P of degree d.  (See helpers TU.)
//   Step 4.  Build φ = θ + P'/P and ω from quadratic.  (See helpers TU.)
//
// Anti-monolith: Step 3 + 4 implementations reside in
// ode_kovacic_case2_helpers.cpp; this TU contains only the dispatcher
// and the Step 1 + Step 2 sub-routines.

#include "ode_kovacic_internal.hpp"
#include "ode_kovacic_case2_helpers.hpp"
#include "ode_kovacic_pf_helpers.hpp"
#include <algorithm>
#include <optional>
#include <vector>
#include <string>

namespace cas::calculus {
namespace kovacic_impl {

namespace {

// ─── PoleData: structured Case 2 pole information ────────────────────────────

struct PoleData {
    ExprPtr  loc;
    bool     is_infinity;
    unsigned order;
    Rational coeff_b;          // For order-2 poles: coefficient `b` of 1/(x-c)².
    bool     coeff_b_set;
};

[[nodiscard]] std::optional<Rational> sqrt_one_plus_4b(const Rational& b) {
    Rational disc = Rational(BigInt(1)) + Rational(BigInt(4)) * b;
    return rational_sqrt(disc);
}

// E_c = {2 + k·√(1+4b) | k = 0, ±2} ∩ Z  (Kovacic §4.1 (c₂) / (∞₂)).
[[nodiscard]] std::optional<std::vector<long long>> build_E_c_order2(
    const Rational& b) {
    auto sq = sqrt_one_plus_4b(b);
    if (!sq) return std::nullopt;
    std::vector<long long> out;
    auto add_if_int = [&](const Rational& cand) {
        auto ival = rational_to_int(cand);
        if (ival && std::find(out.begin(), out.end(), *ival) == out.end()) {
            out.push_back(*ival);
        }
    };
    const Rational two(BigInt(2));
    add_if_int(two);
    add_if_int(two + two * (*sq));
    add_if_int(two - two * (*sq));
    return out;
}

// Extract order-2 coeff b at pole_loc from collected PFPole entries.
[[nodiscard]] std::optional<Rational> extract_order2_coeff(
    const std::vector<PFPole>& pf_poles,
    ExprPtr pole_loc,
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

// Cartesian product enumeration (Step 2).  Caps total combinations by
// ctx.kovacic_case2_max_pole_combinations() to bound explosion.
[[nodiscard]] Result<std::vector<std::vector<std::size_t>>> enumerate_families(
    const std::vector<std::vector<long long>>& E_c_sets,
    const std::vector<long long>& E_inf,
    symbolic::CASContext& ctx) {

    std::size_t total = E_inf.size();
    if (total == 0U) return ok(std::vector<std::vector<std::size_t>>{});
    for (const auto& set : E_c_sets) {
        if (set.empty()) return ok(std::vector<std::vector<std::size_t>>{});
        const std::size_t cap = ctx.kovacic_case2_max_pole_combinations();
        if (total > cap / set.size()) {
            return fail<std::vector<std::vector<std::size_t>>>(kv_unimpl(
                "Kovacic Case 2 Step 2: family enumeration exceeds "
                "ctx.kovacic_case2_max_pole_combinations() = "
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

[[nodiscard]] bool all_even(const std::vector<long long>& vs) {
    for (long long v : vs) if (v % 2 != 0) return false;
    return true;
}

[[nodiscard]] std::optional<long long> compute_d(
    const std::vector<long long>& e_c, long long e_inf) {
    long long sum = 0;
    for (long long e : e_c) sum += e;
    long long diff = e_inf - sum;
    if (diff % 2 != 0) return std::nullopt;
    long long d = diff / 2;
    if (d < 0) return std::nullopt;
    return d;
}

// θ = (1/2) Σ e_c / (x - c)
[[nodiscard]] ExprPtr build_theta(
    const std::vector<ExprPtr>& pole_locs,
    const std::vector<long long>& e_c,
    const Symbol& x, AstArena& a) {
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
    return kv_div(a, sum, kv_int(a, 2));
}

// (Shared PF/∞ helpers moved to ode_kovacic_pf_helpers.{hpp,cpp}.)

} // anonymous namespace

// ─── Public entry point ───────────────────────────────────────────────────────

Result<OmegaPair> case2_omega(
    ExprPtr r, const Symbol& x, symbolic::CASContext& ctx) {

    AstArena& a = ctx.arena();

    auto pf_res = algebra::partial_fractions(r, x, ctx);
    if (pf_res.is_error()) {
        return fail<OmegaPair>(kv_unimpl(
            "Kovacic Case 2: partial_fractions failed: "
            + pf_res.error().message));
    }
    std::vector<PFPole> pf_poles = collect_pf_poles(pf_res.value(), x, ctx);

    // Collect unique poles and their orders.
    std::vector<PoleData> poles;
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
            poles.push_back(PoleData{p.pole, false, p.power, Rational(BigInt(0)), false});
        }
    }
    for (auto& pd : poles) {
        if (pd.order == 2U) {
            auto b_opt = extract_order2_coeff(pf_poles, pd.loc, ctx);
            if (b_opt) {
                pd.coeff_b = *b_opt;
                pd.coeff_b_set = true;
            }
        }
    }

    // Step 1 finite-pole E_c.
    std::vector<std::vector<long long>> E_c_sets;
    std::vector<ExprPtr> pole_locations;
    pole_locations.reserve(poles.size());
    for (const auto& pd : poles) {
        std::vector<long long> e;
        if (pd.order == 1U) {
            e = {4};
        } else if (pd.order == 2U) {
            if (!pd.coeff_b_set) {
                return fail<OmegaPair>(kv_unimpl(
                    "Kovacic Case 2 Step 1 (c₂): could not extract "
                    "order-2 coefficient b at pole."));
            }
            auto E = build_E_c_order2(pd.coeff_b);
            if (!E) {
                return fail<OmegaPair>(kv_unimpl(
                    "Kovacic Case 2 Step 1 (c₂): √(1+4b) ∉ Q at a finite "
                    "order-2 pole; algebraic extension required."));
            }
            e = *E;
        } else {
            e = {static_cast<long long>(pd.order)};
        }
        if (e.empty()) {
            return fail<OmegaPair>(kv_unimpl(
                "Kovacic Case 2 Step 1: empty E_c at pole."));
        }
        E_c_sets.push_back(std::move(e));
        pole_locations.push_back(pd.loc);
    }

    // Step 1 ∞ branch.
    auto inf_opt = compute_infinity_data(r, x, ctx);
    if (!inf_opt) {
        return fail<OmegaPair>(kv_unimpl(
            "Kovacic Case 2 Step 1 (∞): infinity data extraction failed."));
    }
    const long long ord_inf = inf_opt->ord;
    std::vector<long long> E_inf;
    if (ord_inf > 2) {
        E_inf = {0, 2, 4};
    } else if (ord_inf == 2) {
        if (!inf_opt->leading_set) {
            return fail<OmegaPair>(kv_unimpl(
                "Kovacic Case 2 Step 1 (∞₂): leading Laurent coeff at ∞ "
                "extraction failed."));
        }
        auto E = build_E_c_order2(inf_opt->leading_b);
        if (!E) {
            return fail<OmegaPair>(kv_unimpl(
                "Kovacic Case 2 Step 1 (∞₂): √(1+4b) ∉ Q at ∞; algebraic "
                "extension required."));
        }
        E_inf = *E;
    } else {
        E_inf = {ord_inf};
    }
    if (E_inf.empty()) {
        return fail<OmegaPair>(kv_unimpl(
            "Kovacic Case 2 Step 1 (∞): empty E_∞."));
    }

    // Step 2: enumerate families.
    auto fams_res = enumerate_families(E_c_sets, E_inf, ctx);
    if (fams_res.is_error()) return fail<OmegaPair>(fams_res.error());


    // Steps 3-4: iterate, search P, build ω on success.
    for (const auto& idx : fams_res.value()) {
        std::vector<long long> e_c_chosen;
        e_c_chosen.reserve(E_c_sets.size());
        for (std::size_t i = 0U; i < E_c_sets.size(); ++i) {
            e_c_chosen.push_back(E_c_sets[i][idx[i]]);
        }
        const long long e_inf = E_inf[idx.back()];

        std::vector<long long> all_e = e_c_chosen;
        all_e.push_back(e_inf);
        if (all_even(all_e)) continue;
        auto d_opt = compute_d(e_c_chosen, e_inf);
        if (!d_opt) continue;

        ExprPtr theta = build_theta(pole_locations, e_c_chosen, x, a);
        auto P_res = search_polynomial_P_case2(theta, r, *d_opt, x, ctx);
        if (P_res.is_error()) return fail<OmegaPair>(P_res.error());
        if (!P_res.value().has_value()) continue;
        ExprPtr P = *P_res.value();

        auto omega_res = build_omega_from_phi_case2(theta, P, r, x, ctx);
        if (omega_res.is_error()) return fail<OmegaPair>(omega_res.error());
        return omega_res;
    }

    return fail<OmegaPair>(kv_unimpl(
        "Kovacic Case 2: no family yielded a polynomial solution; Case 2 "
        "cannot hold for this DE."));
}

} // namespace kovacic_impl
} // namespace cas::calculus
