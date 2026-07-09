#pragma once

#include <cstddef>
#include <cstdint>

namespace cas::symbolic {

struct CASContextSimplifierParams {
    // ── Simplifier ──────────────────────────────────────────────────────────
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] int max_simplification_depth() const noexcept {
        return max_simplification_depth_;
    }

    // ── Branch cuts (L2-21) ─────────────────────────────────────────────────
    // when true, identities requiring the principal branch are refused without
    // explicit positivity. Default false = historical behaviour.
    void set_strict_branch_cuts(bool strict) noexcept {
        strict_branch_cuts_ = strict;
    }
    [[nodiscard]] bool strict_branch_cuts() const noexcept {
        return strict_branch_cuts_;
    }

    // ── Trig power reduction ─────────────────────────────────────────────────
    void set_max_trig_power_reduction(long long n) noexcept {
        max_trig_power_reduction_ = n;
    }
    [[nodiscard]] long long max_trig_power_reduction() const noexcept {
        return max_trig_power_reduction_;
    }

    // ── Symbolic QR norm-complexity bailout threshold ───────────────────────
    void set_symbolic_qr_max_norm_complexity(long long n) noexcept {
        symbolic_qr_max_norm_complexity_ = n;
    }
    [[nodiscard]] long long symbolic_qr_max_norm_complexity() const noexcept {
        return symbolic_qr_max_norm_complexity_;
    }

    // ── Exact trig evaluation ────────────────────────────────────────────────
    void set_max_trig_exact_denom(int q) noexcept { max_trig_exact_denom_ = q; }
    [[nodiscard]] int max_trig_exact_denom() const noexcept {
        return max_trig_exact_denom_;
    }

    // ── sqrt rational simplification ────────────────────────────────────────
    void set_simplify_sqrt_trial_division_bound(std::size_t n) noexcept {
        simplify_sqrt_trial_division_bound_ = n;
    }
    [[nodiscard]] std::size_t simplify_sqrt_trial_division_bound() const noexcept {
        return simplify_sqrt_trial_division_bound_;
    }

    // F8.0-6.2: opt-in branch-cut aware logexp reduction
    void set_branch_cut_aware_logexp(bool on) noexcept {
        branch_cut_aware_logexp_ = on;
    }
    [[nodiscard]] bool branch_cut_aware_logexp() const noexcept {
        return branch_cut_aware_logexp_;
    }

    // F0.0-A20: Configurable recursion depth bound for simplifier/evaluator
    [[nodiscard]] std::size_t max_recursion_depth() const noexcept {
        return max_recursion_depth_;
    }

    // A30: deterministic per-operation node-visit budget (primary resource
    // gate).  Counted per top-level CASContext operation (ops_count_), so the
    // outcome for a given input+params is machine- and load-independent; the
    // wall-clock timeout stays as an outer anti-hang safety net only.
    // 0 = disabled.  Default 2'000'000 ≈ 17x the heaviest legitimate operation
    // measured across the full quick suite + heavy calculus suites
    // (117k node visits, calibration 2026-07-08); raise it consciously if a
    // legitimate operation ever reports OPS_BUDGET_EXCEEDED.
    // Contract: an explicit set_timeout() call without an explicit
    // set_max_operation_ops() disables this default gate — the caller then
    // owns the operation budget via its wall-clock deadline.
    // (setter declared in CASContext / context_params.cpp)
    [[nodiscard]] std::uint64_t max_operation_ops() const noexcept {
        return max_operation_ops_;
    }

    // HPP-016: number of interning shards (0 = derive from hardware concurrency)
    [[nodiscard]] std::size_t intern_shards() const noexcept {
        return intern_shards_;
    }
    void set_intern_shards(std::size_t n) noexcept {
        intern_shards_ = n;
    }

protected:
    int           max_simplification_depth_{300};
    std::size_t   max_recursion_depth_{256U};
    std::uint64_t max_operation_ops_{2'000'000ULL};
    std::size_t   intern_shards_{0U};
    bool          max_operation_ops_explicit_{false};
    bool          strict_branch_cuts_{false};
    long long     max_trig_power_reduction_{32LL};
    long long     symbolic_qr_max_norm_complexity_{2LL};
    int           max_trig_exact_denom_{100};
    std::size_t   simplify_sqrt_trial_division_bound_{10000U};
    bool          branch_cut_aware_logexp_{false};
};

}  // namespace cas::symbolic

