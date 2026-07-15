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

    // A31 fase 1 (Domain_Conditions_Propagation.md §3.5): bound on the number
    // of distinct DomainCondition entries a single top-level simplify() call
    // may accumulate (SideConditionSet, side_conditions.hpp). Exceeded ->
    // CASContext::emit_side_condition returns a structured Unimplemented
    // (SIDE_CONDITION_BUDGET_EXCEEDED), never a silent drop. Default 256: no
    // legitimate simplify() call needs hundreds of distinct domain
    // assumptions in one pass; mirrors the other small-working-set defaults
    // in this struct (not derived from a formal bound, so kept generous).
    [[nodiscard]] std::size_t max_side_conditions() const noexcept {
        return max_side_conditions_;
    }
    void set_max_side_conditions(std::size_t n) noexcept {
        max_side_conditions_ = n;
    }

    // A7 (Meijer_G_Slater.md §7.3, §8.3): p+q above this -> Unimplemented
    // diagnostic at make_meijerg (never a silent truncation of the parameter
    // list). Default 20: generous over every table §5 entry (p+q <= 5) and
    // every pFq bridge case in §3.1/§10.5 used by this spec.
    [[nodiscard]] std::size_t meijerg_max_param_count() const noexcept {
        return meijerg_max_param_count_;
    }
    void set_meijerg_max_param_count(std::size_t n) noexcept {
        meijerg_max_param_count_ = n;
    }

    // Recursion cap for to_meijerg()'s Product/pFq decomposition (§7.1 of
    // to_meijerg's own future pipeline, Brick 3). Default 8: no elementary
    // conversion in table §5 nests more than 2 levels; generous margin.
    [[nodiscard]] std::size_t meijerg_max_conversion_depth() const noexcept {
        return meijerg_max_conversion_depth_;
    }
    void set_meijerg_max_conversion_depth(std::size_t n) noexcept {
        meijerg_max_conversion_depth_ = n;
    }

    // When true, numeric evaluation of a G-node outside its proven
    // convergence region (§2.1) is refused rather than attempted (Brick 4+).
    // Default true: REGOLA ZERO, never silently evaluate an undefined
    // integral.
    [[nodiscard]] bool meijerg_strict_convergence() const noexcept {
        return meijerg_strict_convergence_;
    }
    void set_meijerg_strict_convergence(bool strict) noexcept {
        meijerg_strict_convergence_ = strict;
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
    std::size_t   max_side_conditions_{256U};
    std::size_t   meijerg_max_param_count_{20U};
    std::size_t   meijerg_max_conversion_depth_{8U};
    bool          meijerg_strict_convergence_{true};
};

}  // namespace cas::symbolic

