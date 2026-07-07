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

protected:
    int           max_simplification_depth_{300};
    std::size_t   max_recursion_depth_{256U};
    bool          strict_branch_cuts_{false};
    long long     max_trig_power_reduction_{32LL};
    long long     symbolic_qr_max_norm_complexity_{2LL};
    int           max_trig_exact_denom_{100};
    std::size_t   simplify_sqrt_trial_division_bound_{10000U};
    bool          branch_cut_aware_logexp_{false};
};

}  // namespace cas::symbolic

