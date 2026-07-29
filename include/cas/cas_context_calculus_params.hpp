#pragma once

#include <cstddef>
#include <cstdint>

namespace cas::symbolic {

struct CASContextCalculusParams {
    // ── Integrator ──────────────────────────────────────────────────────────
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] std::size_t max_integration_depth() const noexcept {
        return max_integration_depth_;
    }

    // ── Bessel recurrence expansion ──────────────────────────────────────────
    // Opt-in: expand BesselJ/Y(n,x) via three-term recurrence. Default false.
    // (setter trivial; kept in CASContext for symmetry with HC-001..003 group)
    [[nodiscard]] bool expand_bessel_recurrence() const noexcept {
        return expand_bessel_recurrence_;
    }

    // ── Integration-by-parts depth ───────────────────────────────────────────
    void set_max_integrate_by_parts_depth(std::size_t n) noexcept {
        max_integrate_by_parts_depth_ = n;
    }
    [[nodiscard]] std::size_t max_integrate_by_parts_depth() const noexcept {
        return max_integrate_by_parts_depth_;
    }

    // ── Risch DE rational ansatz degree (BUG-HANG-001) ──────────────────────
    void set_max_risch_rational_ansatz_degree(std::size_t n) noexcept {
        max_risch_rational_ansatz_degree_ = n;
    }
    [[nodiscard]] std::size_t max_risch_rational_ansatz_degree() const noexcept {
        return max_risch_rational_ansatz_degree_;
    }

    // ── limit log/log fast-path recursion depth (HPP-022 closure) ──────────
    void set_max_log_log_limit_depth(unsigned int d) noexcept {
        max_log_log_limit_depth_ = d;
    }
    [[nodiscard]] unsigned int max_log_log_limit_depth() const noexcept {
        return max_log_log_limit_depth_;
    }

    // ── Adaptive Gauss-Kronrod numerical integration (F6.D) ────────────────
    void set_integration_abs_tol(double t) noexcept {
        integration_abs_tol_ = t;
    }
    [[nodiscard]] double integration_abs_tol() const noexcept {
        return integration_abs_tol_;
    }
    void set_integration_rel_tol(double t) noexcept {
        integration_rel_tol_ = t;
    }
    [[nodiscard]] double integration_rel_tol() const noexcept {
        return integration_rel_tol_;
    }
    void set_integration_max_intervals(std::size_t n) noexcept {
        integration_max_intervals_ = n;
    }
    [[nodiscard]] std::size_t integration_max_intervals() const noexcept {
        return integration_max_intervals_;
    }

    // ── Confluent hypergeometric ₁F₁ series (A7 §6.7 numeric shadow) ───────
    // Kummer series Σ (a)_k/(b)_k · z^k/k!. Entire in z ⇒ converges everywhere,
    // but a large |z| or near-cancelling parameters may need more terms than the
    // default before the relative term ratio drops below tol; both bounds are
    // configurable so a caller can trade cost for reach (Unimplemented otherwise).
    void set_hypergeometric_1f1_max_terms(std::size_t n) noexcept {
        hypergeometric_1f1_max_terms_ = n;
    }
    [[nodiscard]] std::size_t hypergeometric_1f1_max_terms() const noexcept {
        return hypergeometric_1f1_max_terms_;
    }
    void set_hypergeometric_1f1_rel_tol(double t) noexcept {
        hypergeometric_1f1_rel_tol_ = t;
    }
    [[nodiscard]] double hypergeometric_1f1_rel_tol() const noexcept {
        return hypergeometric_1f1_rel_tol_;
    }

    // ── Bessel half-integer recurrence expansion bound (F7.5.E2) ───────────
    void set_max_bessel_half_integer_order(unsigned int n) noexcept {
        max_bessel_half_integer_order_ = n;
    }
    [[nodiscard]] unsigned int max_bessel_half_integer_order() const noexcept {
        return max_bessel_half_integer_order_;
    }

    // ── MRV cycle guards (HPP-F75-AUDIT-CYCLE-GUARD-1/2/3) ───────────────────
    void set_mrv_max_append_depth(unsigned int depth) noexcept {
        mrv_max_append_depth_ = depth;
    }
    [[nodiscard]] unsigned int mrv_max_append_depth() const noexcept {
        return mrv_max_append_depth_;
    }

    void set_diff_field_max_visit_depth(unsigned int depth) noexcept {
        diff_field_max_visit_depth_ = depth;
    }
    [[nodiscard]] unsigned int diff_field_max_visit_depth() const noexcept {
        return diff_field_max_visit_depth_;
    }

    void set_mrv_growth_rank_max_depth(int depth) noexcept {
        mrv_growth_rank_max_depth_ = depth;
    }
    [[nodiscard]] int mrv_growth_rank_max_depth() const noexcept {
        return mrv_growth_rank_max_depth_;
    }

    // HPP-026: max |k| range for trig singularity scan in definite integration.
    void set_integration_singularity_scan_max_k(std::size_t k) noexcept {
        integration_singularity_scan_max_k_ = k;
    }
    [[nodiscard]] std::size_t integration_singularity_scan_max_k() const noexcept {
        return integration_singularity_scan_max_k_;
    }

    // F4.K2 Kovacic Case 2 (Kovacic 1986 §3.2) ────────────────────────────────
    void set_kovacic_case2_max_pole_combinations(std::size_t n) noexcept {
        kovacic_case2_max_pole_combinations_ = n;
    }
    [[nodiscard]] std::size_t kovacic_case2_max_pole_combinations() const noexcept {
        return kovacic_case2_max_pole_combinations_;
    }

    void set_kovacic_case2_max_poly_degree(std::size_t d) noexcept {
        kovacic_case2_max_poly_degree_ = d;
    }
    [[nodiscard]] std::size_t kovacic_case2_max_poly_degree() const noexcept {
        return kovacic_case2_max_poly_degree_;
    }

    // F4.K3 Kovacic Case 3 (Kovacic 1986 §5) ──────────────────────────────────
    void set_kovacic_case3_max_pole_combinations(std::size_t n) noexcept {
        kovacic_case3_max_pole_combinations_ = n;
    }
    [[nodiscard]] std::size_t kovacic_case3_max_pole_combinations() const noexcept {
        return kovacic_case3_max_pole_combinations_;
    }

    void set_kovacic_case3_max_poly_degree(std::size_t d) noexcept {
        kovacic_case3_max_poly_degree_ = d;
    }
    [[nodiscard]] std::size_t kovacic_case3_max_poly_degree() const noexcept {
        return kovacic_case3_max_poly_degree_;
    }

    // Outer wall-clock safety net for the Case 3 §5 machinery: bounds one
    // try_case3_for_n family sweep AND one compute_P_sequence run.  The
    // recurrence is deterministic in PolyExpr form (HC-KV-06 closure); this
    // budget only guards pathological multi-pole inputs.  0 = disabled.
    void set_kovacic_case3_budget_ms(std::size_t ms) noexcept {
        kovacic_case3_budget_ms_ = ms;
    }
    [[nodiscard]] std::size_t kovacic_case3_budget_ms() const noexcept {
        return kovacic_case3_budget_ms_;
    }

    // ── together() polynomial GCD content reduction ──────────────────────────
    void set_together_gcd_enabled(bool enabled) noexcept {
        together_gcd_enabled_ = enabled;
    }
    [[nodiscard]] bool together_gcd_enabled() const noexcept {
        return together_gcd_enabled_;
    }
    void set_together_gcd_max_degree(std::size_t deg) noexcept {
        together_gcd_max_degree_ = deg;
    }
    [[nodiscard]] std::size_t together_gcd_max_degree() const noexcept {
        return together_gcd_max_degree_;
    }
    void set_together_gcd_max_symbols(std::size_t n) noexcept {
        together_gcd_max_symbols_ = n;
    }
    [[nodiscard]] std::size_t together_gcd_max_symbols() const noexcept {
        return together_gcd_max_symbols_;
    }

    // ── u-substitution residual fixpoint passes (A42) ────────────────────────
    // integrate_by_substitution re-runs replace(g→u)+simplify on the
    // residual while it still depends on the integration variable: algebraic
    // reduction (e.g. x^3/(2x) -> x^2/2) can surface a literal occurrence of
    // g only *after* simplify, invisible to the pass that produced it. The
    // loop already self-terminates on a stable fixpoint (structural_equal
    // unchanged); this only bounds a candidate that never stabilizes.
    void set_max_substitution_fixpoint_passes(unsigned int n) noexcept {
        max_substitution_fixpoint_passes_ = n;
    }
    [[nodiscard]] unsigned int max_substitution_fixpoint_passes() const noexcept {
        return max_substitution_fixpoint_passes_;
    }

    // ── Taylor expansion across a removable singularity (A37) ───────────────
    // f = N/D with D(center) = 0 has a Taylor expansion iff val(N) >= val(D).
    // Computing it needs the numerator and denominator series carried
    // val(D) terms deeper than the requested order; val(D) is not known up
    // front, so the search doubles the extra depth (1, 2, 4, ...) until the
    // denominator's leading term becomes visible. This bounds that search:
    // beyond it the expansion is reported Unimplemented (never silently
    // truncated). Default 16 = a denominator vanishing to order 16 at a point,
    // far past anything the truncated-series ring is asked for in practice.
    void set_max_removable_singularity_order(unsigned int n) noexcept {
        max_removable_singularity_order_ = n;
    }
    [[nodiscard]] unsigned int max_removable_singularity_order() const noexcept {
        return max_removable_singularity_order_;
    }

protected:
    std::size_t   max_integration_depth_{16U};
    unsigned int  max_removable_singularity_order_{16U};
    bool          expand_bessel_recurrence_{false};
    unsigned int  max_bessel_half_integer_order_{64U};
    std::size_t   hypergeometric_1f1_max_terms_{4000U};
    double        hypergeometric_1f1_rel_tol_{1e-16};
    double        integration_abs_tol_{1e-10};
    double        integration_rel_tol_{1e-8};
    std::size_t   integration_max_intervals_{4096U};
    std::size_t   max_integrate_by_parts_depth_{8U};
    std::size_t   max_risch_rational_ansatz_degree_{32U};
    unsigned int  max_log_log_limit_depth_{3U};
    unsigned int  max_substitution_fixpoint_passes_{4U};
    std::size_t   integration_singularity_scan_max_k_{1000U};
    std::size_t   kovacic_case2_max_pole_combinations_{1024U};
    std::size_t   kovacic_case2_max_poly_degree_{64U};
    std::size_t   kovacic_case3_max_pole_combinations_{4096U};
    std::size_t   kovacic_case3_max_poly_degree_{32U};
    std::size_t   kovacic_case3_budget_ms_{2000U};
    unsigned int  mrv_max_append_depth_{1024U};
    unsigned int  diff_field_max_visit_depth_{4096U};
    int           mrv_growth_rank_max_depth_{1024};
    bool          together_gcd_enabled_{true};
    std::size_t   together_gcd_max_degree_{64U};
    std::size_t   together_gcd_max_symbols_{8U};
};

}  // namespace cas::symbolic
