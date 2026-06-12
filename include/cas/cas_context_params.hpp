#pragma once

// cas_context_params.hpp — Configurable algorithm parameter storage for CASContext.
//
// Defines CASContextParams: a base struct that carries every tuneable knob
// (thresholds, budgets, quality parameters) exposed by CASContext.
//
// Simple setters (no clamping) are inline here.
// Setters with clamping/validation remain as CASContext out-of-line methods
// in context_core.cpp (they are declared in symbolic.hpp).
//
// CASContext inherits from CASContextParams so all ctx.<getter/setter>() calls
// work without any call-site changes.  Protected field storage is accessible
// to CASContext setters via inheritance.
//
// Anti-monolith: moving field declarations + accessors here keeps symbolic.hpp
// under 500 lines (debt SPLIT-SYMBOLIC-HPP-F2.5 closed).

#include <cstddef>
#include <cstdint>

namespace cas::symbolic {

struct CASContextParams {
    // ── Simplifier ──────────────────────────────────────────────────────────
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] int max_simplification_depth() const noexcept {
        return max_simplification_depth_;
    }

    // ── Integrator ──────────────────────────────────────────────────────────
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] std::size_t max_integration_depth() const noexcept {
        return max_integration_depth_;
    }

    // ── GCD probabilistic error bound ───────────────────────────────────────
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] double gcd_error_probability() const noexcept {
        return gcd_error_probability_;
    }

    // ── Numeric floating-point precision (L3-03) ────────────────────────────
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] unsigned int numeric_precision_digits() const noexcept {
        return numeric_precision_digits_;
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
    // Maximum estimate_complexity() of the column norm N_x = Σ xᵢ² above
    // which symbolic QR (Householder) bails out with Unimplemented when N_x
    // is NOT provably nonnegative.  When N_x is provably nonnegative (via
    // assumptions), the bailout is skipped regardless of this threshold.
    // Default 2 — matches the historical bailout introduced when MGS was the
    // active QR backend.  Bound rationale: complexity ≤ 2 covers numeric
    // entries (sqrt(5), sqrt(2)+1, …) where the simplifier reliably folds
    // sqrt(p)·sqrt(p) → p; symbolic entries with deeper structure require
    // the branch-cut-aware sqrt rule (HC-F8-PENDING-20 residue).
    //   Spec: Householder_Symbolic_Stable.md §HH-3 + HARDCODE_LEDGER.md
    //         HC-F8-QR-HOUSEHOLDER-BAILOUT.
    void set_symbolic_qr_max_norm_complexity(long long n) noexcept {
        symbolic_qr_max_norm_complexity_ = n;
    }
    [[nodiscard]] long long symbolic_qr_max_norm_complexity() const noexcept {
        return symbolic_qr_max_norm_complexity_;
    }

    // ── Exact trig evaluation ────────────────────────────────────────────────
    // Maximum denominator q for exact cos(p·π/q)/sin(p·π/q).
    // Default 100 (Disquisitiones §VII: Gauss constructibility criterion).
    void set_max_trig_exact_denom(int q) noexcept { max_trig_exact_denom_ = q; }
    [[nodiscard]] int max_trig_exact_denom() const noexcept {
        return max_trig_exact_denom_;
    }

    // ── RootOf ──────────────────────────────────────────────────────────────
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] std::size_t max_rootof_explicit_degree() const noexcept {
        return max_rootof_explicit_degree_;
    }

    // ── GCD recursion ────────────────────────────────────────────────────────
    // (setters with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] std::size_t max_gcd_recursion_depth() const noexcept {
        return max_gcd_recursion_depth_;
    }

    [[nodiscard]] std::size_t min_gcd_division_steps() const noexcept {
        return min_gcd_division_steps_;
    }

    // Maximum total recursive calls in the multivariate GCD engine.
    // Derivation: 4096 = 4^6 ≈ (3 offsets × 3 samples)^3 levels.
    // Exceeded → Unimplemented (GCD_MULTIVARIATE_BUDGET_EXCEEDED).
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] std::size_t max_gcd_total_calls() const noexcept {
        return max_gcd_total_calls_;
    }

    // ── Cyclotomic ───────────────────────────────────────────────────────────
    // -1 = auto-derive from polynomial degree (2*(deg+1)).
    // (setter trivial; kept in CASContext for symmetry with other GCD setters)
    [[nodiscard]] int max_cyclotomic_n() const noexcept {
        return max_cyclotomic_n_;
    }

    // ── Half-GCD dispatch threshold ──────────────────────────────────────────
    // For univariate integer GCD: if min(deg(a),deg(b)) > threshold use Half-GCD
    // (Knuth-Schönhage O(M(n)log n)) instead of subresultant PRS. Default 200.
    void set_half_gcd_degree_threshold(std::size_t n) noexcept {
        half_gcd_degree_threshold_ = n;
    }
    [[nodiscard]] std::size_t half_gcd_degree_threshold() const noexcept {
        return half_gcd_degree_threshold_;
    }

    // ── Modular GCD CRT dispatch threshold (B2.1) ───────────────────────────
    // When max coefficient bit-length exceeds this, try CRT GCD before subresultant.
    // Default 48 bits; 0 → always CRT; SIZE_MAX → disable CRT.
    // Derivation: 48 = 3 * 16 (three limbs of 16 bits safe for Mignotte bound).
    void set_modular_gcd_coeff_bits(std::size_t bits) noexcept {
        modular_gcd_coeff_bits_ = bits;
    }
    [[nodiscard]] std::size_t modular_gcd_coeff_bits() const noexcept {
        return modular_gcd_coeff_bits_;
    }

    // ── Berlekamp Q-matrix budget ────────────────────────────────────────────
    // berlekamp_factor_mod_p returns Unimplemented when deg(f)*p > limit.
    // Caller falls back to Cantor-Zassenhaus. Default 1024.
    void set_max_berlekamp_matrix_size(std::size_t n) noexcept {
        max_berlekamp_matrix_size_ = (n < 1U) ? 1U : n;
    }
    [[nodiscard]] std::size_t max_berlekamp_matrix_size() const noexcept {
        return max_berlekamp_matrix_size_;
    }

    // ── HC-001..003 configurable knobs ──────────────────────────────────────
    // (setters with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] std::size_t max_q_alpha_bridge_depth() const noexcept {
        return max_q_alpha_bridge_depth_;
    }

    [[nodiscard]] std::size_t max_gamma_recursion() const noexcept {
        return max_gamma_recursion_;
    }

    [[nodiscard]] std::size_t improper_leading_order_scan() const noexcept {
        return improper_leading_order_scan_;
    }

    // ── Bessel recurrence expansion ──────────────────────────────────────────
    // Opt-in: expand BesselJ/Y(n,x) via three-term recurrence. Default false.
    // (setter trivial; kept in CASContext for symmetry with HC-001..003 group)
    [[nodiscard]] bool expand_bessel_recurrence() const noexcept {
        return expand_bessel_recurrence_;
    }

    // ── Trager tower shift attempts (L3-06) ──────────────────────────────────
    // 0 = auto-derive from discriminant collision bound.
    // (setter trivial; kept in CASContext for symmetry)
    [[nodiscard]] std::size_t max_trager_tower_shift_attempts() const noexcept {
        return max_trager_tower_shift_attempts_;
    }

    // ── Wang multivariate factor good-evaluation-point search radius (F3.2) ──
    // 0 = auto-derive from nterms + main_var degree + 4 (Schwartz-Zippel-style).
    // Exceeded → explicit Unimplemented diagnostic ("no good evaluation point").
    [[nodiscard]] std::size_t max_wang_eval_radius() const noexcept {
        return max_wang_eval_radius_;
    }
    void set_max_wang_eval_radius(std::size_t r) noexcept {
        max_wang_eval_radius_ = r;
    }

    // ── Wang EEZ Hensel lifting attempts (WE-1) ──────────────────────────────
    [[nodiscard]] std::size_t max_hensel_lift_attempts() const noexcept {
        return max_hensel_lift_attempts_;
    }
    void set_max_hensel_lift_attempts(std::size_t attempts) noexcept {
        max_hensel_lift_attempts_ = attempts;
    }

    // ── Kronecker factorization max degree (WE-1) ────────────────────────────
    [[nodiscard]] std::size_t kronecker_max_degree() const noexcept {
        return kronecker_max_degree_;
    }
    void set_kronecker_max_degree(std::size_t degree) noexcept {
        kronecker_max_degree_ = degree;
    }

    // ── F4 Macaulay matrix caps ──────────────────────────────────────────────
    // Hardware guards; exceeded → F4 step falls back to Buchberger.
    void set_f4_max_macaulay_rows(std::size_t n) noexcept {
        f4_max_macaulay_rows_ = n;
    }
    [[nodiscard]] std::size_t f4_max_macaulay_rows() const noexcept {
        return f4_max_macaulay_rows_;
    }

    void set_f4_max_macaulay_monomials(std::size_t n) noexcept {
        f4_max_macaulay_monomials_ = n;
    }
    [[nodiscard]] std::size_t f4_max_macaulay_monomials() const noexcept {
        return f4_max_macaulay_monomials_;
    }

    void set_f4_max_pending_monomials(std::size_t n) noexcept {
        f4_max_pending_monomials_ = n;
    }
    [[nodiscard]] std::size_t f4_max_pending_monomials() const noexcept {
        return f4_max_pending_monomials_;
    }

    // ── FGLM dimension cap ───────────────────────────────────────────────────
    // Maximum standard-monomial dimension D = dim_Q(Q[x]/I) for FGLM order
    // conversion. Exceeded → explicit Unimplemented diagnostic (not silent truncation).
    // Default 512 handles 0-dimensional ideals up to 512 solutions over Q̄.
    void set_fglm_max_dimension(std::size_t n) noexcept {
        fglm_max_dimension_ = n;
    }
    [[nodiscard]] std::size_t fglm_max_dimension() const noexcept {
        return fglm_max_dimension_;
    }

    // ── Integration-by-parts depth ───────────────────────────────────────────
    // Default 8; typical IBP converges in 1-3 levels.
    void set_max_integrate_by_parts_depth(std::size_t n) noexcept {
        max_integrate_by_parts_depth_ = n;
    }
    [[nodiscard]] std::size_t max_integrate_by_parts_depth() const noexcept {
        return max_integrate_by_parts_depth_;
    }

    // ── Risch DE rational ansatz degree (BUG-HANG-001) ──────────────────────
    // Default 32 covers all common elementary integrals.
    void set_max_risch_rational_ansatz_degree(std::size_t n) noexcept {
        max_risch_rational_ansatz_degree_ = n;
    }
    [[nodiscard]] std::size_t max_risch_rational_ansatz_degree() const noexcept {
        return max_risch_rational_ansatz_degree_;
    }

    // ── LLL quality parameter ────────────────────────────────────────────────
    // Default 0.75 (LLL standard, Lenstra-Lenstra-Lovász 1982).
    void set_lll_delta(double delta) noexcept { lll_delta_ = delta; }
    [[nodiscard]] double lll_delta() const noexcept { return lll_delta_; }

    // ── Van Hoeij factorization knobs ────────────────────────────────────────
    // Modular factor count above which van Hoeij LLL replaces subset enumeration.
    // Default 8: 2^8=256 subset budget; van Hoeij O(r^4·n·log(p^a)) for r ≥ 8.
    void set_van_hoeij_threshold(std::size_t t) noexcept {
        van_hoeij_threshold_ = t;
    }
    [[nodiscard]] std::size_t van_hoeij_threshold() const noexcept {
        return van_hoeij_threshold_;
    }

    // Max refinement iterations in van Hoeij LLL knapsack. 0 = auto (4*r).
    void set_max_van_hoeij_iterations(std::size_t n) noexcept {
        max_van_hoeij_iterations_ = n;
    }
    [[nodiscard]] std::size_t max_van_hoeij_iterations() const noexcept {
        return max_van_hoeij_iterations_;
    }

    // Fast-path subset threshold within van Hoeij.
    // Default 10: C(10,5)=252 tractable with Newton-sum Mignotte pruning.
    // r > threshold → mandatory LLL path.
    void set_van_hoeij_lll_threshold(std::size_t t) noexcept {
        van_hoeij_lll_threshold_ = t;
    }
    [[nodiscard]] std::size_t van_hoeij_lll_threshold() const noexcept {
        return van_hoeij_lll_threshold_;
    }

    // ── Timeout check interval ────────────────────────────────────────────────
    // (setter with clamping declared in CASContext / context_core.cpp)
    [[nodiscard]] std::uint64_t timeout_check_interval() const noexcept {
        return timeout_check_interval_;
    }

    // ── F5 signature-based S-pair pruning (F3.3-F5-WIRE) ─────────────────────
    // When true, f4_groebner routes through f5c_groebner which applies the F5
    // criterion (Faugère 2002) + Rewritten criterion to S-pair selection,
    // skipping S-pairs whose signature is divisible by a known syzygy signature.
    // Default false: existing F4/Buchberger pipeline unchanged.
    void set_enable_f5_signature_pruning(bool b) noexcept {
        enable_f5_signature_pruning_ = b;
    }
    [[nodiscard]] bool enable_f5_signature_pruning() const noexcept {
        return enable_f5_signature_pruning_;
    }

    // ── Galois Frobenius / Dedekind prime budget ─────────────────────────────
    // Number of small primes tested for cycle-type evidence in deg-5 Galois
    // identification (Soicher-McKay / Frobenius / Dedekind path).
    // Default 30 = enough to distinguish C5/D5/F20/A5/S5 with high reliability
    // for typical small-coefficient quintics (Chebotarev density argument:
    // expected hit on each cycle type within ~|G|/|class| samples; |G| ≤ 120).
    // Exceeding budget without conclusive evidence → explicit Unimplemented
    // diagnostic (NOT silent wrong group).
    void set_max_galois_frobenius_primes(std::size_t n) noexcept {
        max_galois_frobenius_primes_ = n;
    }
    [[nodiscard]] std::size_t max_galois_frobenius_primes() const noexcept {
        return max_galois_frobenius_primes_;
    }

    // ── Smith Normal Form stabilization ─────────────────────────────────────
    // multiplier used to calculate the stabilization guard (rows + cols) * multiplier.
    // Default 64.
    void set_smith_stabilization_multiplier(std::size_t n) noexcept {
        smith_stabilization_multiplier_ = n;
    }
    [[nodiscard]] std::size_t smith_stabilization_multiplier() const noexcept {
        return smith_stabilization_multiplier_;
    }

    // ── Sparse interpolation (Zippel) retry budget ──────────────────────────
    // Maximum number of distinct prime-offset retries when the candidate
    // skeleton evaluation matrix is singular. Each retry shifts the prime
    // stream by a deterministic offset (no randomness). Exceeded → explicit
    // Unimplemented diagnostic (NOT silent wrong polynomial).
    // Default 5: empirical safe bound for skeletons of moderate size (T ≤ 64).
    void set_sparse_interp_max_retries(std::size_t n) noexcept {
        sparse_interp_max_retries_ = n;
    }
    [[nodiscard]] std::size_t sparse_interp_max_retries() const noexcept {
        return sparse_interp_max_retries_;
    }

    // ── Zippel error probability (ZP-1) ──────────────────────────────────────
    [[nodiscard]] double zippel_error_probability() const noexcept {
        return zippel_error_probability_;
    }
    void set_zippel_error_probability(double prob) noexcept {
        zippel_error_probability_ = prob;
    }

    // ── Zippel density threshold (ZP-3) ──────────────────────────────────────
    [[nodiscard]] double zippel_density_threshold() const noexcept {
        return zippel_density_threshold_;
    }
    void set_zippel_density_threshold(double t) noexcept {
        zippel_density_threshold_ = t;
    }

    // ── sqrt rational simplification ────────────────────────────────────────
    // Trial-division upper bound for perfect-square factor extraction in
    // simplify(sqrt(rational)). Without a bound the O(sqrt(n)) loop is
    // effectively infinite for big numerators (>10²⁰).
    // After trial-division, integer_sqrt provides the perfect-square fallback.
    // Default 10000 (covers all squarefull factors ≤ 10⁸).
    void set_simplify_sqrt_trial_division_bound(std::size_t n) noexcept {
        simplify_sqrt_trial_division_bound_ = n;
    }
    [[nodiscard]] std::size_t simplify_sqrt_trial_division_bound() const noexcept {
        return simplify_sqrt_trial_division_bound_;
    }

    // F8.0-6.2: when true, ln(exp(z)) reduces to z + 2πi·K(z) on z whose
    // realness is not asserted by Assumptions; when false (default), the
    // legacy reduction ln(exp(z)) → z applies regardless of branch. Opt-in
    // because turning this on by default would change the canonical form
    // of every existing log/exp simplification in the codebase.
    void set_branch_cut_aware_logexp(bool on) noexcept {
        branch_cut_aware_logexp_ = on;
    }
    [[nodiscard]] bool branch_cut_aware_logexp() const noexcept {
        return branch_cut_aware_logexp_;
    }

    // ── limit log/log fast-path recursion depth (HPP-022 closure) ──────────
    // Maximum recursion depth for try_log_log_limit (ln(a)/ln(b) fast-path).
    // Bound: max depth ≈ log_2(nesting depth of nested logs in input).
    // Default 3 (sufficient for all known ln/ln cases; never silent wrong
    // answer — exceeding returns nullopt → caller tries other strategies).
    void set_max_log_log_limit_depth(unsigned int d) noexcept {
        max_log_log_limit_depth_ = d;
    }
    [[nodiscard]] unsigned int max_log_log_limit_depth() const noexcept {
        return max_log_log_limit_depth_;
    }

    // ── Adaptive Gauss-Kronrod numerical integration (F6.D) ────────────────
    // Tolerances and resource caps for the priority-queue G7/K15 adaptive
    // integrator (src/numeric/adaptive_integration.cpp). Convergence is
    // declared when  estimated_error ≤ max(abs_tol, rel_tol·|integral|).
    // No silent truncation: when max_intervals is reached without convergence
    // the integrator returns success=false with the partial estimate.
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

    // ── Bessel half-integer recurrence expansion bound (F7.5.E2) ───────────
    // Max |numerator| of rational order p/2 (p odd) for which the half-integer
    // recurrence J/Y/I/K(p/2, x) → elementary closed form is expanded directly.
    // Beyond this bound an Unimplemented diagnostic is returned (no silent
    // truncation). Recurrence base values: ±1/2; iteration linear in |p|.
    // Default 64: covers all practical engineering uses while avoiding very
    // deep nested expressions on adversarial inputs. Reference: spec
    // Bessel_Identities.md §"Mezzi-interi".
    void set_max_bessel_half_integer_order(unsigned int n) noexcept {
        max_bessel_half_integer_order_ = n;
    }
    [[nodiscard]] unsigned int max_bessel_half_integer_order() const noexcept {
        return max_bessel_half_integer_order_;
    }

    // ── Pollard Rho iteration budget (HPP-021 closure) ─────────────────────
    // Maximum iterations per (seed, constant) pair in Pollard Rho.
    // Hardware-safety limit: bound terminates loop with Unimplemented (outer
    // retries with next seed/constant) rather than looping forever.
    // For n with very large factors, rho may need O(n^{1/4}) >> default.
    // Default 4096; configurable in CASContext-aware callers via overload.
    void set_pollard_rho_max_iter(std::size_t n) noexcept {
        pollard_rho_max_iter_ = n;
    }
    [[nodiscard]] std::size_t pollard_rho_max_iter() const noexcept {
        return pollard_rho_max_iter_;
    }

    // ── Special-function integer-argument bit budget (HPP-015 closure) ─────
    // Maximum bit_length of a positive integer argument n accepted by closed-
    // form expansions of Gamma/Digamma/Polygamma/Pochhammer that materialize
    // O(n) AST nodes. Exceeded → explicit Unimplemented diagnostic citing
    // ctx.set_max_special_fn_integer_arg_bits(b). Default 16 (n ≤ 65535).
    // Hardware-safety limit: prevents OOM via runaway AST construction.
    // Increasing past ~24 risks multi-GB allocations on Digamma(2^25).
    void set_max_special_fn_integer_arg_bits(unsigned int b) noexcept {
        max_special_fn_integer_arg_bits_ = b;
    }
    [[nodiscard]] unsigned int max_special_fn_integer_arg_bits() const noexcept {
        return max_special_fn_integer_arg_bits_;
    }

    // ── F5.5 Puiseux multiplicity-decision iteration budget ────────────────
    // Maximum number of repeated-derivative + substitute + simplify rounds
    // used by `puiseux_leading_terms` to decide the multiplicity of a root
    // in the characteristic polynomial.  Exceeded → fall back to the sound
    // lower bound of 1 (NEVER a silently inflated count).
    // Default 32 covers any well-conditioned characteristic; pathological
    // higher orders can be raised here without recompilation.
    void set_max_puiseux_multiplicity_iterations(unsigned int n) noexcept {
        max_puiseux_multiplicity_iterations_ = n;
    }
    [[nodiscard]] unsigned int max_puiseux_multiplicity_iterations() const noexcept {
        return max_puiseux_multiplicity_iterations_;
    }

    // ── F5.6 Residue-theorem numeric driver precision ──────────────────────
    // Working MPFR precision (in decimal digits) of the Aberth root isolator
    // when invoked from `integrate_rational_full_real_line` for irreducible
    // denominator factors of deg ≥ 5 or non-biquadratic quartics.  The
    // residue formula  Res = N(z)/D'(z)  evaluates a sum of complex
    // residues whose imaginary direction suffers catastrophic cancellation
    // when the integrand has real coefficients; the precision must absorb
    // that loss without polluting the to_double() final emission.
    // Default 80 decimal digits (≈ 265 MPFR bits) — empirical safe bound
    // for polynomials of moderate Mahler measure.  Increase for denominators
    // with large root-magnitude spread.
    void set_residue_aberth_precision_digits(unsigned int d) noexcept {
        residue_aberth_precision_digits_ = d;
    }
    [[nodiscard]] unsigned int residue_aberth_precision_digits() const noexcept {
        return residue_aberth_precision_digits_;
    }

    // ── F5.6 Residue-theorem numeric driver iteration budget ───────────────
    // Maximum Aberth main-loop iterations used by the residue driver.
    // Exceeded → explicit Unimplemented diagnostic recommending a precision
    // bump.  Default 500 (≈ cubic Aberth convergence reaches working
    // precision in O(log_2 precision_bits) iterations + safety margin).
    void set_residue_aberth_max_iterations(unsigned int n) noexcept {
        residue_aberth_max_iterations_ = n;
    }
    [[nodiscard]] unsigned int residue_aberth_max_iterations() const noexcept {
        return residue_aberth_max_iterations_;
    }

    // ── Zeilberger creative telescoping (F5.7) ─────────────────────────────
    // Maximum recurrence order J to try (J=1 covers most single-binomial sums;
    // J=2 covers Vandermonde/convolution).
    // Default 2; raising to 3+ enables WZ-provable identities at higher cost.
    void set_max_zeilberger_order(unsigned int j) noexcept {
        max_zeilberger_order_ = (j < 1U) ? 1U : j;
    }
    [[nodiscard]] unsigned int max_zeilberger_order() const noexcept {
        return max_zeilberger_order_;
    }

    // Maximum degree D in n for the operator polynomial coefficients p_i(n).
    // For most proper hypergeometric sums, D=1 suffices.
    // Default 2; D=0 handles constant-coefficient recurrences (e.g., Σ 1 = n).
    void set_max_zeilberger_poly_degree(unsigned int d) noexcept {
        max_zeilberger_poly_degree_ = d;
    }
    [[nodiscard]] unsigned int max_zeilberger_poly_degree() const noexcept {
        return max_zeilberger_poly_degree_;
    }

    // Maximum degree of the Zeilberger certificate R(k) numerator in k.
    // Derivation: deg(R) ≤ J + max(deg(N_r), deg(D_r)) where N_r/D_r = r(k).
    // Default 4; auto-derived from r(k) when 0.
    void set_max_zeilberger_cert_degree(unsigned int d) noexcept {
        max_zeilberger_cert_degree_ = d;
    }
    [[nodiscard]] unsigned int max_zeilberger_cert_degree() const noexcept {
        return max_zeilberger_cert_degree_;
    }

    // ── Bernoulli/Zeta integer-index bit budget (HPP-015 family) ───────────
    // Maximum bit_length of |n| accepted by Zeta(n) closed-form via Bernoulli
    // numbers (positive even or negative odd). Exceeded → explicit Unimplemented
    // diagnostic. Default 30 (|n| ≤ 2^30 ≈ 10^9). Limit: Bernoulli numerator
    // has Θ(n log n) bits, so n ≈ 2^30 gives ~30 GB rationals at upper bound.
    // Increasing risks BigInt OOM in numtheory::bernoulli_number(n).
    void set_max_bernoulli_index_bits(unsigned int b) noexcept {
        max_bernoulli_index_bits_ = b;
    }
    [[nodiscard]] unsigned int max_bernoulli_index_bits() const noexcept {
        return max_bernoulli_index_bits_;
    }

    // ── Swell Guard (F7.0-A3.4) ──────────────────────────────────────────────
    // Maximum number of distinct monomials algebra::expand_power is allowed
    // to produce. Triggers a structured BudgetExceeded-style error (reported
    // as CASErrorKind::Overflow with explicit message) rather than letting
    // the allocator drag the process into OOM-Killer territory.
    //
    // Default 100_000: typical CAS corpus stays well under (a 5-variable
    // poly of degree 20 has C(24,4)=10626 monomials). User can raise for
    // research workloads via set_max_expand_monomials.
    //
    // Estimation upper bound used: n^k for (s_1 + ... + s_n)^k. Conservative
    // — actual distinct monomials = C(n+k-1, n-1).
    void set_max_expand_monomials(std::uint64_t n) noexcept {
        if (n > 0) max_expand_monomials_ = n;
    }
    [[nodiscard]] std::uint64_t max_expand_monomials() const noexcept {
        return max_expand_monomials_;
    }

    // ── solve_inequality (F6.4 / F7.0-A2.1) ──────────────────────────────────
    // Search window half-width used by the Sturm-based univariate inequality
    // solver. Acts as conservative Cauchy bound when the exact bound from
    // poly coefficients exceeds it. Default 10^3 covers typical CAS corpus
    // (poly with rational coefficients of moderate magnitude).
    void set_solve_inequality_search_half_width(long long w) noexcept {
        if (w > 0) solve_inequality_search_half_width_ = w;
    }
    [[nodiscard]] long long solve_inequality_search_half_width() const noexcept {
        return solve_inequality_search_half_width_;
    }
    // Sturm bisection tolerance expressed as inverse denominator: tol = 1/D.
    // Default D = 10^9 (tol = 10^-9) — Mahler-separation safe for deg ≤ 10,
    // H ≤ 100 (Mahler measure root separation lower bound ~10^-12).
    void set_solve_inequality_sturm_tolerance_inv(long long d) noexcept {
        if (d > 0) solve_inequality_sturm_tolerance_inv_ = d;
    }
    [[nodiscard]] long long solve_inequality_sturm_tolerance_inv() const noexcept {
        return solve_inequality_sturm_tolerance_inv_;
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

    void set_fsolve_tolerance_bits(unsigned int bits) noexcept {
        fsolve_tolerance_bits_ = bits;
    }
    [[nodiscard]] unsigned int fsolve_tolerance_bits() const noexcept {
        return fsolve_tolerance_bits_;
    }



protected:
    // All fields with their mathematically-derived or documented defaults.
    // Setters with clamping are CASContext out-of-line methods (context_core.cpp).
    int           max_simplification_depth_{300};
    std::size_t   max_integration_depth_{16U};
    double        gcd_error_probability_{0.001};
    unsigned int  numeric_precision_digits_{15U};
    bool          strict_branch_cuts_{false};
    long long     max_trig_power_reduction_{32LL};
    long long     symbolic_qr_max_norm_complexity_{2LL};
    int           max_trig_exact_denom_{100};
    std::size_t   max_rootof_explicit_degree_{2U};
    std::size_t   max_gcd_recursion_depth_{16U};
    std::size_t   min_gcd_division_steps_{8U};
    std::size_t   max_gcd_total_calls_{4096U};
    int           max_cyclotomic_n_{-1};
    std::size_t   half_gcd_degree_threshold_{200U};
    std::size_t   modular_gcd_coeff_bits_{48U};   // B2.1: CRT GCD dispatch threshold
    std::size_t   max_berlekamp_matrix_size_{1024U};
    std::size_t   max_q_alpha_bridge_depth_{256U};
    std::size_t   max_gamma_recursion_{1024U};
    std::size_t   improper_leading_order_scan_{8U};
    bool          expand_bessel_recurrence_{false};
    unsigned int  max_bessel_half_integer_order_{64U};  // F7.5.E2 half-int recurrence cap
    // F6.D adaptive G7/K15 integrator tolerances (no hardcode allowed in algo).
    double        integration_abs_tol_{1e-10};
    double        integration_rel_tol_{1e-8};
    std::size_t   integration_max_intervals_{4096U};
    std::size_t   max_trager_tower_shift_attempts_{0U};
    std::size_t   max_wang_eval_radius_{0U};  // 0 = auto (nterms + main_deg + 4)
    std::size_t   max_hensel_lift_attempts_{8U};
    std::size_t   kronecker_max_degree_{8U};
    double        zippel_error_probability_{1e-6};
    double        zippel_density_threshold_{5.0};
    std::size_t   f4_max_macaulay_rows_{512U};
    std::size_t   f4_max_macaulay_monomials_{512U};
    std::size_t   f4_max_pending_monomials_{1024U};
    std::size_t   fglm_max_dimension_{512U};
    std::size_t   max_integrate_by_parts_depth_{8U};
    std::size_t   max_risch_rational_ansatz_degree_{32U};
    double        lll_delta_{0.75};
    std::size_t   van_hoeij_threshold_{8U};
    std::size_t   max_van_hoeij_iterations_{0U};  // 0 = auto (4*r)
    std::size_t   van_hoeij_lll_threshold_{10U};  // r > this → mandatory LLL
    std::uint64_t timeout_check_interval_{1024U};
    std::size_t   max_galois_frobenius_primes_{30U};
    std::size_t   smith_stabilization_multiplier_{64U};
    std::size_t   simplify_sqrt_trial_division_bound_{10000U};
    bool          branch_cut_aware_logexp_{false}; // F8.0-6.2 opt-in
    std::size_t   sparse_interp_max_retries_{5U};
    bool          enable_f5_signature_pruning_{false};  // F3.3-F5-WIRE
    unsigned int  max_log_log_limit_depth_{3U};          // HPP-022 closure
    std::size_t   pollard_rho_max_iter_{4096U};          // HPP-021 closure
    unsigned int  max_special_fn_integer_arg_bits_{16U}; // HPP-015 closure
    unsigned int  max_bernoulli_index_bits_{30U};        // HPP-015 closure (Zeta)
    unsigned int  max_puiseux_multiplicity_iterations_{32U};   // F5.5
    unsigned int  residue_aberth_precision_digits_{80U};       // F5.6
    unsigned int  residue_aberth_max_iterations_{500U};        // F5.6
    unsigned int  max_zeilberger_order_{2U};          // F5.7
    unsigned int  max_zeilberger_poly_degree_{2U};    // F5.7
    unsigned int  max_zeilberger_cert_degree_{4U};    // F5.7
    long long     solve_inequality_search_half_width_{1000LL};       // F7.0-A2.1
    long long     solve_inequality_sturm_tolerance_inv_{1000000000LL}; // F7.0-A2.1
    std::uint64_t max_expand_monomials_{100000ULL};                  // F7.0-A3.4 Swell Guard
    unsigned int  mrv_max_append_depth_{1024U};                      // HPP-F75-AUDIT-CYCLE-GUARD-1
    unsigned int  diff_field_max_visit_depth_{4096U};                // HPP-F75-AUDIT-CYCLE-GUARD-2
    int           mrv_growth_rank_max_depth_{1024};                  // HPP-F75-AUDIT-CYCLE-GUARD-3
    unsigned int  fsolve_tolerance_bits_{80U};                       // HPP-006 (Sturm fsolve tolerance)
};

}  // namespace cas::symbolic
