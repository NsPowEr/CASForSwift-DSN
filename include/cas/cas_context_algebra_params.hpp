#pragma once

#include <cstddef>
#include <cstdint>

namespace cas::symbolic {

struct CASContextAlgebraParams {
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

    [[nodiscard]] std::size_t max_gcd_total_calls() const noexcept {
        return max_gcd_total_calls_;
    }

    // ── Cyclotomic ───────────────────────────────────────────────────────────
    // (setter trivial; kept in CASContext for symmetry with other GCD setters)
    [[nodiscard]] int max_cyclotomic_n() const noexcept {
        return max_cyclotomic_n_;
    }

    // ── Half-GCD dispatch threshold ──────────────────────────────────────────
    void set_half_gcd_degree_threshold(std::size_t n) noexcept {
        half_gcd_degree_threshold_ = n;
    }
    [[nodiscard]] std::size_t half_gcd_degree_threshold() const noexcept {
        return half_gcd_degree_threshold_;
    }

    // ── Modular GCD CRT dispatch threshold (B2.1) ───────────────────────────
    void set_modular_gcd_coeff_bits(std::size_t bits) noexcept {
        modular_gcd_coeff_bits_ = bits;
    }
    [[nodiscard]] std::size_t modular_gcd_coeff_bits() const noexcept {
        return modular_gcd_coeff_bits_;
    }

    // ── Berlekamp Q-matrix budget ────────────────────────────────────────────
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

    // ── Trager tower shift attempts (L3-06) ──────────────────────────────────
    [[nodiscard]] std::size_t max_trager_tower_shift_attempts() const noexcept {
        return max_trager_tower_shift_attempts_;
    }

    // ── Wang multivariate factor evaluation search radius ────────────────────
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

    // ── Algebraic tower ring-arithmetic bit-cap (safety belt, HC-F8-PENDING-26) ──
    // Guards Q[y]/(cand_q) Euclidean inverse/divmod against unbounded coefficient
    // growth when cand_q turns out reducible (not a field): returns nullopt past
    // this cap instead of looping on ever-larger BigInt coefficients.
    [[nodiscard]] std::size_t algebraic_tower_eval_max_bits() const noexcept {
        return algebraic_tower_eval_max_bits_;
    }
    void set_algebraic_tower_eval_max_bits(std::size_t bits) noexcept {
        algebraic_tower_eval_max_bits_ = bits;
    }

    // ── F4 Macaulay matrix caps ──────────────────────────────────────────────
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
    void set_fglm_max_dimension(std::size_t n) noexcept {
        fglm_max_dimension_ = n;
    }
    [[nodiscard]] std::size_t fglm_max_dimension() const noexcept {
        return fglm_max_dimension_;
    }

    // ── LLL quality parameter ────────────────────────────────────────────────
    void set_lll_delta(double delta) noexcept { lll_delta_ = delta; }
    [[nodiscard]] double lll_delta() const noexcept { return lll_delta_; }

    // ── Van Hoeij factorization knobs ────────────────────────────────────────
    void set_van_hoeij_threshold(std::size_t t) noexcept {
        van_hoeij_threshold_ = t;
    }
    [[nodiscard]] std::size_t van_hoeij_threshold() const noexcept {
        return van_hoeij_threshold_;
    }

    void set_max_van_hoeij_iterations(std::size_t n) noexcept {
        max_van_hoeij_iterations_ = n;
    }
    [[nodiscard]] std::size_t max_van_hoeij_iterations() const noexcept {
        return max_van_hoeij_iterations_;
    }

    void set_van_hoeij_lll_threshold(std::size_t t) noexcept {
        van_hoeij_lll_threshold_ = t;
    }
    [[nodiscard]] std::size_t van_hoeij_lll_threshold() const noexcept {
        return van_hoeij_lll_threshold_;
    }

    // ── Timeout check interval ────────────────────────────────────────────────
    [[nodiscard]] std::uint64_t timeout_check_interval() const noexcept {
        return timeout_check_interval_;
    }

    // ── F5 signature-based S-pair pruning (F3.3-F5-WIRE) ─────────────────────
    void set_enable_f5_signature_pruning(bool b) noexcept {
        enable_f5_signature_pruning_ = b;
    }
    [[nodiscard]] bool enable_f5_signature_pruning() const noexcept {
        return enable_f5_signature_pruning_;
    }

    // ── Galois Frobenius / Dedekind prime budget ─────────────────────────────
    void set_max_galois_frobenius_primes(std::size_t n) noexcept {
        max_galois_frobenius_primes_ = n;
    }
    [[nodiscard]] std::size_t max_galois_frobenius_primes() const noexcept {
        return max_galois_frobenius_primes_;
    }

    // ── Smith Normal Form stabilization ─────────────────────────────────────
    void set_smith_stabilization_multiplier(std::size_t n) noexcept {
        smith_stabilization_multiplier_ = n;
    }
    [[nodiscard]] std::size_t smith_stabilization_multiplier() const noexcept {
        return smith_stabilization_multiplier_;
    }

    // ── Swell Guard ──────────────────────────────────────────────────────────
    void set_max_expand_monomials(std::uint64_t n) noexcept {
        if (n > 0) max_expand_monomials_ = n;
    }
    [[nodiscard]] std::uint64_t max_expand_monomials() const noexcept {
        return max_expand_monomials_;
    }

    // ── solve_inequality ─────────────────────────────────────────────────────
    void set_solve_inequality_search_half_width(long long w) noexcept {
        if (w > 0) solve_inequality_search_half_width_ = w;
    }
    [[nodiscard]] long long solve_inequality_search_half_width() const noexcept {
        return solve_inequality_search_half_width_;
    }
    void set_solve_inequality_sturm_tolerance_inv(long long d) noexcept {
        if (d > 0) solve_inequality_sturm_tolerance_inv_ = d;
    }
    [[nodiscard]] long long solve_inequality_sturm_tolerance_inv() const noexcept {
        return solve_inequality_sturm_tolerance_inv_;
    }

    // ── fsolve ───────────────────────────────────────────────────────────────
    void set_fsolve_tolerance_bits(unsigned int bits) noexcept {
        fsolve_tolerance_bits_ = bits;
    }
    [[nodiscard]] unsigned int fsolve_tolerance_bits() const noexcept {
        return fsolve_tolerance_bits_;
    }
    void set_fsolve_tolerance(double tol) noexcept {
        if (tol > 0.0) fsolve_tolerance_ = tol;
    }
    [[nodiscard]] double fsolve_tolerance() const noexcept {
        return fsolve_tolerance_;
    }

    // ── Swell Guard & Galois/Zeta Budgets ────────────────────────────────────
    void set_sparse_interp_max_retries(std::size_t n) noexcept {
        sparse_interp_max_retries_ = n;
    }
    [[nodiscard]] std::size_t sparse_interp_max_retries() const noexcept {
        return sparse_interp_max_retries_;
    }

    [[nodiscard]] double zippel_error_probability() const noexcept {
        return zippel_error_probability_;
    }
    void set_zippel_error_probability(double prob) noexcept {
        zippel_error_probability_ = prob;
    }

    [[nodiscard]] double zippel_density_threshold() const noexcept {
        return zippel_density_threshold_;
    }
    void set_zippel_density_threshold(double t) noexcept {
        zippel_density_threshold_ = t;
    }

    void set_pollard_rho_max_iter(std::size_t n) noexcept {
        pollard_rho_max_iter_ = n;
    }
    [[nodiscard]] std::size_t pollard_rho_max_iter() const noexcept {
        return pollard_rho_max_iter_;
    }

    void set_max_special_fn_integer_arg_bits(unsigned int b) noexcept {
        max_special_fn_integer_arg_bits_ = b;
    }
    [[nodiscard]] unsigned int max_special_fn_integer_arg_bits() const noexcept {
        return max_special_fn_integer_arg_bits_;
    }

    void set_max_bernoulli_index_bits(unsigned int b) noexcept {
        max_bernoulli_index_bits_ = b;
    }
    [[nodiscard]] unsigned int max_bernoulli_index_bits() const noexcept {
        return max_bernoulli_index_bits_;
    }

    void set_max_puiseux_multiplicity_iterations(unsigned int n) noexcept {
        max_puiseux_multiplicity_iterations_ = n;
    }
    [[nodiscard]] unsigned int max_puiseux_multiplicity_iterations() const noexcept {
        return max_puiseux_multiplicity_iterations_;
    }

    void set_residue_aberth_precision_digits(unsigned int d) noexcept {
        residue_aberth_precision_digits_ = d;
    }
    [[nodiscard]] unsigned int residue_aberth_precision_digits() const noexcept {
        return residue_aberth_precision_digits_;
    }

    void set_residue_aberth_max_iterations(unsigned int n) noexcept {
        residue_aberth_max_iterations_ = n;
    }
    [[nodiscard]] unsigned int residue_aberth_max_iterations() const noexcept {
        return residue_aberth_max_iterations_;
    }

    // ── Zeilberger creative telescoping (F5.7) ─────────────────────────────
    void set_max_zeilberger_order(unsigned int j) noexcept {
        max_zeilberger_order_ = (j < 1U) ? 1U : j;
    }
    [[nodiscard]] unsigned int max_zeilberger_order() const noexcept {
        return max_zeilberger_order_;
    }

    void set_max_zeilberger_poly_degree(unsigned int d) noexcept {
        max_zeilberger_poly_degree_ = d;
    }
    [[nodiscard]] unsigned int max_zeilberger_poly_degree() const noexcept {
        return max_zeilberger_poly_degree_;
    }

    void set_max_zeilberger_cert_degree(unsigned int d) noexcept {
        max_zeilberger_cert_degree_ = d;
    }
    [[nodiscard]] unsigned int max_zeilberger_cert_degree() const noexcept {
        return max_zeilberger_cert_degree_;
    }

    // Petkovšek Hyper: cap on the number of monic divisors enumerated for the
    // leading/trailing coefficient (bounds the (a,b) candidate search).
    void set_max_hyper_divisors(unsigned int n) noexcept {
        max_hyper_divisors_ = (n < 1U) ? 1U : n;
    }
    [[nodiscard]] unsigned int max_hyper_divisors() const noexcept {
        return max_hyper_divisors_;
    }

protected:
    double        gcd_error_probability_{0.001};
    unsigned int  numeric_precision_digits_{15U};
    std::size_t   max_rootof_explicit_degree_{2U};
    std::size_t   max_gcd_recursion_depth_{16U};
    std::size_t   min_gcd_division_steps_{8U};
    std::size_t   max_gcd_total_calls_{4096U};
    int           max_cyclotomic_n_{-1};
    std::size_t   half_gcd_degree_threshold_{200U};
    std::size_t   modular_gcd_coeff_bits_{48U};
    std::size_t   max_berlekamp_matrix_size_{1024U};
    std::size_t   max_q_alpha_bridge_depth_{256U};
    std::size_t   max_gamma_recursion_{1024U};
    std::size_t   improper_leading_order_scan_{8U};
    std::size_t   max_trager_tower_shift_attempts_{0U};
    std::size_t   max_wang_eval_radius_{0U};
    std::size_t   max_hensel_lift_attempts_{8U};
    std::size_t   kronecker_max_degree_{8U};
    double        zippel_error_probability_{1e-6};
    std::size_t   algebraic_tower_eval_max_bits_{8192U};
    double        zippel_density_threshold_{5.0};
    std::size_t   f4_max_macaulay_rows_{512U};
    std::size_t   f4_max_macaulay_monomials_{512U};
    std::size_t   f4_max_pending_monomials_{1024U};
    std::size_t   fglm_max_dimension_{512U};
    double        lll_delta_{0.75};
    std::size_t   van_hoeij_threshold_{8U};
    std::size_t   max_van_hoeij_iterations_{0U};
    std::size_t   van_hoeij_lll_threshold_{10U};
    std::uint64_t timeout_check_interval_{1024U};
    std::size_t   max_galois_frobenius_primes_{30U};
    std::size_t   smith_stabilization_multiplier_{64U};
    std::size_t   sparse_interp_max_retries_{5U};
    bool          enable_f5_signature_pruning_{false};
    std::size_t   pollard_rho_max_iter_{4096U};
    unsigned int  max_special_fn_integer_arg_bits_{16U};
    unsigned int  max_bernoulli_index_bits_{30U};
    unsigned int  max_puiseux_multiplicity_iterations_{32U};
    unsigned int  residue_aberth_precision_digits_{80U};
    unsigned int  residue_aberth_max_iterations_{500U};
    unsigned int  max_zeilberger_order_{2U};
    unsigned int  max_zeilberger_poly_degree_{2U};
    unsigned int  max_zeilberger_cert_degree_{4U};
    unsigned int  max_hyper_divisors_{64U};
    long long     solve_inequality_search_half_width_{1000LL};
    long long     solve_inequality_sturm_tolerance_inv_{1000000000LL};
    std::uint64_t max_expand_monomials_{100000ULL};
    unsigned int  fsolve_tolerance_bits_{80U};
    double        fsolve_tolerance_{1e-10};
};

}  // namespace cas::symbolic
