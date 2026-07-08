#include "cas/symbolic.hpp"
#include <chrono>
#include <cstdint>

namespace cas::symbolic {

void CASContext::set_timeout_check_interval(std::uint64_t interval) noexcept {
    timeout_check_interval_ = (interval < 64U) ? 64U : interval;
}

void CASContext::set_max_simplification_depth(int depth) noexcept {
    max_simplification_depth_ = (depth < 10) ? 10 : depth;
}

void CASContext::set_max_recursion_depth(std::size_t depth) noexcept {
    max_recursion_depth_ = (depth < 1U) ? 1U : depth;
}

void CASContext::set_max_operation_ops(std::uint64_t ops) noexcept {
    // 0 = disabled (wall-clock safety net only); any positive value accepted.
    max_operation_ops_ = ops;
    max_operation_ops_explicit_ = true;
}

void CASContext::set_max_integration_depth(std::size_t depth) noexcept {
    max_integration_depth_ = (depth < 1U) ? 1U : (depth > 128U) ? 128U : depth;
}

void CASContext::set_gcd_error_probability(double prob) noexcept {
    if (prob < 1e-6) prob = 1e-6;
    if (prob > 0.1)  prob = 0.1;
    gcd_error_probability_ = prob;
}

void CASContext::set_zippel_error_probability(double prob) noexcept {
    if (prob < 1e-15) prob = 1e-15;
    if (prob > 0.1)   prob = 0.1;
    zippel_error_probability_ = prob;
}

void CASContext::set_zippel_density_threshold(double t) noexcept {
    if (t < 0.0) t = 0.0;
    zippel_density_threshold_ = t;
}

void CASContext::set_numeric_precision_digits(unsigned int digits) noexcept {
    // Clamp: minimum 6 digits (≈ 20 bits MPFR), maximum 10000 (~ 33k bits).
    if (digits < 6U) digits = 6U;
    if (digits > 10000U) digits = 10000U;
    numeric_precision_digits_ = digits;
}

void CASContext::set_max_rootof_explicit_degree(std::size_t deg) noexcept {
    max_rootof_explicit_degree_ = (deg < 1U) ? 1U : deg;
}

void CASContext::set_max_gcd_recursion_depth(std::size_t depth) noexcept {
    max_gcd_recursion_depth_ = (depth < 4U) ? 4U : depth;
}

void CASContext::set_min_gcd_division_steps(std::size_t steps) noexcept {
    min_gcd_division_steps_ = (steps < 1U) ? 1U : steps;
}

void CASContext::set_max_gcd_total_calls(std::size_t n) noexcept {
    max_gcd_total_calls_ = (n < 16U) ? 16U : n;
}

void CASContext::set_max_cyclotomic_n(int n) noexcept {
    max_cyclotomic_n_ = n;
}

void CASContext::set_max_q_alpha_bridge_depth(std::size_t depth) noexcept {
    max_q_alpha_bridge_depth_ = (depth < 8U) ? 8U : depth;
}

void CASContext::set_max_gamma_recursion(std::size_t iters) noexcept {
    max_gamma_recursion_ = (iters < 16U) ? 16U : iters;
}

void CASContext::set_improper_leading_order_scan(std::size_t window) noexcept {
    improper_leading_order_scan_ = (window < 1U) ? 1U : window;
}

void CASContext::set_expand_bessel_recurrence(bool enabled) noexcept {
    expand_bessel_recurrence_ = enabled;
}

void CASContext::set_max_trager_tower_shift_attempts(std::size_t attempts) noexcept {
    max_trager_tower_shift_attempts_ = attempts;
}

}  // namespace cas::symbolic
