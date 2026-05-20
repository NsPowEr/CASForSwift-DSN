// CAS-L3-13 — Interval arithmetic (Moore IA-1966) over BigFloat.
//
// Represents a closed real interval [lo, hi] with lo ≤ hi. Arithmetic
// rounds outward to maintain inclusion: for f(I, J) returns smallest
// interval guaranteed to contain {f(x, y) : x ∈ I, y ∈ J}.
//
// All bounds stored as BigFloat (MPFR) → no double precision loss.
// Rounding mode managed via MPFR rounding flags.

#pragma once

#include "cas/bigfloat.hpp"

namespace cas::numeric {

class Interval {
public:
    Interval() noexcept = default;
    Interval(BigFloat lo, BigFloat hi);

    // Singleton interval [v, v].
    explicit Interval(const BigFloat& v);

    // From rational endpoints at default precision.
    [[nodiscard]] static Interval from_rationals(
        const std::string& lo_num, const std::string& lo_den,
        const std::string& hi_num, const std::string& hi_den,
        mpfr_prec_t prec = BigFloat::DEFAULT_PREC);

    [[nodiscard]] const BigFloat& lo() const noexcept { return lo_; }
    [[nodiscard]] const BigFloat& hi() const noexcept { return hi_; }

    // Width = hi - lo (always ≥ 0 by invariant).
    [[nodiscard]] BigFloat width() const;

    // Midpoint (lo + hi) / 2 — exact in BigFloat arithmetic.
    [[nodiscard]] BigFloat midpoint() const;

    // Contains a single value or another interval.
    [[nodiscard]] bool contains(const BigFloat& v) const;
    [[nodiscard]] bool contains(const Interval& other) const;

    // Bound queries (1 if a constant sign holds across interval).
    [[nodiscard]] bool is_positive() const;     // lo > 0
    [[nodiscard]] bool is_negative() const;     // hi < 0
    [[nodiscard]] bool contains_zero() const;   // lo ≤ 0 ≤ hi

    // Arithmetic (rounding outward).
    [[nodiscard]] Interval operator+(const Interval& other) const;
    [[nodiscard]] Interval operator-(const Interval& other) const;
    [[nodiscard]] Interval operator*(const Interval& other) const;
    // Division: requires NOT other.contains_zero(); otherwise returns
    // [-inf, inf] (Moore convention) — caller checks.
    [[nodiscard]] Interval operator/(const Interval& other) const;
    [[nodiscard]] Interval operator-() const;

    // Transcendental functions.
    [[nodiscard]] static Interval sqrt(const Interval& x);
    [[nodiscard]] static Interval exp(const Interval& x);
    [[nodiscard]] static Interval ln(const Interval& x);   // requires x.lo > 0
    [[nodiscard]] static Interval sin(const Interval& x);
    [[nodiscard]] static Interval cos(const Interval& x);

private:
    BigFloat lo_;
    BigFloat hi_;
};

}  // namespace cas::numeric
