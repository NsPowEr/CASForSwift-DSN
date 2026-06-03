// Project-internal complex-BigFloat helper.  Two MPFR BigFloats glued
// together as (real, imag) — sufficient for the Aberth root isolator and the
// residue-theorem numeric driver.  Not a public type; downstream callers see
// the wrapped result through API surfaces such as `ComplexRoot`.
//
// Kept as a thin inline struct on purpose: the algorithms that use it
// (Aberth, residue) perform tight inner loops where the cost of a function
// call boundary would dominate.  Keep operators trivially inlined.

#pragma once

#include "cas/bigfloat.hpp"

#include <utility>

namespace cas::numeric::detail {

struct CBF {
    BigFloat re;
    BigFloat im;

    CBF() = default;
    CBF(BigFloat r, BigFloat i) : re(std::move(r)), im(std::move(i)) {}

    static CBF zero(mpfr_prec_t prec) {
        return CBF{BigFloat(prec), BigFloat(prec)};
    }
    static CBF from_real(BigFloat r) {
        BigFloat zero_im(r.precision_bits());
        return CBF{std::move(r), std::move(zero_im)};
    }
    static CBF from_double_pair(double r, double i, mpfr_prec_t prec) {
        return CBF{BigFloat::from_double(r, prec), BigFloat::from_double(i, prec)};
    }

    [[nodiscard]] bool is_zero() const noexcept {
        return re.is_zero() && im.is_zero();
    }

    [[nodiscard]] CBF operator+(const CBF& o) const { return CBF{re + o.re, im + o.im}; }
    [[nodiscard]] CBF operator-(const CBF& o) const { return CBF{re - o.re, im - o.im}; }
    [[nodiscard]] CBF operator-() const { return CBF{-re, -im}; }
    [[nodiscard]] CBF operator*(const CBF& o) const {
        return CBF{re * o.re - im * o.im, re * o.im + im * o.re};
    }
    [[nodiscard]] CBF operator/(const CBF& o) const {
        BigFloat denom = o.re * o.re + o.im * o.im;
        return CBF{(re * o.re + im * o.im) / denom,
                   (im * o.re - re * o.im) / denom};
    }

    [[nodiscard]] BigFloat abs() const { return BigFloat::sqrt(re * re + im * im); }
    [[nodiscard]] BigFloat abs_sq() const { return re * re + im * im; }
};

}  // namespace cas::numeric::detail
