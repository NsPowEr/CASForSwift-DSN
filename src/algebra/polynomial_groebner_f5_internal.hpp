#ifndef CAS_ALGEBRA_POLYNOMIAL_GROEBNER_F5_INTERNAL_HPP
#define CAS_ALGEBRA_POLYNOMIAL_GROEBNER_F5_INTERNAL_HPP

#include "polynomial_groebner_f4.hpp"
#include "polynomial_groebner_f4_internal.hpp"

#include <algorithm>
#include <cstddef>

// Monomial helpers shared (T-049 anti-monolith split) between the F5C driver
// in polynomial_groebner_f5.cpp and the baseline Buchberger in
// polynomial_groebner_f5_buchberger.cpp.
namespace cas::algebra::f5_detail {

[[nodiscard]] inline bool divides_mon(const Monomial& a, const Monomial& b) {
    for (std::size_t k = 0; k < a.size(); ++k) if (a[k] > b[k]) return false;
    return true;
}

[[nodiscard]] inline Monomial lcm_mon(const Monomial& a, const Monomial& b) {
    Monomial r(a.size());
    for (std::size_t k = 0; k < a.size(); ++k) r[k] = std::max(a[k], b[k]);
    return r;
}

[[nodiscard]] inline unsigned int total_deg(const Monomial& m) {
    unsigned int d = 0; for (unsigned int e : m) d += e; return d;
}

}  // namespace cas::algebra::f5_detail

#endif  // CAS_ALGEBRA_POLYNOMIAL_GROEBNER_F5_INTERNAL_HPP
