#pragma once

#include "polynomial_groebner_f4.hpp"

namespace cas::algebra::detail {

[[nodiscard]] Result<std::vector<PolyF4>> buchberger_groebner(
    std::vector<PolyF4> basis,
    MonomialOrder order);

}  // namespace cas::algebra::detail
