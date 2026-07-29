#pragma once

#include "cas/rational.hpp"
#include <vector>
#include <optional>

namespace cas::algebra {

[[nodiscard]] std::vector<Rational> mat_vec_mul(
    const std::vector<std::vector<Rational>>& M,
    const std::vector<Rational>& v);

[[nodiscard]] std::optional<std::vector<Rational>> linear_dependency(
    const std::vector<std::vector<Rational>>& basis_vecs,
    const std::vector<Rational>& vec);

} // namespace cas::algebra
