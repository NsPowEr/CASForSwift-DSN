#pragma once

#include "cas/algebra.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include <functional>
#include <vector>

namespace cas::algebra {

/**
 * @brief Oracle function that evaluates the hidden polynomial at a given point.
 * The point is a vector of ExprPtr (usually Rational or BigInt).
 */
using InterpolationOracle = std::function<Result<ExprPtr>(const std::vector<ExprPtr>&)>;

/**
 * @brief Implements Zippel's sparse interpolation algorithm.
 * 
 * @param oracle Black-box that returns P(x1, ..., xn) for given values.
 * @param variables The symbols to interpolate.
 * @param degree_bounds Upper bounds for degrees of each variable.
 * @param ctx CAS context for arena and operations.
 * @return Result<ExprPtr> The interpolated polynomial as an expression.
 */
[[nodiscard]] Result<ExprPtr> sparse_interpolate(
    const InterpolationOracle& oracle,
    const std::vector<Symbol>& variables,
    const std::vector<std::size_t>& degree_bounds,
    symbolic::CASContext& ctx);

} // namespace cas::algebra
