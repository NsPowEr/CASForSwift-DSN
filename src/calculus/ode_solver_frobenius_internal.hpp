// L2-01: Frobenius series solution - private declarations and helpers.
#pragma once

#include "cas/ast.hpp"
#include "cas/bigint.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "cas/error.hpp"

#include <vector>
#include <optional>
#include <string>

namespace cas::calculus {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message);

[[nodiscard]] inline ExprPtr make_int(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] BigInt factorial_big(unsigned int n);

[[nodiscard]] bool is_literal_zero(ExprPtr e);

[[nodiscard]] bool contains_undefined_constant(ExprPtr e);

[[nodiscard]] Result<ExprPtr> taylor_coefficient(
    ExprPtr expr,
    const Symbol& x,
    unsigned int k,
    symbolic::CASContext& ctx);

[[nodiscard]] ExprPtr make_x_to_r(ExprPtr r, const Symbol& x, AstArena& arena);

[[nodiscard]] Result<ExprPtr> indicial_value(
    ExprPtr p0,
    ExprPtr q0,
    ExprPtr r_arg,
    symbolic::CASContext& ctx);

[[nodiscard]] std::optional<unsigned int> integer_gap(
    ExprPtr r1, ExprPtr r2, symbolic::CASContext& ctx);

// Series generation functions implemented in ode_solver_frobenius_series.cpp
[[nodiscard]] Result<std::vector<ExprPtr>> compute_recurrence(
    ExprPtr root_r,
    const std::vector<ExprPtr>& p_coeffs,
    const std::vector<ExprPtr>& q_coeffs,
    ExprPtr p0,
    ExprPtr q0,
    unsigned int num_terms,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> build_log_branch(
    ExprPtr r1,
    ExprPtr r2,
    unsigned int N,
    const std::vector<ExprPtr>& a_coeffs,
    const std::vector<ExprPtr>& p_coeffs,
    const std::vector<ExprPtr>& q_coeffs,
    ExprPtr p0,
    ExprPtr q0,
    unsigned int num_terms,
    ExprPtr y_1_series,
    const Symbol& x,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> build_series(
    ExprPtr root_r,
    const std::vector<ExprPtr>& c,
    const Symbol& x,
    symbolic::CASContext& ctx);

}  // namespace cas::calculus
