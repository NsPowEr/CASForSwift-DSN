#pragma once

#include "integrate_definite_patterns.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/error.hpp"
#include "cas/extended_real.hpp"

#include <optional>
#include <vector>

namespace cas::calculus {

struct PolyProductMatch {
    BigInt m;
    BigInt n;
    std::vector<ExprPtr> other_factors;
};

void flatten_mul_factors(ExprPtr expr, std::vector<ExprPtr>& out);

[[nodiscard]] bool depends_on_var(ExprPtr expr, const Symbol& var);

[[nodiscard]] std::optional<BigInt> match_poly_call(
    ExprPtr expr, BuiltinOp op, const Symbol& var);

[[nodiscard]] bool is_literal_rational(ExprPtr expr, long long num, long long den);

[[nodiscard]] bool is_literal_neg_one(ExprPtr expr);

[[nodiscard]] std::optional<PolyProductMatch> match_two_poly_product(
    ExprPtr integrand_normalized,
    BuiltinOp op,
    const Symbol& var);

[[nodiscard]] Result<std::optional<ExprPtr>> apply_constant_factors(
    ExprPtr base,
    const std::vector<ExprPtr>& others,
    const Symbol& var,
    symbolic::CASContext& ctx);

}  // namespace cas::calculus
