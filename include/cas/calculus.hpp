#pragma once

#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <vector>

namespace cas::calculus {

struct TaylorExpansion {
    ExprPtr polynomial;
    ExprPtr remainder;
    unsigned int computed_order;
};

[[nodiscard]] Result<ExprPtr> diff(
    ExprPtr expr,
    const Symbol& var,
    unsigned int order,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> partial_diff(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> implicit_diff(
    ExprPtr relation,
    const Symbol& dependent,
    const Symbol& independent,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::vector<ExprPtr>> gradient(
    ExprPtr expr,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> jacobian(
    const std::vector<ExprPtr>& funcs,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> hessian(
    ExprPtr expr,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> integrate(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> definite_integral(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> limit(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    LimitDirection dir,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> sum(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> residue(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr pole,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<TaylorExpansion> taylor_series(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr point,
    unsigned int order,
    symbolic::CASContext& ctx);

}  // namespace cas::calculus
