// algebraic_tower_bridge_internal.hpp — Private declarations for tower bridge.
#pragma once

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_internal.hpp"

namespace cas {
namespace algebra {

class RootOfExplicitDegreeGuard {
public:
    RootOfExplicitDegreeGuard(symbolic::CASContext& ctx, std::size_t temp)
        : ctx_(ctx), saved_(ctx.max_rootof_explicit_degree()) {
        ctx_.set_max_rootof_explicit_degree(temp);
    }

    ~RootOfExplicitDegreeGuard() {
        ctx_.set_max_rootof_explicit_degree(saved_);
    }

private:
    symbolic::CASContext& ctx_;
    std::size_t saved_;
};

[[nodiscard]] Result<ExprPtr> clone_expr_raw(ExprPtr expr, symbolic::CASContext& ctx);

[[nodiscard]] Result<ExprPtr> canonicalize_root_expr(ExprPtr expr, symbolic::CASContext& ctx);

[[nodiscard]] Result<PolyExpr> make_constant_poly_raw(ExprPtr coefficient, symbolic::CASContext& ctx);

[[nodiscard]] Result<PolyExpr> poly_pow_raw(PolyExpr base, std::size_t exponent, symbolic::CASContext& ctx);

[[nodiscard]] Result<PolyExpr> parse_polynomial_raw(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& ctx);

}  // namespace algebra
}  // namespace cas
