#pragma once

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/result.hpp"
#include "calculus_internal.hpp"

#include <optional>

namespace cas::calculus {

class LimitEngine {
public:
    explicit LimitEngine(symbolic::CASContext& context) noexcept : context_(context), arena_(context.arena()) {}

    [[nodiscard]] Result<ExprPtr> compute(ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir);

    [[nodiscard]] Result<ExprPtr> substitute_and_simplify(ExprPtr expr, const Symbol& var, ExprPtr point);

    [[nodiscard]] Result<ExprPtr> compute_recursive(
        ExprPtr expr,
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth);

    [[nodiscard]] std::optional<Result<ExprPtr>> try_log_log_limit(
        const QuotientView& quotient, const Symbol& var, ExprPtr point,
        LimitDirection dir, unsigned int depth);

    [[nodiscard]] Result<ExprPtr> compute_quotient_limit(
        ExprPtr expr,
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth,
        const QuotientView& quotient,
        const Result<ExprPtr>& direct);

    [[nodiscard]] Result<ExprPtr> compute_lhopital_or_taylor(
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir,
        unsigned int depth,
        const QuotientView& quotient);

    // Signed infinity of a pole numerator/denominator at a finite point, recovered
    // from the reciprocal denominator/numerator → 0 (general; nullopt when the
    // sign cannot be decided exactly).
    [[nodiscard]] std::optional<Result<ExprPtr>> try_signed_pole_via_reciprocal(
        ExprPtr numerator,
        ExprPtr denominator,
        const Symbol& var,
        ExprPtr point,
        LimitDirection dir);

private:
    symbolic::CASContext& context_;
    AstArena& arena_;
    unsigned int max_depth_budget_{16U};
};

}  // namespace cas::calculus
