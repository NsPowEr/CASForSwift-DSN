#pragma once

#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/result.hpp"
#include "calculus_internal.hpp"

#include <optional>

namespace cas::calculus {

// Result of the cut-edge analysis of one branch-cut argument
// (limit_branch_cut.cpp): the approach side of Im(arg) for each direction
// (+1 = top edge / principal, −1 = bottom edge) and −Re(arg(point)) > 0.
struct CutEdgeAnalysis {
    int side_right{0};
    int side_left{0};
    ExprPtr neg_re0{nullptr};
};

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

    // HC-F8-PENDING-20 §3.2 — direction-limit table at branch-cut edges
    // (limit_branch_cut.cpp). nullopt → no cut involvement, legacy path.
    [[nodiscard]] std::optional<Result<ExprPtr>> try_branch_cut_directional(
        ExprPtr expr, const Symbol& var, ExprPtr point, LimitDirection dir,
        unsigned int depth);

    [[nodiscard]] std::optional<CutEdgeAnalysis> analyze_cut_edge(
        ExprPtr arg, const Symbol& var, ExprPtr point);

    [[nodiscard]] Result<ExprPtr> cut_edge_value(
        BuiltinOp op, ExprPtr neg_re0, int side);

private:
    symbolic::CASContext& context_;
    AstArena& arena_;
    unsigned int max_depth_budget_{16U};
};

}  // namespace cas::calculus
