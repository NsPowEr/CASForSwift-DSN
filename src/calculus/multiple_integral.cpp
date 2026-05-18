#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

#include <string>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}

// Returns true iff `expr` contains any free reference to `var`.
// Uses symbolic substitute: replace var → fresh_sentinel, check if result differs.
[[nodiscard]] bool contains_var(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    // Simple structural check: walk the AST looking for the symbol.
    // We delegate to a helper visitor via substitute with a sentinel value.
    // Sentinel: replace var with var+1; if the expression changes, var appeared.
    auto sentinel = ctx.arena().make<Binary>(
        BinaryOp::Add,
        ctx.arena().make<Symbol>(var.name),
        ctx.arena().make<IntegerLit>(BigInt(1)));
    auto substituted = symbolic::substitute(expr, var, sentinel, ctx);
    if (substituted.is_error()) return true; // conservative
    // If expr is structurally identical after substitution, var was not present.
    return substituted.value() != expr;
}

} // namespace

Result<ExprPtr> multiple_integral(
    ExprPtr integrand,
    const std::vector<IntegralSpec>& specs,
    symbolic::CASContext& ctx) {

    if (specs.empty()) {
        return ok(integrand);
    }
    if (!integrand) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "multiple_integral: null integrand"));
    }

    ExprPtr current = integrand;
    for (const IntegralSpec& spec : specs) {
        if (!spec.lower || !spec.upper) {
            return fail<ExprPtr>(make_error(
                CASErrorKind::InvalidArgument,
                "multiple_integral: null bound for variable " + spec.var.name));
        }
        auto result = definite_integral(current, spec.var, spec.lower, spec.upper, ctx);
        if (result.is_error()) return result;
        current = result.value();
    }
    return ok(current);
}

Result<ExprPtr> fubini_swap(
    ExprPtr integrand,
    const Symbol& x, ExprPtr ax, ExprPtr bx,
    const Symbol& y, ExprPtr ay, ExprPtr by,
    symbolic::CASContext& ctx) {

    // Guard: rectangular domain requires bounds to be independent of the other variable.
    if (contains_var(ax, y, ctx) || contains_var(bx, y, ctx)) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::Unimplemented,
            "fubini_swap: bounds of x-integral depend on y — non-rectangular domain not supported"));
    }
    if (contains_var(ay, x, ctx) || contains_var(by, x, ctx)) {
        return fail<ExprPtr>(make_error(
            CASErrorKind::Unimplemented,
            "fubini_swap: bounds of y-integral depend on x — non-rectangular domain not supported"));
    }

    // Compute ∫_ay^by ∫_ax^bx f dx dy  (x inner, y outer)
    auto inner_x = definite_integral(integrand, x, ax, bx, ctx);
    if (inner_x.is_error()) return inner_x;
    auto outer_y = definite_integral(inner_x.value(), y, ay, by, ctx);
    if (outer_y.is_error()) return outer_y;

    // Simplify and return
    auto simplified = ctx.simplify(outer_y.value());
    if (simplified.is_error()) return outer_y;
    return ok(simplified.value());
}

} // namespace cas::calculus
