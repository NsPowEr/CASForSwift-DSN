#include "cas/algebra.hpp"
#include "cas/numeric.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace cas::algebra {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}

// Convert equation to f(x)=0 form: if Binary(=, L, R) → L-R; else identity.
[[nodiscard]] ExprPtr to_zero_form(ExprPtr eq, symbolic::CASContext& ctx) {
    if (const auto* b = expr_cast<Binary>(eq)) {
        if (b->op == BinaryOp::Equal) {
            return ctx.arena().make<Binary>(BinaryOp::Sub, b->left, b->right);
        }
    }
    return eq;
}

// Pack a list of double roots into an (n×1) Matrix of DecimalLit nodes.
[[nodiscard]] ExprPtr roots_to_matrix(const std::vector<double>& roots, AstArena& arena) {
    std::vector<ExprPtr> elems;
    elems.reserve(roots.size());
    for (double r : roots) {
        // Represent each root as a DecimalLit (rounded to 10 significant digits)
        elems.push_back(arena.make<DecimalLit>(r));
    }
    return arena.make<Matrix>(static_cast<int>(roots.size()), 1, std::move(elems));
}

} // namespace

Result<ExprPtr> fsolve(
    ExprPtr equation,
    const Symbol& var,
    symbolic::CASContext& ctx,
    double search_low,
    double search_high) {

    if (!equation) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "fsolve: null equation"));
    }

    ExprPtr f = to_zero_form(equation, ctx);

    // Phase 1: try symbolic polynomial solve
    auto sym_res = solve_polynomial(f, var, ctx);
    if (sym_res.is_ok() && !sym_res.value().empty()) {
        // Wrap symbolic solutions in a matrix
        auto& solns = sym_res.value();
        // Filter to real-valued solutions (exclude those still containing imaginary units)
        std::vector<ExprPtr> real_solns;
        for (ExprPtr s : solns) {
            // Quick check: if the solution string doesn't contain 'i' literal it's likely real.
            // More robust: try numeric evaluation; if imaginary part is near zero, keep it.
            // Here we conservatively keep all symbolic solutions and let the user decide.
            real_solns.push_back(s);
        }
        std::vector<ExprPtr> elems(real_solns.begin(), real_solns.end());
        return ok(ctx.arena().make<Matrix>(static_cast<int>(elems.size()), 1, std::move(elems)));
    }

    // Phase 2a: if `f` parses as a univariate polynomial over Q, use the
    // exact Sturm-sequence root isolation (Sturm 1829). This counts the
    // number of distinct real roots in [low, high] exactly and refines
    // each within `tolerance`. Newton polish drives the final precision.
    constexpr double kTolerance = 1e-10;
    {
        auto sturm_res = numeric::find_polynomial_roots_sturm(
            f, var.name, ctx, search_low, search_high, kTolerance);
        if (sturm_res.is_ok()) {
            const auto& roots = sturm_res.value();
            if (roots.empty()) {
                return ok(ctx.arena().make<Matrix>(0, 1, std::vector<ExprPtr>{}));
            }
            return ok(roots_to_matrix(roots, ctx.arena()));
        }
        // fall through to grid scan if expr is not a rational polynomial
    }

    // Phase 2b: transcendental fallback — grid scan + bisection + Newton.
    // The 400-sample grid here is a *bootstrap*, not a precision bound:
    // every detected sign change is refined to `tolerance` by bisection
    // and then polished by Newton. This path is replaced by Lipschitz
    // dyadic refinement in Step 6 (see Roadmap WAVE A.2).
    numeric::MultiRootOptions opts;
    opts.low  = search_low;
    opts.high = search_high;
    opts.num_samples = 400;
    opts.dedup_tolerance = 1e-6;
    opts.root_opts.tolerance = kTolerance;
    opts.root_opts.max_iterations = 100;

    auto roots_res = numeric::find_roots_on_interval(f, var.name, ctx, opts);
    if (roots_res.is_error()) return fail<ExprPtr>(roots_res.error());

    const auto& roots = roots_res.value();
    if (roots.empty()) {
        // Return empty matrix — no roots found in search range
        return ok(ctx.arena().make<Matrix>(0, 1, std::vector<ExprPtr>{}));
    }

    return ok(roots_to_matrix(roots, ctx.arena()));
}

} // namespace cas::algebra
