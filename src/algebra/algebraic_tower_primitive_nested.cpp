// algebraic_tower_primitive_nested.cpp — F3.4-DEBT-01 closure.
//
// Provides nested_lift_min_poly: for an outer RootOf α whose defining
// polynomial f(x) has coefficients in Q(β) for some other collected
// RootOf β with rational min-poly g, computes the absolute polynomial
// R(x) = Res_y(g(y), f(x, y)) ∈ Q[x] (Cohen "A Course in Computational
// Algebraic Number Theory" §3.6.1; Trager 1976).
//
// Output: monic R if squarefree; otherwise explicit Unimplemented
// diagnostic (no silent failure — closes Cat 4 ban from CLAUDE.md).
//
// This file is the dedicated split target to keep
// algebraic_tower_primitive.cpp under the 500-LOC anti-monolith limit.

#include "cas/algebraic_tower_bridge.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include "polynomial_internal.hpp"
#include "algebraic_tower_primitive_internal.hpp"
#include "algebraic_tower_primitive_nested.hpp"

#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

namespace cas {
namespace algebra {
namespace primitive_internal {

namespace {

// Recursive walker: substitute occurrences of inner_expr (by structural
// equality) with replacement inside expr.
struct SubstWalker {
    ExprPtr inner_expr;
    ExprPtr replacement;
    AstArena& arena;
    ExprPtr walk(ExprPtr e) {
        if (!e) return e;
        if (structural_equal(e, inner_expr)) return replacement;
        ExprPtr result = e;
        visit_expr(e, [&](const auto& node) {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<Node, Unary>) {
                ExprPtr new_op = walk(node.operand);
                if (new_op != node.operand)
                    result = arena.make<Unary>(node.op, new_op);
            } else if constexpr (std::is_same_v<Node, Binary>) {
                ExprPtr nl = walk(node.left);
                ExprPtr nr = walk(node.right);
                if (nl != node.left || nr != node.right)
                    result = arena.make<Binary>(node.op, nl, nr);
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> na;
                na.reserve(node.args.size());
                bool changed = false;
                for (ExprPtr a : node.args) {
                    ExprPtr w = walk(a);
                    if (w != a) changed = true;
                    na.push_back(w);
                }
                if (changed) result = arena.make<FuncCall>(node.name, std::move(na));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> nt;
                nt.reserve(node.terms.size());
                bool changed = false;
                for (ExprPtr t : node.terms) {
                    ExprPtr w = walk(t);
                    if (w != t) changed = true;
                    nt.push_back(w);
                }
                if (changed) result = arena.make<Sum>(std::move(nt));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> nf;
                nf.reserve(node.factors.size());
                bool changed = false;
                for (ExprPtr f : node.factors) {
                    ExprPtr w = walk(f);
                    if (w != f) changed = true;
                    nf.push_back(w);
                }
                if (changed) result = arena.make<Product>(std::move(nf));
            }
        });
        return result;
    }
};

}  // namespace

Result<std::optional<AlgebraicNumber::CoeffVec>> try_nested_lift_min_poly(
    const RootOf& outer,
    ExprPtr beta_canon,
    const AlgebraicNumber::CoeffVec& beta_min_poly,
    symbolic::CASContext& ctx,
    const Deadline& deadline) {

    // Substitute β → fresh symbol y in outer's polynomial.
    const Symbol y_sym = ctx.make_fresh_symbol("__nested_y");
    ExprPtr y_expr = ctx.arena().make<Symbol>(y_sym);
    SubstWalker walker{beta_canon, y_expr, ctx.arena()};
    ExprPtr substituted = walker.walk(outer.polynomial);
    if (substituted == outer.polynomial) {
        // β not present in outer polynomial.
        return ok(std::optional<AlgebraicNumber::CoeffVec>{});
    }

    // Parse outer in its variable x ⇒ PolyExpr with ExprPtr coefficients in y.
    auto parsed_x = parse_polynomial(substituted, outer.variable, ctx);
    if (parsed_x.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
    const auto& poly_x = parsed_x.value();

    // Convert each coefficient (in y) to RatPoly.
    std::vector<RatPoly> f_xy;
    f_xy.reserve(poly_x.size());
    for (std::size_t k = 0U; k < poly_x.size(); ++k) {
        auto parsed_y = parse_polynomial(poly_x[k], y_sym, ctx);
        if (parsed_y.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
        auto rp = poly_to_rational_poly(parsed_y.value());
        if (rp.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
        f_xy.push_back(std::move(rp.value()));
    }

    // Inner min-poly g(y) as RatPoly.
    RatPoly g_y = vec_to_ratpoly(vec_make_monic(beta_min_poly));

    auto R_res = compute_absolute_resultant_xy(g_y, f_xy, deadline);
    if (R_res.is_error()) {
        return fail<std::optional<AlgebraicNumber::CoeffVec>>(R_res.error());
    }
    RatPoly R = std::move(R_res.value());
    if (R.is_zero() || R.degree() == 0U) {
        return ok(std::optional<AlgebraicNumber::CoeffVec>{});
    }

    // Squarefree check: in char 0, R is squarefree iff f(x, β) has no
    // repeated absolute root. Non-squarefree → irreducible-factor selection
    // over Q(β)[x] would be required, not yet implemented.
    if (!ratpoly_is_squarefree(R)) {
        return fail<std::optional<AlgebraicNumber::CoeffVec>>(CASError{
            CASErrorKind::Unimplemented,
            "try_nested_lift_min_poly (F3.4-DEBT-01): absolute resultant "
            "Res_y(g(y), f(x, y)) is not squarefree. Irreducible-factor "
            "selection over Q(β)[x] not yet implemented.",
            std::nullopt});
    }
    R = ratpoly_make_monic(std::move(R));
    return ok(std::optional<AlgebraicNumber::CoeffVec>(R.coefficients()));
}

}  // namespace primitive_internal

// ── detect_tower_n_level (relocated from algebraic_tower_primitive.cpp) ──────
//
// Detection: collects all distinct RootOf subexpressions in expr.  Each is
// canonicalized (RootOf-form preserved) and given a rational min-poly.  When
// a min-poly cannot be obtained because its coefficients depend on another
// collected RootOf (1-level nesting, F3.4-DEBT-01), the absolute resultant
// lift in try_nested_lift_min_poly resolves it.  Failure modes:
//
//   - RootOf depending on a non-collected algebraic element → Unimplemented.
//   - Absolute resultant not squarefree → Unimplemented (factor selection
//     not yet implemented).
//
// On success the collapsed primitive element is computed by Trager
// (compute_primitive_element).

Result<std::optional<PrimitiveElementResult>> detect_tower_n_level(
    ExprPtr expr,
    symbolic::CASContext& ctx) {
    if (!expr) return ok(std::optional<PrimitiveElementResult>{});

    // Collect all distinct RootOf subexpressions (up to 16).
    std::vector<ExprPtr> roots;
    struct Collector {
        std::vector<ExprPtr>& out;
        void collect(ExprPtr e) {
            if (!e) return;
            if (expr_is<RootOf>(e)) {
                bool seen = false;
                for (ExprPtr ex : out) {
                    if (structural_equal(ex, e)) { seen = true; break; }
                }
                if (!seen && out.size() < 16U) out.push_back(e);
                collect(expr_ref<RootOf>(e).polynomial);
                return;
            }
            visit_expr(e, [&](const auto& node) {
                using Node = std::decay_t<decltype(node)>;
                if constexpr (std::is_same_v<Node, Unary>) {
                    collect(node.operand);
                } else if constexpr (std::is_same_v<Node, Binary>) {
                    collect(node.left); collect(node.right);
                } else if constexpr (std::is_same_v<Node, FuncCall>) {
                    for (ExprPtr arg : node.args) collect(arg);
                } else if constexpr (std::is_same_v<Node, Sum>) {
                    for (ExprPtr term : node.terms) collect(term);
                } else if constexpr (std::is_same_v<Node, Product>) {
                    for (ExprPtr factor : node.factors) collect(factor);
                } else if constexpr (std::is_same_v<Node, Matrix>) {
                    for (ExprPtr item : node.elements) collect(item);
                }
            });
        }
    };
    Collector col{roots};
    col.collect(expr);

    if (roots.size() <= 1U) return ok(std::optional<PrimitiveElementResult>{});

    struct RootEntry {
        ExprPtr canon;
        const RootOf* node{nullptr};
        std::optional<AlgebraicNumber::CoeffVec> mp;
    };

    std::vector<RootEntry> entries;
    entries.reserve(roots.size());
    for (ExprPtr root_expr : roots) {
        const std::size_t saved = ctx.max_rootof_explicit_degree();
        ctx.set_max_rootof_explicit_degree(1U);
        auto canon = ctx.simplify(root_expr);
        ctx.set_max_rootof_explicit_degree(saved);
        if (canon.is_error()) return ok(std::optional<PrimitiveElementResult>{});
        const auto* root_node = expr_cast<RootOf>(canon.value());
        if (!root_node) return ok(std::optional<PrimitiveElementResult>{});
        RootEntry entry{canon.value(), root_node, std::nullopt};
        auto mp_res = rootof_min_poly(*root_node, ctx);
        if (mp_res.is_ok()) entry.mp = mp_res.value();
        entries.push_back(std::move(entry));
    }

    // Second pass: lift nested entries via absolute resultant (Cohen §3.6.1).
    const primitive_internal::Deadline lift_deadline =
        std::chrono::steady_clock::now() + ctx.timeout();

    for (std::size_t i = 0U; i < entries.size(); ++i) {
        if (entries[i].mp.has_value()) continue;
        bool lifted = false;
        for (std::size_t j = 0U; j < entries.size() && !lifted; ++j) {
            if (j == i || !entries[j].mp.has_value()) continue;
            auto lift = primitive_internal::try_nested_lift_min_poly(
                *entries[i].node, entries[j].canon, entries[j].mp.value(),
                ctx, lift_deadline);
            if (lift.is_error()) {
                return fail<std::optional<PrimitiveElementResult>>(lift.error());
            }
            if (lift.value().has_value()) {
                entries[i].mp = std::move(lift.value().value());
                lifted = true;
            }
        }
        if (!lifted) {
            return fail<std::optional<PrimitiveElementResult>>(CASError{
                CASErrorKind::Unimplemented,
                "detect_tower_n_level (F3.4-DEBT-01): RootOf at generator "
                "index " + std::to_string(i) + " has non-rational min-poly "
                "coefficients and no single rational-RootOf β was found "
                "covering all non-rational terms. Multi-β iterated lift "
                "(Cohen §3.6.4) not yet implemented.",
                std::nullopt});
        }
    }

    std::vector<ExprPtr> alphas;
    std::vector<AlgebraicNumber::CoeffVec> min_polys_vec;
    alphas.reserve(entries.size());
    min_polys_vec.reserve(entries.size());
    for (auto& e : entries) {
        alphas.push_back(e.canon);
        min_polys_vec.push_back(std::move(e.mp.value()));
    }

    auto result = compute_primitive_element(alphas, min_polys_vec, ctx);
    if (result.is_error()) {
        if (result.error().kind == CASErrorKind::Unimplemented)
            return ok(std::optional<PrimitiveElementResult>{});
        return fail<std::optional<PrimitiveElementResult>>(result.error());
    }
    return ok(std::optional<PrimitiveElementResult>(std::move(result.value())));
}

}  // namespace algebra
}  // namespace cas
