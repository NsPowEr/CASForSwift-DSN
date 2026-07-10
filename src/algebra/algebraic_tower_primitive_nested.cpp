// algebraic_tower_primitive_nested.cpp — F3.4-DEBT-01 and A9 closure.
//
// Provides nested_lift_min_poly and detect_tower_n_level for arbitrary
// N-level algebraic extensions Q(α₁, ..., α_n) with multi-β nesting
// and non-squarefree/reducible lifted minimal polynomials.
//
// Solves Cohen "A Course in Computational Algebraic Number Theory" §3.6.1-3.6.4;
// Trager 1976.

#include "cas/algebraic_tower_bridge.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include "polynomial_internal.hpp"
#include "algebraic_tower_primitive_internal.hpp"
#include "algebraic_tower_primitive_nested.hpp"

#include "algebraic_tower_primitive_chain.hpp"

#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

namespace cas {
namespace algebra {
namespace primitive_internal {

namespace {

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

struct PresenceWalker {
    ExprPtr target;
    bool found{false};
    void check(ExprPtr e) {
        if (!e || found) return;
        if (structural_equal(e, target)) {
            found = true;
            return;
        }
        visit_expr(e, [&](const auto& node) {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<Node, Unary>) {
                check(node.operand);
            } else if constexpr (std::is_same_v<Node, Binary>) {
                check(node.left); check(node.right);
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                for (ExprPtr a : node.args) check(a);
            } else if constexpr (std::is_same_v<Node, Sum>) {
                for (ExprPtr t : node.terms) check(t);
            } else if constexpr (std::is_same_v<Node, Product>) {
                for (ExprPtr f : node.factors) check(f);
            } else if constexpr (std::is_same_v<Node, Matrix>) {
                for (ExprPtr item : node.elements) check(item);
            }
        });
    }
};

inline ExprPtr poly_to_expr(const std::vector<Rational>& coeffs, ExprPtr var_expr, AstArena& arena) {
    if (coeffs.empty()) return arena.make<IntegerLit>(BigInt(0));
    ExprPtr result = arena.make<RationalLit>(coeffs.back().numerator(), coeffs.back().denominator());
    for (std::size_t i = coeffs.size() - 1U; i > 0U; --i) {
        ExprPtr mul = arena.make<Binary>(BinaryOp::Mul, result, var_expr);
        ExprPtr term = arena.make<RationalLit>(coeffs[i - 1U].numerator(), coeffs[i - 1U].denominator());
        result = arena.make<Binary>(BinaryOp::Add, mul, term);
    }
    return result;
}

}  // namespace

Result<std::optional<AlgebraicNumber::CoeffVec>> try_nested_lift_min_poly_multi(
    const RootOf& outer,
    ExprPtr outer_canon,
    const std::vector<ResolvedGen>& resolved_gens,
    symbolic::CASContext& ctx,
    const Deadline& deadline) {

    std::vector<ResolvedGen> matched;
    for (const auto& rg : resolved_gens) {
        if (!rg.canon) continue;
        PresenceWalker pw{rg.canon, false};
        pw.check(outer.polynomial);
        if (pw.found) {
            matched.push_back(rg);
        }
    }

    if (matched.empty()) {
        return ok(std::optional<AlgebraicNumber::CoeffVec>{});
    }

    const Symbol y_sym = ctx.make_fresh_symbol("__nested_y");
    ExprPtr y_expr = ctx.arena().make<Symbol>(y_sym);

    std::vector<RatPoly> f_xy;
    RatPoly g_y;

    if (matched.size() == 1U) {
        SubstWalker walker{matched[0].canon, y_expr, ctx.arena()};
        ExprPtr substituted = walker.walk(outer.polynomial);
        if (substituted == outer.polynomial) {
            return ok(std::optional<AlgebraicNumber::CoeffVec>{});
        }
        auto parsed_x = parse_polynomial(substituted, outer.variable, ctx);
        if (parsed_x.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
        const auto& poly_x = parsed_x.value();

        f_xy.reserve(poly_x.size());
        for (std::size_t k = 0U; k < poly_x.size(); ++k) {
            auto parsed_y = parse_polynomial(poly_x[k], y_sym, ctx);
            if (parsed_y.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
            auto rp = poly_to_rational_poly(parsed_y.value());
            if (rp.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
            f_xy.push_back(std::move(rp.value()));
        }
        g_y = vec_to_ratpoly(vec_make_monic(matched[0].mp));
    } else {
        std::vector<ExprPtr> alphas;
        std::vector<AlgebraicNumber::CoeffVec> min_polys;
        alphas.reserve(matched.size());
        min_polys.reserve(matched.size());
        for (const auto& m : matched) {
            alphas.push_back(m.canon);
            min_polys.push_back(m.mp);
        }
        auto prim_res = compute_primitive_element(alphas, min_polys, ctx);
        if (prim_res.is_error()) {
            return fail<std::optional<AlgebraicNumber::CoeffVec>>(prim_res.error());
        }
        const auto& prim = prim_res.value();

        ExprPtr substituted = outer.polynomial;
        for (std::size_t idx = 0U; idx < matched.size(); ++idx) {
            ExprPtr theta_poly_expr = poly_to_expr(prim.alphas_in_theta[idx], y_expr, ctx.arena());
            SubstWalker walker{matched[idx].canon, theta_poly_expr, ctx.arena()};
            substituted = walker.walk(substituted);
        }

        auto parsed_x = parse_polynomial(substituted, outer.variable, ctx);
        if (parsed_x.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
        const auto& poly_x = parsed_x.value();

        f_xy.reserve(poly_x.size());
        for (std::size_t k = 0U; k < poly_x.size(); ++k) {
            auto parsed_y = parse_polynomial(poly_x[k], y_sym, ctx);
            if (parsed_y.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
            auto rp = poly_to_rational_poly(parsed_y.value());
            if (rp.is_error()) return ok(std::optional<AlgebraicNumber::CoeffVec>{});
            f_xy.push_back(std::move(rp.value()));
        }
        g_y = vec_to_ratpoly(vec_make_monic(prim.min_poly_theta));
    }

    auto R_res = compute_absolute_resultant_xy(g_y, f_xy, deadline);
    if (R_res.is_error()) {
        return fail<std::optional<AlgebraicNumber::CoeffVec>>(R_res.error());
    }
    RatPoly R = std::move(R_res.value());
    if (R.is_zero() || R.degree() == 0U) {
        return ok(std::optional<AlgebraicNumber::CoeffVec>{});
    }
    R = ratpoly_make_monic(std::move(R));

    auto all_factors = collect_irred_factors_over_q(R, ctx);
    if (all_factors.size() == 1U) {
        // Single candidate: either R is irreducible, or R = f^k (min-poly f),
        // or the factorization fell back to {R} on an internal error.  In the
        // fallback case R may be non-squarefree — a non-squarefree "min-poly"
        // is never valid, so guard it (REGOLA ZERO, no silent-wrong).
        if (!ratpoly_is_squarefree(all_factors[0])) {
            return make_unimplemented<std::optional<AlgebraicNumber::CoeffVec>>(UnimplementedInfo{
                .module      = "algebra",
                .function    = "try_nested_lift_min_poly_multi",
                .input_shape = "non-squarefree absolute resultant, factorization unavailable",
                .reason      = cas::error::reason_codes::POLY_FACTOR_EXTENSION,
                .suggestion  = "Squarefree decomposition of the absolute resultant failed",
                .ticket      = "A9"
            });
        }
        return ok(std::optional<AlgebraicNumber::CoeffVec>(all_factors[0].coefficients()));
    }

    std::vector<RatPoly> matching;
    matching.reserve(all_factors.size());
    for (const auto& f_cand : all_factors) {
        if (cand_vanishes_at_theta_expr(f_cand, outer_canon, ctx)) {
            matching.push_back(f_cand);
        }
    }
    // Exactly one certified factor: two distinct monic irreducible factors
    // cannot share a root, so a proven-vanishing factor is THE min-poly of α.
    // No certificate (the simplifier cannot reduce cand(α) to a literal zero):
    // picking an arbitrary factor would assign α a polynomial it may not
    // satisfy — silent-wrong downstream.  Bail explicitly (REGOLA ZERO).
    if (matching.size() == 1U) {
        return ok(std::optional<AlgebraicNumber::CoeffVec>(matching[0].coefficients()));
    }
    return make_unimplemented<std::optional<AlgebraicNumber::CoeffVec>>(UnimplementedInfo{
        .module      = "algebra",
        .function    = "try_nested_lift_min_poly_multi",
        .input_shape = "reducible absolute resultant: " +
                       std::to_string(all_factors.size()) + " irreducible factors, " +
                       std::to_string(matching.size()) + " certified vanishing at the outer root",
        .reason      = cas::error::reason_codes::POLY_FACTOR_EXTENSION,
        .suggestion  = "Factor selection needs an exact vanishing certificate "
                       "(e.g. isolating-interval arithmetic on the nested RootOf)",
        .ticket      = "A9"
    });
}

Result<std::optional<AlgebraicNumber::CoeffVec>> try_nested_lift_min_poly(
    const RootOf& outer,
    ExprPtr beta_canon,
    const AlgebraicNumber::CoeffVec& beta_min_poly,
    symbolic::CASContext& ctx,
    const Deadline& deadline) {
    std::vector<ResolvedGen> gens;
    gens.push_back(ResolvedGen{beta_canon, nullptr, beta_min_poly});
    ExprPtr outer_canon = ctx.arena().make<RootOf>(outer.polynomial, outer.variable, outer.root_index);
    return try_nested_lift_min_poly_multi(outer, outer_canon, gens, ctx, deadline);
}


}  // namespace primitive_internal

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

    // Iterative resolution pass for multi-level / multi-beta nested extensions.
    const primitive_internal::Deadline lift_deadline =
        std::chrono::steady_clock::now() + ctx.timeout();

    std::vector<primitive_internal::ResolvedGen> resolved;
    resolved.reserve(entries.size());
    for (const auto& e : entries) {
        if (e.mp.has_value()) {
            resolved.push_back(primitive_internal::ResolvedGen{e.canon, e.node, e.mp.value()});
        }
    }

    bool progress = true;
    while (progress && resolved.size() < entries.size()) {
        progress = false;
        for (std::size_t i = 0U; i < entries.size(); ++i) {
            if (entries[i].mp.has_value()) continue;

            auto lift = primitive_internal::try_nested_lift_min_poly_multi(
                *entries[i].node, entries[i].canon, resolved, ctx, lift_deadline);
            if (lift.is_error()) {
                return fail<std::optional<PrimitiveElementResult>>(lift.error());
            }
            if (lift.value().has_value()) {
                entries[i].mp = std::move(lift.value().value());
                resolved.push_back(primitive_internal::ResolvedGen{entries[i].canon, entries[i].node, entries[i].mp.value()});
                progress = true;
            }
        }
    }

    if (resolved.size() < entries.size()) {
        return make_unimplemented<std::optional<PrimitiveElementResult>>(UnimplementedInfo{
            .module      = "algebra",
            .function    = "detect_tower_n_level",
            .input_shape = "RootOf with non-rational min-poly coefficients",
            .reason      = cas::error::reason_codes::POLY_FACTOR_EXTENSION,
            .suggestion  = "Check for cyclic or uncollected algebraic dependencies",
            .ticket      = "F3.4-DEBT-01"
        });
    }

    std::vector<ExprPtr> alphas;
    std::vector<const RootOf*> nodes;
    std::vector<AlgebraicNumber::CoeffVec> min_polys_vec;
    alphas.reserve(entries.size());
    nodes.reserve(entries.size());
    min_polys_vec.reserve(entries.size());
    for (auto& e : entries) {
        alphas.push_back(e.canon);
        nodes.push_back(e.node);
        min_polys_vec.push_back(std::move(e.mp.value()));
    }

    // A9: a sequentially nested tower (each level linear in the previous one)
    // is flattened without a single shift resultant. Returns nullopt — never a
    // wrong answer — when the structure is not a chain, and we fall through to
    // the generic Trager merge below.
    auto chain = primitive_internal::try_primitive_element_from_chain(
        alphas, nodes, min_polys_vec, ctx);
    if (chain.is_error()) {
        return fail<std::optional<PrimitiveElementResult>>(chain.error());
    }
    if (chain.value().has_value()) {
        return ok(std::optional<PrimitiveElementResult>(
            std::move(chain.value().value())));
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
