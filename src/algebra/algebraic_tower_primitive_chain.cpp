// A9 / F3.4-DEBT-01 — Primitive element of a *sequentially nested* tower.
//
// The generic Trager merge in compute_primitive_element assumes the generators
// are algebraically independent: at step k it builds the shift resultant
// Res_y(q(y), m_k(z − s·y)), whose degree is deg(q)·deg(m_k). For a nested
// radical chain that product wildly overshoots the true field degree. Measured
// on β₁ = √3, β₂ = √(2+√3), α = √(1+√(2+√3)) — where Q(β₁,β₂,α) = Q(α) has
// degree 8 — the merge of α (deg 8) with β₂ (deg 4) builds a degree-32
// resultant, factors it, and recovers the correct degree-8 factor after 105 s;
// the next merge (expected degree 16) does not finish at all.
//
// But such a tower needs no resultant. Each level is defined by a polynomial
// that is *linear in its parent*:
//
//     β₂ = RootOf(z² − 2 − β₁)   ⇒   β₁ = β₂² − 2
//     α  = RootOf(z² − 1 − β₂)   ⇒   β₂ = α²  − 1
//
// so every generator is a polynomial in the deepest one. In general, if the
// defining polynomial of a child c is p(z) = A(z) + B(z)·g with A, B ∈ Q[z] —
// i.e. degree 1 in the parent g — then p(c) = 0 gives
//
//     g = −A(c) · B(c)⁻¹        (B(c) ≠ 0, invertible since M is irreducible)
//
// Taking θ = the generator of largest degree D and M = its Q-minimal
// polynomial, we resolve parents outward from θ, each as an element of
// Q[y]/(M). No resultant is ever formed; the whole tower costs a handful of
// multiplications mod M.
//
// Soundness. The relation g = −A(θ)/B(θ) is an identity in Q(θ) derived from
// p(θ) = 0, which holds for the actual root θ (not merely some conjugate), so
// the representation of every parent is the actual parent, with no root-index
// ambiguity: the child determines the parent uniquely. We nevertheless verify
// the result rather than trust the derivation — each generator's own absolute
// minimal polynomial must annihilate its computed representative mod M. A
// failed certificate, an unresolved generator, or a level that is not linear
// in its parent all return nullopt, and the caller falls back to the generic
// merge. Nothing here can produce a wrong answer silently.

#include "algebraic_tower_primitive_chain.hpp"

#include "algebraic_tower_primitive_internal.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {
namespace primitive_internal {

namespace {

// Rebuild `e` with every subexpression structurally equal to `target`
// replaced by `replacement`.
struct SubExprReplacer {
    ExprPtr target;
    ExprPtr replacement;
    AstArena& arena;

    ExprPtr walk(ExprPtr e) {  // NOLINT(misc-no-recursion)
        if (!e) return e;
        if (structural_equal(e, target)) return replacement;
        ExprPtr result = e;
        visit_expr(e, [&](const auto& node) {
            using Node = std::decay_t<decltype(node)>;
            if constexpr (std::is_same_v<Node, Unary>) {
                ExprPtr op = walk(node.operand);
                if (op != node.operand) result = arena.make<Unary>(node.op, op);
            } else if constexpr (std::is_same_v<Node, Binary>) {
                ExprPtr l = walk(node.left);
                ExprPtr r = walk(node.right);
                if (l != node.left || r != node.right)
                    result = arena.make<Binary>(node.op, l, r);
            } else if constexpr (std::is_same_v<Node, FuncCall>) {
                std::vector<ExprPtr> args;
                args.reserve(node.args.size());
                bool changed = false;
                for (ExprPtr a : node.args) {
                    ExprPtr w = walk(a);
                    changed = changed || (w != a);
                    args.push_back(w);
                }
                if (changed) result = arena.make<FuncCall>(node.func_id, std::move(args));
            } else if constexpr (std::is_same_v<Node, Sum>) {
                std::vector<ExprPtr> terms;
                terms.reserve(node.terms.size());
                bool changed = false;
                for (ExprPtr t : node.terms) {
                    ExprPtr w = walk(t);
                    changed = changed || (w != t);
                    terms.push_back(w);
                }
                if (changed) result = arena.make<Sum>(std::move(terms));
            } else if constexpr (std::is_same_v<Node, Product>) {
                std::vector<ExprPtr> factors;
                factors.reserve(node.factors.size());
                bool changed = false;
                for (ExprPtr f : node.factors) {
                    ExprPtr w = walk(f);
                    changed = changed || (w != f);
                    factors.push_back(w);
                }
                if (changed) result = arena.make<Product>(std::move(factors));
            }
        });
        return result;
    }
};

// Horner evaluation of p ∈ Q[z] at a ring element `x` ∈ Q[y]/(m).
[[nodiscard]] RatPoly eval_at_ring_element(const RatPoly& p, const RatPoly& x,
                                           const RatPoly& m) {
    RatPoly acc;  // zero
    for (std::size_t i = p.size(); i-- > 0U;) {
        acc = ratpoly_mulmod(acc, x, m);
        RatPoly c(std::vector<Rational>{p[i]});
        c.normalize([](const Rational& r) { return r.numerator().is_zero(); });
        acc = ratpoly_mod(add_rational_poly(acc, c), m);
    }
    return acc;
}

// Inverse of `a` in Q[y]/(m); nullopt when gcd(a, m) is not a unit (only
// possible if a ≡ 0, since m is irreducible).
[[nodiscard]] std::optional<RatPoly> ring_inverse(const RatPoly& a,
                                                  const RatPoly& m) {
    if (a.is_zero()) return std::nullopt;
    auto [g, s, t] = extended_gcd_rational_poly(a, m);
    (void)t;
    if (g.is_zero() || g.degree() != 0U) return std::nullopt;
    const Rational inv_lc = Rational(BigInt(1)) / g[0];
    RatPoly out = s;
    for (auto& c : out.coefficients()) c = c * inv_lc;
    return ratpoly_mod(out, m);
}

[[nodiscard]] std::vector<Rational> pad_to(const RatPoly& p, std::size_t deg) {
    std::vector<Rational> out(deg, Rational(BigInt(0)));
    for (std::size_t i = 0U; i < p.size() && i < deg; ++i) out[i] = p[i];
    return out;
}

}  // namespace

Result<std::optional<PrimitiveElementResult>> try_primitive_element_from_chain(
    const std::vector<ExprPtr>& alphas,
    const std::vector<const RootOf*>& nodes,
    const std::vector<AlgebraicNumber::CoeffVec>& min_polys,
    symbolic::CASContext& ctx) {
    const std::size_t n = alphas.size();
    if (n < 2U || nodes.size() != n || min_polys.size() != n) {
        return ok(std::optional<PrimitiveElementResult>{});
    }

    // θ = the generator of largest degree; every other degree must divide it
    // (a necessary condition for Q(α_i) ⊆ Q(θ)).
    std::size_t theta_idx = 0U;
    for (std::size_t i = 1U; i < n; ++i) {
        if (min_polys[i].size() > min_polys[theta_idx].size()) theta_idx = i;
    }
    const std::size_t deg_theta = min_polys[theta_idx].size() - 1U;
    if (deg_theta < 2U) return ok(std::optional<PrimitiveElementResult>{});
    for (std::size_t i = 0U; i < n; ++i) {
        const std::size_t d = min_polys[i].size() - 1U;
        if (d == 0U || (deg_theta % d) != 0U) {
            return ok(std::optional<PrimitiveElementResult>{});
        }
    }

    const RatPoly m_theta = vec_to_ratpoly(vec_make_monic(min_polys[theta_idx]));

    // reps[i] = α_i as an element of Q[y]/(m_theta); θ itself is y.
    std::vector<std::optional<RatPoly>> reps(n);
    reps[theta_idx] =
        RatPoly(std::vector<Rational>{Rational(BigInt(0)), Rational(BigInt(1))});

    const Symbol w_sym = ctx.make_fresh_symbol("__chain_w");
    ExprPtr w_expr = ctx.arena().make<Symbol>(w_sym);

    // Resolve parents outward from θ: whenever a resolved child's defining
    // polynomial is linear in an unresolved generator, that generator is a
    // rational function of the child, hence an element of Q[y]/(m_theta).
    bool progress = true;
    while (progress) {
        progress = false;
        for (std::size_t c = 0U; c < n; ++c) {
            if (!reps[c].has_value()) continue;
            for (std::size_t g = 0U; g < n; ++g) {
                if (reps[g].has_value()) continue;

                SubExprReplacer repl{alphas[g], w_expr, ctx.arena()};
                ExprPtr subst = repl.walk(nodes[c]->polynomial);
                if (subst == nodes[c]->polynomial) continue;  // g not present

                // Must be degree 1 in the parent g.
                auto in_w = parse_polynomial(subst, w_sym, ctx);
                if (in_w.is_error() || poly_degree(in_w.value()) != 1U) continue;

                // A(z), B(z) must have rational coefficients: no third
                // generator may hide inside them.
                auto a_pe = parse_polynomial(in_w.value()[0], nodes[c]->variable, ctx);
                auto b_pe = parse_polynomial(in_w.value()[1], nodes[c]->variable, ctx);
                if (a_pe.is_error() || b_pe.is_error()) continue;
                auto a_rp = poly_to_rational_poly(a_pe.value());
                auto b_rp = poly_to_rational_poly(b_pe.value());
                if (a_rp.is_error() || b_rp.is_error()) continue;

                const RatPoly a_val = eval_at_ring_element(a_rp.value(), *reps[c], m_theta);
                const RatPoly b_val = eval_at_ring_element(b_rp.value(), *reps[c], m_theta);
                auto b_inv = ring_inverse(b_val, m_theta);
                if (!b_inv.has_value()) continue;

                RatPoly neg_a = a_val;
                for (auto& coeff : neg_a.coefficients()) {
                    coeff = Rational(BigInt(0)) - coeff;
                }
                reps[g] = ratpoly_mulmod(neg_a, *b_inv, m_theta);
                progress = true;
            }
        }
    }

    for (std::size_t i = 0U; i < n; ++i) {
        if (!reps[i].has_value()) return ok(std::optional<PrimitiveElementResult>{});
    }

    // Certificate: each generator's own absolute minimal polynomial must
    // annihilate its representative in Q[y]/(m_theta). Derivation aside, an
    // unverified representative is never returned.
    for (std::size_t i = 0U; i < n; ++i) {
        const RatPoly m_i = vec_to_ratpoly(vec_make_monic(min_polys[i]));
        const RatPoly value = eval_at_ring_element(m_i, *reps[i], m_theta);
        if (!value.is_zero()) return ok(std::optional<PrimitiveElementResult>{});
    }

    PrimitiveElementResult res;
    res.theta_expr = alphas[theta_idx];
    res.shifts = {};
    res.min_poly_theta = vec_make_monic(min_polys[theta_idx]);
    res.alphas_in_theta.reserve(n);
    for (std::size_t i = 0U; i < n; ++i) {
        res.alphas_in_theta.push_back(pad_to(*reps[i], deg_theta));
    }
    return ok(std::optional<PrimitiveElementResult>(std::move(res)));
}

}  // namespace primitive_internal
}  // namespace algebra
}  // namespace cas
