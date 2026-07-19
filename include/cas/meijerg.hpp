#pragma once

// A7 Brick 2 — Meijer G-function factory, validation, canonicalization and
// typed accessor. The node itself is FuncCall{BuiltinOp::MeijerG, args}
// (no new NodeType — see Meijer_G_Slater.md §7); this header is the ONLY
// sanctioned way to build or decode one.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Meijer_G_Slater.md §7.

#include "cas/ast.hpp"
#include "cas/result.hpp"

#include <cstddef>
#include <vector>

namespace cas::symbolic {

class CASContext;

// Typed view of a validated MeijerG FuncCall's flat argument layout
// (args = [m, n, p, q, a_1..a_p, b_1..b_q, z], Meijer_G_Slater.md §7.1).
struct MeijerGView {
    std::size_t m;
    std::size_t n;
    std::size_t p;
    std::size_t q;
    std::vector<ExprPtr> a;  // size p, upper parameters
    std::vector<ExprPtr> b;  // size q, lower parameters
    ExprPtr z;
};

// Builds G^{m,n}_{p,q}(z | a; b) as a validated, canonicalized FuncCall.
// p = a.size(), q = b.size() (derived, not separate parameters).
// - REQUIRE 0 <= m <= q, 0 <= n <= p (else InvalidArgument, structured).
// - REQUIRE p + q <= ctx.meijerg_max_param_count() (else Unimplemented,
//   structured, diagnostic -- never silent truncation).
// - Pole-overlap check (§2.2): for j in 1..n, k in 1..m, ctx.simplify(a[j-1]
//   - b[k-1]) decided to be a positive IntegerLit -> PoleOverlap error. Only
//   applied when decidable; symbolic differences that don't reduce to a
//   literal are accepted (the guard is exact-when-decidable, not exhaustive
//   -- matches DLMF's own "unless" phrasing, never a false rejection).
// - Canonicalizes each of the four sub-groups (a_1..a_n / a_{n+1}..a_p /
//   b_1..b_m / b_{m+1}..b_q) independently via canonical_compare -- the
//   G-value is invariant under permuting WITHIN a sub-group (each is a
//   product of Gammas), so this gives structural sharing + pointer equality
//   between equivalent constructions (Regola 2). The order BETWEEN the four
//   groups is semantically rigid and is never touched.
[[nodiscard]] Result<ExprPtr> make_meijerg(
    CASContext& ctx,
    std::size_t m, std::size_t n,
    std::vector<ExprPtr> a, std::vector<ExprPtr> b,
    ExprPtr z);

// Decodes a FuncCall{BuiltinOp::MeijerG, ...}'s flat args back into a typed
// MeijerGView. Fails structured (InternalError) if `call` is not a
// well-formed MeijerG node (wrong func_id, size mismatch, or m/n/p/q not
// IntegerLit) -- this should only happen if a node bypassed make_meijerg.
[[nodiscard]] Result<MeijerGView> view_meijerg(const FuncCall& call);

// A7 step 3 (Meijer_G_Slater.md §5, §3.1) -- elementary/hypergeometric to
// Meijer G. Returns an expression MATHEMATICALLY EQUAL to `expr` whose
// transcendental core is a MeijerG node (typically prefactor * G); the input
// is returned untouched only never -- unsupported shapes fail with a
// structured Unimplemented (never a silent passthrough). Table entries §5.1-
// 5.8 are the fast path; hypergeometric nodes route through the GENERAL
// DLMF 16.18.1 bridge below (CLAUDE.md cat. 8: table = fast path only).
// Not yet representable (missing engine prerequisite, see TASKLIST A7):
// incomplete gamma §5.9 (no BuiltinOp for Gamma(a,z)/gamma(a,z) exists).
[[nodiscard]] Result<ExprPtr> to_meijerg(CASContext& ctx, ExprPtr expr);

// A7 step 4 (Meijer_G_Slater.md §3.2) -- Slater expansion, the GENERAL
// inverse path (DLMF 16.17.2): G^{m,n}_{p,q}(z) = sum_{k=1..m} A_k * pFq_k.
// Preconditions enforced: p <= q; no two of b_1..b_m may differ by a
// decidable integer (confluent poles need logarithmic terms -- structured
// Unimplemented, never a wrong series). Each summand's pF(q-1) is emitted
// as a closed form for (p,q-1) in {(0,0): exp, (1,0): binomial} or as a
// Hypergeometric0F1/1F1/2F1 node; other arities have no engine node yet ->
// structured Unimplemented (the G node simply stays, spec §9.4).
[[nodiscard]] Result<ExprPtr> slater_expand(CASContext& ctx, const FuncCall& g);

// A7 step 4 -- single-node inverse conversion: table §5 inverse fast path
// (shape-gated on the canonical argument forms produced by to_meijerg,
// recovering sin/cos/sinh/cosh/arctan/arcsin/erf/BesselJ/ln exactly) with
// slater_expand as the general fallback (CLAUDE.md cat. 8). Fails
// structured when neither path applies -- the caller may legitimately keep
// the G node (it is a first-class node, not an error state).
[[nodiscard]] Result<ExprPtr> from_meijerg(CASContext& ctx, const FuncCall& g);

// Maps from_meijerg over every MeijerG node in `expr`, leaving nodes with
// no known expansion INTACT (spec §9.4: staying in G form is a valid
// outcome, not an error). Returns the original pointer when nothing
// changed (Regola 2).
[[nodiscard]] Result<ExprPtr> expand_meijerg_nodes(CASContext& ctx, ExprPtr expr);

// A7 step 5 (Meijer_G_Slater.md §6) -- rewrite identities. Every rebuild
// goes back through make_meijerg (the §2.2 guards re-run; an overlap is a
// structured refusal). Formulas verified numerically in the spec.
// §6.2 (DLMF 16.19.2): z^mu G(z|a;b) = G(z|a+mu;b+mu) -- returns the
// shifted node (the caller owns the z^mu bookkeeping).
[[nodiscard]] Result<ExprPtr> meijerg_power_shift(
    CASContext& ctx, const FuncCall& g, ExprPtr mu);
// §6.3 (DLMF 16.19.3): cancels every decidably-equal pair (upper n-group
// param, lower outside-m param) to fixpoint; no pair -> identical rebuild.
[[nodiscard]] Result<ExprPtr> meijerg_cancel_common_param(
    CASContext& ctx, const FuncCall& g);
// §6.1 (DLMF 16.19.1): builds G^{n,m}_{q,p}(1/z | 1-b ; 1-a) from g.
[[nodiscard]] Result<ExprPtr> meijerg_invert_argument(
    CASContext& ctx, const FuncCall& g);
// §6.6 (h=-1): int G(t|a;b) dt = z * G^{m,n+1}_{p+1,q+1}(z | 0,a ; b,-1).
// The argument of g must BE the integration variable's expression; callers
// with G(c*x^r) substitute first (see integrate_meijerg.cpp).
[[nodiscard]] Result<ExprPtr> meijerg_antiderivative(
    CASContext& ctx, const FuncCall& g);
// §6.5 (h=+1, DLMF 16.19.5 family; numerically certified 2026-07-19,
// mpmath 3 shapes x 4 irrational points, err < 1e-41):
//   z * d/dz G^{m,n}_{p,q}(z|a;b) = G^{m,n+1}_{p+1,q+1}(z | 0,a ; b, +1).
// Returns the SHIFTED NODE only; the caller owns the (dz/dx)/z chain-rule
// prefactor (see differentiate_rules.cpp).
[[nodiscard]] Result<ExprPtr> meijerg_derivative_shift(
    CASContext& ctx, const FuncCall& g);

// General pFq -> G bridge (DLMF 16.18.1, spec §3.1), valid for p <= q+1:
//   pFq(alpha; beta; z) = [prod Gamma(beta) / prod Gamma(alpha)] *
//     G^{1,p}_{p,q+1}(-z | 1-alpha_1..1-alpha_p ; 0, 1-beta_1..1-beta_q).
// Fails structured when p > q+1, or when any alpha_k / beta_k is decidably a
// nonpositive integer (Gamma pole: the pFq degenerates to a polynomial /
// is undefined -- the bridge formula does not apply there).
[[nodiscard]] Result<ExprPtr> pfq_to_meijerg(
    CASContext& ctx,
    std::vector<ExprPtr> alpha, std::vector<ExprPtr> beta, ExprPtr z);

}  // namespace cas::symbolic
