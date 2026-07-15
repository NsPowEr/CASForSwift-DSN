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

}  // namespace cas::symbolic
