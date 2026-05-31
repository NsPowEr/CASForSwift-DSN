#pragma once

#include "cas/algebra.hpp"
#include "cas/algebraic_number.hpp"
#include "cas/algebraic_tower.hpp"
#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas {
namespace algebra {

struct TowerGenerators {
    ExprPtr alpha_1;
    AlgebraicNumber::CoeffVec min_poly_1;

    ExprPtr alpha_2;
    std::vector<AlgebraicNumber> min_poly_2;
};

[[nodiscard]] Result<std::optional<TowerGenerators>> detect_two_level_tower(
    ExprPtr expr,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::optional<AlgebraicTowerTwoLevel>> try_express_in_tower_two_level(
    ExprPtr expr,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx);

[[nodiscard]] ExprPtr tower_to_expr(
    const AlgebraicTowerTwoLevel& value,
    const TowerGenerators& gens,
    AstArena& arena);

// L3-06: factorization of p(x) in Q[x] over the two-level tower Q(α_1, α_2)
// via composite Trager shift (s_1, s_2) and iterated absolute resultant.
[[nodiscard]] Result<Factorization> factor_polynomial_tower(
    ExprPtr poly,
    const Symbol& var,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx);

// ── F3.4 Primitive Element Theorem (Trager 1976; Cohen §3.6.2; GCL §8.7) ──
//
// PrimitiveElementResult represents the collapse of Q(α₁, …, α_n) → Q(θ).
//
// Fields:
//   theta_expr        — ExprPtr for θ built from the α_i using integer shifts
//                       θ = α₁ + s₂·α₂ + s₃·α₃ + … + s_n·α_n.
//   shifts            — shifts[i] = s_{i+1} ∈ ℤ (0-indexed: shifts[0] = s₂,
//                       …).  shifts.size() == n-1.
//   min_poly_theta    — Monic min-poly of θ over Q (ascending coeff, ≥2 entries).
//                       degree = deg(m₁) · … · deg(m_n).
//   alphas_in_theta   — alphas_in_theta[i] = coefficients of αᵢ as a Q-poly
//                       in θ of degree < deg(min_poly_theta).
//                       alphas_in_theta.size() == n (one entry per generator).

struct PrimitiveElementResult {
    ExprPtr theta_expr;
    std::vector<BigInt> shifts;
    AlgebraicNumber::CoeffVec min_poly_theta;
    std::vector<std::vector<Rational>> alphas_in_theta;
};

// Compute primitive element θ for Q(α₁, …, α_n).
//
// Each α_i is given as an ExprPtr (used only for constructing theta_expr and
// the verification output) together with its monic rational minimal polynomial
// min_polys[i] over Q.
//
// The incremental Trager shift search proceeds as:
//   K₁ = Q(α₁),  θ₁ = α₁.
//   For k = 2…n:
//     Try θ_k = θ_{k-1} + s·α_k for s = 1, 2, …
//     Candidate min-poly R_s(y) = Res_x( m_k(x),  Q_{k-1}(y − s·x) )
//     where Q_{k-1} is the min-poly of θ_{k-1}.
//     Accept iff R_s is squarefree over Q (gcd(R_s, R_s') == constant).
//     Express α_k and θ_{k-1} in Q(θ_k) via extended-GCD linear algebra.
//
// Returns Unimplemented (with attempt count in message) if shift search
// exhausts ctx.max_trager_tower_shift_attempts() for any step.
//
// ZERO DEBT: all limits derived from ctx; no silent truncation.
[[nodiscard]] Result<PrimitiveElementResult> compute_primitive_element(
    const std::vector<ExprPtr>& alphas,
    const std::vector<AlgebraicNumber::CoeffVec>& min_polys,
    symbolic::CASContext& ctx);

// ── F3.5 n-level tower factorisation (Trager via primitive element) ──
//
// Generic input form: a list of generators (α₁, …, α_n) with their rational
// minimal polynomials.  Pipeline:
//   1. Reduce Q(α₁, …, α_n) → Q(θ) via compute_primitive_element (F3.4).
//   2. Single-extension Trager: find integer s such that
//        N(x) = Res_y( q_θ(y), f(x − s·y) ) ∈ Q[x] is square-free.
//   3. Factor N over Q via factor_polynomial; each Q-factor g_i lifts to
//        f_i(x) = gcd_{Q(θ)[x]}( f(x), g_i(x + s·θ) ).
//   4. Trager invariant: Σ deg(f_i) == deg(f) — internal_error if violated.
//
// This is dramatically faster than the 2-level iterated resultant path
// (one resultant of degree D·deg(f) instead of nested resultants of total
// degree D₁·D₂·deg(f)).  It also handles n ≥ 3 and irreducible inputs
// (N becomes a perfect power of an irreducible — square-free fails and
// we increment s, but the small single resultant is computed quickly).
//
// Shift budget: ctx.max_trager_tower_shift_attempts() (shared with F3.4).
// Cap exceeded → explicit Unimplemented (NO silent truncation).
//
// Returned factors are expressed in the original generators (re-rendered
// via alphas_in_theta back-substitution at the call site).
struct TowerGeneratorsN {
    std::vector<ExprPtr> alphas;
    std::vector<AlgebraicNumber::CoeffVec> min_polys;
};

[[nodiscard]] Result<Factorization> factor_polynomial_tower_n(
    ExprPtr poly,
    const Symbol& var,
    const TowerGeneratorsN& gens,
    symbolic::CASContext& ctx);

// Detect an n-level tower Q(α₁, …, α_n) from a nested-radical ExprPtr
// (e.g. sqrt(2) + sqrt(3) + sqrt(5)) and compute its primitive element.
//
// Detection collects all distinct RootOf subexpressions in `expr`.  For n == 0
// or n == 1 returns std::nullopt (no tower needed).  For n == 2 this is the
// fast path delegating to detect_two_level_tower; for n ≥ 3 uses
// compute_primitive_element.
//
// F3.4-DEBT-01 (2026-05-30): 1-level nesting now supported via absolute
// resultant lift (Cohen §3.6.1).  When a RootOf coefficient depends on
// another collected RootOf β with rational min-poly g(y), the outer
// polynomial f(x,β) is lifted to R(x) = Res_y(g(y), f(x, y)) ∈ Q[x] and
// used as effective min-poly if squarefree.  Residual OPEN cases (R
// reducible, multi-β iterated nesting, γ not in pool) return explicit
// Unimplemented diagnostics.
[[nodiscard]] Result<std::optional<PrimitiveElementResult>> detect_tower_n_level(
    ExprPtr expr,
    symbolic::CASContext& ctx);

}  // namespace algebra
}  // namespace cas
