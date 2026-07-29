// differential_field_integration.cpp
// Rational-function integration over a differential field (Hermite
// reduction + Rothstein-Trager), split out of differential_field.cpp
// (T-046 anti-monolith, F5.1 field-ops vs integration).

#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "cas/error_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"
#include "calculus_internal.hpp"

#include <algorithm>
#include <string>

namespace cas::calculus {

Result<ExprPtr> compute_resultant(ExprPtr A, ExprPtr B, const Symbol& x, symbolic::CASContext& ctx) {
    return algebra::polynomial_resultant(A, B, x, ctx);
}

// Horowitz-Ostrogradsky one-step reduction for ∫ P/V^n (n >= 2, V squarefree).
// Finds A of deg < deg(V) such that (n-1)*A*V' ≡ -P (mod V) via extended GCD.
// Returns updated P and adds -A/((n-1)*V^(n-1)) to rational_part.
static bool ho_reduce_step(
    algebra::PolyExpr& P_poly,     // in/out: numerator poly (updated)
    const algebra::PolyExpr& V_poly,
    unsigned int n,
    ExprPtr V,
    const Symbol& x,
    ExprPtr& rational_part,        // accumulates rational contributions
    const DifferentialField& field,
    symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    // V' as poly
    auto V_expr_res = algebra::polynomial_to_expr(V_poly, x, ctx);
    if (V_expr_res.is_error()) return false;
    auto Vd_res = field.derive_in_generators(V_expr_res.value(), ctx);
    if (Vd_res.is_error()) return false;
    auto Vd_poly_res = algebra::parse_polynomial(Vd_res.value(), x, ctx);
    if (Vd_poly_res.is_error()) return false;
    const auto& Vd_poly = Vd_poly_res.value();

    // (n-1)*V' as poly
    algebra::PolyExpr nmo_poly;
    nmo_poly.push_back(algebra::poly_make_integer(arena, static_cast<long long>(n - 1)));
    auto nmo_Vd_res = algebra::poly_multiply(nmo_poly, Vd_poly, ctx);
    if (nmo_Vd_res.is_error()) return false;

    // Extended GCD: s*(n-1)*V' + t*V = g  (g = constant since gcd(V,V')=1 for squarefree V)
    auto xgcd_res = algebra::poly_extended_gcd(nmo_Vd_res.value(), V_poly, ctx);
    if (xgcd_res.is_error()) return false;
    const auto& xgcd = xgcd_res.value();

    // P mod V
    auto P_rem_res = algebra::divide_poly_with_remainder(P_poly, V_poly, ctx);
    if (P_rem_res.is_error()) return false;
    const auto& P_rem = P_rem_res.value().remainder;

    // A_raw = -P_rem * s
    auto neg_P_res = algebra::poly_negate(P_rem, ctx);
    if (neg_P_res.is_error()) return false;
    auto A_raw_res = algebra::poly_multiply(neg_P_res.value(), xgcd.s, ctx);
    if (A_raw_res.is_error()) return false;

    // A_poly = A_raw mod V  (deg A < deg V)
    auto A_rem_res = algebra::divide_poly_with_remainder(A_raw_res.value(), V_poly, ctx);
    if (A_rem_res.is_error()) return false;
    algebra::PolyExpr A_poly = A_rem_res.value().remainder;

    // Divide by gcd (should be a nonzero constant)
    if (!algebra::is_zero_poly(xgcd.gcd) && algebra::poly_degree(xgcd.gcd) == 0) {
        ExprPtr g_val = algebra::leading_coefficient(xgcd.gcd);
        auto A_scaled = algebra::poly_divide_by_scalar(A_poly, g_val, ctx);
        if (A_scaled.is_error()) return false;
        A_poly = A_scaled.value();
    }

    // rational_part += A / V^(n-1)
    // From H-O derivation: d/dx[A/V^(n-1)] = (A'V - (n-1)AV')/V^n, so P/V^n = d/dx[A/V^(n-1)] + B/V^(n-1)
    auto A_expr_res = algebra::polynomial_to_expr(A_poly, x, ctx);
    if (A_expr_res.is_error()) return false;
    ExprPtr A_expr = A_expr_res.value();
    ExprPtr Vn1 = arena.make<Binary>(BinaryOp::Pow, V, arena.make<IntegerLit>(BigInt(n - 1)));
    ExprPtr contrib = arena.make<Binary>(BinaryOp::Div, A_expr, Vn1);
    rational_part = arena.make<Binary>(BinaryOp::Add, rational_part, contrib);

    // Compute A' (polynomial derivative of A_poly via field derivation)
    auto A_prime_res = field.derive_in_generators(A_expr, ctx);
    if (A_prime_res.is_error()) return false;
    auto A_prime_poly_res = algebra::parse_polynomial(A_prime_res.value(), x, ctx);
    if (A_prime_poly_res.is_error()) return false;
    algebra::PolyExpr A_prime_poly = A_prime_poly_res.value();

    // New P = (P - A'*V + (n-1)*A*V') / V  (H-O standard formula)
    auto A_Vd_res = algebra::poly_multiply(A_poly, Vd_poly, ctx);
    if (A_Vd_res.is_error()) return false;
    auto nmo_A_Vd_res = algebra::poly_multiply(nmo_poly, A_Vd_res.value(), ctx);
    if (nmo_A_Vd_res.is_error()) return false;
    auto new_P_sum_res = algebra::poly_add(P_poly, nmo_A_Vd_res.value(), ctx);
    if (new_P_sum_res.is_error()) return false;
    algebra::PolyExpr new_P_full = new_P_sum_res.value();
    if (!algebra::is_zero_poly(A_prime_poly)) {
        auto A_prime_V_res = algebra::poly_multiply(A_prime_poly, V_poly, ctx);
        if (A_prime_V_res.is_ok()) {
            auto sub_res = algebra::poly_subtract(new_P_full, A_prime_V_res.value(), ctx);
            if (sub_res.is_ok()) new_P_full = sub_res.value();
        }
    }
    auto new_P_div_res = algebra::divide_poly_with_remainder(new_P_full, V_poly, ctx);
    if (new_P_div_res.is_error()) return false;
    P_poly = new_P_div_res.value().quotient;

    return true;
}

/// Performs Hermite reduction on P/Q with respect to generator t_var
Result<HermiteReduction> hermite_reduce(
    ExprPtr P, ExprPtr Q, const Symbol& t_var, const DifferentialField& field, symbolic::CASContext& ctx) {

    // If we are in the base field (no extensions), use the exact RatPoly-based reduction.
    // This is more robust for rational functions over Q(parameters).
    if (field.extensions().empty() && t_var.name == field.base_var().name) {
        return hermite_reduction_exact(P, Q, t_var, ctx);
    }

    AstArena& arena = ctx.arena();
    const Symbol& x = field.base_var();

    // Square-free factorization of Q
    auto sqf_res = algebra::square_free_factorization(Q, x, ctx);
    if (sqf_res.is_error()) {
        return ok(HermiteReduction{
            .rational_part = arena.make<IntegerLit>(BigInt(0)),
            .remaining_P = P,
            .remaining_Q = Q
        });
    }
    const auto& sqf = sqf_res.value();

    bool has_repeated = false;
    for (const auto& f : sqf.factors) {
        if (f.multiplicity >= 2) { has_repeated = true; break; }
    }
    if (!has_repeated) {
        return ok(HermiteReduction{
            .rational_part = arena.make<IntegerLit>(BigInt(0)),
            .remaining_P = P,
            .remaining_Q = Q
        });
    }

    // Parse current numerator as polynomial
    auto P_poly_res = algebra::parse_polynomial(P, x, ctx);
    if (P_poly_res.is_error()) {
        return ok(HermiteReduction{
            .rational_part = arena.make<IntegerLit>(BigInt(0)),
            .remaining_P = P,
            .remaining_Q = Q
        });
    }
    algebra::PolyExpr P_poly = P_poly_res.value();
    ExprPtr rational_part = arena.make<IntegerLit>(BigInt(0));
    ExprPtr current_Q = Q;

    for (const auto& factor : sqf.factors) {
        if (factor.multiplicity < 2) continue;
        ExprPtr V = factor.factor;
        auto V_poly_res = algebra::parse_polynomial(V, x, ctx);
        if (V_poly_res.is_error()) continue;
        const auto& V_poly = V_poly_res.value();

        for (unsigned int n = factor.multiplicity; n >= 2; --n) {
            if (!ho_reduce_step(P_poly, V_poly, n, V, x, rational_part, field, ctx)) break;

            // Remove one power of V from the denominator
            auto Q_poly_res = algebra::parse_polynomial(current_Q, x, ctx);
            if (Q_poly_res.is_error()) break;
            auto Q_div = algebra::divide_poly_with_remainder(Q_poly_res.value(), V_poly, ctx);
            if (Q_div.is_error()) break;
            auto new_Q_res = algebra::polynomial_to_expr(Q_div.value().quotient, x, ctx);
            if (new_Q_res.is_error()) break;
            current_Q = new_Q_res.value();
        }
    }

    auto P_expr_res = algebra::polynomial_to_expr(P_poly, x, ctx);
    ExprPtr remaining_P = P_expr_res.is_ok() ? P_expr_res.value() : P;

    {
        auto s = ctx.simplify(rational_part);
        if (s.is_ok()) rational_part = s.value();
    }
    {
        auto s = ctx.simplify(remaining_P);
        if (s.is_ok()) remaining_P = s.value();
    }
    {
        auto s = ctx.simplify(current_Q);
        if (s.is_ok()) current_Q = s.value();
    }

    return ok(HermiteReduction{
        .rational_part = rational_part,
        .remaining_P = remaining_P,
        .remaining_Q = current_Q
    });
}

Result<ExprPtr> integrate_rothstein_trager(
    ExprPtr P, ExprPtr Q, const Symbol& t_var, const DifferentialField& field, symbolic::CASContext& ctx) {

    AstArena& arena = ctx.arena();

    // A27: the derivative must stay in the generator representation.  A plain
    // derive() of a generator symbol returns D(t_i) in ORIGINAL form (e.g.
    // 1/(x·ln x)), and resultants/GCDs over the resulting mixed expression
    // silently lose factors — this produced the wrong antiderivative 0 for
    // ∫ 1/(x·ln x·ln ln x) dx (root accepted with gcd v = 1, term ln|1| = 0).
    auto dQ_res = field.derive_in_generators(Q, ctx);
    if (dQ_res.is_error()) return fail<ExprPtr>(dQ_res.error());
    ExprPtr dQ = dQ_res.value();

    // t must be a fresh variable not appearing in P or Q.
    const Symbol t_fresh = ctx.make_fresh_symbol("rt_t");
    ExprPtr t_sym = arena.make<Symbol>(t_fresh.name);
    ExprPtr t_dQ = arena.make<Product>(std::vector<ExprPtr>{t_sym, dQ});
    ExprPtr A = arena.make<Binary>(BinaryOp::Sub, P, t_dQ);

    auto R_res = compute_resultant(A, Q, t_var, ctx);
    if (R_res.is_error()) return fail<ExprPtr>(R_res.error());
    ExprPtr R = R_res.value();

    // Check P == 0 early (integral of 0 is 0)
    {
        bool P_zero = false;
        if (const auto* il = expr_cast<IntegerLit>(P)) P_zero = il->value.is_zero();
        if (P_zero) return ok(arena.make<IntegerLit>(BigInt(0)));
    }

    auto roots_res = algebra::solve_polynomial(R, t_fresh, ctx);
    if (roots_res.is_error()) return fail<ExprPtr>(roots_res.error());

    // Deduplicate roots (structural equality)
    std::vector<ExprPtr> unique_roots;
    for (ExprPtr r : roots_res.value()) {
        bool dup = false;
        for (ExprPtr u : unique_roots) { if (structural_equal(r, u)) { dup = true; break; } }
        if (!dup) unique_roots.push_back(r);
    }

    if (unique_roots.empty()) {
        // No rational roots found but P != 0 → integral has non-rational parts (arctan etc.)
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Rothstein-Trager: resultant has no rational roots; integral requires arctan or algebraic extensions",
        });
    }

    std::vector<ExprPtr> integral_terms;
    // (root, v) pairs of the emitted log terms — needed for the
    // hyperexponential residue computation below (Bronstein §5.9).
    std::vector<std::pair<ExprPtr, ExprPtr>> log_args;
    for (ExprPtr root_raw : unique_roots) {
        // A27: Rothstein-Trager is only valid for CONSTANT roots (the c_i in
        // Σ c_i·ln v_i live in Const(K)).  A root depending on the integration
        // variable or on a tower generator is not a valid residue; emitting a
        // term for it produces a silently wrong antiderivative.
        ExprPtr root = root_raw;
        if (auto s = ctx.simplify(root); s.is_ok()) root = s.value();
        if (depends_on(root, field.base_var())) continue;
        bool depends_on_generator = false;
        for (const auto& ext : field.extensions()) {
            if (depends_on(root, ext.t_var)) { depends_on_generator = true; break; }
        }
        if (depends_on_generator) continue;

        ExprPtr root_dQ = arena.make<Product>(std::vector<ExprPtr>{root, dQ});
        ExprPtr A_root = arena.make<Binary>(BinaryOp::Sub, P, root_dQ);
        { auto s = ctx.simplify(A_root); if (s.is_ok()) A_root = s.value(); }

        auto v_res = algebra::polynomial_gcd(A_root, Q, t_var, ctx);
        if (v_res.is_error()) continue;
        ExprPtr v = v_res.value();

        // A27: a gcd that is trivial in t (degree 0) contributes ln(const),
        // i.e. nothing — such a root is spurious for the log part.  Emitting
        // ln|1| terms here masked the "no valid terms" diagnostic below.
        if (!depends_on(v, t_var)) continue;

        // Use ln(|v|) for real-domain correctness
        ExprPtr abs_v = arena.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{v});
        ExprPtr ln_v = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{abs_v});
        ExprPtr term = arena.make<Product>(std::vector<ExprPtr>{root, ln_v});
        integral_terms.push_back(term);
        log_args.emplace_back(root, v);
    }

    if (integral_terms.empty()) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Rothstein-Trager: no valid GCD terms computed (roots may be algebraic/RootOf)",
        });
    }

    // Bronstein §5.9 (IntegrateHyperexponential): when the topmost monomial t
    // is hyperexponential, deg_t(Dv) = deg_t(v), so D(ln v) = Dv/v carries a
    // polynomial part and the residue h − Σ cᵢ·Dvᵢ/vᵢ is in general a NONZERO
    // element of k that must be integrated recursively over k.  Dropping it
    // silently returned −ln|eˣ+1| for ∫ 1/(1+eˣ) dx (correct: x − ln(eˣ+1)).
    // For primitive (log) towers and the rational base case the residue is 0
    // by Theorem 5.6.1, so this block is exp-only by theorem, not shortcut.
    const bool topmost_is_exp = !field.extensions().empty()
        && field.extensions().back().type == ExtensionType::Exponential
        && field.extensions().back().t_var.name == t_var.name;
    if (topmost_is_exp) {
        std::vector<ExprPtr> parts;
        parts.push_back(arena.make<Binary>(BinaryOp::Div, P, Q));
        for (const auto& [root, v] : log_args) {
            auto dv_res = field.derive_in_generators(v, ctx);
            if (dv_res.is_error()) return fail<ExprPtr>(dv_res.error());
            parts.push_back(arena.make<Product>(std::vector<ExprPtr>{
                arena.make<IntegerLit>(BigInt(-1)), root, dv_res.value(),
                arena.make<Binary>(BinaryOp::Pow, v,
                    arena.make<IntegerLit>(BigInt(-1)))}));
        }
        ExprPtr residual = arena.make<Sum>(std::move(parts));
        if (auto tog = algebra::together(residual, ctx); tog.is_ok())
            residual = tog.value();
        if (auto s = ctx.simplify(residual); s.is_ok()) residual = s.value();
        if (depends_on(residual, t_var)) {
            // The RT invariant h − Σ cᵢ·Dvᵢ/vᵢ ∈ k could not be certified
            // (dropped root or simplify shortfall) — bail, never guess.
            return make_unimplemented<ExprPtr>(
                "calculus", "integrate_rothstein_trager",
                "hyperexponential residue still depends on the topmost "
                "generator; log part alone would be silently wrong",
                cas::error::reason_codes::RISCH_EXPONENTIAL_DE,
                "Strengthen together/simplify over the tower or handle the "
                "dropped resultant roots (Bronstein §5.9)",
                "T1-2026-07-16");
        }
        bool residual_zero = false;
        if (const auto* il = expr_cast<IntegerLit>(residual))
            residual_zero = il->value.is_zero();
        if (const auto* rl = expr_cast<RationalLit>(residual))
            residual_zero = rl->numerator.is_zero();
        if (!residual_zero) {
            auto orig_res = field.from_field_generators(residual, ctx);
            if (orig_res.is_error()) return orig_res;
            // residual ∈ k: one extension fewer — recursion terminates.
            auto rec = integrate(orig_res.value(), field.base_var(), ctx);
            if (rec.is_error()) return rec;  // clean bail, never silent
            integral_terms.push_back(rec.value());
        }
    }

    if (integral_terms.size() == 1) {
        return ok(integral_terms[0]);
    }
    return ok(arena.make<Sum>(std::move(integral_terms)));
}

} // namespace cas::calculus
