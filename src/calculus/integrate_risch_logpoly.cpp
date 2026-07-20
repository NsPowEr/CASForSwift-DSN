// integrate_risch_logpoly.cpp — ∫ poly-in-ln(u) dx (Bronstein §5.10).
// Anti-monolith split from integrate_risch_hermite.cpp (zero logic changes
// beyond the split itself): this TU holds integrate_log_polynomial_part();
// integrate_risch_poly_and_rational_part() and its Hermite/Rothstein-Trager
// machinery stay in integrate_risch_hermite.cpp.

#include "integrate_risch_internal.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::calculus {

// Integrate the polynomial-in-t part of a single logarithmic extension
// tower:  ∫ Σ_{k=0..n} a_k(x) * t^k dx,   where t = ln(u(x)),  Dt = u'/u.
//
// Standard ansatz (Bronstein §5.10, IntegratePrimitivePolynomial): the
// antiderivative is again a polynomial in t, of degree ONE MORE than the
// integrand,
//   B(t) = Σ_{k=0..n+1} b_k(x) * t^k,   with
//
//   d/dx B(t)
//     = Σ b_k'(x) * t^k + Σ k * b_k(x) * (u'/u) * t^{k-1}
//     = Σ_{k} [ b_k' + (k+1) * b_{k+1} * (u'/u) ] * t^k.
//
// Matching coefficients with a_k * t^k gives, for every k,
//
//   b_k' + (k+1) * b_{k+1} * η  =  a_k,        η = Dt = u'/u ∈ k.       (§5.10)
//
// The t^{n+1} coefficient forces b_{n+1}' = 0, i.e. b_{n+1} ∈ Const(k): the top
// coefficient of B is a FREE CONSTANT, not zero.  Level k is therefore not a
// plain integration but the LIMITED INTEGRATION PROBLEM (Bronstein §7.2, eq.
// 7.30) in the lower field k:
//
//   a_k - (k+1)*b_{k+1}^{(0)}*η  =  D(b_k) + c*η                        (7.30)
//
// whose solution (b_k, c) also fixes the residual freedom of the level above,
// b_{k+1} += c/(k+1) — adding a constant to b_{k+1} leaves its own equation
// (which only sees b_{k+1}') untouched, so the descent stays consistent.
// This is the A38 wiring: limited_integrate_field() reaches the parametric
// Risch DE tower solver (A1 + A26), which until now no integrate() path called.
//
// The k = 0 level has no constant to determine (an additive constant of B is
// the integration constant), so it is a plain integration in k.
//
// When the limited-integration route is unavailable for a given level (the
// parametric machinery reports Unimplemented, e.g. the deep-tower ConstantSystem
// residue of A1), we fall back to the previous behaviour — plain ∫rhs in the
// lower field, which is the c = 0 special case of (7.30).
// Failures propagate as Unimplemented.
Result<ExprPtr> integrate_log_polynomial_part(
    const algebra::PolyExpr& quot,
    ExprPtr u_arg,
    const Symbol& t_top,
    const Symbol& var,
    const DifferentialField& lower_field,
    symbolic::CASContext& context) {
    AstArena& arena = context.arena();

    if (quot.empty()) {
        return ok(arena.make<IntegerLit>(BigInt(0)));
    }
    const std::size_t deg = quot.size() - 1U;

    // η = u'/u (simplified once up front)
    auto du_res = diff(u_arg, var, 1U, context);
    if (du_res.is_error()) return fail<ExprPtr>(du_res.error());
    ExprPtr du_over_u = arena.make<Binary>(BinaryOp::Div, du_res.value(), u_arg);
    if (auto s = context.simplify(du_over_u); s.is_ok()) du_over_u = s.value();
    // The coefficients a_k live in generator form (gen_expr), so η must too,
    // otherwise the tower solver compares mixed representations (A27).
    ExprPtr eta_gen = du_over_u;
    if (auto g = lower_field.to_field_generators(du_over_u, context); g.is_ok())
        eta_gen = g.value();

    // deg+2 slots: b_{deg+1} is the free constant of the §5.10 ansatz.
    std::vector<ExprPtr> b(deg + 2U, ExprPtr{});

    // Invariant: every slot of b[] is kept in GENERATOR form (like a_k/eta_gen),
    // so corrections never mix representations (A27).  The plain-integrate
    // fallback is the only place that needs the ORIGINAL form (integrate()
    // recognises FuncCall(Ln/Exp), not bare generator symbols); it round-trips
    // through lower_field there and converts the result straight back.
    for (std::ptrdiff_t k = static_cast<std::ptrdiff_t>(deg); k >= 0; --k) {
        const std::size_t kz = static_cast<std::size_t>(k);
        ExprPtr a_k = (kz < quot.size()) ? quot[kz] : ExprPtr{};
        if (!a_k) a_k = arena.make<IntegerLit>(BigInt(0));

        // rhs = a_k  -  (k+1) * b_{k+1} * η   (all terms in generator form)
        ExprPtr rhs = a_k;
        if (b[kz + 1U]) {
            ExprPtr kp1 = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(kz + 1U)));
            ExprPtr correction = arena.make<Product>(std::vector<ExprPtr>{kp1, b[kz + 1U], eta_gen});
            rhs = arena.make<Binary>(BinaryOp::Sub, rhs, correction);
        }
        if (auto s = context.simplify(rhs); s.is_ok()) rhs = s.value();

        if (kz >= 1U) {
            // Limited integration in the lower field: rhs = D(b_k) + c·η.
            auto lim = limited_integrate_field(rhs, {eta_gen}, lower_field, context);
            if (lim.is_ok()) {
                b[kz] = lim.value().v;
                const Rational& c = lim.value().c.front();
                if (!c.numerator().is_zero()) {
                    const Rational share =
                        c / Rational(BigInt(static_cast<std::int64_t>(kz + 1U)));
                    ExprPtr share_e = (share.denominator() == BigInt(1))
                        ? static_cast<ExprPtr>(arena.make<IntegerLit>(share.numerator()))
                        : static_cast<ExprPtr>(arena.make<RationalLit>(
                              share.numerator(), share.denominator()));
                    ExprPtr updated = b[kz + 1U]
                        ? static_cast<ExprPtr>(arena.make<Binary>(
                              BinaryOp::Add, b[kz + 1U], share_e))
                        : share_e;
                    if (auto s = context.simplify(updated); s.is_ok()) updated = s.value();
                    b[kz + 1U] = updated;
                }
                continue;
            }
            // Fall through: c = 0 special case (previous behaviour).
        }

        // Fallback when limited integration is not available for this level
        // (our §7.2 solver is incomplete, so "no solution" here does not yet
        // prove non-integrability the way Bronstein's LimitedIntegrate does).
        //
        // TERMINATION (HC-A38-01, Bronstein §5.2/§5.7): the recursion must move
        // along one of exactly two axes — eliminate a monomial (recurse on
        // k = C(t_1..t_{n-1})), or lower deg_t(p) at fixed tower.  A root-restart
        // on the generic integrate() does neither: it REBUILDS the full tower
        // from the expression, so the sub-problem can be structurally identical
        // to the current one (observed empirically on a nested log-in-log tower:
        // kz=1 → kz=0 → back to kz=1 with the same state, forever).  Dispatching
        // on `lower_field` instead keeps the tower height strictly decreasing,
        // which is the measure the algorithm's termination proof relies on.
        auto compute_b_k = [&]() -> Result<ExprPtr> {
            if (lower_field.extensions().empty()) {
                // k = Q(x): no generators left for a tower rebuild to recreate,
                // so the generic integrator is the base case of the recursion
                // and provably cannot re-enter this function.  It also needs
                // FuncCall(Ln/Exp) form rather than opaque generator symbols.
                ExprPtr rhs_original = rhs;
                if (auto orig = lower_field.from_field_generators(rhs, context); orig.is_ok())
                    rhs_original = orig.value();
                return integrate(rhs_original, var, context);
            }
            // k ⊋ Q(x): stay inside the field-aware pipeline one level down.
            // rhs is already an element of k in generator form, which is exactly
            // the `gen_expr` contract of integrate_risch_poly_and_rational_part.
            ExprPtr zero = arena.make<IntegerLit>(BigInt(0));
            return integrate_risch_poly_and_rational_part(
                rhs, rhs, var, lower_field, zero, context);
        };
        auto b_k_res = compute_b_k();
        if (b_k_res.is_error()) {
            return b_k_res;
        }
        ExprPtr b_k = b_k_res.value();
        if (auto g = lower_field.to_field_generators(b_k, context); g.is_ok())
            b_k = g.value();
        b[kz] = b_k;
    }

    // Build B(t) = Σ b_k * t^k.
    std::vector<ExprPtr> terms;
    terms.reserve(b.size());
    ExprPtr t_sym = arena.make<Symbol>(t_top.name);
    for (std::size_t k = 0; k < b.size(); ++k) {
        if (!b[k]) continue;
        if (const auto* il = expr_cast<IntegerLit>(b[k]); il && il->value.is_zero()) continue;
        ExprPtr term;
        if (k == 0U) {
            term = b[k];
        } else if (k == 1U) {
            term = arena.make<Binary>(BinaryOp::Mul, b[k], t_sym);
        } else {
            ExprPtr t_pow = arena.make<Binary>(
                BinaryOp::Pow,
                t_sym,
                arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k))));
            term = arena.make<Binary>(BinaryOp::Mul, b[k], t_pow);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return ok(arena.make<IntegerLit>(BigInt(0)));
    if (terms.size() == 1U) return ok(terms.front());
    ExprPtr raw = arena.make<Sum>(std::move(terms));
    if (auto s = context.simplify(raw); s.is_ok()) return ok(s.value());
    return ok(raw);
}

} // namespace cas::calculus
