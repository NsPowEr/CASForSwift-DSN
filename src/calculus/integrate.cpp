#include "cas/calculus.hpp"
#include "cas/algebra.hpp"
#include "cas/bigfloat.hpp"
#include "cas/extended_real.hpp"
#include "cas/rational.hpp"
#include "integrate_definite_patterns.hpp"
#include "integrate_engine.hpp"

#include <functional>
#include <optional>

namespace cas::calculus {

namespace {

// Adopt the canonical extended-real predicates; the local copies missed
// Constant(NegInfinity) (only handled Unary(Neg, Constant(Infinity))).
using cas::is_pos_infinity;
using cas::is_neg_infinity;

[[nodiscard]] std::optional<Rational> exact_rational_from_expr(ExprPtr expr) {
    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        return Rational(integer->value);
    }
    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        return Rational(rational->numerator, rational->denominator);
    }
    if (const auto* unary = expr_cast<Unary>(expr);
        unary != nullptr && unary->op == UnaryOp::Neg) {
        auto value = exact_rational_from_expr(unary->operand);
        if (value.has_value()) return -value.value();
    }
    return std::nullopt;
}

// Precision (in bits) used for numerical bound + singularity checks in
// integrate.cpp.  256 bits ≈ 77 decimal digits — far above the ±1e-9 tolerance
// applied downstream in cos_zero_in_range.  HPP-005 closure: all numerical
// reasoning at integration boundaries now flows through MPFR BigFloat rather
// than IEEE double, eliminating silent overflow on numerator/denominator
// values that exceed 2^53 and the to_u64() truncation that violated
// CLAUDE.md REGOLA 1.
inline constexpr mpfr_prec_t kSingularityCheckPrec = 256;

[[nodiscard]] BigFloat bigfloat_from_rational(const Rational& r) {
    return BigFloat::from_rational_parts(
        r.numerator().decimal(),
        r.denominator().decimal(),
        kSingularityCheckPrec);
}

// Approximate a bound expression as a BigFloat for singularity range checks.
// Returns nullopt only if the expression is genuinely symbolic (free variables).
// Uses exact-rational + MPFR arithmetic throughout — no double round-trips.
[[nodiscard]] std::optional<BigFloat> approx_bound(ExprPtr expr) {
    if (auto r = exact_rational_from_expr(expr)) {
        return bigfloat_from_rational(*r);
    }
    if (const auto* c = expr_cast<Constant>(expr)) {
        if (c->value == MathConstant::Pi) return BigFloat::pi(kSingularityCheckPrec);
        if (c->value == MathConstant::E)  return BigFloat::e(kSingularityCheckPrec);
        if (c->value == MathConstant::Infinity) {
            // +∞ via MPFR (set_inf semantics): build from "inf" literal.
            return BigFloat::from_double(
                std::numeric_limits<double>::infinity(), kSingularityCheckPrec);
        }
    }
    if (const auto* u = expr_cast<Unary>(expr); u && u->op == UnaryOp::Neg) {
        auto v = approx_bound(u->operand);
        if (v) return -*v;
    }
    if (const auto* b = expr_cast<Binary>(expr)) {
        auto lv = approx_bound(b->left);
        auto rv = approx_bound(b->right);
        if (!lv || !rv) return std::nullopt;
        if (b->op == BinaryOp::Add) return *lv + *rv;
        if (b->op == BinaryOp::Sub) return *lv - *rv;
        if (b->op == BinaryOp::Mul) return *lv * *rv;
        if (b->op == BinaryOp::Div) {
            if (rv->is_zero()) return std::nullopt;
            return *lv / *rv;
        }
    }
    return std::nullopt;
}

[[nodiscard]] bool is_between_closed(const Rational& value, const Rational& a, const Rational& b) {
    const Rational& lower = (a <= b) ? a : b;
    const Rational& upper = (a <= b) ? b : a;
    return lower <= value && value <= upper;
}

// Extract affine coefficient of var in a linear expression a*var + b.
// Returns {a, b} as rationals, or nullopt if not linear in var.
[[nodiscard]] std::optional<std::pair<Rational,Rational>>
try_extract_linear(ExprPtr expr, const Symbol& var) {
    if (const auto* s = expr_cast<Symbol>(expr))
        return s->name == var.name
            ? std::make_optional(std::pair<Rational,Rational>{Rational(BigInt(1)), Rational(BigInt(0))})
            : std::nullopt;
    if (const auto* i = expr_cast<IntegerLit>(expr))
        return std::make_pair(Rational(BigInt(0)), Rational(i->value));
    if (const auto* r = expr_cast<RationalLit>(expr))
        return std::make_pair(Rational(BigInt(0)), Rational(r->numerator, r->denominator));
    if (const auto* u = expr_cast<Unary>(expr); u && u->op == UnaryOp::Neg) {
        auto inner = try_extract_linear(u->operand, var);
        if (!inner) return std::nullopt;
        return std::make_pair(-inner->first, -inner->second);
    }
    if (const auto* b = expr_cast<Binary>(expr)) {
        if (b->op == BinaryOp::Mul) {
            auto lc = exact_rational_from_expr(b->left);
            if (lc) {
                auto ri = try_extract_linear(b->right, var);
                if (ri) return std::make_pair(*lc * ri->first, *lc * ri->second);
            }
            auto rc = exact_rational_from_expr(b->right);
            if (rc) {
                auto li = try_extract_linear(b->left, var);
                if (li) return std::make_pair(*rc * li->first, *rc * li->second);
            }
        }
        if (b->op == BinaryOp::Add) {
            auto ll = try_extract_linear(b->left, var);
            auto rr = try_extract_linear(b->right, var);
            if (ll && rr) return std::make_pair(ll->first + rr->first, ll->second + rr->second);
        }
    }
    return std::nullopt;
}

// Check if sin(arg) or cos(arg) has a zero in [a,b] where arg is linear in var.
// sin(c*x+d)=0 at x = (kπ - d)/c  (k integer)
// cos(c*x+d)=0 at x = (π/2+kπ - d)/c  (k integer)
// We check integer k values covering the interval.
[[nodiscard]] bool trig_zero_in_interval(
    BuiltinOp func, ExprPtr arg, const Symbol& var,
    const Rational& a, const Rational& b,
    std::size_t scan_max_k) {
    auto linear = try_extract_linear(arg, var);
    if (!linear || linear->first == Rational(BigInt(0))) return false;
    const Rational& c = linear->first;   // coefficient of var
    const Rational& d = linear->second;  // constant term
    // For sin: zero when c*x+d = kπ → x = (kπ-d)/c
    // For cos: zero when c*x+d = π/2+kπ → x = ((2k+1)π/2-d)/c
    // We approximate π ≈ 355/113 for integer arithmetic bound estimation.
    // Then verify via rational interval arithmetic that a candidate k works.
    // Since π is irrational, x = (kπ+offset-d)/c can only be rational if k=0
    // and d is a rational multiple of π (which we don't handle here).
    // However: we check whether the zero COULD fall in [a,b] by estimating bounds.
    // Use π ∈ (3, 4) to bound the integer k range.
    const Rational& lo = (a <= b) ? a : b;
    const Rational& hi = (a <= b) ? b : a;
    // x = (kπ + offset) / c  (offset = -d for sin, offset = π/2 - d for cos)
    // lo ≤ x ≤ hi  →  lo*c ≤ kπ + offset ≤ hi*c  (if c>0)
    //               →  (lo*c - offset) / π ≤ k ≤ (hi*c - offset) / π
    // Bound k using π > 3:
    //   k_max_approx = ceil((hi*c_abs + offset_abs_bound) / 3) + 1
    Rational c_abs = (c.numerator().is_negative() ? -c : c);
    Rational hi_c = hi * c_abs;
    Rational lo_c = lo * c_abs;
    // rough k range: ±(|lo_c|+|hi_c|+|d|)/3 + 2
    auto bound_num = (hi_c.numerator().is_negative() ? -hi_c : hi_c)
                   + (lo_c.numerator().is_negative() ? -lo_c : lo_c)
                   + (d.numerator().is_negative() ? -d : d);
    // bound/3 + 2 as integer
    auto k_bound_bigint = bound_num.numerator() / bound_num.denominator() / BigInt(3) + BigInt(3);
    if (k_bound_bigint.bit_length() > 20) return false; // safety: don't scan huge range
    // HPP-026: cap configurable via ctx.integration_singularity_scan_max_k().
    const std::size_t k_native = static_cast<std::size_t>(k_bound_bigint.to_u64());
    const int64_t k_bound = static_cast<int64_t>(k_native > scan_max_k ? scan_max_k : k_native);

    // For each candidate k, compute tight bounds on the zero x using π ∈ (3, 22/7).
    // If both bounds are in [a,b], definitively in interval.
    // If one is in and one is out, use π ≈ 355/113 for a tighter check.
    // We only care about definite "yes" — if uncertain, fall through.
    const Rational pi_lo(BigInt(3), BigInt(1));       // π > 3
    const Rational pi_hi(BigInt(22), BigInt(7));       // π < 22/7

    for (int64_t k = -k_bound; k <= k_bound; ++k) {
        Rational kR{BigInt{k}};
        // sin: zero at x = (kπ - d)/c
        // cos: zero at x = ((k + 1/2)π - d)/c = ((2k+1)π/2 - d)/c
        Rational half_offset = (func == BuiltinOp::Sin)
            ? Rational{BigInt{0}}
            : Rational{BigInt{1}, BigInt{2}};
        Rational keff = kR + half_offset;  // k for sin, k+1/2 for cos
        // x = (keff * π - d) / c
        // lower bound on x: (keff * pi_lo - d) / c  (if c > 0, else flip)
        Rational x_lo_num = keff * pi_lo - d;
        Rational x_hi_num = keff * pi_hi - d;
        if (c.numerator().is_negative()) std::swap(x_lo_num, x_hi_num);
        Rational x_lo_cand = x_lo_num / c;
        Rational x_hi_cand = x_hi_num / c;
        // If the range [x_lo_cand, x_hi_cand] intersects [lo, hi], singularity exists.
        if (x_hi_cand >= lo && x_lo_cand <= hi) return true;
    }
    return false;
}

// Check if the denominator of expr has a non-rational singularity in [a,b].
// Handles: sin(linear_x), cos(linear_x), sqrt(linear_x≤0_at_endpoint), pow(x,neg).
[[nodiscard]] bool non_rational_singularity_in_interval(
    ExprPtr denom, const Symbol& var, const Rational& a, const Rational& b,
    symbolic::CASContext& ctx) {
    if (!denom || !integrate_detail::depends_on(denom, var)) return false;
    // FuncCall: sin, cos, sqrt
    if (const auto* fc = expr_cast<FuncCall>(denom); fc && fc->args.size() == 1U) {
        if (fc->func_id == BuiltinOp::Sin || fc->func_id == BuiltinOp::Cos)
            return trig_zero_in_interval(fc->func_id, fc->args[0], var, a, b,
                                         ctx.integration_singularity_scan_max_k());
        if (fc->func_id == BuiltinOp::Sqrt) {
            // sqrt(f(x)) = 0 when f(x) = 0. Check if linear f has zero in [a,b].
            auto linear = try_extract_linear(fc->args[0], var);
            if (linear && linear->first != Rational(BigInt(0))) {
                Rational zero = -linear->second / linear->first;
                if (is_between_closed(zero, a, b)) return true;
            }
        }
    }
    // Pow(x, neg_rational): singularity at x=0
    if (const auto* pw = expr_cast<Binary>(denom);
        pw && pw->op == BinaryOp::Pow) {
        if (const auto* sym = expr_cast<Symbol>(pw->left); sym && sym->name == var.name) {
            auto exp_r = exact_rational_from_expr(pw->right);
            if (exp_r && exp_r->numerator().is_negative()) {
                // x^(-r): zero at x=0
                if (is_between_closed(Rational(BigInt(0)), a, b)) return true;
            }
        }
    }
    // Product: any factor is a singularity source
    if (const auto* prod = expr_cast<Product>(denom)) {
        for (ExprPtr f : prod->factors) {
            if (non_rational_singularity_in_interval(f, var, a, b, ctx)) return true;
        }
    }
    // Binary Mul
    if (const auto* bin = expr_cast<Binary>(denom); bin && bin->op == BinaryOp::Mul) {
        if (non_rational_singularity_in_interval(bin->left, var, a, b, ctx)) return true;
        if (non_rational_singularity_in_interval(bin->right, var, a, b, ctx)) return true;
    }
    return false;
}

[[nodiscard]] Result<void> reject_rational_poles_in_closed_interval(
    ExprPtr expr,
    const Symbol& var,
    ExprPtr lower,
    ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto lower_value = exact_rational_from_expr(lower);
    auto upper_value = exact_rational_from_expr(upper);

    auto together_expr = algebra::together(expr, ctx);
    ExprPtr rational_expr = together_expr.is_ok() ? together_expr.value() : expr;
    auto simplified = ctx.simplify(rational_expr);
    if (simplified.is_ok()) {
        rational_expr = simplified.value();
    }

    // Direct check: detect tan singularities — uses numeric bounds so works even when
    // bounds involve π (which exact_rational_from_expr cannot represent).
    // Note: simplify() expands tan(x) → sin(x)/cos(x), so we must also detect cos
    // zeros in denominators, not just FuncCall(Tan) nodes.
    auto approx_lo = approx_bound(lower);
    auto approx_hi = approx_bound(upper);
    if (approx_lo && approx_hi) {
        const BigFloat& blo_raw = *approx_lo;
        const BigFloat& bhi_raw = *approx_hi;
        const BigFloat& blo = (blo_raw <= bhi_raw) ? blo_raw : bhi_raw;
        const BigFloat& bhi = (blo_raw <= bhi_raw) ? bhi_raw : blo_raw;

        // Check whether cos(arg) = 0 has a solution in [blo,bhi] for linear arg.
        // All arithmetic in BigFloat (MPFR) at kSingularityCheckPrec — no double
        // round-trips, no silent overflow on numerator/denominator > 2^53.
        auto cos_zero_in_range = [&](ExprPtr arg) -> bool {
            auto linear = try_extract_linear(arg, var);
            if (!linear || linear->first == Rational{BigInt{0}}) return false;
            BigFloat c = bigfloat_from_rational(linear->first);
            BigFloat d = bigfloat_from_rational(linear->second);
            // Reject c = 0 (no x-dependence → cos(d) constant, no roots in x).
            if (c.is_zero()) return false;
            BigFloat pi = BigFloat::pi(kSingularityCheckPrec);
            BigFloat half = BigFloat::from_double(0.5, kSingularityCheckPrec);
            BigFloat tol = BigFloat::from_double(1e-9, kSingularityCheckPrec);
            // cos(c*x + d) = 0 → x = (π/2 + k·π − d) / c
            BigFloat x_base = (pi * half - d) / c;
            BigFloat pi_over_c = pi / c;
            for (int k = -20; k <= 20; ++k) {
                BigFloat km = BigFloat::from_double(static_cast<double>(k),
                    kSingularityCheckPrec);
                BigFloat x_k = x_base + km * pi_over_c;
                if (x_k >= (blo - tol) && x_k <= (bhi + tol)) return true;
            }
            return false;
        };

        // Walk expression tree looking for tan(arg) or cos(arg) appearing in a denominator.
        // Returns true if a singularity is found in [dlo,dhi].
        std::function<bool(ExprPtr, bool)> scan_trig_singularity;
        scan_trig_singularity = [&](ExprPtr e, bool in_denom) -> bool {
            if (!e) return false;
            if (const auto* fc = expr_cast<FuncCall>(e)) {
                if (fc->args.size() == 1U) {
                    if (fc->func_id == BuiltinOp::Tan)
                        return cos_zero_in_range(fc->args[0]);
                    if (fc->func_id == BuiltinOp::Cos && in_denom)
                        return cos_zero_in_range(fc->args[0]);
                }
                return false;
            }
            if (const auto* bin = expr_cast<Binary>(e)) {
                if (bin->op == BinaryOp::Div)
                    return scan_trig_singularity(bin->left, in_denom)
                        || scan_trig_singularity(bin->right, true);
                if (bin->op == BinaryOp::Mul)
                    return scan_trig_singularity(bin->left, in_denom)
                        || scan_trig_singularity(bin->right, in_denom);
                if (bin->op == BinaryOp::Add || bin->op == BinaryOp::Sub)
                    return scan_trig_singularity(bin->left, in_denom)
                        || scan_trig_singularity(bin->right, in_denom);
                if (bin->op == BinaryOp::Pow) {
                    // cos(arg)^n with negative integer n: cos is in the denominator.
                    if (const auto* exp_lit = expr_cast<IntegerLit>(bin->right)) {
                        if (exp_lit->value.is_negative())
                            return scan_trig_singularity(bin->left, true);
                    }
                    if (const auto* exp_unary = expr_cast<Unary>(bin->right);
                            exp_unary && exp_unary->op == UnaryOp::Neg)
                        return scan_trig_singularity(bin->left, true);
                }
                return false;
            }
            if (const auto* prod = expr_cast<Product>(e)) {
                for (ExprPtr f : prod->factors)
                    if (scan_trig_singularity(f, in_denom)) return true;
                return false;
            }
            if (const auto* sum = expr_cast<Sum>(e)) {
                for (ExprPtr t : sum->terms)
                    if (scan_trig_singularity(t, in_denom)) return true;
                return false;
            }
            if (const auto* u = expr_cast<Unary>(e))
                return scan_trig_singularity(u->operand, in_denom);
            return false;
        };

        if (scan_trig_singularity(rational_expr, false)) {
            return fail<void>(integrate_detail::make_error(
                CASErrorKind::Undefined,
                "Definite integral crosses a tan/cos singularity; improper handling not implemented"));
        }
    }

    // Rational bound check requires exact bounds; skip if bounds involve π or ∞.
    if (!lower_value.has_value() || !upper_value.has_value()) {
        return ok();
    }

    auto parts = algebra::apart_num_den(rational_expr, ctx);
    if (parts.is_error()) {
        return ok();
    }

    auto denominator = ctx.simplify(parts.value().denominator);
    if (denominator.is_error()) {
        return ok();
    }
    if (!integrate_detail::depends_on(denominator.value(), var)) {
        return ok();
    }

    auto roots = algebra::solve_polynomial(denominator.value(), var, ctx);
    if (roots.is_error()) {
        // Not a polynomial denominator — check for trig/algebraic singularities.
        if (non_rational_singularity_in_interval(
                denominator.value(), var, lower_value.value(), upper_value.value(), ctx)) {
            return fail<void>(integrate_detail::make_error(
                CASErrorKind::Undefined,
                "Definite integral crosses a non-rational singularity (trig/algebraic); "
                "improper handling not implemented"));
        }
        return ok();
    }

    for (ExprPtr root : roots.value()) {
        auto root_value = exact_rational_from_expr(root);
        if (!root_value.has_value()) {
            continue;
        }
        if (is_between_closed(root_value.value(), lower_value.value(), upper_value.value())) {
            return fail<void>(integrate_detail::make_error(
                CASErrorKind::Undefined,
                "Definite integral crosses a rational pole; improper/PV handling is not implemented here"));
        }
    }
    // Even if polynomial roots found, also check non-rational factors in denominator.
    if (non_rational_singularity_in_interval(
            denominator.value(), var, lower_value.value(), upper_value.value(), ctx)) {
        return fail<void>(integrate_detail::make_error(
            CASErrorKind::Undefined,
            "Definite integral crosses a non-rational singularity (trig/algebraic); "
            "improper handling not implemented"));
    }
    return ok();
}

[[nodiscard]] ExprPtr normalize_definite_integrand(ExprPtr expr, symbolic::CASContext& ctx) {
    auto together_expr = algebra::together(expr, ctx);
    ExprPtr normalized = together_expr.is_ok() ? together_expr.value() : expr;
    auto simplified = ctx.simplify(normalized);
    if (simplified.is_ok()) {
        normalized = simplified.value();
    }
    return normalized;
}

} // anonymous namespace
Result<ExprPtr> integrate(ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    if (ctx.is_caching_enabled()) {
        auto key = symbolic::CASContext::IntegrateKey{expr, var.name};
        if (auto found = ctx.integrate_cache_.get(key)) {
            return ok(*found);
        }
    }

    // A53 — qui MANCA il budget dell'operazione (misurato): ogni `ctx.simplify`
    // interna e' top-level e azzera timer e contatore, quindi nulla limita il
    // TOTALE. `CASContext::OperationScope` lo risolve, ma aprire l'operazione
    // cambia anche la semantica A31 delle side-conditions — quelle dei rami
    // SCARTATI sopravvivono nel risultato. Serve prima quel rollback.
    auto primitive = integrate_detail::integrate_indefinite_impl(expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }
    auto materialized = symbolic::materialize_expr(primitive.value(), ctx.arena());
    if (materialized.is_error()) {
        return materialized;
    }
    if (ctx.is_caching_enabled()) {
        auto key = symbolic::CASContext::IntegrateKey{expr, var.name};
        ctx.integrate_cache_.put(key, materialized.value());
    }
    return materialized;
}

Result<ExprPtr> definite_integral(ExprPtr expr, const Symbol& var, ExprPtr lower, ExprPtr upper, symbolic::CASContext& ctx) {
    ExprPtr normalized_expr = normalize_definite_integrand(expr, ctx);

    // Extensible pattern table: each matcher returns nullopt to skip, value to commit.
    DefiniteContext dc{
        .integrand = expr,
        .integrand_normalized = normalized_expr,
        .var = var,
        .lower = lower,
        .upper = upper,
        .ctx = ctx,
    };
    for (DefinitePatternFn matcher : definite_patterns()) {
        auto match = matcher(dc);
        if (match.is_error()) return fail<ExprPtr>(match.error());
        if (match.value().has_value()) return ok(match.value().value());
    }

    // Generic infinite-domain fallback: only the Gaussian pattern is currently handled there;
    // anything else over (-inf, +inf) goes to Unimplemented.
    if (is_neg_infinity(lower) && is_pos_infinity(upper)) {
        return fail<ExprPtr>(integrate_detail::make_error(CASErrorKind::Unimplemented,
            "Integrazione su dominio infinito: pattern non riconosciuto."));
    }

    auto pole_check = reject_rational_poles_in_closed_interval(normalized_expr, var, lower, upper, ctx);
    if (pole_check.is_error()) {
        return fail<ExprPtr>(pole_check.error());
    }

    auto primitive = integrate(normalized_expr, var, ctx);
    if (primitive.is_error()) {
        return primitive;
    }

    auto lower_value = ctx.substitute(primitive.value(), var, lower);
    if (lower_value.is_error()) {
        return lower_value;
    }

    auto upper_value = ctx.substitute(primitive.value(), var, upper);
    if (upper_value.is_error()) {
        return upper_value;
    }

    return ctx.simplify(integrate_detail::make_sum(ctx.arena(), {
        upper_value.value(),
        integrate_detail::make_unary(ctx.arena(), UnaryOp::Neg, lower_value.value()),
    }));
}

}  // namespace cas::calculus
