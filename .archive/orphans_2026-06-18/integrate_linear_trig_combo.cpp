// W3.1 — Closed form for ∫ dx / (a + b·sin(c·x+d) + c'·cos(c·x+d)).
//
// Derivation via Weierstrass t = tan((c·x+d)/2):
//   sin(c·x+d) = 2t/(1+t²),  cos(c·x+d) = (1-t²)/(1+t²),
//   d(c·x+d) = c·dx, so dx = (2/(c·(1+t²))) dt.
// Integrand · dx = 2/(c · [(a+c') + 2b·t + (a-c')·t²]) dt.
//
// Discriminant of the quadratic in t:
//   Δ = 4·b² − 4·(a-c')·(a+c') = 4·(b² − a² + c'²) = −4·(a² − b² − c'²).
//
// Case Δ < 0 (i.e. a² > b² + c'²) — denominator irreducible over Q[t]:
//   ∫ 2 dt/(α·t² + β·t + γ) = (2/√(4αγ − β²)) · arctan((2α·t + β)/√(4αγ − β²)),
//   with α = a−c', β = 2b, γ = a+c'.  After multiplying by 1/c:
//   ∫ dx/(a+b·sin+c'·cos) = (2/(c·√(a²−b²−c'²)))
//                          · arctan(((a−c')·tan((c·x+d)/2) + b) / √(a²−b²−c'²)).
//
// Case Δ ≥ 0 (denominator factors): handled by the standard rational pipeline
// after Weierstrass substitution — this fast-path bails to caller.

#include "integrate_engine.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas::calculus::integrate_detail {

namespace {

struct TrigPart {
    Rational coeff{BigInt(0)};
    ExprPtr arg{nullptr};
};

// Try to pull a single sin/cos(linear_arg) factor times a Rational coefficient
// out of a term.  Returns the coefficient and the arg expression; nullopt if
// the term is not of that shape.
struct ExtractedTrig {
    BuiltinOp kind;
    Rational coeff;
    ExprPtr arg;
};

[[nodiscard]] std::optional<ExtractedTrig> extract_signed_sincos(ExprPtr term) {
    Rational sign(BigInt(1));
    while (true) {
        if (const auto* un = expr_cast<Unary>(term); un && un->op == UnaryOp::Neg) {
            sign = sign * Rational(BigInt(-1));
            term = un->operand;
            continue;
        }
        break;
    }
    Rational coeff(BigInt(1));
    ExprPtr trig_factor = term;
    if (const auto* p = expr_cast<Product>(term)) {
        coeff = Rational(BigInt(1));
        ExprPtr trig = nullptr;
        for (ExprPtr f : p->factors) {
            if (auto s = exact_scalar_from_expr(f); s.has_value()) {
                coeff = coeff * s.value();
                continue;
            }
            if (trig != nullptr) return std::nullopt;
            trig = f;
        }
        if (!trig) return std::nullopt;
        trig_factor = trig;
    } else {
        if (auto s = exact_scalar_from_expr(term); s.has_value()) {
            // Pure constant — not a trig term.
            return std::nullopt;
        }
    }
    const auto* fc = expr_cast<FuncCall>(trig_factor);
    if (!fc || fc->args.size() != 1U) return std::nullopt;
    if (fc->func_id != BuiltinOp::Sin && fc->func_id != BuiltinOp::Cos) return std::nullopt;
    return ExtractedTrig{fc->func_id, sign * coeff, fc->args[0]};
}

}  // namespace

[[nodiscard]] Result<ExprPtr> try_integrate_linear_trig_combo(
    ExprPtr numerator, ExprPtr denominator, const Symbol& var, AstArena& arena,
    symbolic::CASContext& ctx) {

    // Numerator must be a non-zero constant (independent of var).
    auto num_scalar = exact_scalar_from_expr(numerator);
    if (!num_scalar || num_scalar->numerator().is_zero()
        || depends_on(numerator, var))
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_linear_trig_combo",
            "numerator is not a non-zero rational constant",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");

    // Denominator must be a Sum.
    const auto* sum = expr_cast<Sum>(denominator);
    if (!sum)
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_linear_trig_combo",
            "denominator is not a Sum",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");

    Rational a(BigInt(0));      // constant term
    Rational b_sin(BigInt(0));  // coeff of sin(arg)
    Rational c_cos(BigInt(0));  // coeff of cos(arg)
    ExprPtr trig_arg = nullptr;

    for (ExprPtr t : sum->terms) {
        if (!depends_on(t, var)) {
            auto s = exact_scalar_from_expr(t);
            if (!s) return ::cas::make_unimplemented<ExprPtr>(
                "calculus", "try_integrate_linear_trig_combo",
                "non-rational constant term in Sum",
                cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");
            a = a + s.value();
            continue;
        }
        auto trig = extract_signed_sincos(t);
        if (!trig)
            return ::cas::make_unimplemented<ExprPtr>(
                "calculus", "try_integrate_linear_trig_combo",
                "term is not const·sin(linear) or const·cos(linear)",
                cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");
        if (!trig_arg) {
            trig_arg = trig->arg;
        } else if (!structural_equal(trig_arg, trig->arg)) {
            return ::cas::make_unimplemented<ExprPtr>(
                "calculus", "try_integrate_linear_trig_combo",
                "mixed trig arguments",
                cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");
        }
        if (trig->kind == BuiltinOp::Sin) b_sin = b_sin + trig->coeff;
        else                              c_cos = c_cos + trig->coeff;
    }
    if (!trig_arg)
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_linear_trig_combo",
            "Sum has no trig term",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");

    // arg = c · var + d (affine in var), c ≠ 0.
    auto aff = extract_affine_argument(trig_arg, var);
    if (!aff || aff->coefficient.numerator().is_zero())
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_linear_trig_combo",
            "trig arg is not affine in var",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");

    // Compute a² − b² − c'².
    const Rational a_sq = a * a;
    const Rational disc = a_sq - b_sin * b_sin - c_cos * c_cos;
    if (disc.numerator().is_zero() || disc.numerator().is_negative())
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_linear_trig_combo",
            "discriminant ≤ 0 — denominator factors, defer to Weierstrass+rational",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");

    // Emit: (2 · num_scalar / (c_coeff · √disc)) ·
    //       arctan( ((a − c') · tan((c·x+d)/2) + b) / √disc )
    const Rational c_coeff = aff->coefficient;
    ExprPtr disc_expr = make_rational(arena, disc);
    ExprPtr sqrt_disc = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{disc_expr});

    // Emit tan((c·x+d)/2) via the identity tan(θ/2) = sin(θ)/(1+cos(θ)) so
    // the output uses the same full-angle trig functions as the Maxima
    // reference (avoids half-angle mismatch in the cert).
    ExprPtr sin_arg = arena.make<FuncCall>(BuiltinOp::Sin, std::vector<ExprPtr>{trig_arg});
    ExprPtr cos_arg = arena.make<FuncCall>(BuiltinOp::Cos, std::vector<ExprPtr>{trig_arg});
    ExprPtr one_plus_cos = arena.make<Sum>(std::vector<ExprPtr>{make_integer(arena, 1), cos_arg});
    ExprPtr tan_half = arena.make<Binary>(BinaryOp::Div, sin_arg, one_plus_cos);

    // Inner numerator: (a − c') · tan(half) + b.
    const Rational a_minus_c = a - c_cos;
    std::vector<ExprPtr> inner_terms;
    if (!a_minus_c.numerator().is_zero()) {
        inner_terms.push_back(arena.make<Product>(std::vector<ExprPtr>{
            make_rational(arena, a_minus_c), tan_half}));
    }
    if (!b_sin.numerator().is_zero()) {
        inner_terms.push_back(make_rational(arena, b_sin));
    }
    ExprPtr inner_num = inner_terms.empty() ? make_integer(arena, 0)
        : (inner_terms.size() == 1U ? inner_terms[0]
           : arena.make<Sum>(std::move(inner_terms)));
    ExprPtr arctan_arg = arena.make<Binary>(BinaryOp::Div, inner_num, sqrt_disc);
    ExprPtr arctan_expr = arena.make<FuncCall>(BuiltinOp::Atan,
        std::vector<ExprPtr>{arctan_arg});

    // Outer scale: 2 · num / (c_coeff · √disc).
    const Rational two_num = Rational(BigInt(2)) * num_scalar.value() / c_coeff;
    ExprPtr scale_rat = make_rational(arena, two_num);
    ExprPtr scale_full = arena.make<Binary>(BinaryOp::Div, scale_rat, sqrt_disc);

    ExprPtr primitive = arena.make<Product>(std::vector<ExprPtr>{scale_full, arctan_expr});
    if (auto simp = ctx.simplify(primitive); simp.is_ok()) primitive = simp.value();
    return ok(primitive);
}

}  // namespace cas::calculus::integrate_detail
