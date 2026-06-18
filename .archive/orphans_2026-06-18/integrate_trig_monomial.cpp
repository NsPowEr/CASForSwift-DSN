// W3.1 — Trig-monomial integrator.
//
// Generic algorithm for ∫ sin^a(c·x+d) · cos^b(c·x+d) dx with a, b ∈ Z
// (positive or negative integer exponents).
//
// Dispatched from `Integrator::integrate_product` BEFORE the Weierstrass
// substitution so that polynomial trig integrands like sin^3·cos^2 do not
// produce the enormous tan(x/2)-form output that the Weierstrass pipeline
// emits for non-rational cases.
//
// Two structural cases (no per-input lookup — general algorithm):
//
//   1. If `a` is a positive odd integer, substitute u = cos(c·x+d).
//      sin^(a-1) = (1-u^2)^((a-1)/2) is a polynomial in u of degree a-1.
//      ∫ sin^a · cos^b dx = (-1/c) ∫ (1-u^2)^((a-1)/2) · u^b du.
//      The right-hand side is a polynomial (or Laurent polynomial when b<0)
//      that is handed back to the standard integrator on the fresh symbol u.
//
//   2. Symmetric: `b` positive odd → u = sin(c·x+d).
//
//   3. Otherwise (both even, including negative), substitute t = tan(c·x+d).
//      sin^2 = t^2/(1+t^2),  cos^2 = 1/(1+t^2),  dx = dt / (c·(1+t^2)).
//      Integrand becomes t^a / (c · (1+t^2)^((a+b)/2 + 1)) — a rational
//      function of t — and is delegated to the standard integrator.
//
// All three cases produce a rational(u) or rational(t) integrand that the
// existing partial-fractions / power dispatch already handles correctly.

#include "integrate_engine.hpp"
#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas::calculus::integrate_detail {

namespace {

// Identify a single Product factor as one of:
//   sin(arg)          → kind=Sin, power=+1
//   cos(arg)          → kind=Cos, power=+1
//   sin(arg)^k        → kind=Sin, power=k  (k integer, can be negative)
//   cos(arg)^k        → kind=Cos, power=k
// Anything else returns std::nullopt and the caller falls back to constants.
struct TrigFactor {
    enum Kind { Sin, Cos };
    Kind kind;
    long long power;
    ExprPtr arg;
};

[[nodiscard]] std::optional<TrigFactor> classify_trig_factor(ExprPtr expr) {
    // Direct sin(arg) or cos(arg).
    if (const auto* fc = expr_cast<FuncCall>(expr); fc && fc->args.size() == 1U) {
        if (fc->func_id == BuiltinOp::Sin) return TrigFactor{TrigFactor::Sin, 1, fc->args[0]};
        if (fc->func_id == BuiltinOp::Cos) return TrigFactor{TrigFactor::Cos, 1, fc->args[0]};
    }
    // Pow(sin(arg), k) or Pow(cos(arg), k).
    if (const auto* pw = expr_cast<Binary>(expr); pw && pw->op == BinaryOp::Pow) {
        const auto* il = expr_cast<IntegerLit>(pw->right);
        if (!il) return std::nullopt;
        // Cap absolute exponent at 64 — beyond this the binomial expansion is
        // not worth the cost; let the Weierstrass / Risch paths take over.
        const BigInt abs_v = il->value.abs();
        if (abs_v > BigInt(64)) return std::nullopt;
        const long long k = il->value.is_negative()
            ? -static_cast<long long>(abs_v.to_u64())
            : static_cast<long long>(abs_v.to_u64());
        if (k == 0) return std::nullopt;
        if (const auto* fc = expr_cast<FuncCall>(pw->left); fc && fc->args.size() == 1U) {
            if (fc->func_id == BuiltinOp::Sin) return TrigFactor{TrigFactor::Sin, k, fc->args[0]};
            if (fc->func_id == BuiltinOp::Cos) return TrigFactor{TrigFactor::Cos, k, fc->args[0]};
        }
    }
    return std::nullopt;
}

// (1 - u^2)^m as an explicit polynomial Sum of Pow(u, 2k) terms, m ≥ 0.
[[nodiscard]] ExprPtr expand_one_minus_u_sq_pow(long long m, ExprPtr u_sym, AstArena& arena) {
    // (1 - u^2)^m = Σ_{k=0..m} C(m,k) · (-1)^k · u^(2k).
    auto binom = [](long long n, long long k) -> long long {
        if (k < 0 || k > n) return 0;
        if (k > n - k) k = n - k;
        long long r = 1;
        for (long long i = 0; i < k; ++i) {
            r *= (n - i);
            r /= (i + 1);
        }
        return r;
    };
    std::vector<ExprPtr> terms;
    terms.reserve(static_cast<std::size_t>(m + 1));
    for (long long k = 0; k <= m; ++k) {
        const long long c = ((k & 1LL) ? -1LL : 1LL) * binom(m, k);
        if (c == 0) continue;
        ExprPtr u_pow;
        if (k == 0) {
            u_pow = make_integer(arena, 1);
        } else {
            u_pow = arena.make<Binary>(BinaryOp::Pow, u_sym,
                make_integer(arena, 2 * k));
        }
        if (c == 1) {
            terms.push_back(u_pow);
        } else if (c == -1) {
            terms.push_back(arena.make<Unary>(UnaryOp::Neg, u_pow));
        } else {
            terms.push_back(arena.make<Product>(std::vector<ExprPtr>{
                make_integer(arena, c), u_pow}));
        }
    }
    if (terms.empty()) return make_integer(arena, 0);
    if (terms.size() == 1U) return terms[0];
    return arena.make<Sum>(std::move(terms));
}

// Build a · u^k where k is an integer (positive, zero, or negative).
[[nodiscard]] ExprPtr u_power(ExprPtr u_sym, long long k, AstArena& arena) {
    if (k == 0) return make_integer(arena, 1);
    if (k == 1) return u_sym;
    return arena.make<Binary>(BinaryOp::Pow, u_sym, make_integer(arena, k));
}

// Build (1 + u_sym^2) symbolically.
[[nodiscard]] ExprPtr one_plus_u_sq(ExprPtr u_sym, AstArena& arena) {
    return arena.make<Sum>(std::vector<ExprPtr>{
        make_integer(arena, 1),
        arena.make<Binary>(BinaryOp::Pow, u_sym, make_integer(arena, 2))});
}

// Multiply an arbitrary expression by 1/c using a Rational coefficient.
[[nodiscard]] ExprPtr scale_by_inverse(ExprPtr expr, const Rational& c, AstArena& arena) {
    if (c.numerator() == BigInt(1) && c.denominator() == BigInt(1)) return expr;
    const Rational inv(c.denominator(), c.numerator());
    return arena.make<Product>(std::vector<ExprPtr>{
        make_rational(arena, inv), expr});
}

}  // namespace

// Public entry. Returns ok(...) when the product is a sin/cos monomial with
// matching linear argument and we successfully integrate; returns
// Unimplemented otherwise so the caller can fall through.
[[nodiscard]] Result<ExprPtr> try_integrate_trig_monomial(
    const Product& product, const Symbol& var, AstArena& arena,
    symbolic::CASContext& ctx) {

    long long sin_pow = 0;
    long long cos_pow = 0;
    ExprPtr trig_arg = nullptr;
    std::vector<ExprPtr> consts;

    for (ExprPtr f : product.factors) {
        if (!depends_on(f, var)) { consts.push_back(f); continue; }
        auto cls = classify_trig_factor(f);
        if (!cls) {
            // Any var-dependent non-trig factor: not a pure trig monomial.
            return ::cas::make_unimplemented<ExprPtr>(
                "calculus", "try_integrate_trig_monomial",
                "product contains var-dependent non-sin/cos factor",
                cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
                "Use a different strategy",
                "W3.1");
        }
        if (trig_arg == nullptr) {
            trig_arg = cls->arg;
        } else {
            // All trig args must be the same linear function of var.
            if (!structural_equal(trig_arg, cls->arg)) {
                return ::cas::make_unimplemented<ExprPtr>(
                    "calculus", "try_integrate_trig_monomial",
                    "trig factors have different arguments",
                    cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
                    "Try simplification first",
                    "W3.1");
            }
        }
        if (cls->kind == TrigFactor::Sin) sin_pow += cls->power;
        else                              cos_pow += cls->power;
    }

    if (trig_arg == nullptr) {
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_trig_monomial",
            "no trig factor present",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
            "Use a different strategy",
            "W3.1");
    }
    if (sin_pow == 0 && cos_pow == 0) {
        // sin^0·cos^0 = 1 — fall through (constant times 1).
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_trig_monomial",
            "degenerate sin^0·cos^0",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
            "",
            "W3.1");
    }

    // Argument must be affine in var: c·var + d, c ≠ 0.
    auto aff = extract_affine_argument(trig_arg, var);
    if (!aff || aff->coefficient.numerator().is_zero()) {
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_trig_monomial",
            "trig argument is not affine in var",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
            "Try u-substitution",
            "W3.1");
    }
    const Rational c_coeff = aff->coefficient;

    Symbol u = ctx.make_fresh_symbol("u");
    ExprPtr u_sym = arena.make<Symbol>(u.name);

    ExprPtr integrand_in_u = nullptr;

    // Case 1: sin_pow odd positive → u = cos(arg). du = -c·sin·dx.
    //   ∫ sin^a · cos^b dx = (-1/c) ∫ (1-u²)^((a-1)/2) · u^b du.
    if (sin_pow > 0 && (sin_pow % 2) == 1) {
        const long long m = (sin_pow - 1) / 2;
        ExprPtr poly = expand_one_minus_u_sq_pow(m, u_sym, arena);
        ExprPtr u_b = u_power(u_sym, cos_pow, arena);
        ExprPtr prod = arena.make<Product>(std::vector<ExprPtr>{poly, u_b});
        integrand_in_u = arena.make<Unary>(UnaryOp::Neg,
            scale_by_inverse(prod, c_coeff, arena));
    }
    // Case 2: cos_pow odd positive → u = sin(arg). du = c·cos·dx.
    //   ∫ sin^a · cos^b dx = (1/c) ∫ u^a · (1-u²)^((b-1)/2) du.
    else if (cos_pow > 0 && (cos_pow % 2) == 1) {
        const long long m = (cos_pow - 1) / 2;
        ExprPtr poly = expand_one_minus_u_sq_pow(m, u_sym, arena);
        ExprPtr u_a = u_power(u_sym, sin_pow, arena);
        ExprPtr prod = arena.make<Product>(std::vector<ExprPtr>{u_a, poly});
        integrand_in_u = scale_by_inverse(prod, c_coeff, arena);
    }
    // Case 4: both positive even ≥0 → power-reduction Fourier form.
    //   sin²θ = (1-cos(2θ))/2, cos²θ = (1+cos(2θ))/2.
    //   sin^a · cos^b = ((1-c)/2)^(a/2) · ((1+c)/2)^(b/2) with c = cos(2·arg).
    //   Expand the polynomial in c, then reduce each c^n to a Fourier sum
    //   ∑_m β_{n,m}·cos(m·2·arg) via the product-to-sum recurrence
    //   cos(y)·cos(m·y) = (cos((m+1)y) + cos((m-1)y))/2 (with cos(0)=1).
    //   Integrate term by term: ∫cos(m·2·arg) dx = sin(m·2·arg)/(2m·c_coeff).
    else if (sin_pow >= 0 && cos_pow >= 0
             && (sin_pow % 2) == 0 && (cos_pow % 2) == 0
             && (sin_pow + cos_pow) > 0) {
        const long long j = sin_pow / 2;
        const long long k = cos_pow / 2;
        // Polynomial in c = cos(2·arg): coefs[i] = integer coefficient of c^i.
        std::vector<BigInt> coefs{BigInt(1)};
        for (long long ii = 0; ii < j; ++ii) {
            std::vector<BigInt> r(coefs.size() + 1U, BigInt(0));
            for (std::size_t i = 0; i < coefs.size(); ++i) {
                r[i] = r[i] + coefs[i];
                r[i + 1U] = r[i + 1U] - coefs[i];
            }
            coefs = std::move(r);
        }
        for (long long ii = 0; ii < k; ++ii) {
            std::vector<BigInt> r(coefs.size() + 1U, BigInt(0));
            for (std::size_t i = 0; i < coefs.size(); ++i) {
                r[i] = r[i] + coefs[i];
                r[i + 1U] = r[i + 1U] + coefs[i];
            }
            coefs = std::move(r);
        }
        // Build Fourier representation of (1/2^(j+k)) · sum coefs[n]·c^n by
        // iteratively computing c^n in Fourier form and accumulating.
        // fourier[m] is the rational coefficient of cos(m·2·arg).
        Rational denom_scale(BigInt(1));
        for (long long ii = 0; ii < j + k; ++ii) denom_scale = denom_scale * Rational(BigInt(1), BigInt(2));
        std::vector<Rational> total_fourier;  // index = m
        std::vector<Rational> cosn{Rational(BigInt(1))};  // c^0 = cos(0·y) = 1
        const std::size_t deg = coefs.size() - 1U;
        for (std::size_t n = 0; n <= deg; ++n) {
            if (!coefs[n].is_zero()) {
                const Rational w = Rational(coefs[n]) * denom_scale;
                for (std::size_t m = 0; m < cosn.size(); ++m) {
                    if (m >= total_fourier.size()) total_fourier.resize(m + 1U, Rational(BigInt(0)));
                    total_fourier[m] = total_fourier[m] + w * cosn[m];
                }
            }
            if (n < deg) {
                std::vector<Rational> next(cosn.size() + 1U, Rational(BigInt(0)));
                // m=0: cos(y)·1 = cos(y) → next[1] += cosn[0].
                if (!cosn.empty()) next[1] = next[1] + cosn[0];
                // m≥1: cos(y)·cos(my) = (cos((m+1)y)+cos((m-1)y))/2.
                for (std::size_t m = 1; m < cosn.size(); ++m) {
                    Rational half = cosn[m] * Rational(BigInt(1), BigInt(2));
                    next[m + 1U] = next[m + 1U] + half;
                    next[m - 1U] = next[m - 1U] + half;
                }
                cosn = std::move(next);
            }
        }
        // Emit: constant term → coeff · var ; cos(m·2·arg) term → coeff·sin(m·2·arg)/(2m·c_coeff).
        std::vector<ExprPtr> terms;
        if (!total_fourier.empty() && !total_fourier[0].numerator().is_zero()) {
            ExprPtr var_expr = arena.make<Symbol>(var.name);
            terms.push_back(arena.make<Product>(std::vector<ExprPtr>{
                make_rational(arena, total_fourier[0]), var_expr}));
        }
        for (std::size_t m = 1; m < total_fourier.size(); ++m) {
            if (total_fourier[m].numerator().is_zero()) continue;
            ExprPtr m2 = make_integer(arena, 2LL * static_cast<long long>(m));
            ExprPtr arg_scaled = arena.make<Product>(std::vector<ExprPtr>{m2, trig_arg});
            ExprPtr sin_term = arena.make<FuncCall>(BuiltinOp::Sin, std::vector<ExprPtr>{arg_scaled});
            const Rational two_m_c = Rational(BigInt(2LL * static_cast<long long>(m))) * c_coeff;
            const Rational coef = total_fourier[m] / two_m_c;
            terms.push_back(arena.make<Product>(std::vector<ExprPtr>{
                make_rational(arena, coef), sin_term}));
        }
        ExprPtr primitive = terms.empty() ? make_integer(arena, 0)
            : (terms.size() == 1U ? terms[0] : arena.make<Sum>(std::move(terms)));
        if (!consts.empty()) {
            consts.push_back(primitive);
            primitive = arena.make<Product>(std::move(consts));
        }
        if (auto simp = ctx.simplify(primitive); simp.is_ok()) primitive = simp.value();
        return ok(primitive);
    }
    // Case 3 (+ degenerate mixed odd-negative): t = tan(arg).
    //   sin² = t²/(1+t²), cos² = 1/(1+t²), dx = dt / (c · (1+t²)).
    //   integrand · dx/dt = t^a · (1+t²)^(-(a+b)/2 - 1) / c.
    // Restrict to inputs with at least one NEGATIVE exponent (csc / sec mix):
    // for both-positive-even inputs, the Weierstrass + simplify_trig path
    // already produces a canonical answer (sin²cos², sin²cos⁴, sin⁴cos⁴) and
    // the t=tan substitution generates intermediate forms the cert cannot
    // reduce. So we bail here and let the downstream pipeline take over.
    else if ((sin_pow < 0 || cos_pow < 0)
             && ((sin_pow + cos_pow) % 2 == 0)) {
        const long long pow_exp = -((sin_pow + cos_pow) / 2 + 1);
        ExprPtr u_a = u_power(u_sym, sin_pow, arena);
        ExprPtr one_plus = one_plus_u_sq(u_sym, arena);
        ExprPtr one_plus_pow = (pow_exp == 0)
            ? make_integer(arena, 1)
            : arena.make<Binary>(BinaryOp::Pow, one_plus, make_integer(arena, pow_exp));
        ExprPtr prod = arena.make<Product>(std::vector<ExprPtr>{u_a, one_plus_pow});
        integrand_in_u = scale_by_inverse(prod, c_coeff, arena);
    }
    else {
        return ::cas::make_unimplemented<ExprPtr>(
            "calculus", "try_integrate_trig_monomial",
            "both-positive-even trig powers handled by downstream Weierstrass + simplify_trig",
            cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
            "Use downstream pipeline", "W3.1");
    }

    // Simplify the integrand in u BEFORE handing off — collapses
    // Product[Pow(u,-2), (1+u²)] to 1/u²+1 so the rational integrator sees a
    // clean polynomial-or-rational shape, not a residual nested Pow/Product.
    if (auto simp_in = ctx.simplify(integrand_in_u); simp_in.is_ok())
        integrand_in_u = simp_in.value();
    // Rewriting via together() folds nested rational factors (relevant when
    // u-substitution produces t^a · (1+t²)^(-k) with k>1 — the integrator
    // expects a clean Div(N, D) for partial-fractions dispatch).
    if (auto tg = algebra::together(integrand_in_u, ctx); tg.is_ok())
        integrand_in_u = tg.value();

    // Integrate the polynomial- (or Laurent- or rational-) integrand in u.
    auto inner = integrate_indefinite_impl(integrand_in_u, u, ctx);
    if (inner.is_error()) return fail<ExprPtr>(inner.error());

    // Back-substitute: u = cos(arg) (case 1), u = sin(arg) (case 2),
    // u = tan(arg) (cases 3/4).
    BuiltinOp back_op;
    if (sin_pow > 0 && (sin_pow % 2) == 1)      back_op = BuiltinOp::Cos;
    else if (cos_pow > 0 && (cos_pow % 2) == 1) back_op = BuiltinOp::Sin;
    else                                         back_op = BuiltinOp::Tan;
    ExprPtr back_sub = arena.make<FuncCall>(back_op, std::vector<ExprPtr>{trig_arg});
    auto subst = symbolic::substitute(inner.value(), u, back_sub, ctx);
    if (subst.is_error()) return fail<ExprPtr>(subst.error());

    ExprPtr primitive = subst.value();
    if (!consts.empty()) {
        consts.push_back(primitive);
        primitive = arena.make<Product>(std::move(consts));
    }
    auto simp = ctx.simplify(primitive);
    if (simp.is_ok()) primitive = simp.value();
    return ok(primitive);
}

// Wrapper for Binary(Pow, Product[trig_factors], integer_exp). Distributes
// the exponent over the inner Product and dispatches to the main entry.
// Catches integrands like 1/(sin²x·cos²x) which parse as
// Pow(Product[sin²,cos²], -1) at the AST root.
[[nodiscard]] Result<ExprPtr> try_integrate_trig_monomial_pow(
    const Binary& power, const Symbol& var, AstArena& arena,
    symbolic::CASContext& ctx) {
    const auto* base_prod = expr_cast<Product>(power.left);
    if (!base_prod) return ::cas::make_unimplemented<ExprPtr>(
        "calculus", "try_integrate_trig_monomial_pow",
        "Pow base is not a Product", cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
        "", "W3.1");
    const auto* il = expr_cast<IntegerLit>(power.right);
    if (!il || il->value.is_zero()) return ::cas::make_unimplemented<ExprPtr>(
        "calculus", "try_integrate_trig_monomial_pow",
        "Pow exponent is not a non-zero integer literal",
        cas::error::reason_codes::INTEGRATE_NO_STRATEGY, "", "W3.1");
    std::vector<ExprPtr> new_factors;
    new_factors.reserve(base_prod->factors.size());
    for (ExprPtr f : base_prod->factors) {
        new_factors.push_back(arena.make<Binary>(BinaryOp::Pow, f,
            arena.make<IntegerLit>(il->value)));
    }
    Product distributed{std::move(new_factors)};
    return try_integrate_trig_monomial(distributed, var, arena, ctx);
}

}  // namespace cas::calculus::integrate_detail
