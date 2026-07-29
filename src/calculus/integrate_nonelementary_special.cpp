// A43 §5 — Famiglie non esponenziali del fallback non elementare:
//   gaussiana  ∫ e^{A x² + B x + C} dx        → erf / erfi
//   trig/lin   ∫ trig(a x + b)/(x − r) dx     → Si, Ci, Shi, Chi
//   log-int    ∫ x^s/ln x dx, ∫ c/ln(αx+β) dx → li
//   dilog      ∫ ln(α x + β)/(x − r) dx       → Li₂
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Nonelementary_Antiderivatives.md
// Dispatcher e famiglia Ei: `integrate_nonelementary.cpp`.
//
// Ogni formula qui sotto è scritta in modo che la sua derivata sia verificabile
// termine a termine — i test di `test_a43_nonelementary_integrate.cpp` fanno
// esattamente quello (D(F) = f), che è l'unico criterio di accettazione dopo
// A42: un'antiderivata non derivabile non è verificabile né componibile.

#include "integrate_nonelementary.hpp"

#include "cas/error_helpers.hpp"

#include <utility>
#include <vector>

namespace cas::calculus::nonelementary {

namespace {

using integrate_detail::make_binary;
using integrate_detail::make_function;
using integrate_detail::make_product;
using integrate_detail::make_sum;
using integrate_detail::make_unary;

// Vedi la nota gemella in `integrate_nonelementary.cpp`: denominatore 1 →
// `IntegerLit`, altrimenti la costante sopravvive come `RationalLit(n,1)`
// dentro `exp`/`ln` e l'antiderivata non è confrontabile strutturalmente.
[[nodiscard]] ExprPtr make_rational(AstArena& arena, const Rational& value) {
    if (value.is_integer()) return arena.make<IntegerLit>(value.numerator());
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] CASError not_this_family(const char* what) {
    return integrate_detail::make_error(CASErrorKind::Unimplemented, what);
}

[[nodiscard]] bool is_unit_exponent(const Factor& factor) {
    return factor.exponent == Rational(BigInt(1));
}

// Il singolo fattore `(x − r)^{-1}` di una vista, con il suo fattore di scala.
// Restituisce `false` se la vista non contiene esattamente un polo semplice
// lineare più un altro fattore.
[[nodiscard]] bool split_simple_pole(
    const ProductView& view, const Symbol& var,
    LinearPole& pole, ExprPtr& other) {
    bool pole_found = false;
    bool other_found = false;
    for (const Factor& factor : view.factors) {
        if (auto candidate = as_linear_pole(factor, var, /*max_order=*/1U);
            candidate.has_value() && candidate->order == 1) {
            if (pole_found) return false;
            pole = candidate.value();
            pole_found = true;
            continue;
        }
        if (!is_unit_exponent(factor) || other_found) return false;
        other = factor.base;
        other_found = true;
    }
    return pole_found && other_found;
}

// √v per v razionale positivo — lasciato simbolico, la costant-folding del
// simplifier riduce i quadrati perfetti.
[[nodiscard]] ExprPtr make_sqrt(AstArena& arena, const Rational& value) {
    return make_function(arena, "sqrt", {make_rational(arena, value)});
}

}  // namespace

// ∫ e^{A x² + B x + C} dx, A ≠ 0. Completamento del quadrato:
//   A x² + B x + C = A·(x + B/2A)² + (C − B²/4A)
// e con t = x + B/2A:
//   A > 0 : ∫e^{A t²}dt = ½·√(π/A)·erfi(√A·t)
//   A < 0 : ∫e^{−a t²}dt = ½·√(π/a)·erf(√a·t),  a = −A
// Le due derivate: d/dt erfi(u) = (2/√π)e^{u²}·u′ e d/dt erf(u) = (2/√π)e^{−u²}·u′
// riproducono esattamente l'integranda (spec §3, verificate mpmath).
Result<ExprPtr> try_gaussian(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (view.factors.size() != 1U || !is_unit_exponent(view.factors.front())) {
        return fail<ExprPtr>(not_this_family("gaussian: integrand is not a single exponential"));
    }
    ExprPtr argument{};
    if (!is_exponential(view.factors.front().base, argument)) {
        return fail<ExprPtr>(not_this_family("gaussian: factor is not an exponential"));
    }
    auto quadratic = integrate_detail::extract_quadratic_argument(argument, var);
    if (!quadratic.has_value() || quadratic->quadratic.numerator().is_zero()) {
        return fail<ExprPtr>(not_this_family("gaussian: exponent is not a genuine quadratic"));
    }

    const Rational A = quadratic->quadratic;
    const Rational B = quadratic->linear;
    const Rational C = quadratic->constant;
    const Rational two(BigInt(2));
    const Rational four(BigInt(4));
    const Rational shift = B / (two * A);              // x + shift = t
    const Rational constant_term = C - (B * B) / (four * A);
    const bool positive_leading = !A.numerator().is_negative();
    const Rational magnitude = positive_leading ? A : -A;

    ExprPtr t = make_sum(arena, {arena.make<Symbol>(var), make_rational(arena, shift)});
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    // ½·√(π/|A|)
    ExprPtr prefactor = make_product(arena, {
        make_rational(arena, Rational(BigInt(1), BigInt(2))),
        make_function(arena, "sqrt",
            {make_binary(arena, BinaryOp::Div, pi, make_rational(arena, magnitude))}),
    });
    ExprPtr inner = make_product(arena, {make_sqrt(arena, magnitude), t});

    return ctx.simplify(make_product(arena, {
        make_rational(arena, view.coefficient),
        make_function(arena, "exp", {make_rational(arena, constant_term)}),
        prefactor,
        make_function(arena, positive_leading ? "erfi" : "erf", {inner}),
    }));
}

// ∫ trig(a x + b)/(x − r) dx. Con t = x − r l'argomento diventa a t + φ,
// φ = a r + b, e le formule di addizione danno una combinazione delle quattro
// primitive non elementari (spec §5; le identità di addizione sono esatte, non
// una tabella per ogni b):
//   sin :  cos φ · Si(a t) + sin φ · Ci(a t)
//   cos :  cos φ · Ci(a t) − sin φ · Si(a t)
//   sinh: cosh φ · Shi(a t) + sinh φ · Chi(a t)
//   cosh: cosh φ · Chi(a t) + sinh φ · Shi(a t)
Result<ExprPtr> try_trig_over_linear(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    LinearPole pole{Rational(), Rational(), 0};
    ExprPtr other{};
    if (!split_simple_pole(view, var, pole, other)) {
        return fail<ExprPtr>(not_this_family("trig/linear: no simple pole plus one factor"));
    }
    const auto* call = expr_cast<FuncCall>(other);
    if (call == nullptr || call->args.size() != 1U) {
        return fail<ExprPtr>(not_this_family("trig/linear: numerator is not a function call"));
    }

    const char* aligned = nullptr;   // primitiva del termine in fase con l'argomento
    const char* shifted = nullptr;   // primitiva del termine in quadratura
    bool hyperbolic = false;
    bool negate_shifted = false;
    switch (call->func_id) {
        case BuiltinOp::Sin:  aligned = "Si";  shifted = "Ci";  break;
        case BuiltinOp::Cos:  aligned = "Ci";  shifted = "Si";  negate_shifted = true; break;
        case BuiltinOp::Sinh: aligned = "Shi"; shifted = "Chi"; hyperbolic = true; break;
        case BuiltinOp::Cosh: aligned = "Chi"; shifted = "Shi"; hyperbolic = true; break;
        default:
            return fail<ExprPtr>(not_this_family("trig/linear: not sin/cos/sinh/cosh"));
    }

    auto affine = integrate_detail::extract_affine_argument(call->args.front(), var);
    if (!affine.has_value() || affine->coefficient.numerator().is_zero()) {
        return fail<ExprPtr>(not_this_family("trig/linear: argument is not affine in the variable"));
    }

    const Rational a = affine->coefficient;
    const Rational phase = a * pole.root + affine->constant;
    ExprPtr t = make_sum(arena, {arena.make<Symbol>(var), make_rational(arena, -pole.root)});
    ExprPtr scaled = make_product(arena, {make_rational(arena, a), t});
    ExprPtr prefactor = make_rational(arena, view.coefficient * pole.scale);

    if (phase.numerator().is_zero()) {
        // Caso in fase: nessun termine in quadratura, la forma resta quella
        // canonica della spec §5 (∫sin x/x = Si x, ∫cosh x/x = Chi x, …).
        return ctx.simplify(make_product(arena, {
            prefactor, make_function(arena, aligned, {scaled})}));
    }

    ExprPtr phase_expr = make_rational(arena, phase);
    ExprPtr aligned_weight = make_function(arena, hyperbolic ? "cosh" : "cos", {phase_expr});
    ExprPtr shifted_weight = make_function(arena, hyperbolic ? "sinh" : "sin", {phase_expr});
    ExprPtr shifted_term = make_product(arena, {shifted_weight, make_function(arena, shifted, {scaled})});
    if (negate_shifted) {
        shifted_term = make_unary(arena, UnaryOp::Neg, shifted_term);
    }
    return ctx.simplify(make_product(arena, {
        prefactor,
        make_sum(arena, {
            make_product(arena, {aligned_weight, make_function(arena, aligned, {scaled})}),
            shifted_term,
        }),
    }));
}

// Integrale logaritmico. Due letture della stessa identità ∫ g′/ln g = li(g):
//   (L1) ∫ c/ln(α x + β) dx = (c/α)·li(α x + β)      — g = α x + β, g′ = α
//   (L2) ∫ x^s/ln x dx      = li(x^{s+1}),  s ≠ −1   — g = x^{s+1}, poiché
//        ln(x^{s+1}) = (s+1)·ln x e (x^{s+1})′ = (s+1)x^s
// (L2) contiene (L1) per α = 1, β = 0; entrambe sono verificate derivando.
Result<ExprPtr> try_log_integral(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr log_argument{};
    ExprPtr numerator{};
    for (const Factor& factor : view.factors) {
        const auto* call = expr_cast<FuncCall>(factor.base);
        const bool is_log = call != nullptr && call->args.size() == 1U
            && (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log);
        if (is_log && factor.exponent == -Rational(BigInt(1))) {
            if (log_argument) return fail<ExprPtr>(not_this_family("log-integral: two logarithms"));
            log_argument = call->args.front();
            continue;
        }
        if (numerator) return fail<ExprPtr>(not_this_family("log-integral: too many factors"));
        numerator = make_binary(arena, BinaryOp::Pow, factor.base,
                                make_rational(arena, factor.exponent));
    }
    if (!log_argument) {
        return fail<ExprPtr>(not_this_family("log-integral: no logarithm in the denominator"));
    }

    ExprPtr coefficient = make_rational(arena, view.coefficient);
    if (!numerator) {
        // (L1): il numeratore è la sola costante.
        auto affine = integrate_detail::extract_affine_argument(log_argument, var);
        if (!affine.has_value() || affine->coefficient.numerator().is_zero()) {
            return fail<ExprPtr>(not_this_family("log-integral: logarithm argument is not affine"));
        }
        return ctx.simplify(make_product(arena, {
            make_rational(arena, view.coefficient / affine->coefficient),
            make_function(arena, "li", {log_argument}),
        }));
    }

    // (L2): il logaritmo deve essere ln(x) e il numeratore una potenza di x.
    if (!integrate_detail::is_same_symbol(log_argument, var)) {
        return fail<ExprPtr>(not_this_family("log-integral: power case needs ln(x) exactly"));
    }
    Rational power(BigInt(0));
    bool power_found = false;
    for (const Factor& factor : view.factors) {
        if (integrate_detail::is_same_symbol(factor.base, var)) {
            power = factor.exponent;
            power_found = true;
        }
    }
    const Rational minus_one = -Rational(BigInt(1));
    if (!power_found || power == minus_one) {
        return fail<ExprPtr>(not_this_family("log-integral: numerator is not x^s with s != -1"));
    }
    ExprPtr shifted_power = make_binary(arena, BinaryOp::Pow, arena.make<Symbol>(var),
                                        make_rational(arena, power + Rational(BigInt(1))));
    return ctx.simplify(make_product(arena, {
        coefficient, make_function(arena, "li", {shifted_power})}));
}

// ∫ ln(α x + β)/(x − r) dx. Con t = x − r e c = α r + β (≠ 0):
//   ln(α t + c) = ln c + ln(1 + (α/c) t)
//   ∫ ln(1 + u)/t dt con u = (α/c)t  →  −Li₂(−(α/c)·t)   (spec §2, §5)
// quindi il risultato è ln(c)·ln(t) − Li₂(−(α/c)·t).
// Per c = 0 l'integranda è elementare (ln(αt)/t) e non appartiene a questa
// famiglia: il match non scatta e il controllo prosegue.
Result<ExprPtr> try_dilogarithm(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    LinearPole pole{Rational(), Rational(), 0};
    ExprPtr other{};
    if (!split_simple_pole(view, var, pole, other)) {
        return fail<ExprPtr>(not_this_family("dilog: no simple pole plus one factor"));
    }
    const auto* call = expr_cast<FuncCall>(other);
    if (call == nullptr || call->args.size() != 1U
        || (call->func_id != BuiltinOp::Ln && call->func_id != BuiltinOp::Log)) {
        return fail<ExprPtr>(not_this_family("dilog: numerator is not a natural logarithm"));
    }
    auto affine = integrate_detail::extract_affine_argument(call->args.front(), var);
    if (!affine.has_value() || affine->coefficient.numerator().is_zero()) {
        return fail<ExprPtr>(not_this_family("dilog: logarithm argument is not affine"));
    }

    const Rational alpha = affine->coefficient;
    const Rational c = alpha * pole.root + affine->constant;
    if (c.numerator().is_zero()) {
        return fail<ExprPtr>(not_this_family("dilog: shifted logarithm is homogeneous (elementary)"));
    }

    ExprPtr t = make_sum(arena, {arena.make<Symbol>(var), make_rational(arena, -pole.root)});
    ExprPtr dilog_argument = make_product(arena, {make_rational(arena, -alpha / c), t});
    std::vector<ExprPtr> terms;
    if (c != Rational(BigInt(1))) {
        terms.push_back(make_product(arena, {
            make_function(arena, "ln", {make_rational(arena, c)}),
            make_function(arena, "ln", {t}),
        }));
    }
    terms.push_back(make_unary(arena, UnaryOp::Neg,
        make_function(arena, "dilog", {dilog_argument})));

    return ctx.simplify(make_product(arena, {
        make_rational(arena, view.coefficient * pole.scale),
        make_sum(arena, std::move(terms)),
    }));
}

}  // namespace cas::calculus::nonelementary
