// A43 §5 — Fallback per le antiderivate NON elementari: dispatcher, vista
// moltiplicativa piatta, e la famiglia dell'integrale esponenziale (Ei).
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Nonelementary_Antiderivatives.md
// (primitive verificate con mpmath.quad su intervalli definiti, §7).
//
// Vincolo di collocazione (spec §5, ultimo capoverso): questo fallback va
// tentato DOPO Risch. Risch è la procedura di DECISIONE — se risponde "non
// elementare", allora e solo allora si cerca una forma in questa famiglia
// estesa. Anticiparlo maschererebbe i casi elementari. È montato in
// `integrate_core.cpp` fra Weierstrass e il fallback Meijer G: prima di
// Meijer perché su questa classe Ei/Si/erfi è la forma canonica che gli
// oracoli emettono, mentre la via Meijer restituisce un ₁F₁ equivalente ma
// meno leggibile (misurato: ∫e^{x²} usciva come `x·₁F₁(1/2,3/2,x²)`).
//
// Confine dichiarato (spec §6): denominatori che non si spezzano in fattori
// lineari su Q richiederebbero `Ei` ad argomento algebrico/complesso, che la
// spec lascia fuori scope — quelle forme semplicemente non fanno match e il
// controllo prosegue verso Meijer G e poi verso l'errore NO_STRATEGY già
// esistente. Nessun risultato silenziosamente sbagliato, nessun bail-out nuovo.

#include "integrate_nonelementary.hpp"

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"

#include <utility>
#include <vector>

namespace cas::calculus {

namespace nonelementary {

namespace {

using integrate_detail::make_binary;
using integrate_detail::make_function;
using integrate_detail::make_integer;
using integrate_detail::make_product;
using integrate_detail::make_sum;
using integrate_detail::make_unary;

// Un razionale a denominatore 1 deve uscire come `IntegerLit`, non come
// `RationalLit(n,1)`: la seconda forma sopravvive dentro `exp(...)` e rende
// l'antiderivata non confrontabile strutturalmente (misurato: l'antiderivata
// di e^x/(x+1) usciva come `exp(-1/1)·Ei(x+1)`). Il `make_rational` condiviso
// dell'integratore non normalizza, e cambiarlo lì avrebbe raggio d'azione su
// tutto il modulo — la normalizzazione resta locale ad A43.
[[nodiscard]] ExprPtr make_rational(AstArena& arena, const Rational& value) {
    if (value.is_integer()) return arena.make<IntegerLit>(value.numerator());
    return arena.make<RationalLit>(value.numerator(), value.denominator());
}

[[nodiscard]] std::optional<Rational> numeric_value(ExprPtr expr) {
    return integrate_detail::exact_scalar_from_expr(expr);
}

// Ricorsione di `flatten_product`: `sign` vale +1 al numeratore, −1 sotto una
// divisione, e si propaga agli esponenti.
void collect_factors(ExprPtr expr, const Rational& sign, ProductView& out) {
    if (!expr) return;

    if (auto value = numeric_value(expr); value.has_value()) {
        // Il fold nel coefficiente è esatto solo per esponente ±1: una potenza
        // razionale di un numero (2^{1/2}) non è un razionale, e 1/0 non è
        // rappresentabile. Negli altri casi il fattore resta opaco, così il
        // riconoscimento fallisce invece di produrre un valore falso.
        const Rational one(BigInt(1));
        if (sign == one) {
            out.coefficient *= value.value();
            return;
        }
        if (sign == -one && !value->numerator().is_zero()) {
            out.coefficient /= value.value();
            return;
        }
        out.factors.push_back(Factor{expr, sign});
        return;
    }

    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) {
        out.coefficient *= Rational(BigInt(-1));
        collect_factors(unary->operand, sign, out);
        return;
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        for (ExprPtr factor : product->factors) {
            collect_factors(factor, sign, out);
        }
        return;
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Mul) {
            collect_factors(binary->left, sign, out);
            collect_factors(binary->right, sign, out);
            return;
        }
        if (binary->op == BinaryOp::Div) {
            collect_factors(binary->left, sign, out);
            collect_factors(binary->right, -sign, out);
            return;
        }
        if (binary->op == BinaryOp::Pow) {
            if (auto exponent = numeric_value(binary->right); exponent.has_value()) {
                collect_factors(binary->left, sign * exponent.value(), out);
                return;
            }
        }
    }

    out.factors.push_back(Factor{expr, sign});
}

}  // namespace

ProductView flatten_product(ExprPtr expr) {
    ProductView view;
    collect_factors(expr, Rational(BigInt(1)), view);
    return view;
}

bool is_exponential(ExprPtr expr, ExprPtr& argument) {
    if (const auto* call = expr_cast<FuncCall>(expr);
        call != nullptr && call->func_id == BuiltinOp::Exp && call->args.size() == 1U) {
        argument = call->args.front();
        return true;
    }
    if (const auto* binary = expr_cast<Binary>(expr); binary != nullptr && binary->op == BinaryOp::Pow) {
        if (const auto* base = expr_cast<Constant>(binary->left);
            base != nullptr && base->value == MathConstant::E) {
            argument = binary->right;
            return true;
        }
    }
    return false;
}

std::optional<long long> bounded_degree(const BigInt& value, std::size_t max_order) {
    if (value.is_negative()) return std::nullopt;
    if (value > BigInt(static_cast<long long>(max_order))) return std::nullopt;
    return static_cast<long long>(value.to_u64());
}

std::optional<LinearPole> as_linear_pole(
    const Factor& factor, const Symbol& var, std::size_t max_order) {
    if (!factor.exponent.is_integer() || !factor.exponent.numerator().is_negative()) {
        return std::nullopt;
    }
    // Un polo di ordine k costa k−1 passi di riduzione per parti: oltre il
    // budget di ricorsione dell'integratore il match semplicemente non scatta.
    const auto order = bounded_degree(-factor.exponent.numerator(), max_order);
    if (!order.has_value() || order.value() < 1) return std::nullopt;

    auto affine = integrate_detail::extract_affine_argument(factor.base, var);
    if (!affine.has_value() || affine->coefficient.numerator().is_zero()) {
        return std::nullopt;
    }
    // (α x + β)^{-k} = α^{-k} (x − r)^{-k}, r = −β/α.
    Rational scale(BigInt(1));
    for (long long i = 0; i < order.value(); ++i) {
        scale /= affine->coefficient;
    }
    return LinearPole{
        .root = -affine->constant / affine->coefficient,
        .scale = scale,
        .order = order.value(),
    };
}

Result<ExprPtr> exp_over_linear_power(
    const Rational& a, const Rational& b, const Rational& r, long long k,
    const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (auto irq = ctx.check_interrupt(); irq.is_error()) return fail<ExprPtr>(irq.error());
    if (a.numerator().is_zero() || k < 1) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::InvalidArgument, "exp_over_linear_power requires a != 0 and k >= 1"));
    }

    ExprPtr x = arena.make<Symbol>(var);
    // t = x − r
    ExprPtr shifted = make_sum(arena, {x, make_rational(arena, -r)});
    // e^{a·x+b}
    ExprPtr exponential = make_function(arena, "exp",
        {make_sum(arena, {make_product(arena, {make_rational(arena, a), x}), make_rational(arena, b)})});

    if (k == 1) {
        // ∫ e^{a x + b}/(x − r) dx = e^{a r + b}·Ei(a·(x − r))
        // (sostituzione t = a(x−r): l'integranda diventa e^{ar+b}·e^t/t).
        return ok(make_product(arena, {
            make_function(arena, "exp", {make_rational(arena, a * r + b)}),
            make_function(arena, "Ei", {make_product(arena, {make_rational(arena, a), shifted})}),
        }));
    }

    // Riduzione per parti (u = e^{a x + b}, dv = (x−r)^{-k} dx):
    //   ∫ e^{a x+b}(x−r)^{-k} dx
    //     = −e^{a x+b}/((k−1)(x−r)^{k−1}) + a/(k−1)·∫ e^{a x+b}(x−r)^{-(k−1)} dx
    // Ricorsione fino a k = 1: è la riduzione richiesta dalla spec §5, non una
    // riga di tabella per ogni k.
    auto reduced = exp_over_linear_power(a, b, r, k - 1, var, ctx);
    if (reduced.is_error()) return reduced;

    const Rational one_over_k_minus_one(BigInt(1), BigInt(k - 1));
    ExprPtr boundary = make_unary(arena, UnaryOp::Neg, make_product(arena, {
        make_rational(arena, one_over_k_minus_one),
        exponential,
        make_binary(arena, BinaryOp::Pow, shifted, make_integer(arena, -(k - 1))),
    }));
    return ok(make_sum(arena, {
        boundary,
        make_product(arena, {make_rational(arena, a * one_over_k_minus_one), reduced.value()}),
    }));
}

Result<ExprPtr> exp_times_monomial(
    const Rational& a, const Rational& b, long long n,
    const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    if (auto irq = ctx.check_interrupt(); irq.is_error()) return fail<ExprPtr>(irq.error());
    if (a.numerator().is_zero() || n < 0) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::InvalidArgument, "exp_times_monomial requires a != 0 and n >= 0"));
    }

    ExprPtr x = arena.make<Symbol>(var);
    ExprPtr exponential = make_function(arena, "exp",
        {make_sum(arena, {make_product(arena, {make_rational(arena, a), x}), make_rational(arena, b)})});
    const Rational inverse_a = Rational(BigInt(1)) / a;

    if (n == 0) {
        return ok(make_product(arena, {make_rational(arena, inverse_a), exponential}));
    }

    // ∫ e^{a x+b} x^n dx = x^n e^{a x+b}/a − (n/a)·∫ e^{a x+b} x^{n−1} dx.
    // Parte elementare: la copre la stessa ricorsione per parti, così la
    // decomposizione dell'integranda razionale non deve rientrare nel
    // dispatcher (che l'ha già rifiutata prima di arrivare qui).
    auto reduced = exp_times_monomial(a, b, n - 1, var, ctx);
    if (reduced.is_error()) return reduced;
    return ok(make_sum(arena, {
        make_product(arena, {
            make_rational(arena, inverse_a),
            make_binary(arena, BinaryOp::Pow, x, make_integer(arena, n)),
            exponential,
        }),
        make_unary(arena, UnaryOp::Neg, make_product(arena, {
            make_rational(arena, Rational(BigInt(n)) * inverse_a),
            reduced.value(),
        })),
    }));
}

namespace {

// Estrae dalla vista l'unico fattore esponenziale con argomento affine.
// Restituisce false se gli esponenziali sono zero o più d'uno (in quel caso
// l'integranda non appartiene alla classe ∫e^{L}·R).
[[nodiscard]] bool split_exponential(
    const ProductView& view, const Symbol& var,
    integrate_detail::AffineArgument& exponent, std::vector<Factor>& rest) {
    bool found = false;
    for (const Factor& factor : view.factors) {
        ExprPtr argument{};
        if (factor.exponent == Rational(BigInt(1)) && is_exponential(factor.base, argument)) {
            auto affine = integrate_detail::extract_affine_argument(argument, var);
            if (affine.has_value() && !affine->coefficient.numerator().is_zero()) {
                if (found) return false;
                exponent = affine.value();
                found = true;
                continue;
            }
        }
        rest.push_back(factor);
    }
    return found;
}

// ∫ e^{a x + b}·R(x) dx con R razionale: decomposizione in fratti semplici,
// poi riduzione di ogni pezzo (spec §5 — "implementare la riduzione, non le
// righe"). I poli semplici producono Ei, quelli multipli vi si riducono per
// parti, la parte polinomiale resta elementare.
[[nodiscard]] Result<ExprPtr> try_exponential_rational(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    integrate_detail::AffineArgument exponent{Rational(), Rational()};
    std::vector<Factor> rest;
    if (!split_exponential(view, var, exponent, rest)) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::Unimplemented, "no single exponential factor with affine argument"));
    }

    // Ricostruisci la parte razionale e decomponila.
    std::vector<ExprPtr> rational_factors;
    rational_factors.push_back(make_rational(arena, view.coefficient));
    for (const Factor& factor : rest) {
        rational_factors.push_back(
            make_binary(arena, BinaryOp::Pow, factor.base, make_rational(arena, factor.exponent)));
    }
    auto rational_part = ctx.simplify(make_product(arena, std::move(rational_factors)));
    if (rational_part.is_error()) return fail<ExprPtr>(rational_part.error());

    // La decomposizione in fratti semplici serve solo se la parte razionale NON
    // e' gia' un singolo termine della decomposizione. Chiamarla comunque su
    // `(x−r)^{-k}` la fa passare per espansione e rifattorizzazione del
    // denominatore — costo misurato: `∫e^x/(x−1)²` 4.0 s contro 0.58 s di
    // `∫e^x/x²`, e `∫e^{2x}/(x−1)³` oltre il minuto. Riconoscere il termine
    // gia' decomposto non e' una scorciatoia: e' evitare di decomporre due volte.
    std::vector<ExprPtr> terms;
    {
        const ProductView direct = flatten_product(rational_part.value());
        const bool already_decomposed =
            direct.factors.empty()
            || (direct.factors.size() == 1U
                && (as_linear_pole(direct.factors.front(), var, ctx.max_integration_depth()).has_value()
                    || (integrate_detail::is_same_symbol(direct.factors.front().base, var)
                        && direct.factors.front().exponent.is_integer())));
        if (already_decomposed) {
            terms.push_back(rational_part.value());
        } else if (auto decomposed = algebra::partial_fractions(rational_part.value(), var, ctx);
                   decomposed.is_ok() && !decomposed.value().empty()) {
            terms = decomposed.value();
        } else {
            terms.push_back(rational_part.value());
        }
    }

    std::vector<ExprPtr> pieces;

    for (ExprPtr term : terms) {
        auto simplified = ctx.simplify(term);
        if (simplified.is_error()) return fail<ExprPtr>(simplified.error());
        const ProductView piece = flatten_product(simplified.value());

        // Un pezzo ammissibile è `c · (x−r)^{-k}` oppure `c · x^n`; qualunque
        // altra forma (denominatore irriducibile di grado ≥ 2, radici non
        // razionali) è fuori dal confine dichiarato dalla spec §6.
        if (piece.factors.empty()) {
            auto monomial = exp_times_monomial(exponent.coefficient, exponent.constant, 0, var, ctx);
            if (monomial.is_error()) return monomial;
            pieces.push_back(make_product(arena, {make_rational(arena, piece.coefficient), monomial.value()}));
            continue;
        }
        if (piece.factors.size() != 1U) {
            return fail<ExprPtr>(integrate_detail::make_error(
                CASErrorKind::Unimplemented, "partial-fraction term is not c·(x−r)^k"));
        }

        const Factor& factor = piece.factors.front();
        if (auto pole = as_linear_pole(factor, var, ctx.max_integration_depth()); pole.has_value()) {
            auto reduced = exp_over_linear_power(
                exponent.coefficient, exponent.constant, pole->root, pole->order, var, ctx);
            if (reduced.is_error()) return reduced;
            pieces.push_back(make_product(arena, {
                make_rational(arena, piece.coefficient * pole->scale), reduced.value()}));
            continue;
        }
        if (integrate_detail::is_same_symbol(factor.base, var) && factor.exponent.is_integer()) {
            const auto degree = bounded_degree(
                factor.exponent.numerator(), ctx.max_integration_depth());
            if (!degree.has_value()) {
                return fail<ExprPtr>(integrate_detail::make_error(
                    CASErrorKind::Unimplemented, "monomial degree outside the reduction budget"));
            }
            auto monomial = exp_times_monomial(
                exponent.coefficient, exponent.constant, degree.value(), var, ctx);
            if (monomial.is_error()) return monomial;
            pieces.push_back(make_product(arena, {
                make_rational(arena, piece.coefficient), monomial.value()}));
            continue;
        }
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::Unimplemented, "partial-fraction term outside the Ei reduction class"));
    }

    if (pieces.empty()) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::Unimplemented, "empty decomposition"));
    }
    return ctx.simplify(make_sum(arena, std::move(pieces)));
}

// ∫ e^{a x + b}·ln(α x + β) dx — per parti UNA volta, poi la coda cade nella
// famiglia Ei (spec §5: "le ultime tre si ottengono per parti dalle prime"):
//   = e^{a x+b}·ln(α x+β)/a − (1/a)·e^{a r + b}·Ei(a·(x − r)),  r = −β/α.
[[nodiscard]] Result<ExprPtr> try_exponential_log(
    const ProductView& view, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    integrate_detail::AffineArgument exponent{Rational(), Rational()};
    std::vector<Factor> rest;
    if (!split_exponential(view, var, exponent, rest) || rest.size() != 1U) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::Unimplemented, "not exp(affine)·log(affine)"));
    }
    const Factor& log_factor = rest.front();
    if (log_factor.exponent != Rational(BigInt(1))) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::Unimplemented, "logarithm is not to the first power"));
    }
    const auto* call = expr_cast<FuncCall>(log_factor.base);
    if (call == nullptr || call->args.size() != 1U
        || (call->func_id != BuiltinOp::Ln && call->func_id != BuiltinOp::Log)) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::Unimplemented, "second factor is not a natural logarithm"));
    }
    auto argument = integrate_detail::extract_affine_argument(call->args.front(), var);
    if (!argument.has_value() || argument->coefficient.numerator().is_zero()) {
        return fail<ExprPtr>(integrate_detail::make_error(
            CASErrorKind::Unimplemented, "logarithm argument is not affine in the variable"));
    }

    const Rational a = exponent.coefficient;
    const Rational b = exponent.constant;
    const Rational root = -argument->constant / argument->coefficient;
    const Rational inverse_a = Rational(BigInt(1)) / a;

    ExprPtr x = arena.make<Symbol>(var);
    ExprPtr exponential = make_function(arena, "exp",
        {make_sum(arena, {make_product(arena, {make_rational(arena, a), x}), make_rational(arena, b)})});
    ExprPtr tail = make_product(arena, {
        make_function(arena, "exp", {make_rational(arena, a * root + b)}),
        make_function(arena, "Ei", {make_product(arena, {
            make_rational(arena, a),
            make_sum(arena, {x, make_rational(arena, -root)}),
        })}),
    });

    return ctx.simplify(make_product(arena, {
        make_rational(arena, view.coefficient),
        make_sum(arena, {
            make_product(arena, {make_rational(arena, inverse_a), exponential, log_factor.base}),
            make_unary(arena, UnaryOp::Neg,
                make_product(arena, {make_rational(arena, inverse_a), tail})),
        }),
    }));
}

}  // namespace

}  // namespace nonelementary

Result<ExprPtr> integrate_nonelementary_fallback(
    ExprPtr expr, const Symbol& var, symbolic::CASContext& ctx) {
    if (auto irq = ctx.check_interrupt(); irq.is_error()) return fail<ExprPtr>(irq.error());

    const nonelementary::ProductView view = nonelementary::flatten_product(expr);

    // Le famiglie sono disgiunte per forma (l'esponenziale gaussiano ha
    // argomento di grado 2, quello di Ei grado 1, le altre non hanno
    // esponenziale affatto), quindi l'ordine serve solo a fermarsi al primo
    // match: nessuna famiglia può mascherarne un'altra.
    if (auto result = nonelementary::try_gaussian(view, var, ctx); result.is_ok()) {
        return result;
    }
    if (auto result = nonelementary::try_trig_over_linear(view, var, ctx); result.is_ok()) {
        return result;
    }
    if (auto result = nonelementary::try_log_integral(view, var, ctx); result.is_ok()) {
        return result;
    }
    if (auto result = nonelementary::try_dilogarithm(view, var, ctx); result.is_ok()) {
        return result;
    }
    if (auto result = nonelementary::try_exponential_rational(view, var, ctx); result.is_ok()) {
        return result;
    }
    if (auto result = nonelementary::try_exponential_log(view, var, ctx); result.is_ok()) {
        return result;
    }

    return make_unimplemented<ExprPtr>(
        "calculus", "integrate_nonelementary_fallback",
        "integrand outside the Ei/Si/Ci/Shi/Chi/li/Li2/erfi families",
        cas::error::reason_codes::INTEGRATE_NO_STRATEGY,
        "Extend Nonelementary_Antiderivatives.md §5 with this family, or see "
        "§6 for the declared boundary (complex/algebraic argument)",
        "A43");
}

}  // namespace cas::calculus
