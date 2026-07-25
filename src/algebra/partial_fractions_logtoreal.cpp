// A46 — LogToReal: resa in forma reale CHIUSA della parte logaritmica quando il
// resultant di Rothstein-Trager e' irriducibile di grado > 2.
//
// Spec: Symbolic_Integration_I.md (Bronstein) — LogToAtan righe 2188-2200,
// LogToReal righe 2381-2395, Teorema 2.8.4 (Rioboo) righe 2320-2340.
//
// La somma formale sui residui (`RootSum`) e' matematicamente corretta ma non
// derivabile: per renderla in forma chiusa servono i residui, e i residui di un
// quartico generico non sono esprimibili per radicali reali. Sono esprimibili
// esattamente quando vivono in un'estensione QUADRATICA di Q, cioe' quando il
// cubico risolvente del quartico ha una radice razionale: in quel caso il
// quartico si spezza in due quadratici reali e ciascuno passa per la forma
// chiusa di Rioboo grado-2 che il motore ha gia'. Fuori da quella classe la
// funzione dichiara di non saper fare e il chiamante torna al `RootSum`.

#include "partial_fractions_rioboo.hpp"

#include "cas/algebra.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"

#include <utility>
#include <vector>

namespace cas::algebra {

namespace {

[[nodiscard]] Rational rat(long long value) { return Rational(BigInt(value)); }

[[nodiscard]] Result<PolyExpr> scale_poly(const PolyExpr& poly, ExprPtr scalar,
                                          symbolic::CASContext& ctx) {
    return poly_multiply(poly, PolyExpr({scalar}), ctx);
}

// LogToAtan(A, B) — spec righe 2188-2200. Restituisce una somma di arctangenti
// di POLINOMI f tale che df/dx = d/dx [ i·log((A+iB)/(A−iB)) ].
//
// Perche' la ricorsione e non il solo primo ramo: fermarsi a 2·arctan(A/B)
// quando B non divide A da' un argomento che e' una funzione razionale, con due
// conseguenze misurate — l'antiderivata ha salti dove il denominatore si annulla
// (Rioboo la costruisce continua) e derivarla per verificarla costa ordini di
// grandezza in piu' (certificato D(F)=f in timeout con coefficienti radicali).
//
// Terminazione: max(deg C, deg D) < max(deg A, deg B) a ogni passo, quindi il
// bound e' derivato dall'input (nessuna costante magica).
[[nodiscard]] Result<ExprPtr> log_to_atan(
    PolyExpr A, PolyExpr B, const Symbol& var, symbolic::CASContext& ctx, std::size_t budget) {
    AstArena& arena = ctx.arena();
    normalize_poly(A);
    normalize_poly(B);

    if (is_zero_poly(B)) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "LogToAtan richiede B != 0 (spec Symbolic_Integration_I.md:2192)",
            .hint = std::nullopt,
        });
    }
    if (budget == 0U) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "LogToAtan: ricorsione oltre il bound derivato dai gradi di A e B",
            .hint = std::nullopt,
        });
    }

    ExprPtr two = poly_make_integer(arena, 2);
    auto atan_of = [&](ExprPtr argument) -> Result<ExprPtr> {
        auto simplified = ctx.simplify(argument);
        ExprPtr arg = simplified.is_ok() ? simplified.value() : argument;
        return ok(rioboo_mul(arena, two,
            arena.make<FuncCall>(BuiltinOp::Atan, std::vector<ExprPtr>{arg})));
    };

    // if B | A then return 2·arctan(A/B)
    auto division = divide_poly_with_remainder(A, B, ctx);
    if (division.is_ok() && is_zero_poly(division.value().remainder)) {
        auto quotient = polynomial_to_expr(division.value().quotient, var, ctx);
        if (quotient.is_error()) return fail<ExprPtr>(quotient.error());
        return atan_of(quotient.value());
    }

    // if deg(A) < deg(B) then return LogToAtan(−B, A)
    if (poly_degree(A) < poly_degree(B)) {
        auto negated = poly_negate(B, ctx);
        if (negated.is_error()) return fail<ExprPtr>(negated.error());
        return log_to_atan(negated.value(), std::move(A), var, ctx, budget - 1U);
    }

    // (D, C, G) <- ExtendedEuclidean(B, −A), cioe' B·D − A·C = G
    auto neg_A = poly_negate(A, ctx);
    if (neg_A.is_error()) return fail<ExprPtr>(neg_A.error());
    auto bezout = poly_extended_gcd(B, neg_A.value(), ctx);
    if (bezout.is_error()) return fail<ExprPtr>(bezout.error());
    const PolyExpr D = bezout.value().s;
    const PolyExpr C = bezout.value().t;
    const PolyExpr G = bezout.value().gcd;
    if (is_zero_poly(G) || is_zero_poly(C)) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "LogToAtan: Bezout degenere (C o G nullo), fuori dalle ipotesi del "
                       "Teorema 2.8.1",
            .hint = std::nullopt,
        });
    }

    auto AD = poly_multiply(A, D, ctx);
    auto BC = poly_multiply(B, C, ctx);
    if (AD.is_error()) return fail<ExprPtr>(AD.error());
    if (BC.is_error()) return fail<ExprPtr>(BC.error());
    auto numerator = poly_add(AD.value(), BC.value(), ctx);
    if (numerator.is_error()) return fail<ExprPtr>(numerator.error());

    // (A·D + B·C)/G e' un polinomio per il Teorema 2.8.1; se la divisione non e'
    // esatta si tiene il quoziente razionale invece di mentire sul risultato.
    ExprPtr argument;
    auto exact = divide_poly_with_remainder(numerator.value(), G, ctx);
    if (exact.is_ok() && is_zero_poly(exact.value().remainder)) {
        auto quotient = polynomial_to_expr(exact.value().quotient, var, ctx);
        if (quotient.is_error()) return fail<ExprPtr>(quotient.error());
        argument = quotient.value();
    } else {
        auto num_expr = polynomial_to_expr(numerator.value(), var, ctx);
        auto den_expr = polynomial_to_expr(G, var, ctx);
        if (num_expr.is_error()) return fail<ExprPtr>(num_expr.error());
        if (den_expr.is_error()) return fail<ExprPtr>(den_expr.error());
        argument = arena.make<Binary>(BinaryOp::Div, num_expr.value(), den_expr.value());
    }

    auto head = atan_of(argument);
    if (head.is_error()) return fail<ExprPtr>(head.error());
    auto tail = log_to_atan(D, C, var, ctx, budget - 1U);
    if (tail.is_error()) return fail<ExprPtr>(tail.error());
    return ok(arena.make<Sum>(std::vector<ExprPtr>{head.value(), tail.value()}));
}

}  // namespace

ExprPtr rioboo_mul(AstArena& arena, ExprPtr a, ExprPtr b) {
    if (!a) return b;
    if (!b) return a;
    if (is_zero_poly(PolyExpr({a})) || is_zero_poly(PolyExpr({b}))) {
        return poly_make_integer(arena, 0);
    }
    return arena.make<Binary>(BinaryOp::Mul, a, b);
}


Result<ExprPtr> rioboo_quadratic_real_form(
    ExprPtr a,
    ExprPtr b,
    const PolyExpr& G_z_x,
    const Symbol& var,
    const Symbol& z_var,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr one = poly_make_integer(arena, 1);
    ExprPtr two = poly_make_integer(arena, 2);
    PolyExpr monic_r({b, a, one});

    // Q₁·z + Q₀ = S(z,x) mod (z² + a·z + b), coefficiente per coefficiente in x.
    PolyExpr Q1;
    PolyExpr Q0;
    Q1.resize(G_z_x.size(), poly_make_integer(arena, 0));
    Q0.resize(G_z_x.size(), poly_make_integer(arena, 0));
    for (std::size_t i = 0; i < G_z_x.size(); ++i) {
        if (!G_z_x[i]) continue;
        auto g_i_z = parse_polynomial(G_z_x[i], z_var, ctx);
        if (g_i_z.is_error()) continue;
        auto division = divide_poly_with_remainder(g_i_z.value(), monic_r, ctx);
        if (division.is_error()) return fail<ExprPtr>(division.error());
        const PolyExpr& remainder = division.value().remainder;
        if (remainder.size() > 0U && remainder[0]) Q0[i] = remainder[0];
        if (remainder.size() > 1U && remainder[1]) Q1[i] = remainder[1];
    }
    normalize_poly(Q1);
    normalize_poly(Q0);

    auto Q1_res = polynomial_to_expr(Q1, var, ctx);
    ExprPtr Q1_x = Q1_res.is_ok() ? Q1_res.value() : poly_make_integer(arena, 0);
    auto Q0_res = polynomial_to_expr(Q0, var, ctx);
    ExprPtr Q0_x = Q0_res.is_ok() ? Q0_res.value() : poly_make_integer(arena, 0);

    // disc = 4b − a² > 0 per una coppia di residui coniugati.
    auto disc_res = ctx.simplify(arena.make<Binary>(BinaryOp::Sub,
        rioboo_mul(arena, poly_make_integer(arena, 4), b),
        rioboo_mul(arena, a, a)));
    if (disc_res.is_error()) return fail<ExprPtr>(disc_res.error());
    ExprPtr sqrt_disc = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{disc_res.value()});
    auto sqrt_simplified = ctx.simplify(sqrt_disc);
    if (sqrt_simplified.is_ok()) sqrt_disc = sqrt_simplified.value();

    // A46: `expand` prima di `simplify` non e' cosmetica. Con a, b radicali il
    // norm arriva come (1/16)·(4x−4√2)² + … e `simplify` da sola non distribuisce
    // il quadrato: il logaritmo resta in una forma che a valle fa esplodere
    // `together`/`diff` (misurato: certificato D(F)=f in timeout). Distribuito,
    // lo stesso norm e' x² − √2·x + 1.
    auto normalize = [&ctx](ExprPtr expr) -> Result<ExprPtr> {
        auto expanded = expand(expr, ctx);
        if (expanded.is_error()) return ctx.simplify(expr);
        return ctx.simplify(expanded.value());
    };

    auto norm_res = normalize(arena.make<Sum>(std::vector<ExprPtr>{
        rioboo_mul(arena, Q0_x, Q0_x),
        rioboo_mul(arena, arena.make<Unary>(UnaryOp::Neg, a), rioboo_mul(arena, Q1_x, Q0_x)),
        rioboo_mul(arena, b, rioboo_mul(arena, Q1_x, Q1_x))}));
    if (norm_res.is_error()) return fail<ExprPtr>(norm_res.error());

    auto half_neg_a = ctx.simplify(arena.make<Binary>(BinaryOp::Div,
        arena.make<Unary>(UnaryOp::Neg, a), two));
    if (half_neg_a.is_error()) return fail<ExprPtr>(half_neg_a.error());

    // log(c·u) = log(u) + costante: dentro un'antiderivata il fattore costante e'
    // libero, e tenerlo costa. Nella pipeline via frazioni parziali il norm arriva
    // con contenuto 1/729 e ogni coefficiente del logaritmo se lo porta dietro,
    // rendendo il risultato (e la sua verifica) molto piu' grosso del necessario.
    ExprPtr norm_expr = norm_res.value();
    if (auto norm_poly = parse_polynomial(norm_expr, var, ctx);
        norm_poly.is_ok() && poly_degree(norm_poly.value()) > 0U) {
        if (auto monic = normalize_poly_monic(norm_poly.value(), ctx); monic.is_ok()) {
            if (auto monic_expr = polynomial_to_expr(monic.value(), var, ctx); monic_expr.is_ok()) {
                norm_expr = monic_expr.value();
            }
        }
    }

    ExprPtr ln_term = rioboo_mul(arena, half_neg_a.value(),
        arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{norm_expr}));

    // I due residui sono −a/2 ± i·√disc/2, quindi S(α,x) = Q₁·α + Q₀ si separa in
    //   A = Q₀ − (a/2)·Q₁      B = (√disc/2)·Q₁
    // e la parte reale della coppia coniugata e' (spec righe 2391-2393)
    //   (−a/2)·log(A² + B²) + (√disc/2)·LogToAtan(A, B),
    // con A² + B² = Q₀² − a·Q₀·Q₁ + b·Q₁² = norm calcolato sopra.
    auto half_sqrt_disc = ctx.simplify(arena.make<Binary>(BinaryOp::Div, sqrt_disc, two));
    if (half_sqrt_disc.is_error()) return fail<ExprPtr>(half_sqrt_disc.error());

    // LogToAtan vale su un campo K con √−1 ∉ K (spec riga 2191): serve che i due
    // residui siano una coppia CONIUGATA, cioe' disc = 4b − a² > 0. Con disc < 0
    // i residui sono reali, B = (√disc/2)·Q₁ diventa immaginario e la ricorsione
    // di Rioboo esce dalle ipotesi del Teorema 2.8.1 (misurato: 670s per
    // ∫(x²−1)/(x⁴+1) dx, il cui resultant 8z²−1 ha disc = −1/2). In quel caso si
    // usa la forma chiusa precedente, che il simplifier riporta ad atanh.
    //
    // Il segno di disc: se e' razionale la decisione e' esatta qui. Se non lo e',
    // l'unico produttore di a, b non razionali e' la via quartica, che ha GIA'
    // stabilito con Sturm che R non ha radici reali — quindi entrambi i fattori
    // quadratici hanno residui coniugati e disc > 0.
    bool residues_are_conjugate = true;
    if (auto disc_rational = poly_to_rational_poly(PolyExpr({disc_res.value()}));
        disc_rational.is_ok() && disc_rational.value().size() <= 1U) {
        residues_are_conjugate = rat(0) < disc_rational.value().constant_term();
    }

    ExprPtr atan_term;
    auto half_a = ctx.simplify(arena.make<Binary>(BinaryOp::Div, a, two));
    auto scaled_Q1 = half_a.is_ok() ? scale_poly(Q1, half_a.value(), ctx)
                                    : fail<PolyExpr>(half_a.error());
    auto A_poly = scaled_Q1.is_ok() ? poly_subtract(Q0, scaled_Q1.value(), ctx)
                                    : fail<PolyExpr>(scaled_Q1.error());
    auto B_poly = scale_poly(Q1, half_sqrt_disc.value(), ctx);
    if (residues_are_conjugate && A_poly.is_ok() && B_poly.is_ok()) {
        const std::size_t budget = poly_degree(A_poly.value()) + poly_degree(B_poly.value()) + 2U;
        auto atan_sum = log_to_atan(A_poly.value(), B_poly.value(), var, ctx, budget);
        if (atan_sum.is_ok()) {
            atan_term = rioboo_mul(arena, half_sqrt_disc.value(), atan_sum.value());
        }
    }

    if (!atan_term) {
        // Ripiego: 2·arctan(A/B) troncato al primo ramo di LogToAtan. Resta
        // un'antiderivata corretta (differisce per una costante a tratti), ma
        // con argomento razionale invece che polinomiale.
        auto atan_numerator = normalize(rioboo_mul(arena, sqrt_disc, Q1_x));
        auto atan_denominator = normalize(arena.make<Binary>(BinaryOp::Sub,
            rioboo_mul(arena, two, Q0_x),
            rioboo_mul(arena, a, Q1_x)));
        if (atan_numerator.is_error()) return fail<ExprPtr>(atan_numerator.error());
        if (atan_denominator.is_error()) return fail<ExprPtr>(atan_denominator.error());

        ExprPtr atan_arg = arena.make<Binary>(BinaryOp::Div,
            atan_numerator.value(), atan_denominator.value());
        auto atan_arg_simplified = ctx.simplify(atan_arg);
        if (atan_arg_simplified.is_ok()) atan_arg = atan_arg_simplified.value();

        atan_term = rioboo_mul(arena, arena.make<Unary>(UnaryOp::Neg, sqrt_disc),
            arena.make<FuncCall>(BuiltinOp::Atan, std::vector<ExprPtr>{atan_arg}));
    }

    return ok(arena.make<Sum>(std::vector<ExprPtr>{ln_term, atan_term}));
}


}  // namespace cas::algebra
