// A46 — LogToReal / Rioboo (Bronstein, Symbolic_Integration_I.md §2.8): resa in
// forma reale CHIUSA della parte logaritmica quando il resultant di
// Rothstein-Trager e' irriducibile di grado > 2.
//
// Spec: LogToAtan righe 2188-2200, LogToReal righe 2381-2395, Teorema 2.8.4
// (Rioboo) righe 2320-2340.
//
// Nessun valore atteso hardcoded: i test verificano il certificato D(F) = f,
// proprieta' strutturali (assenza di RootSum, fattorizzazione che si
// ri-moltiplica al polinomio di partenza) e il conteggio esatto di radici reali.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/ast_debug.hpp"
#include "cas/symbolic.hpp"

#include "algebra/polynomial_internal.hpp"
#include "algebra/partial_fractions_rioboo.hpp"

#include <chrono>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra {
namespace {

class LogToRealA46Test : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& src) {
        Lexer lexer(src);
        auto tokens = lexer.tokenize();
        EXPECT_TRUE(tokens.is_ok()) << "lex failed for: " << src;
        if (tokens.is_error()) return ExprPtr{};
        Parser parser(tokens.value(), ctx.arena());
        auto parsed = parser.parse();
        EXPECT_TRUE(parsed.is_ok()) << "parse failed for: " << src;
        return parsed.is_ok() ? parsed.value() : ExprPtr{};
    }

    [[nodiscard]] static bool contains_root_sum(ExprPtr expr) {
        if (!expr) return false;
        if (const auto* call = expr_cast<FuncCall>(expr)) {
            if (call->func_id == BuiltinOp::RootSum) return true;
            for (auto arg : call->args) if (contains_root_sum(arg)) return true;
            return false;
        }
        if (const auto* sum = expr_cast<Sum>(expr)) {
            for (auto term : sum->terms) if (contains_root_sum(term)) return true;
            return false;
        }
        if (const auto* product = expr_cast<Product>(expr)) {
            for (auto factor : product->factors) if (contains_root_sum(factor)) return true;
            return false;
        }
        if (const auto* binary = expr_cast<Binary>(expr)) {
            return contains_root_sum(binary->left) || contains_root_sum(binary->right);
        }
        if (const auto* unary = expr_cast<Unary>(expr)) return contains_root_sum(unary->operand);
        return false;
    }

    // Certificato multi-punto ESATTO: D(F) − f valutata in k punti razionali.
    // Non e' un ripiego numerico — la differenza e' una funzione razionale di
    // grado noto (num e den limitati da deg(Q) + deg(F)), quindi annullarsi in
    // piu' punti del suo grado la costringe a essere identicamente nulla. Serve
    // perche' ridurre D(F) − f a una sola frazione in x, sulla forma che produce
    // la pipeline via frazioni parziali, sfora qualunque budget ragionevole
    // (misurato: >120s), mentre valutare in un punto e' immediato.
    [[nodiscard]] bool antiderivative_is_correct_pointwise(
        ExprPtr F, const std::string& p_src, const std::string& q_src, const Symbol& var,
        const std::vector<int> points = {1, 2, 3, -1, -2, 5, 7, -4, 9, -6, 11, 13}) {
        auto dF = calculus::diff(F, var, 1U, ctx);
        if (dF.is_error()) {
            ADD_FAILURE() << "diff failed: " << dF.error().message;
            return false;
        }
        ExprPtr integrand = ctx.arena().make<Binary>(BinaryOp::Div, parse(p_src), parse(q_src));
        ExprPtr difference = ctx.arena().make<Binary>(BinaryOp::Sub, dF.value(), integrand);
        for (int point : points) {
            ExprPtr value = ctx.arena().make<IntegerLit>(BigInt(point));
            auto at_point = ctx.substitute(difference, var, value);
            if (at_point.is_error()) {
                ADD_FAILURE() << "substitute failed in x=" << point << ": "
                              << at_point.error().message;
                return false;
            }
            auto combined = together(at_point.value(), ctx);
            ExprPtr candidate = combined.is_ok() ? combined.value() : at_point.value();
            auto simplified = ctx.simplify(candidate);
            if (simplified.is_error()) {
                ADD_FAILURE() << "simplify failed in x=" << point << ": "
                              << simplified.error().message;
                return false;
            }
            if (!is_zero_poly(PolyExpr({simplified.value()}))) {
                ADD_FAILURE() << "D(F) - f != 0 in x=" << point << ": "
                              << debug_print(simplified.value()) << "\n  F = " << debug_print(F);
                return false;
            }
        }
        return true;
    }

};

// Il conteggio di Sturm a ±∞ deve essere esatto: e' la decisione che stabilisce
// se un quartico ha residui coniugati (nessuna radice reale) o no.
TEST_F(LogToRealA46Test, SturmRealRootCountIsExact) {
    const std::vector<std::pair<std::string, std::size_t>> cases = {
        {"z^2 + 1", 0U},          // coppia coniugata
        {"z^2 - 1", 2U},
        {"z^4 + 1", 0U},          // il caso di ∫dx/(x^4+1)
        {"256*z^4 + 1", 0U},      // il resultant stesso
        {"z^4 - 2", 2U},          // irriducibile su Q ma con due radici reali
        {"z^3 - z", 3U},
        {"z^2 - 2*z + 1", 1U},    // radice doppia: Sturm conta le DISTINTE
        {"1296*z^4 - 36*z^2 + 1", 0U},
    };
    Symbol z("z");
    for (const auto& [src, expected] : cases) {
        auto poly = parse_polynomial(parse(src), z, ctx);
        ASSERT_TRUE(poly.is_ok()) << src;
        auto rational = poly_to_rational_poly(poly.value());
        ASSERT_TRUE(rational.is_ok()) << src;
        EXPECT_EQ(count_real_roots_rational(rational.value()), expected) << src;
    }
}

// La fattorizzazione in due quadratici reali non viene mai assunta: se torna,
// deve ri-moltiplicarsi esattamente al quartico monico di partenza.
TEST_F(LogToRealA46Test, QuarticSplitReproducesTheInputWhenItSucceeds) {
    Symbol z("z");
    const std::vector<std::string> in_class = {"256*z^4 + 1", "1296*z^4 - 36*z^2 + 1", "z^4 + 1"};
    for (const auto& src : in_class) {
        auto poly = parse_polynomial(parse(src), z, ctx);
        ASSERT_TRUE(poly.is_ok()) << src;
        auto split = real_quadratic_factors_of_quartic(poly.value(), z, ctx);
        ASSERT_TRUE(split.is_ok()) << src;
        ASSERT_EQ(split.value().size(), 2U) << "atteso split in due quadratici reali per " << src;

        ExprPtr one = poly_make_integer(ctx.arena(), 1);
        PolyExpr first({split.value()[0].second, split.value()[0].first, one});
        PolyExpr second({split.value()[1].second, split.value()[1].first, one});
        auto product = poly_multiply(first, second, ctx);
        ASSERT_TRUE(product.is_ok()) << src;
        auto product_expr = polynomial_to_expr(product.value(), z, ctx);
        ASSERT_TRUE(product_expr.is_ok()) << src;

        auto rational = poly_to_rational_poly(poly.value());
        ASSERT_TRUE(rational.is_ok()) << src;
        auto monic = normalize_poly_monic(poly.value(), ctx);
        ASSERT_TRUE(monic.is_ok()) << src;
        auto monic_expr = polynomial_to_expr(monic.value(), z, ctx);
        ASSERT_TRUE(monic_expr.is_ok()) << src;

        auto difference = expand(ctx.arena().make<Binary>(BinaryOp::Sub,
            product_expr.value(), monic_expr.value()), ctx);
        ASSERT_TRUE(difference.is_ok()) << src;
        auto simplified = ctx.simplify(difference.value());
        ASSERT_TRUE(simplified.is_ok()) << src;
        EXPECT_TRUE(is_zero_poly(PolyExpr({simplified.value()})))
            << src << ": prodotto dei quadratici != quartico, resto "
            << debug_print(simplified.value());
    }
}

// Fuori dalla classe supportata la funzione deve dire NO (vettore vuoto), non
// tirare a indovinare: un quartico con radici reali non ha due fattori
// quadratici a residui coniugati.
TEST_F(LogToRealA46Test, QuarticSplitDeclinesOutOfClassInputs) {
    Symbol z("z");
    const std::vector<std::string> out_of_class = {"z^4 - 2", "z^4 - 5*z^2 + 1"};
    for (const auto& src : out_of_class) {
        auto poly = parse_polynomial(parse(src), z, ctx);
        ASSERT_TRUE(poly.is_ok()) << src;
        auto split = real_quadratic_factors_of_quartic(poly.value(), z, ctx);
        ASSERT_TRUE(split.is_ok()) << src;
        EXPECT_TRUE(split.value().empty())
            << src << ": atteso rifiuto (radici reali presenti)";
    }
}

// Il cuore di A46: resultant irriducibile di grado 4 -> forma chiusa reale,
// non piu' RootSum, e il certificato D(F) = f deve valere.
TEST_F(LogToRealA46Test, QuarticResultantYieldsVerifiableClosedForm) {
    Symbol x("x");
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"1", "x^4 + 1"},
        {"x^2", "x^4 + 1"},
        {"1", "x^6 + 1"},
    };
    for (const auto& [p_src, q_src] : cases) {
        auto result = integrate_rational_lrt(parse(p_src), parse(q_src), x, ctx);
        ASSERT_TRUE(result.is_ok()) << "LRT failed on (" << p_src << ")/(" << q_src
                                    << "): " << result.error().message;
        EXPECT_FALSE(contains_root_sum(result.value()))
            << "(" << p_src << ")/(" << q_src << ") deve essere in forma chiusa, got "
            << debug_print(result.value());
        EXPECT_TRUE(antiderivative_is_correct_pointwise(result.value(), p_src, q_src, x))
            << "D(F) != f per (" << p_src << ")/(" << q_src << ")";
    }
}

// La pipeline completa (integrate) non deve piu' rifiutare questi integrali:
// il rifiuto di A45 scattava sul RootSum, che ora non viene piu' prodotto.
TEST_F(LogToRealA46Test, IntegratePipelineReturnsClosedFormForQuarticDenominators) {
    Symbol x("x");
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"1", "x^4 + 1"},
        {"x^2", "x^4 + 1"},
        {"1", "x^6 + 1"},
    };
    for (const auto& [p_src, q_src] : cases) {
        ExprPtr integrand = ctx.arena().make<Binary>(BinaryOp::Div, parse(p_src), parse(q_src));
        auto result = calculus::integrate(integrand, x, ctx);
        ASSERT_TRUE(result.is_ok()) << "integrate failed on (" << p_src << ")/(" << q_src
                                    << "): " << result.error().message;
        EXPECT_FALSE(contains_root_sum(result.value()))
            << "(" << p_src << ")/(" << q_src << ") deve essere in forma chiusa";
        EXPECT_TRUE(antiderivative_is_correct_pointwise(result.value(), p_src, q_src, x))
            << "D(F) != f per (" << p_src << ")/(" << q_src << ")";
    }
}

// La forma sintattica dell'integranda non deve decidere il costo: `Pow(Q,−1)` e
// `Product(N, Pow(D,−1))` devono entrare nel pipeline razionale come `Div(N,D)`.
// Bound larghi (non misure di riferimento): servono a cogliere il ritorno del
// blow-up misurato prima di A46 — 40s per 1/(x⁴+1) scritto come Pow, 41s per
// 1/(x⁶+1), 157s e poi FALLIMENTO sulla seconda frazione parziale.
TEST_F(LogToRealA46Test, RationalShapeDoesNotDecideCost) {
    Symbol x("x");
    const std::vector<std::string> denominators = {"x^4 + 1", "x^6 + 1", "x^4 - x^2 + 1"};
    for (const auto& den : denominators) {
        symbolic::CASContext fresh;   // a caldo la seconda misura e' 0ms (cache)
        auto parse_fresh = [&fresh](const std::string& src) {
            Lexer lexer(src);
            auto tokens = lexer.tokenize();
            Parser parser(tokens.value(), fresh.arena());
            return parser.parse().value();
        };
        // stessa integranda, tre forme sintattiche diverse
        ExprPtr as_div = fresh.arena().make<Binary>(BinaryOp::Div,
            parse_fresh("1"), parse_fresh(den));
        ExprPtr as_pow = fresh.arena().make<Binary>(BinaryOp::Pow,
            parse_fresh(den), fresh.arena().make<IntegerLit>(BigInt(-1)));
        ExprPtr as_product = fresh.arena().make<Product>(std::vector<ExprPtr>{
            parse_fresh("x^2"),
            fresh.arena().make<Binary>(BinaryOp::Pow, parse_fresh(den),
                fresh.arena().make<IntegerLit>(BigInt(-1)))});

        for (auto [label, integrand] : std::vector<std::pair<const char*, ExprPtr>>{
                {"Div", as_div}, {"Pow^-1", as_pow}, {"Product", as_product}}) {
            const auto started = std::chrono::steady_clock::now();
            auto result = calculus::integrate(integrand, x, fresh);
            const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - started);
            EXPECT_TRUE(result.is_ok())
                << "integrate ha fallito su 1/(" << den << ") in forma " << label
                << ": " << (result.is_error() ? result.error().message : "");
            EXPECT_LT(elapsed.count(), 30)
                << "1/(" << den << ") in forma " << label << " ha impiegato "
                << elapsed.count() << "s: la forma sta di nuovo decidendo il costo";
        }
    }
}

}  // namespace
}  // namespace cas::algebra
