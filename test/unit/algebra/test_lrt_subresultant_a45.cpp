// A45 — Lazard-Rioboo-Trager: correttezza della catena subresultante e
// dell'antiderivata razionale che ne deriva.
//
// Spec: Symbolic_Integration_I.md (Bronstein) — SubResultant righe 1044-1071,
// IntRationalLogPart righe 1815-1833, correzione Mulders righe 1829-1830.
//
// Il bug chiuso da A45: `resultant(A,B) = R_k` (ultimo elemento della PRS di
// grado 0) vale SOLO se `deg(R_{k-1}) = 1`; quando l'ultimo passo e' difettoso
// serve il fattore tau_k del Teorema 1.5.3, e senza di esso il resultant esce
// di grado troppo basso. Nessun valore atteso e' hardcoded qui: i test
// verificano proprieta' (D(F) = f, struttura della catena), non stringhe.

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/ast_debug.hpp"
#include "cas/symbolic.hpp"

#include "algebra/polynomial_internal.hpp"

#include <chrono>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace cas::algebra {
namespace {

class LrtSubresultantA45Test : public ::testing::Test {
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

    // Certificato: d/dx F(x) - P/Q deve ridursi a zero. Confronto su forma
    // razionale comune (together), l'idioma gia' usato altrove nel motore per
    // evitare falsi negativi da forme non confluenti.
    [[nodiscard]] bool antiderivative_is_correct(ExprPtr F, const std::string& p_src,
                                                 const std::string& q_src, const Symbol& var) {
        auto dF = calculus::diff(F, var, 1U, ctx);
        if (dF.is_error()) {
            ADD_FAILURE() << "diff failed: " << dF.error().message << "\n  F = " << debug_print(F);
            return false;
        }
        ExprPtr integrand = ctx.arena().make<Binary>(BinaryOp::Div, parse(p_src), parse(q_src));
        ExprPtr diff_node = ctx.arena().make<Binary>(BinaryOp::Sub, dF.value(), integrand);
        auto together_res = together(diff_node, ctx);
        if (together_res.is_error()) {
            ADD_FAILURE() << "together failed: " << together_res.error().message;
            return false;
        }
        auto simplified = ctx.simplify(together_res.value());
        if (simplified.is_error()) {
            ADD_FAILURE() << "simplify failed: " << simplified.error().message;
            return false;
        }
        const bool zero = is_zero_poly(PolyExpr({simplified.value()}));
        if (!zero) {
            ADD_FAILURE() << "D(F) - f = " << debug_print(simplified.value())
                          << "\n  F = " << debug_print(F);
        }
        return zero;
    }

    // Un risultato in forma RootSum e' la somma formale sulle radici di
    // Bronstein (spec riga 1831): corretto, ma non derivabile dal motore
    // (diff non implementa RootSum), quindi il certificato D(F)=f non e'
    // applicabile. Serve solo a separare i due gruppi di casi.
    [[nodiscard]] static bool is_root_sum(ExprPtr F) {
        const auto* call = expr_cast<FuncCall>(F);
        return call != nullptr && call->func_id == BuiltinOp::RootSum;
    }
};

// Sequenze PRS con ultimo passo difettoso (deg(R_{k-1}) > 1): e' il caso in cui
// serve il fattore tau_k. Due dei quattro integrali golden difettosi
// (integrate/basic.jsonl:117-120) hanno resultant 4(8z^2±1)^2, quindi fattore
// quadratico e forma chiusa verificabile; gli altri due hanno resultant
// irriducibile di grado 4 e finiscono in RootSum (vedi test successivo).
TEST_F(LrtSubresultantA45Test, DefectiveLastStepIntegralsAreCorrect) {
    Symbol x("x");
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"x^2 + 1", "x^4 + 1"},
        {"x^2 - 1", "x^4 + 1"},
        {"x^3", "x^4 + 1"},
    };
    for (const auto& [p_src, q_src] : cases) {
        auto result = integrate_rational_lrt(parse(p_src), parse(q_src), x, ctx);
        ASSERT_TRUE(result.is_ok()) << "LRT failed on (" << p_src << ")/(" << q_src
                                    << "): " << result.error().message;
        ASSERT_FALSE(is_root_sum(result.value()))
            << "(" << p_src << ")/(" << q_src << ") must reduce to closed form";
        EXPECT_TRUE(antiderivative_is_correct(result.value(), p_src, q_src, x))
            << "D(F) != f for (" << p_src << ")/(" << q_src << ")";
    }
}

// Controllo di non-regressione sulle sequenze regolari.
TEST_F(LrtSubresultantA45Test, RegularSequenceIntegralsStayCorrect) {
    Symbol x("x");
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"1", "x^2 + 1"},
        {"1", "x^2 - 1"},
        {"x", "x^4 + 1"},
        {"1", "x^3 - x"},
        {"1", "x^2 + x + 1"},
        {"x", "x^2 + 1"},
    };
    for (const auto& [p_src, q_src] : cases) {
        auto result = integrate_rational_lrt(parse(p_src), parse(q_src), x, ctx);
        ASSERT_TRUE(result.is_ok()) << "LRT failed on (" << p_src << ")/(" << q_src
                                    << "): " << result.error().message;
        if (is_root_sum(result.value())) continue;
        EXPECT_TRUE(antiderivative_is_correct(result.value(), p_src, q_src, x))
            << "D(F) != f for (" << p_src << ")/(" << q_src << ")";
    }
}

// Quando il resultant e' irriducibile di grado > 2 su Q, la spec riga 1831
// prescrive la somma formale sulle radici e il motore emette RootSum.
//
// A46 ha spostato il confine: i quartici i cui residui vivono in un'estensione
// QUADRATICA di Q (cubico risolvente con radice razionale) ora escono in forma
// chiusa reale — 256z^4+1 di ∫dx/(x^4+1) e' uno di quelli, e il suo caso vive in
// test_logtoreal_a46.cpp. Resta RootSum quando i residui non sono esprimibili
// per radicali reali: qui il resultant di ∫dx/(x^8+1) (grado 8, Phi_16
// irriducibile). Il test fissa il confine, non un valore.
TEST_F(LrtSubresultantA45Test, ResidueOutsideRadicalsStillYieldsFormalRootSum) {
    Symbol x("x");
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"1", "x^8 + 1"},
    };
    for (const auto& [p_src, q_src] : cases) {
        auto result = integrate_rational_lrt(parse(p_src), parse(q_src), x, ctx);
        ASSERT_TRUE(result.is_ok()) << "LRT failed on (" << p_src << ")/(" << q_src
                                    << "): " << result.error().message;
        EXPECT_TRUE(is_root_sum(result.value()))
            << "expected the formal RootSum of spec line 1831 for (" << p_src << ")/(" << q_src
            << "), got " << debug_print(result.value());
    }
}

// Proprieta' strutturali della catena esposta da resultant_generic:
// i gradi devono essere strettamente decrescenti (definizione di PRS) e
// l'ultimo elemento deve avere grado 0 quando il resultant e' non nullo.
TEST_F(LrtSubresultantA45Test, ChainOutIsAStrictlyDecreasingPrs) {
    Symbol x("x");
    Symbol z_var = ctx.make_fresh_symbol("lrt_t");
    ExprPtr z_sym = ctx.arena().make<Symbol>(z_var.name);

    const std::vector<std::pair<std::string, std::string>> cases = {
        {"x^2 + 1", "x^4 + 1"},  // ultimo passo difettoso
        {"x^2", "x^4 + 1"},      // regolare
        {"x", "x^5 + x + 1"},    // difettoso a meta' catena
    };
    for (const auto& [p_src, q_src] : cases) {
        auto q_poly = parse_polynomial(parse(q_src), x, ctx);
        ASSERT_TRUE(q_poly.is_ok());
        auto q_prime = calculus::diff(parse(q_src), x, 1U, ctx);
        ASSERT_TRUE(q_prime.is_ok());
        ExprPtr target = ctx.arena().make<Binary>(
            BinaryOp::Sub, parse(p_src),
            ctx.arena().make<Binary>(BinaryOp::Mul, z_sym, q_prime.value()));
        auto target_poly = parse_polynomial(target, x, ctx);
        ASSERT_TRUE(target_poly.is_ok());

        std::vector<std::vector<ExprPtr>> chain;
        auto res = resultant_generic<ExprPtr>(q_poly.value().coefficients(),
                                              target_poly.value().coefficients(), &ctx,
                                              std::nullopt, &chain);
        ASSERT_TRUE(res.is_ok()) << "resultant failed on (" << p_src << ")/(" << q_src << ")";
        ASSERT_GE(chain.size(), 2U) << "chain must contain at least R_0 and R_1";

        std::size_t previous = chain.front().size();
        for (std::size_t i = 1; i < chain.size(); ++i) {
            EXPECT_LT(chain[i].size(), previous)
                << "PRS degrees must strictly decrease, at index " << i << " for (" << p_src
                << ")/(" << q_src << ")";
            previous = chain[i].size();
        }
        EXPECT_EQ(chain.back().size(), 1U)
            << "last PRS element must have degree 0 when the resultant is nonzero";
    }
}

// La regola della spec (riga 1828) e' applicabile: per ogni molteplicita' i
// della fattorizzazione squarefree del resultant esiste un elemento di catena
// di grado x pari a i (oppure i = deg(Q), gestito dalla riga 1826).
TEST_F(LrtSubresultantA45Test, SubresultantOfDegreeMatchingMultiplicityExists) {
    Symbol x("x");
    Symbol z_var = ctx.make_fresh_symbol("lrt_t");
    ExprPtr z_sym = ctx.arena().make<Symbol>(z_var.name);

    const std::vector<std::pair<std::string, std::string>> cases = {
        {"x^2 + 1", "x^4 + 1"},
        {"x^2 - 1", "x^4 + 1"},
        {"1", "x^4 + 1"},
        {"1", "x^3 - x"},
    };
    for (const auto& [p_src, q_src] : cases) {
        auto q_poly = parse_polynomial(parse(q_src), x, ctx);
        ASSERT_TRUE(q_poly.is_ok());
        auto q_prime = calculus::diff(parse(q_src), x, 1U, ctx);
        ASSERT_TRUE(q_prime.is_ok());
        ExprPtr target = ctx.arena().make<Binary>(
            BinaryOp::Sub, parse(p_src),
            ctx.arena().make<Binary>(BinaryOp::Mul, z_sym, q_prime.value()));
        auto target_poly = parse_polynomial(target, x, ctx);
        ASSERT_TRUE(target_poly.is_ok());

        std::vector<std::vector<ExprPtr>> chain;
        auto res = resultant_generic<ExprPtr>(q_poly.value().coefficients(),
                                              target_poly.value().coefficients(), &ctx,
                                              std::nullopt, &chain);
        ASSERT_TRUE(res.is_ok());
        auto sf = square_free_factorization(res.value(), z_var, ctx);
        ASSERT_TRUE(sf.is_ok()) << "squarefree factorization of the resultant failed";

        const std::size_t deg_q = poly_degree(q_poly.value());
        for (const auto& factor : sf.value().factors) {
            auto factor_poly = parse_polynomial(factor.factor, z_var, ctx);
            if (factor_poly.is_error() || poly_degree(factor_poly.value()) == 0U) continue;
            const std::size_t multiplicity = static_cast<std::size_t>(factor.multiplicity);
            if (multiplicity == deg_q) continue;  // spec riga 1826: S_i = D
            bool found = false;
            for (const auto& step : chain) {
                PolyExpr candidate{step};
                normalize_poly(candidate);
                if (!is_zero_poly(candidate) && poly_degree(candidate) == multiplicity) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << "no PRS element of degree " << multiplicity << " for (" << p_src
                               << ")/(" << q_src << ")";
        }
    }
}


// A45 fase 3 — regressione di performance. Il debito originale ("il resultant
// di x^2/(x^4+1) impiega 8.3s e sfora l'ops-budget") era in realta' un effetto
// del bug di correttezza: la ricorrenza sbagliata produceva divisioni non
// esatte e quindi coefficienti razionali che gonfiavano a ogni passo. Con la
// ricorrenza corretta il costo e' sceso di due ordini di grandezza e
// l'aritmetica esatta RatPoly non serve. Il test fissa il risultato misurato
// cosi' che una regressione futura sia visibile subito; il bound e' largo
// (10x il misurato: 197/633/406 ms su M1 Pro) per non essere flaky sotto
// carico — vedi memoria "quick-suite flaky load-dependent".
TEST_F(LrtSubresultantA45Test, HigherDegreeDenominatorsStayFast) {
    Symbol x("x");
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"1", "x^6 + 1"},
        {"1", "x^5 + x + 1"},
        {"1", "x^8 + 1"},
    };
    for (const auto& [p_src, q_src] : cases) {
        const auto started = std::chrono::steady_clock::now();
        auto result = integrate_rational_lrt(parse(p_src), parse(q_src), x, ctx);
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - started);
        EXPECT_TRUE(result.is_ok()) << "LRT failed on (" << p_src << ")/(" << q_src
                                    << "): " << result.error().message;
        EXPECT_LT(elapsed.count(), 8000)
            << "LRT on (" << p_src << ")/(" << q_src << ") took " << elapsed.count()
            << " ms — the coefficient swell of A45 may be back";
    }
}

}  // namespace
}  // namespace cas::algebra
