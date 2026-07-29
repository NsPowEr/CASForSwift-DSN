// A50 (identità 1/3) — generalizzazione di `sqrt(a)*sqrt(b) -> sqrt(a*b)`
// oltre i letterali razionali (`collapse_sqrt_pairs`, Phase 2,
// `simplify_arithmetic_chain_sqrt.cpp`).
//
// Prima del fix: `get_rat_pos` accettava SOLO `sqrt(IntegerLit)` o
// `sqrt(RationalLit)` come argomento — `sqrt(pi/2)` (argomento
// `Product([Pi, RationalLit(1,2)])`, non un letterale nudo) non veniva mai
// riconosciuto, quindi `sqrt(pi/2)*sqrt(2)` restava due fattori separati e
// non collassava a `sqrt(pi)`. Bloccava `D(F)=f` sulla gaussiana generale
// `∫e^{2x²}` (A > 0, |A| non quadrato perfetto).
//
// Fix: `extract_monomial` (la stessa decomposizione coefficiente-razionale +
// fattori-simbolici che Step 4 usa per la raccolta dei like-term) scompone
// l'argomento di ciascun sqrt, e la non-negatività si dimostra
// componendo `known_nonneg` (già usato da Phase 1 per `sqrt(a)*sqrt(a)->a`)
// su coefficiente e fattori — non un caso speciale per `pi`.

#include <gtest/gtest.h>

#include <string>

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A50SqrtProductGeneralTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    [[nodiscard]] ExprPtr simplified(const std::string& s) {
        auto r = ctx.simplify(parse(s));
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

// Il repro esatto che sblocca A50: sqrt(pi/2)*sqrt(2) -> sqrt(pi).
TEST_F(A50SqrtProductGeneralTest, RationalCoefficientTimesPiCollapses) {
    ExprPtr result = simplified("sqrt(pi/2)*sqrt(2)");
    ExprPtr expected = simplified("sqrt(pi)");
    EXPECT_TRUE(structural_equal(result, expected))
        << "sqrt(pi/2)*sqrt(2) deve collassare a sqrt(pi)";
}

// Simmetrico: sqrt(A)*sqrt(pi/A) per A intero > 1.
TEST_F(A50SqrtProductGeneralTest, IntegerTimesRationalOverPiCollapses) {
    ExprPtr result = simplified("sqrt(3)*sqrt(pi/3)");
    ExprPtr expected = simplified("sqrt(pi)");
    EXPECT_TRUE(structural_equal(result, expected));
}

// Il caso puramente letterale (comportamento storico) resta invariato — la
// generalizzazione lo contiene come caso con resto simbolico vuoto.
TEST_F(A50SqrtProductGeneralTest, PureLiteralCaseUnchanged) {
    ExprPtr result = simplified("sqrt(2)*sqrt(8)");
    // sqrt(16) = 4.
    const auto* il = expr_cast<IntegerLit>(result);
    ASSERT_NE(il, nullptr) << "sqrt(2)*sqrt(8) deve ridursi a un intero";
    EXPECT_EQ(il->value, BigInt(4));
}

// Controllo negativo — due costanti positive SENZA relazione (coefficiente
// banale su entrambe, nessuna base condivisa) NON devono fondersi: sarebbe
// la semplice concatenazione che `builtin_rewrite_algebraic.cpp` disfa nella
// direzione opposta (canonicalizza `sqrt(prodotto non-negativo)` nella forma
// SEPARATA) — le due regole oscillerebbero. La fusione generalizzata serve
// SOLO quando produce una riduzione vera (guardia in Phase 2).
TEST_F(A50SqrtProductGeneralTest, UnrelatedPositiveConstantsDoNotCollapse) {
    ExprPtr result = simplified("sqrt(pi)*sqrt(e)");
    ExprPtr would_be_wrong = simplified("sqrt(pi*e)");
    EXPECT_FALSE(structural_equal(result, would_be_wrong));
}

// Simbolo con assunzione di positività esplicita: sqrt(2*a)*sqrt(a) deve
// collassare (il fattore `a` compare a esponente dispari in entrambi, quindi
// richiede `known_nonneg(a)`, soddisfatta dall'assunzione).
TEST_F(A50SqrtProductGeneralTest, PositiveAssumedSymbolCollapses) {
    ctx.assumptions().assume_positive(Symbol("a"));
    ExprPtr result = simplified("sqrt(2*a)*sqrt(a)");
    ExprPtr expected = simplified("a*sqrt(2)");
    EXPECT_TRUE(structural_equal(result, expected))
        << "sqrt(2*a)*sqrt(a) con a>0 deve collassare a a*sqrt(2)";
}

// Controllo negativo — SOUNDNESS: senza alcuna informazione di segno su x, y,
// sqrt(x)*sqrt(y) NON deve collassare a sqrt(x*y) (se uno dei due fosse
// negativo, sqrt(x)*sqrt(y) e sqrt(x*y) divergono su C — l'identità vale solo
// per operandi non-negativi).
TEST_F(A50SqrtProductGeneralTest, UnconstrainedSymbolsDoNotCollapse) {
    ExprPtr result = simplified("sqrt(x)*sqrt(y)");
    ExprPtr would_be_wrong = simplified("sqrt(x*y)");
    EXPECT_FALSE(structural_equal(result, would_be_wrong))
        << "sqrt(x)*sqrt(y) non deve collassare senza prova di non-negatività";
}

// Controllo negativo — un solo fattore con segno provato non basta se l'ALTRO
// fattore del prodotto misto non lo è: sqrt(2*x)*sqrt(3) con x libero (segno
// ignoto) non deve collassare a sqrt(6*x).
TEST_F(A50SqrtProductGeneralTest, MixedUnknownFactorDoesNotCollapse) {
    ExprPtr result = simplified("sqrt(2*x)*sqrt(3)");
    ExprPtr would_be_wrong = simplified("sqrt(6*x)");
    EXPECT_FALSE(structural_equal(result, would_be_wrong));
}

}  // namespace
