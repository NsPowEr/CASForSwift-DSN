// A47 — costo di `together`: lavoro ripetuto in `apart_num_den`.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Together_Polynomial_GCD_Reduction.md
// (contratto di `together`: totale, output matematicamente equale, fallback
// silenzioso quando la riduzione GCD non si applica).
//
// Difetto corretto, dentro l'aggregazione a coppie: `multiply_exprs`/`add_exprs`
// semplificano ognuna il proprio risultato, ma l'espressione costruita finisce
// SEMPRE in `normalize_rational_parts`, che semplifica di nuovo prima di
// qualunque test di forma — sei `simplify` per combinazione invece di due,
// pagate a ogni passo sull'accumulato che cresce.
//
// Misura (suite WeierstrassSubstitutionTest): 31.1 s → 22.2 s complessivi;
// `∫dx/cos²` 17.7 s → 11.7 s (−34%); `∫dx/(2+sin 2x)` 12.7 s → 9.9 s (−22%).
//
// NB: la parte di A47 su «la verifica D(F)=f brucia l'intero budget» era già
// stata chiusa da A53 — misurato prima di toccare qualunque codice qui:
// `verify:simplify` da 60'118 ms (cap) a 70 ms, `verify:together` da 45'880 ms
// a 2'877 ms.

#include <gtest/gtest.h>

#include <string>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A47TogetherCostTest : public ::testing::Test {
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
};

// Denominatori uguali: `N₁/D + N₂/D = (N₁+N₂)/D`. Il denominatore NON deve
// essere quadrato dalla cross-moltiplicazione — è il lavoro che veniva creato e
// subito disfatto dalla riduzione GCD. Verifica strutturale sul risultato, non
// sul tempo: il tempo misura la macchina (lezione A51).
TEST_F(A47TogetherCostTest, EqualDenominatorsAreNotCrossMultiplied) {
    auto res = algebra::together(parse("a/(y+1) + b/(y+1)"), ctx);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // (a+b)/(y+1): il denominatore compare a esponente 1. Se fosse stato
    // cross-moltiplicato senza riduzione avremmo (y+1)^2 al denominatore.
    auto expected = ctx.simplify(parse("(a+b)/(y+1)"));
    ASSERT_TRUE(expected.is_ok());
    auto got = ctx.simplify(res.value());
    ASSERT_TRUE(got.is_ok());
    auto eq = symbolic::mathematically_equal(got.value(), expected.value(), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// Il contratto di `together` non cambia: l'output resta matematicamente uguale
// all'input su forme miste, incluse quelle che la riduzione GCD tocca davvero.
TEST_F(A47TogetherCostTest, TogetherStaysEquivalentOnMixedForms) {
    const char* inputs[] = {
        "1/x + 1/(x+1)",
        "a/(y+1) + b/(y+1)^2",
        "x/(x^2+1) - 1/(x^2+1)",
        "1/(x*y) + 1/(x*y) + 2/(x*y)",
        "(x^2-1)/(x-1)",
    };
    for (const char* in : inputs) {
        ExprPtr e = parse(in);
        auto tog = algebra::together(e, ctx);
        ASSERT_TRUE(tog.is_ok()) << in << ": " << tog.error().message;
        auto eq = symbolic::mathematically_equal(tog.value(), e, ctx);
        ASSERT_TRUE(eq.is_ok()) << in;
        EXPECT_TRUE(eq.value()) << "together ha cambiato il valore di " << in;
    }
}

// Zero e denominatore unitario restano i casi limite della spec §5: `together`
// è totale e non deve inventare errori su di essi.
TEST_F(A47TogetherCostTest, BoundaryCasesFromTheSpec) {
    for (const char* in : {"0/x + 0/(x+1)", "x/1 + y/1", "1/x - 1/x"}) {
        auto tog = algebra::together(parse(in), ctx);
        ASSERT_TRUE(tog.is_ok()) << in << ": " << tog.error().message;
        auto eq = symbolic::mathematically_equal(tog.value(), parse(in), ctx);
        ASSERT_TRUE(eq.is_ok()) << in;
        EXPECT_TRUE(eq.value()) << in;
    }
}

// NB — due difetti PREESISTENTI trovati misurando qui, nessuno introdotto da
// questo lavoro (verificato rieseguendo la misura con la modifica in stash):
//   * A54: la scorciatoia esatta `N₁/D + N₂/D = (N₁+N₂)/D` cambia la forma che
//     arriva a Risch e gli fa produrre un'antiderivata SBAGLIATA su
//     `∫(1/x + 1/(x·ln x))·ln(ln x)`. Per questo la scorciatoia — che varrebbe
//     un altro ~25% — resta disattivata: un guadagno di costo non si paga con
//     un silent-wrong.
//   * A55: `together` non riduce le potenze del denominatore comune ed espande
//     i binomi. Su `a/(y+1) + … + f/(y+1)` produce
//     `a·y⁵/(y+1)⁶ + 5a·y⁴/(y+1)⁶ + …`, cioe' `a·(y+1)⁵/(y+1)⁶` mai collassato
//     a `a/(y+1)`. Corretto ma non ridotto — ed e' anche parte del costo.
// Finche' A55 vive, una somma lunga a denominatore comune non ha forma attesa
// stabile da asserire, quindi qui non c'e' un test su di essa.

}  // namespace
