// A56 — rientro della scorciatoia `N₁/D + N₂/D = (N₁+N₂)/D` in `add_parts`/
// `subtract_parts` (`src/algebra/factorization_num_den.cpp`), guardata da
// `structural_equal(lhs.denominator, rhs.denominator)`.
//
// Perche' non basta testare `together()`: la funzione pubblica applica SEMPRE
// `reduce_rational_by_gcd` dopo `apart_num_den` (spec
// .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Together_Polynomial_GCD_Reduction.md),
// quindi anche senza la scorciatoia un denominatore cross-moltiplicato
// (y+1)^2 verrebbe comunque ridotto a (y+1) a valle — il test su `together()`
// in test_a47_together_cost.cpp (EqualDenominatorsAreNotCrossMultiplied) non
// distingue le due strade. Qui si chiama `apart_num_den` direttamente, PRIMA
// della riduzione GCD, per verificare che la scorciatoia sia davvero quella a
// evitare il quadrato del denominatore.
//
// Rientrata dopo A54 (falso negativo del verificatore IBP, ora chiuso): senza
// A54, questa stessa scorciatoia faceva fallire
// `RischLimitedIntegrate.EndToEnd_NestedLogTower_Solved` — quel test resta la
// prova di non-regressione end-to-end (rieseguito nel gate, verde).

#include <gtest/gtest.h>

#include <string>

#include "cas/algebra.hpp"
#include "cas/ast_nodes.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A56TogetherShortcutTest : public ::testing::Test {
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

// `a/(y+1) + b/(y+1)`: denominatori uguali. Il denominatore di apart_num_den
// (PRIMA della riduzione GCD di `together`) deve restare `(y+1)`, non
// `(y+1)^2` — la prova diretta che `add_parts` ha preso la scorciatoia.
TEST_F(A56TogetherShortcutTest, AddPartsEqualDenominatorsSkipsCrossMultiply) {
    auto parts = algebra::apart_num_den(parse("a/(y+1) + b/(y+1)"), ctx);
    ASSERT_TRUE(parts.is_ok()) << parts.error().message;

    ExprPtr expected_den = simplified("y+1");
    ExprPtr expected_num = simplified("a+b");
    EXPECT_TRUE(structural_equal(parts.value().denominator, expected_den))
        << "denominatore atteso (y+1), non cross-moltiplicato";
    EXPECT_TRUE(structural_equal(parts.value().numerator, expected_num));
}

// `a/(y+1) - b/(y+1)`: stesso argomento, ma e' il trigger empirico misurato
// (il grosso del guadagno viene da `subtract_parts`, non da `add_parts`).
TEST_F(A56TogetherShortcutTest, SubtractPartsEqualDenominatorsSkipsCrossMultiply) {
    auto parts = algebra::apart_num_den(parse("a/(y+1) - b/(y+1)"), ctx);
    ASSERT_TRUE(parts.is_ok()) << parts.error().message;

    ExprPtr expected_den = simplified("y+1");
    ExprPtr expected_num = simplified("a-b");
    EXPECT_TRUE(structural_equal(parts.value().denominator, expected_den))
        << "denominatore atteso (y+1), non cross-moltiplicato";
    EXPECT_TRUE(structural_equal(parts.value().numerator, expected_num));
}

// Controllo negativo: denominatori DIVERSI devono continuare a cross-
// moltiplicare — la guardia `structural_equal` non deve scattare per errore
// su forme solo simili.
TEST_F(A56TogetherShortcutTest, DifferentDenominatorsStillCrossMultiply) {
    auto parts = algebra::apart_num_den(parse("a/(y+1) + b/(y+2)"), ctx);
    ASSERT_TRUE(parts.is_ok()) << parts.error().message;

    ExprPtr single_factor_1 = simplified("y+1");
    ExprPtr single_factor_2 = simplified("y+2");
    EXPECT_FALSE(structural_equal(parts.value().denominator, single_factor_1));
    EXPECT_FALSE(structural_equal(parts.value().denominator, single_factor_2));
}

// Il contratto matematico resta quello di sempre — la scorciatoia non cambia
// il valore, solo la forma intermedia prima della riduzione GCD.
TEST_F(A56TogetherShortcutTest, ShortcutPreservesValue) {
    const char* inputs[] = {
        "a/(y+1) + b/(y+1)",
        "a/(y+1) - b/(y+1)",
        "1/(x*y) + 1/(x*y) + 2/(x*y)",
        "x/(x^2+1) - 1/(x^2+1)",
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

}  // namespace
