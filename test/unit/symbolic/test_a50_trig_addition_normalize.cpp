// A50 (identità 2) — formule di addizione trig/iperboliche in
// `mathematically_equal` (`src/algebra/trig_addition_normalize.cpp`).
//
// `trig_addition_normalize` è privato (`algebra_internal.hpp`, come
// `hyperbolic_normalize`/`nonelementary_normalize` che segue come pattern):
// testato solo tramite l'API pubblica `mathematically_equal`, mai incluso
// direttamente — stessa convenzione del resto del file.
//
// Sblocca `D(F)=f` su `∫sin(x+c)/x` (fase non nulla) e sulla famiglia a polo
// traslato `∫sin(x)/(x-c)` (verificato end-to-end in
// test_a43_nonelementary_integrate.cpp, ShiftedTrigUsesAdditionFormulas).

#include <gtest/gtest.h>

#include <string>

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A50TrigAdditionNormalizeTest : public ::testing::Test {
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

    void expect_eq(const std::string& lhs, const std::string& rhs, bool expected = true) {
        auto eq = symbolic::mathematically_equal(parse(lhs), parse(rhs), ctx);
        ASSERT_TRUE(eq.is_ok()) << lhs << " vs " << rhs;
        EXPECT_EQ(eq.value(), expected) << lhs << (expected ? " deve" : " NON deve")
                                         << " uguagliare " << rhs;
    }
};

// sin(P+C) = sin(P)cos(C) + cos(P)sin(C), fase positiva.
TEST_F(A50TrigAdditionNormalizeTest, SinAdditionPositivePhase) {
    expect_eq("sin(x+1)", "sin(x)*cos(1) + cos(x)*sin(1)");
}

// cos(P+C) = cos(P)cos(C) - sin(P)sin(C), fase positiva.
TEST_F(A50TrigAdditionNormalizeTest, CosAdditionPositivePhase) {
    expect_eq("cos(x+1)", "cos(x)*cos(1) - sin(x)*sin(1)");
}

// Fase NEGATIVA (il caso che ha esposto il bug del letterale negativo non
// unificato con Unary(Neg,...) dalle regole di parità esistenti):
// sin(x-2) = sin(x)cos(2) - cos(x)sin(2)  [sin dispari: sin(-2)=-sin(2)]
// cos(x-2) = cos(x)cos(2) + sin(x)sin(2)  [cos pari: cos(-2)=cos(2)]
TEST_F(A50TrigAdditionNormalizeTest, SinCosAdditionNegativePhase) {
    expect_eq("sin(x-2)", "sin(x)*cos(2) - cos(x)*sin(2)");
    expect_eq("cos(x-2)", "cos(x)*cos(2) + sin(x)*sin(2)");
}

// Iperboliche, fase positiva e negativa (sinh dispari, cosh pari).
TEST_F(A50TrigAdditionNormalizeTest, HyperbolicAddition) {
    expect_eq("sinh(x+1)", "sinh(x)*cosh(1) + cosh(x)*sinh(1)");
    expect_eq("cosh(x+1)", "cosh(x)*cosh(1) + sinh(x)*sinh(1)");
    expect_eq("sinh(x-2)", "sinh(x)*cosh(2) - cosh(x)*sinh(2)");
    expect_eq("cosh(x-2)", "cosh(x)*cosh(2) - sinh(x)*sinh(2)");
}

// La famiglia a polo traslato (il repro reale di A50): la derivata di
// cos(c)*Si(x-c)+sin(c)*Ci(x-c) deve tornare a sin(x)/(x-c) — verificato qui
// direttamente sull'identità di addizione che la rende possibile, non solo
// end-to-end via integrate (quel percorso è in test_a43).
TEST_F(A50TrigAdditionNormalizeTest, ShiftedPoleCrossTermsCancel) {
    // cos(2)*sin(x-2) + sin(2)*cos(x-2) deve ridursi a sin(x) — richiede sia
    // l'espansione dell'addizione sia l'identità pitagorica cos²+sin²=1.
    expect_eq("cos(2)*sin(x-2) + sin(2)*cos(x-2)", "sin(x)");
}

// Controllo negativo — la formula di addizione non deve provare uguaglianze
// FALSE: un'espansione errata (es. segno scambiato) sarebbe un silent-wrong.
TEST_F(A50TrigAdditionNormalizeTest, WrongSignIsNotProvenEqual) {
    expect_eq("sin(x+1)", "sin(x)*cos(1) - cos(x)*sin(1)", /*expected=*/false);
    expect_eq("cos(x-2)", "cos(x)*cos(2) - sin(x)*sin(2)", /*expected=*/false);
}

}  // namespace
