// A49 — `∫e^{P(x)}·R(x)` con R dotato di polo genuino non entra nella catena IBP.
//
// Nessuna spec formale pertinente (verificato in MISSING_FEATURES_SPECS): il
// criterio è un teorema, non una feature da specificare. Per Liouville
// quell'integrale non ha primitiva elementare, quindi l'integrazione per parti
// non può chiuderlo: la sua forma chiusa vive nella famiglia Ei (fallback A43).
//
// Il difetto misurato prima del fix: la catena IBP partiva lo stesso, e ogni suo
// livello rilanciava l'INTERA pipeline sul sotto-integrale — Risch compreso (56%
// dei campioni, ricorsione profonda 13). Non erano né i passi IBP (18 contro 18
// col polo in 0) né la dimensione delle espressioni (6 contro 8 nodi) né le
// chiamate ad expand/parse_polynomial (+10%): era il costo per livello.
//
//   ∫e^x/(x−1)²    316'681 ops / 15.0 s  →  6'235 ops / 0.29 s
//   ∫e^x/(x+5)²    non terminava (>300 s) →  6'297 ops / 0.27 s
//   ∫e^{2x}/(x−1)³ 500'094 ops (tetto A53) / 35 s → 12'627 ops / 0.57 s

#include <gtest/gtest.h>

#include <string>

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A49ExpPoleTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

// Il costo non deve più dipendere dalla TRASLAZIONE del polo: è quella
// dipendenza — e non l'ordine del polo — che A49 ha isolato come difetto.
// La soglia è sul BUDGET DETERMINISTICO (ops), non sul tempo: un test a
// wall-clock misurerebbe la macchina, non il codice (lezione A51).
TEST_F(A49ExpPoleTest, TranslatedPoleCostsLikeThePoleAtZero) {
    const auto ops_for = [this](const std::string& integrand) {
        ctx.reset_ops_high_water();
        auto F = calculus::integrate(parse(integrand), x, ctx);
        EXPECT_TRUE(F.is_ok()) << integrand << ": " << F.error().message;
        return ctx.ops_high_water();
    };

    const std::uint64_t at_zero = ops_for("exp(x)/x^2");
    EXPECT_GT(at_zero, 0U);
    for (const char* translated : {"exp(x)/(x-1)^2", "exp(x)/(x+5)^2"}) {
        symbolic::CASContext fresh;  // niente cache condivisa fra i due
        Symbol fx{"x"};
        fresh.reset_ops_high_water();
        auto tokens = Lexer(translated).tokenize();
        ASSERT_TRUE(tokens.is_ok());
        Parser parser(tokens.value(), fresh.arena());
        auto parsed = parser.parse();
        ASSERT_TRUE(parsed.is_ok());
        auto F = calculus::integrate(parsed.value(), fx, fresh);
        ASSERT_TRUE(F.is_ok()) << translated << ": " << F.error().message;
        // Prima del fix il rapporto misurato era 5.6x (e per `x+5` divergeva);
        // 4x lascia margine al rumore di misura senza tollerare il difetto.
        EXPECT_LT(fresh.ops_high_water(), at_zero * 4U)
            << translated << " costa " << fresh.ops_high_water()
            << " ops contro " << at_zero << " del polo in 0";
    }
}

// Il guard non deve toccare ciò che l'integrazione per parti chiude davvero:
// senza polo, IBP resta la strada e il risultato non cambia.
TEST_F(A49ExpPoleTest, IntegrationByPartsStillHandlesPoleFreeProducts) {
    for (const char* integrand : {"x^2*exp(x)", "x*exp(x)", "x*log(x)"}) {
        auto F = calculus::integrate(parse(integrand), x, ctx);
        ASSERT_TRUE(F.is_ok()) << integrand << ": " << F.error().message;
        auto D = calculus::diff(F.value(), x, 1U, ctx);
        ASSERT_TRUE(D.is_ok()) << integrand;
        auto eq = symbolic::mathematically_equal(D.value(), parse(integrand), ctx);
        ASSERT_TRUE(eq.is_ok()) << integrand;
        EXPECT_TRUE(eq.value()) << "D(F) != f per " << integrand;
    }
}

// Un `Sum` che CONTIENE un polo non è un polo dell'integranda: il guard guarda i
// fattori, non il sottoalbero, e questo caso deve restare integrabile.
TEST_F(A49ExpPoleTest, SumContainingAPoleIsNotTreatedAsAPoleFactor) {
    auto F = calculus::integrate(parse("exp(x)*(x + 1/x)"), x, ctx);
    // Non si pretende la forma chiusa (il termine e^x/x è Ei): si pretende che
    // il motore non scambi questa forma per un fattore-polo e la rifiuti a
    // priori. Qualunque esito diverso da un errore di dispatch va bene.
    if (F.is_error()) {
        EXPECT_EQ(F.error().message.find("Liouville"), std::string::npos)
            << "il guard A49 ha catturato una forma che non è un fattore-polo";
    }
}

}  // namespace
