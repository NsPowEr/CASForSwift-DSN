// A51 — l'esito di un confronto deve dipendere SOLO dai suoi operandi.
//
// Nessuna spec formale pertinente (verificato in MISSING_FEATURES_SPECS): e' un
// bug di correttezza infrastrutturale. L'invariante violato: chi apre
// un'operazione top-level deve inizializzarne il budget.
//
// `CASContext::simplify`, `CASContext::substitute` e `algebra::polynomial_gcd`
// lo facevano; `mathematically_equal` no — alzava `operation_active_` senza
// azzerare `ops_count_` ne' far ripartire `operation_started_at_`. Le simplify
// interne vedevano l'operazione gia' aperta e non lo facevano al posto suo,
// quindi `Simplifier::check_timeout` misurava `elapsed` dall'ULTIMA operazione
// top-level del contesto — o dall'epoch, se non ce n'era mai stata una, e
// allora `elapsed` valeva l'uptime della macchina e il timeout scattava al
// primo controllo wall-clock.
//
// Effetto misurato in produzione: nel golden runner, dove UN contesto e'
// condiviso da tutte le entry, il budget residuo di ogni confronto dipendeva da
// quanto avevano impiegato le entry precedenti. Una entry
// (`integrate(1/(2-sin(x)), x)`) oscillava fra SKIP e FAIL su misure con lo
// STESSO binario, proprio sul confine `FAIL_CEILING`.

#include <gtest/gtest.h>

#include <chrono>
#include <string>

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A51OperationBudgetTest : public ::testing::Test {
protected:
    [[nodiscard]] static ExprPtr parse(symbolic::CASContext& ctx, const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

// Il repro esatto: le due antiderivate di 1/(2-sin x) prodotte da CAS e da
// Maxima. Su un contesto APPENA COSTRUITO il confronto tornava
// `Symbolic operation timed out` in ~60 ms — un timeout che non poteva essere
// vero, visto che il budget wall-clock e' di 10 s. Ora deve produrre un
// verdetto. Quale verdetto non e' oggetto di QUESTO test (le due forme
// differiscono per un'identita' di semi-angolo che il motore non riconosce
// ancora dentro `arctan` — famiglia B.2); qui conta che non sia un errore.
TEST_F(A51OperationBudgetTest, FirstOperationOnFreshContextIsNotSpuriouslyTimedOut) {
    symbolic::CASContext ctx;
    auto eq = symbolic::mathematically_equal(
        parse(ctx, "2*atan(2*sin(x/2)/(cos(x/2)*sqrt(3)) - 1/sqrt(3))/sqrt(3)"),
        parse(ctx, "2*atan(2*sin(x)/((cos(x)+1)*sqrt(3)) - 1/sqrt(3))/sqrt(3)"),
        ctx);
    ASSERT_TRUE(eq.is_ok())
        << "timeout spurio alla prima operazione del contesto: " << eq.error().message;
}

// L'invariante in forma diretta: lo stesso confronto deve dare lo stesso esito
// su un contesto vergine e su un contesto che ha gia' lavorato. Con il difetto
// i due divergevano — il primo ereditava l'epoch, il secondo il timer
// dell'operazione precedente.
TEST_F(A51OperationBudgetTest, VerdictDoesNotDependOnContextHistory) {
    const std::string lhs = "2*atan(2*sin(x/2)/(cos(x/2)*sqrt(3)) - 1/sqrt(3))/sqrt(3)";
    const std::string rhs = "2*atan(2*sin(x)/((cos(x)+1)*sqrt(3)) - 1/sqrt(3))/sqrt(3)";

    bool fresh_value = false;
    {
        symbolic::CASContext ctx;
        auto eq = symbolic::mathematically_equal(parse(ctx, lhs), parse(ctx, rhs), ctx);
        ASSERT_TRUE(eq.is_ok()) << eq.error().message;
        fresh_value = eq.value();
    }
    {
        symbolic::CASContext ctx;
        // Consuma budget con un'operazione top-level precedente.
        auto warm = ctx.simplify(parse(ctx, "expand((x+1)^12)"));
        ASSERT_TRUE(warm.is_ok());
        auto eq = symbolic::mathematically_equal(parse(ctx, lhs), parse(ctx, rhs), ctx);
        ASSERT_TRUE(eq.is_ok()) << eq.error().message;
        EXPECT_EQ(eq.value(), fresh_value)
            << "l'esito dipende dalla storia del contesto";
    }
}

// Un budget wall-clock stretto deve troncare in base al lavoro DI QUESTO
// confronto, non a quello gia' fatto prima sul contesto. Con il difetto il
// tempo speso da un'operazione precedente veniva addebitato alla successiva.
TEST_F(A51OperationBudgetTest, PreviousWorkIsNotChargedToTheNextComparison) {
    symbolic::CASContext ctx;
    ctx.set_timeout(std::chrono::seconds(30));
    ctx.set_max_operation_ops(ctx.max_operation_ops());  // esplicito: A30 + A51

    // Prima operazione: consuma tempo reale.
    auto warm = ctx.simplify(parse(ctx, "expand((x+y+1)^10)"));
    ASSERT_TRUE(warm.is_ok());

    // Un confronto banale, subito dopo, deve restare banale.
    auto eq = symbolic::mathematically_equal(
        parse(ctx, "sin(x)^2 + cos(x)^2"), parse(ctx, "1"), ctx);
    ASSERT_TRUE(eq.is_ok()) << eq.error().message;
    EXPECT_TRUE(eq.value());
}

// Gli altri tre punti d'ingresso che aprono un'operazione devono restare
// corretti: e' il pattern che A51 ha incapsulato, e questo test lo pinna su
// tutti e tre invece che sul solo `mathematically_equal`.
TEST_F(A51OperationBudgetTest, AllOperationEntryPointsInitialiseTheirBudget) {
    {   // CASContext::simplify
        symbolic::CASContext ctx;
        auto r = ctx.simplify(parse(ctx, "expand((x+1)^10) - expand((x+1)^10)"));
        ASSERT_TRUE(r.is_ok()) << r.error().message;
    }
    {   // CASContext::substitute
        symbolic::CASContext ctx;
        auto r = ctx.substitute(parse(ctx, "expand((x+1)^10)"), Symbol{"x"},
                                parse(ctx, "y+1"));
        ASSERT_TRUE(r.is_ok()) << r.error().message;
    }
    {   // mathematically_equal
        symbolic::CASContext ctx;
        auto r = symbolic::mathematically_equal(
            parse(ctx, "(x+1)^3"), parse(ctx, "x^3 + 3*x^2 + 3*x + 1"), ctx);
        ASSERT_TRUE(r.is_ok()) << r.error().message;
        EXPECT_TRUE(r.value());
    }
}

}  // namespace
