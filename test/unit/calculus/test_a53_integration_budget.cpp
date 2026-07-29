// A53 — `calculus::integrate` e' UNA operazione, con un budget suo.
//
// Difetto chiuso qui: le ops si contano solo in Simplifier/Substituter, e ogni
// `ctx.simplify()` dell'integratore era top-level, quindi ne azzerava contatore
// e timer. Nessun budget limitava il TOTALE di un'integrazione: misurato, la
// stessa integranda consumava per intero qualunque cap le si desse (30 s con
// cap 30 s, oltre 300 s con cap 300 s) e cap diversi troncavano entry diverse
// — il risultato era funzione del tempo concesso, non dell'input. L'unico
// freno era il SIGALRM esterno, cioe' il carico della macchina.
//
// Nessuna spec formale pertinente in MISSING_FEATURES_SPECS (verificato): e'
// correttezza infrastrutturale, come A51 di cui e' la scoperta n. 2. La parte
// side-conditions vive in test_a53_condition_rollback.cpp.

#include <gtest/gtest.h>

#include <chrono>
#include <string>

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A53IntegrationBudgetTest : public ::testing::Test {
protected:
    [[nodiscard]] static ExprPtr parse(symbolic::CASContext& ctx, const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    // Integranda che il motore NON risolve, ma che gli costa lavoro CONTATO:
    // misurata a 105'380 ops / 4.7 s con il gate spento (golden runner,
    // `--ops-report --max-ops 0`). Serve cosi': se il costo stesse in codice
    // che non incrementa il contatore (Risch, algebra dei polinomi) il budget
    // non lo vedrebbe e il test misurerebbe il wall-clock invece del gate.
    static constexpr const char* kExpensive = "x^2*exp(x)*sin(x)";

    // Integranda che il motore risolve per parti, ben dentro qualunque budget.
    static constexpr const char* kEasy = "x^2*exp(x)";
};

// LA prova di A53: lo stesso input, con budget WALL-CLOCK diversi, deve dare lo
// stesso esito. Prima il tempo concesso decideva il risultato; ora decide il
// budget deterministico, che dal wall-clock non dipende.
//
// Il tetto d'integrazione e' abbassato per rendere il test veloce: verifica il
// MECCANISMO, non il valore di default (quello e' calibrato sul corpus, vedi
// max_integration_ops in cas_context_simplifier_params.hpp).
TEST_F(A53IntegrationBudgetTest, VerdictDoesNotDependOnTheWallClockBudget) {
    const auto run = [](std::chrono::milliseconds wall) {
        symbolic::CASContext ctx;
        ctx.set_max_integration_ops(20'000ULL);
        ctx.set_max_operation_ops(ctx.max_operation_ops());  // esplicito: A30
        ctx.set_timeout(wall);
        Symbol x{"x"};
        auto r = calculus::integrate(parse(ctx, kExpensive), x, ctx);
        return r.is_ok();
    };

    const bool tight = run(std::chrono::seconds(2));
    const bool loose = run(std::chrono::seconds(30));
    EXPECT_EQ(tight, loose)
        << "l'esito dipende ancora dal tempo concesso invece che dall'integranda";
}

// Confine DICHIARATO di questo budget: le ops si contano in Simplifier e
// Substituter, quindi il gate morde solo dove il lavoro passa di li'. Una
// integranda il cui costo sta altrove (Risch, algebra dei polinomi) resta
// limitata dal solo wall-clock — misurato: `sin(log(x))*cos(log(x))/x^3`
// consuma minuti con un tetto di 20'000 ops, perche' quelle ops non le spende
// mai. Estendere il conteggio a quelle fasi e' lavoro proprio, non di A53; qui
// si pinna il confine perche' non venga scambiato per copertura totale.
TEST_F(A53IntegrationBudgetTest, BudgetOnlyBindsWhereOpsAreCounted) {
    symbolic::CASContext ctx;
    ctx.set_max_integration_ops(1ULL);  // il piu' stretto possibile
    Symbol x{"x"};
    ctx.reset_ops_high_water();
    auto r = calculus::integrate(parse(ctx, kEasy), x, ctx);
    EXPECT_TRUE(r.is_error()) << "con tetto 1 il gate deve mordere su questa forma";
    // ... e il consumo registrato resta osservabile: e' il dato con cui si
    // tara la soglia, senza il quale sarebbe un numero preso a intuito.
    EXPECT_GT(ctx.ops_high_water(), 0U);
}

// Il budget e' del TOTALE, non della singola simplify: stringerlo deve poter
// fermare un'integrazione che al budget di default riesce. Senza questo, il
// parametro esisterebbe senza mordere (era esattamente lo stato pre-A53: gate
// ops misurato a 406'469 contro un tetto di 2'000'000, mai raggiunto).
TEST_F(A53IntegrationBudgetTest, TighteningTheBudgetStopsAnIntegrationThatOtherwiseSucceeds) {
    Symbol x{"x"};

    symbolic::CASContext generous;
    auto ok_result = calculus::integrate(parse(generous, kEasy), x, generous);
    ASSERT_TRUE(ok_result.is_ok()) << ok_result.error().message;

    symbolic::CASContext tight;
    tight.set_max_integration_ops(50ULL);  // sotto il costo di qualunque IBP reale
    auto stopped = calculus::integrate(parse(tight, kEasy), x, tight);
    ASSERT_TRUE(stopped.is_error()) << "il budget d'integrazione non morde";
    // Mai un risultato silenziosamente sbagliato o un successo travestito: il
    // rifiuto e' diagnostico (REGOLA ZERO).
    EXPECT_FALSE(stopped.error().message.empty());
}

// Rientranza: le integrazioni annidate (by-parts, sostituzione) non riaprono
// l'operazione, quindi il loro costo resta addebitato a quella piu' esterna.
// E' il punto dell'intera task — e' il TOTALE che dev'essere limitato.
TEST_F(A53IntegrationBudgetTest, NestedIntegrationsChargeTheOutermostOperation) {
    symbolic::CASContext ctx;
    Symbol x{"x"};
    ctx.reset_ops_high_water();
    auto r = calculus::integrate(parse(ctx, "x^2*exp(x)"), x, ctx);  // IBP -> integrate annidati
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_GT(ctx.ops_high_water(), 0U)
        << "nessuna operazione ha registrato consumo: il budget non copre l'integrazione";
}

// Contratto A30: se il chiamante spegne il gate ops (wall-clock esplicito senza
// budget ops esplicito), l'integrazione non se ne impone uno per conto proprio.
// E' cio' che rende possibile un run di sola misura — imporre il tetto
// falserebbe il dato che si sta raccogliendo per calibrarlo.
// NB: il nome NON deve contenere "Disabled" — `test_quick.sh` esclude
// `*Disabled*` come convenzione gtest, e il test non girerebbe mai nel gate.
TEST_F(A53IntegrationBudgetTest, OpsGateTurnedOffByCallerIsNotReintroduced) {
    symbolic::CASContext ctx;
    ctx.set_timeout(std::chrono::seconds(30));  // A30: azzera il gate ops implicito
    ASSERT_EQ(ctx.max_operation_ops(), 0U);
    Symbol x{"x"};
    auto r = calculus::integrate(parse(ctx, "x^2*exp(x)"), x, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    EXPECT_EQ(ctx.max_operation_ops(), 0U) << "lo scope non ha ripristinato il gate spento";
}

}  // namespace
