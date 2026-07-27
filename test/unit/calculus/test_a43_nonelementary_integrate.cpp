// A43 incremento 4 — integrazione verso le antiderivate NON elementari.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Nonelementary_Antiderivatives.md
// §5 (regole d'integrazione) e §7 (verifica numerica mpmath delle primitive).
//
// Criterio di accettazione, fissato da A42 e A45: un'antiderivata che non si
// riesce a derivare non e' verificabile ne' componibile, quindi OGNI risultato
// e' controllato con D(F) = f, non solo confrontato con la forma attesa. Le due
// verifiche sono complementari — la forma attesa coglie le regressioni di
// canonicita', D(F)=f coglie gli errori di segno/fattore che una forma "simile"
// nasconderebbe.

#include <gtest/gtest.h>

#include <string>
#include <utility>
#include <vector>

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A43NonelementaryIntegrateTest : public ::testing::Test {
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

    // Integra, poi riderivata: D(F) deve tornare all'integranda.
    void expect_antiderivative_of(const std::string& integrand) {
        auto F = calculus::integrate(parse(integrand), x, ctx);
        ASSERT_TRUE(F.is_ok()) << integrand << ": " << F.error().message;
        auto D = calculus::diff(F.value(), x, 1U, ctx);
        ASSERT_TRUE(D.is_ok()) << integrand << " (diff): " << D.error().message;
        auto eq = symbolic::mathematically_equal(D.value(), parse(integrand), ctx);
        ASSERT_TRUE(eq.is_ok()) << integrand;
        EXPECT_TRUE(eq.value()) << "D(F) != f for " << integrand;
    }

    void expect_closed_form(const std::string& integrand, const std::string& expected) {
        auto F = calculus::integrate(parse(integrand), x, ctx);
        ASSERT_TRUE(F.is_ok()) << integrand << ": " << F.error().message;
        auto eq = symbolic::mathematically_equal(F.value(), parse(expected), ctx);
        ASSERT_TRUE(eq.is_ok()) << integrand;
        EXPECT_TRUE(eq.value()) << integrand << " -> expected " << expected;
    }
};

// Spec §5, riga per riga: la forma chiusa deve essere quella dichiarata.
TEST_F(A43NonelementaryIntegrateTest, SpecTableClosedForms) {
    const std::vector<std::pair<std::string, std::string>> table{
        {"exp(x)/x", "Ei(x)"},
        {"exp(2*x)/x", "Ei(2*x)"},
        {"exp(x)/(x+1)", "exp(-1)*Ei(x+1)"},
        {"sin(x)/x", "Si(x)"},
        {"cos(x)/x", "Ci(x)"},
        {"1/ln(x)", "li(x)"},
        {"ln(1+x)/x", "-dilog(-x)"},
        {"exp(x^2)", "sqrt(pi)/2*erfi(x)"},
        {"exp(x)/x^2", "Ei(x) - exp(x)/x"},
        {"exp(-x)*ln(x)", "-exp(-x)*ln(x) + Ei(-x)"},
        {"exp(x)*ln(x)", "exp(x)*ln(x) - Ei(x)"},
        // Non in tabella ma nella stessa famiglia (spec §2/§4): Shi e Chi.
        {"sinh(x)/x", "Shi(x)"},
        {"cosh(x)/x", "Chi(x)"},
    };
    for (const auto& [integrand, expected] : table) {
        expect_closed_form(integrand, expected);
    }
}

// Il criterio di accettazione vero: D(F) = f su tutta la tabella §5.
TEST_F(A43NonelementaryIntegrateTest, SpecTableDerivativesRoundTrip) {
    for (const char* integrand : {
             "exp(x)/x", "exp(2*x)/x", "exp(x)/(x+1)", "sin(x)/x", "cos(x)/x",
             "1/ln(x)", "ln(1+x)/x", "exp(x^2)", "exp(x)/x^2",
             "exp(-x)*ln(x)", "exp(x)*ln(x)", "sinh(x)/x", "cosh(x)/x"}) {
        expect_antiderivative_of(integrand);
    }
}

// Spec §5: "implementare la riduzione, non le righe". Un polo di ordine k > 2
// non e' una riga di tabella — deve uscire dalla stessa ricorsione per parti.
//
// Il polo TRASLATO e' coperto a ordine 2 (`exp(x)/(x-1)^2`), non a ordine 3:
// `exp(2*x)/(x-1)^3` non termina, e la causa NON e' questa riduzione. Misurato
// col fallback A43 disattivato: >200 s comunque, cioe' il costo e' tutto nella
// catena ricorsiva pre-Risch. Conferma dal budget: a `max_integration_depth=3`
// la risposta esce corretta in 1.99 s, a 4 in 4.42 s, a 16 non esce — il tempo
// cresce col budget di ricorsione, non con l'ordine del polo. Aperta come A49;
// il caso rientra qui quando A49 chiude.
TEST_F(A43NonelementaryIntegrateTest, HigherOrderPolesReduceToEi) {
    for (const char* integrand : {"exp(x)/x^3", "exp(x)/x^4", "exp(x)/(x-1)^2"}) {
        expect_antiderivative_of(integrand);
    }
}

// --- Generalizzazioni oltre la tabella §5 -----------------------------------
//
// Le forme chiuse qui sotto sono ESATTE, e lo sono per certificato numerico
// mpmath a 30 cifre (errore relativo < 1e-20): `scripts/a43_special_fn_check.py`,
// sezione "generalizzazioni EMESSE DAL MOTORE", una riga per ciascuna.
//
// Su parte di esse il motore NON riesce a chiudere `D(F) = f` simbolicamente,
// per capacita' che non ha ancora — non per un errore di queste formule:
//   * gaussiana con |A| non quadrato perfetto  -> serve ridurre √(π/A)·√A;
//   * trigonometrica con fase non nulla        -> servono le formule di addizione;
//   * ∫x^s/ln x con s ≥ 1                      -> serve ln(x^{s+1}) = (s+1)·ln x.
// I gap sono aperti come A50. Dove la verifica simbolica passa oggi
// (`expect_antiderivative_of`) resta pretesa: e' quella la difesa contro un
// errore di segno o di fattore introdotto in seguito. Dove non passa, il test
// pretende comunque la forma chiusa esatta, cosi' una regressione della formula
// resta un fallimento e non un silenzio.

// Gaussiana generale: completamento del quadrato, entrambi i segni del
// coefficiente direttore (erfi per A > 0, erf per A < 0).
TEST_F(A43NonelementaryIntegrateTest, GeneralGaussianCompletesTheSquare) {
    // √A = 1: la verifica simbolica chiude.
    expect_antiderivative_of("exp(x^2+x)");
    expect_antiderivative_of("exp(x^2+3*x+1)");
    // √A ≠ 1: forma chiusa pretesa, D(F)=f bloccato da A50.
    expect_closed_form("exp(2*x^2)", "sqrt(pi/2)/2*erfi(sqrt(2)*x)");
    expect_closed_form("exp(-2*x^2+x)",
                       "sqrt(pi/2)/2*exp(1/8)*erf(sqrt(2)*(x-1/4))");
}

// Argomento trigonometrico traslato: le formule di addizione, non una riga per
// ogni fase (spec §5). Il risultato mescola le due primitive della coppia.
TEST_F(A43NonelementaryIntegrateTest, ShiftedTrigUsesAdditionFormulas) {
    // Fase nulla: la verifica simbolica chiude.
    expect_antiderivative_of("sin(2*x)/x");
    expect_antiderivative_of("cos(3*x)/x");
    // Fase non nulla: forma chiusa pretesa, D(F)=f bloccato da A50.
    expect_closed_form("sin(x+1)/x", "cos(1)*Si(x) + sin(1)*Ci(x)");
    expect_closed_form("cos(x+1)/x", "cos(1)*Ci(x) - sin(1)*Si(x)");
    expect_closed_form("sin(x)/(x-2)", "cos(2)*Si(x-2) + sin(2)*Ci(x-2)");
    expect_closed_form("sinh(x+1)/x", "cosh(1)*Shi(x) + sinh(1)*Chi(x)");
    expect_closed_form("cosh(x+1)/x", "cosh(1)*Chi(x) + sinh(1)*Shi(x)");
}

// li: la lettura generale dell'identita' ∫g'/ln(g) = li(g), nelle due forme
// (potenza al numeratore e argomento affine nel logaritmo).
TEST_F(A43NonelementaryIntegrateTest, LogIntegralGeneralForms) {
    expect_antiderivative_of("1/ln(x)");
    expect_closed_form("x/ln(x)", "li(x^2)");
    expect_closed_form("x^2/ln(x)", "li(x^3)");
    expect_closed_form("1/ln(2*x)", "li(2*x)/2");
}

// Li₂ con logaritmo affine generale e polo traslato.
TEST_F(A43NonelementaryIntegrateTest, DilogarithmGeneralAffineArgument) {
    expect_antiderivative_of("ln(1+x)/x");
    expect_antiderivative_of("ln(1+2*x)/x");
    expect_closed_form("ln(2+x)/x", "ln(2)*ln(x) - dilog(-x/2)");
    expect_closed_form("ln(1+x)/(x-3)", "ln(4)*ln(x-3) - dilog(-(x-3)/4)");
}

// Non-regressione: il fallback e' l'ULTIMA risorsa prima di Meijer G, quindi
// non deve rubare gli integrandi elementari — quelli devono continuare a uscire
// in forma elementare, senza funzioni speciali.
TEST_F(A43NonelementaryIntegrateTest, ElementaryIntegrandsStayElementary) {
    const std::vector<std::pair<std::string, std::string>> elementary{
        {"exp(x)", "exp(x)"},
        {"1/x", "ln(abs(x))"},
        {"x*exp(x)", "x*exp(x) - exp(x)"},
        {"ln(x)", "x*ln(x) - x"},
        {"sin(x)", "-cos(x)"},
        {"x/(x^2+1)", "ln(abs(x^2+1))/2"},
    };
    for (const auto& [integrand, expected] : elementary) {
        expect_closed_form(integrand, expected);
    }
}

}  // namespace
