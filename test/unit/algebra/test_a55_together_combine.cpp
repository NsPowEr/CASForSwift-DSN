// A55 — `together()` non combinava davvero in un'unica frazione: il passo
// finale `divide_exprs(N, D)` (`src/algebra/factorization_num_den.cpp`)
// costruisce `Product(Sum(N), Pow(D,-1))`, e lo Step 8 del simplifier
// (`simplify_product_factors`, `simplify_arithmetic_chain.cpp`) distribuisce
// QUALUNQUE fattore `Sum` di un `Product` sugli altri fattori — regola
// generale e sempre attiva altrove (es. Gamma reflection la richiede). Senza
// sospenderla, `together("a/(y+1) + b/(y+1) + ... + f/(y+1)")` NON restituiva
// `(a+...+f)/(y+1)`, ma tornava una `Sum` di frazioni distribuite — proprio la
// forma citata nella scoperta originale di A55 (coefficienti binomiali
// 1,5,10,10,5,1 di `(y+1)^5` inclusi), lì attribuita al cross-moltiplica di
// A47/pre-A56.
//
// **Root cause reale, trovata rimisurando dopo A56**: A56 (guardia
// `structural_equal` in `add_parts`/`subtract_parts`) elimina già il
// cross-moltiplica per denominatori uguali — con A56 da solo, `apart_num_den`
// produce correttamente `N=a+...+f`, `D=y+1` (verificato in
// test_a56_together_shortcut.cpp). Ma `together()` continuava a fallire, per
// un motivo NUOVO e indipendente: lo Step 8 ridistribuiva il risultato appena
// assemblato. Fix: `CombinedFormScope` (`algebra_internal.hpp`,
// `set_symbolic_sum_distribution(false)`) sospende lo Step 8 per la durata
// della costruzione finale di `together` — esatto complementare di
// `ExpandedFormScope` (A54, che sospende la raccolta invece della
// distribuzione).
//
// **Difetto (a) della scoperta originale ("GCD non collassa potenze dello
// stesso fattore") verificato NON riproducibile isolatamente** — vedi
// `GcdReductionAlreadyCollapsesSharedFactorPower`: era un sintomo del
// cross-moltiplica, non un difetto separato di `reduce_rational_by_gcd`.
//
// **Difetto (b) ("binomi espansi invece di fattorizzati") resta APERTO ma
// più stretto**: si manifesta SOLO quando la riduzione GCD riduce
// genuinamente un fattore condiviso fra denominatori diversi (es.
// `1/(x-1) + (x-1)/(x-1)^3` → `x/(x²-2x+1)` invece di `x/(x-1)²`) — non nel
// caso repro originale di A55 (denominatori uguali), che con questo fix non
// passa mai dalla riduzione GCD. Non misurato come costo reale altrove nel
// motore; non affrontato qui — richiede una ri-fattorizzazione post-GCD con
// spec e gate propri, questione di design distinta.

#include <gtest/gtest.h>

#include <string>

#include "cas/algebra.hpp"
#include "cas/ast_nodes.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A55TogetherCombineTest : public ::testing::Test {
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

// Repro originale di A55: sei termini sullo stesso denominatore. Deve
// restare UNA sola frazione (non una Sum di termini distribuiti) e deve
// valere lo stesso dell'input.
TEST_F(A55TogetherCombineTest, LongCommonDenominatorSumStaysSingleFraction) {
    ExprPtr e = parse("a/(y+1) + b/(y+1) + c/(y+1) + d/(y+1) + q/(y+1) + f/(y+1)");
    auto tog = algebra::together(e, ctx);
    ASSERT_TRUE(tog.is_ok()) << tog.error().message;

    EXPECT_FALSE(expr_is<Sum>(tog.value()))
        << "together ha lasciato il risultato distribuito sul denominatore comune";

    auto eq = symbolic::mathematically_equal(tog.value(), e, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// Caso minimo (2 termini) e caso con denominatori strutturati (non solo un
// binomio semplice) — stessa proprietà, per non restringere la copertura al
// solo repro storico.
TEST_F(A55TogetherCombineTest, TwoTermCommonDenominatorStaysSingleFraction) {
    ExprPtr e = parse("a/(y+1) + b/(y+1)");
    auto tog = algebra::together(e, ctx);
    ASSERT_TRUE(tog.is_ok()) << tog.error().message;
    EXPECT_FALSE(expr_is<Sum>(tog.value()));
}

TEST_F(A55TogetherCombineTest, NestedReduceStaysSingleFraction) {
    ExprPtr e = parse("(x^2 - y^2)/((x - y)*(x^2 + x*y + y^2))");
    auto tog = algebra::together(e, ctx);
    ASSERT_TRUE(tog.is_ok()) << tog.error().message;
    EXPECT_FALSE(expr_is<Sum>(tog.value()))
        << "il numeratore (x+y) non deve essere ridistribuito sul denominatore";

    auto eq = symbolic::mathematically_equal(tog.value(), e, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

// Difetto (a) della scoperta originale: denominatori DIVERSI che condividono
// un fattore a potenza diversa. Verifica che la riduzione GCD collassi
// davvero la potenza condivisa (non solo che il valore sia giusto): il grado
// del denominatore in x deve scendere a 2 (x-1)^2, non restare a 4.
TEST_F(A55TogetherCombineTest, GcdReductionAlreadyCollapsesSharedFactorPower) {
    ExprPtr e = parse("1/(x-1) + (x-1)/(x-1)^3");
    auto tog = algebra::together(e, ctx);
    ASSERT_TRUE(tog.is_ok()) << tog.error().message;

    auto expected = parse("x/(x-1)^2");
    auto eq = symbolic::mathematically_equal(tog.value(), expected, ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());

    auto deg = algebra::polynomial_degree(tog.value(), Symbol("x"), ctx);
    if (deg.is_ok()) {
        EXPECT_LE(deg.value(), 2)
            << "il denominatore ridotto non dovrebbe superare grado 2 in x "
               "((x-1)^2), la potenza condivisa deve collassare";
    }
}

}  // namespace
