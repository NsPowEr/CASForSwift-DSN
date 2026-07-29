// A48 — `simplify` deve raggiungere il punto fisso su un Sum che una
// distribuzione ha annidato sotto un Product.
//
// Invariante violato (nessuna spec formale: e' un bug di correttezza su codice
// esistente, non una feature mancante): `simplify(A - A)` deve valere 0 per
// OGNI A. Falliva per ogni A = (polinomio in x)·(funzione trascendente), con
// le due A strutturalmente IDENTICHE — quindi non era un problema di forme
// diverse ma di punto fisso mancato.
//
// RESIDUO DICHIARATO (A52) — confine misurato, non supposto.
// Restano fuori i casi in cui il polinomio ha un monomio a coefficiente
// UNITARIO e c'e' UN SOLO fattore trascendente:
//     (x-2)·e^{2x}   (x+1)·sin x   (x^2-2)·cos x   (x^2+x-1)·cos x   -> BUG
//     (2x-2)·ln x    (4x-6)·ln x   (2x^2-3x)·cos x                   -> OK
//     (x-1)·e^x·sin x  -> OK  (coefficiente unitario ma DUE trascendenti)
// Causa radice diversa da quelle chiuse qui: i termini superstiti stanno al
// livello del Sum, quindi competono al collettore primario (Step 4), non alla
// raccolta dei coefficienti introdotta in A48. Fallivano anche PRIMA di A48 —
// non e' una regressione. Tre tentativi di diagnosi falliti -> fermato per
// protocollo anti-loop invece di patchare alla cieca. Rientrano qui con A52.
//
// Due cause distinte, entrambe necessarie:
//   1. `try_merge_symbolic_like_terms` costruiva il Sum dei coefficienti
//      grezzo (splice + sort, nessuna raccolta): (-2 + 2) e (-2x + 2x)
//      restavano affiancati.
//   2. `decompose_term` non ricorreva sul ramo Neg: `-(2·x)` dava
//      coeff=-1/fattori=[Product(2,x)] mentre `2·x` dava coeff=2/fattori=[x],
//      due chiavi diverse per termini opposti.
// Fissare solo la (1) lasciava `(-(2x) + 2x)·e^{2x}` — misurato.

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A48SumFixpointTest : public ::testing::Test {
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

    [[nodiscard]] bool simplifies_to_zero(const std::string& src) {
        ExprPtr a = ctx.simplify(parse(src)).value();
        auto d = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Sub, a, a));
        if (!d.is_ok()) return false;
        const auto* lit = expr_cast<IntegerLit>(d.value());
        return lit != nullptr && lit->value.is_zero();
    }
};

// La tabella del repro, misurata prima del fix: le prime due righe passavano
// gia' (polinomio puro), tutte le altre no.
TEST_F(A48SumFixpointTest, SelfDifferenceIsZeroForPolynomialTimesTranscendental) {
    for (const char* src : {
             "(2*x-2)",              // polinomio puro — passava gia'
             "(2*x-2)*x",            // polinomio puro — passava gia'
             "(2*x-2)*exp(2*x)",     // il repro originale
             "(2*x-2)*exp(x)",       // non e' l'argomento dell'esponenziale
             "(2*x+2)*exp(2*x)",     // non e' il segno (nessun Neg nel Sum)
             "(3*x-2)*exp(2*x)",     // non e' il coefficiente
             "(2*x-2)*sin(x)",       // non e' `exp`: qualunque trascendente
             "(2*x^2-3*x)*cos(x)",   // grado > 1
             "(2*x-2)*ln(x)",        // logaritmo
             "(4*x-6)*ln(x)",        // coefficienti non unitari, entrambi pari
             "(x-1)*exp(x)*sin(x)",  // due fattori trascendenti
         }) {
        EXPECT_TRUE(simplifies_to_zero(src)) << "simplify(A - A) != 0 per A = " << src;
    }
}

// Il fix non deve trasformare il confronto in un timbro: differenze GENUINE
// devono restare diverse da zero.
TEST_F(A48SumFixpointTest, GenuineDifferencesDoNotCollapse) {
    const std::vector<std::pair<std::string, std::string>> pairs{
        {"(2*x-2)*exp(2*x)", "(2*x-3)*exp(2*x)"},
        {"(2*x-2)*exp(2*x)", "(2*x-2)*exp(3*x)"},
        {"(2*x-2)*sin(x)", "(2*x-2)*cos(x)"},
        {"x*exp(x)", "exp(x)"},
    };
    for (const auto& [lhs, rhs] : pairs) {
        ExprPtr a = ctx.simplify(parse(lhs)).value();
        ExprPtr b = ctx.simplify(parse(rhs)).value();
        auto d = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Sub, a, b));
        ASSERT_TRUE(d.is_ok()) << lhs << " - " << rhs;
        const auto* lit = expr_cast<IntegerLit>(d.value());
        EXPECT_FALSE(lit != nullptr && lit->value.is_zero())
            << lhs << " - " << rhs << " non deve annullarsi";
    }
}

// La raccolta introdotta non deve rompere l'idempotenza: e' esattamente la
// proprieta' che il commento al sito chiamante proteggeva vietando di
// ri-simplificare il prodotto fuso.
TEST_F(A48SumFixpointTest, MergeStaysIdempotent) {
    for (const char* src : {
             "(2*x-2)*exp(2*x) + (x+1)*exp(2*x)",
             "3*exp(x) + 5*exp(x)",
             "x^3 + x",                       // guardia residuo: non fattorizzare
             "(a+b)*sin(x) + (a-b)*sin(x)",   // coefficienti simbolici
         }) {
        auto once = ctx.simplify(parse(src));
        ASSERT_TRUE(once.is_ok()) << src;
        auto twice = ctx.simplify(once.value());
        ASSERT_TRUE(twice.is_ok()) << src;
        EXPECT_TRUE(structural_equal(once.value(), twice.value()))
            << src << ": simplify deve essere idempotente";
    }
}

// Il caso che ha fatto emergere A48: due antiderivate equivalenti la cui
// differenza incrociata ha esattamente questa forma (A43).
TEST_F(A48SumFixpointTest, CrossMultipliedAntiderivativesCompareEqual) {
    auto eq = symbolic::mathematically_equal(
        parse("2*(2*x-2)^-1*exp(2*x)"), parse("(x-1)^-1*exp(2*x)"), ctx);
    ASSERT_TRUE(eq.is_ok());
    EXPECT_TRUE(eq.value());
}

}  // namespace
