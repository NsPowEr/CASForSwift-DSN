// A54 — `expand` deve garantire la forma espansa, qualunque sia il TIPO dei
// fattori.
//
// Difetto chiuso qui: `expand` costruisce ogni Sum tramite `simplify`, che
// contiene la raccolta simbolica F1.4 (`a·u + b·u → (a+b)·u`) — l'inversa
// esatta della distribuzione.  Il guard di F1.4 salta solo il caso puramente
// polinomiale, quindi con almeno una base condivisa non-Symbol la raccolta
// scattava e `expand` restituiva un Product contenente un Sum:
//
//   expand(x·(ln x + 1)·z)         → x·z·ln x + x·z          (espanso)
//   expand(x·(ln x + 1)·ln(ln x))  → x·(ln x + 1)·ln(ln x)   (NON espanso)
//
// Conseguenza a valle: una differenza identicamente nulla come
// `x·(L+1)·M − x·L·M − x·M` non si annullava più, perché il termine che
// contiene il Sum ha chiave-monomio opaca per il collettore di Step 4.  Ogni
// verifica costruita su `expand` (fra cui `D(F) = f` dopo un'integrazione)
// dava quindi un FALSO NEGATIVO: antiderivate corrette dichiarate sbagliate.
// Il caso d'origine — `∫(1/x + 1/(x·ln x))·ln(ln x)` — è verificato qui
// end-to-end nelle due forme equivalenti dell'integranda.
//
// Verifica SOUND: mai toString(); cancellazione a zero, post-condizione
// strutturale, e D(∫f) − f ≡ 0.

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>

namespace cas::calculus {
namespace {

ExprPtr parse_expr(const std::string& s, AstArena& arena) {
    auto t = Lexer(s).tokenize();
    EXPECT_TRUE(t.is_ok()) << s << ": " << t.error().message;
    Parser p(t.value(), arena);
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << s << ": " << r.error().message;
    return r.value();
}

[[nodiscard]] bool is_zero_expr(ExprPtr e, symbolic::CASContext& ctx) {
    auto tog = algebra::together(e, ctx);
    ExprPtr x = tog.is_ok() ? tog.value() : e;
    auto s = ctx.simplify(x);
    ExprPtr z = s.is_ok() ? s.value() : x;
    if (const auto* il = expr_cast<IntegerLit>(z)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(z)) return rl->numerator.is_zero();
    auto parts = algebra::apart_num_den(z, ctx);
    if (parts.is_error()) return false;
    ExprPtr num = parts.value().numerator;
    if (auto ex = algebra::expand(num, ctx); ex.is_ok()) num = ex.value();
    if (auto sn = ctx.simplify(num); sn.is_ok()) num = sn.value();
    if (const auto* il = expr_cast<IntegerLit>(num)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(num)) return rl->numerator.is_zero();
    return false;
}

// Post-condizione di expand: nessun Sum può restare sotto un Product.
[[nodiscard]] bool has_sum_under_product(ExprPtr e, bool inside_product = false) {
    if (!e) return false;
    if (inside_product && expr_is<Sum>(e)) return true;
    if (const auto* product = expr_cast<Product>(e)) {
        for (ExprPtr f : product->factors) {
            if (has_sum_under_product(f, true)) return true;
        }
        return false;
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        for (ExprPtr t : sum->terms) {
            if (has_sum_under_product(t, false)) return true;
        }
        return false;
    }
    if (const auto* unary = expr_cast<Unary>(e)) {
        return has_sum_under_product(unary->operand, inside_product);
    }
    return false;
}

class A54FormDependence : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    AstArena& arena() { return ctx.arena(); }
};

// Una differenza identicamente nulla deve annullarsi in tutti e tre i casi: il
// numero di fattori trascendenti nel Product non può cambiare il verdetto.
TEST_F(A54FormDependence, ExpandCancelsNullDifference) {
    struct Case { const char* label; const char* src; };
    const Case cases[] = {
        {"simboli nudi",     "x*(y+1)*z - x*y*z - x*z"},
        {"un trascendente",  "x*(ln(x)+1)*z - x*ln(x)*z - x*z"},
        {"due trascendenti", "x*(ln(x)+1)*ln(ln(x)) - x*ln(x)*ln(ln(x)) - x*ln(ln(x))"},
    };
    for (const auto& c : cases) {
        ExprPtr e = parse_expr(c.src, arena());
        auto expanded = algebra::expand(e, ctx);
        ASSERT_TRUE(expanded.is_ok()) << c.label << ": " << expanded.error().message;
        auto simplified = ctx.simplify(expanded.value());
        ASSERT_TRUE(simplified.is_ok()) << c.label << ": " << simplified.error().message;
        const auto* il = expr_cast<IntegerLit>(simplified.value());
        EXPECT_TRUE(il != nullptr && il->value.is_zero())
            << c.label << ": " << c.src << " non si annulla";
    }
}

// La post-condizione strutturale, indipendente dalla forma dei fattori.
TEST_F(A54FormDependence, ExpandLeavesNoSumUnderProduct) {
    const char* srcs[] = {
        "x*(y+1)",
        "x*(ln(x)+1)",
        "x*(ln(x)+1)*z",
        "x*(ln(x)+1)*ln(ln(x))",
        "(ln(x)+1)*(ln(ln(x))+2)*exp(x)",
    };
    for (const char* src : srcs) {
        ExprPtr e = parse_expr(src, arena());
        auto expanded = algebra::expand(e, ctx);
        ASSERT_TRUE(expanded.is_ok()) << src << ": " << expanded.error().message;
        EXPECT_FALSE(has_sum_under_product(expanded.value()))
            << src << ": expand ha lasciato un Sum sotto un Product";
    }
}

// End-to-end: la stessa integranda scritta come somma di frazioni o come
// frazione unica deve dare lo stesso verdetto.  L'antiderivata attesa è
// ln(ln x)²/2 + ln x·ln(ln x) − ln x (certificata numericamente, mpmath 40
// cifre, su 1.3 ≤ x ≤ 100).
TEST_F(A54FormDependence, IntegrateVerdictDoesNotDependOnTheForm) {
    const char* forms[] = {
        "(1/x + 1/(x*ln(x)))*ln(ln(x))",
        "((ln(x) + 1)/(x*ln(x)))*ln(ln(x))",
    };
    for (const char* src : forms) {
        Symbol x("x");
        ExprPtr f = parse_expr(src, arena());
        auto F = integrate(f, x, ctx);
        ASSERT_TRUE(F.is_ok()) << src << ": " << F.error().message;
        auto dF = diff(F.value(), x, 1U, ctx);
        ASSERT_TRUE(dF.is_ok()) << src << ": " << dF.error().message;
        ExprPtr delta = arena().make<Binary>(BinaryOp::Sub, dF.value(), f);
        EXPECT_TRUE(is_zero_expr(delta, ctx))
            << "D(∫f) ≠ f — antiderivata corretta dichiarata sbagliata su: " << src;
    }
}

}  // namespace
}  // namespace cas::calculus
