// A53 fase 1 — rollback delle side-conditions (A31) emesse da un tentativo
// abbandonato. Il meccanismo: CASContext::SideConditionRollback.
//
// Motivo: dare un budget d'operazione a calculus::integrate (OperationScope)
// ferma l'azzeramento che le simplify top-level facevano dell'accumulatore,
// quindi le condizioni di un ramo scartato sopravvivevano nel risultato di un
// altro — misurato su `∫e^{-x²}`, che usciva con `x>0` emessa dal fallback
// Meijer/Mellin di una sotto-integrazione poi scartata.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Domain_Conditions_Propagation.md

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/side_conditions.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>

namespace cas::symbolic {
namespace {

[[nodiscard]] ExprPtr parse_expr(const std::string& input, CASContext& ctx) {
    auto tokens = Lexer(input).tokenize();
    EXPECT_TRUE(tokens.is_ok()) << input;
    Parser parser(tokens.value(), ctx.arena());
    auto result = parser.parse();
    EXPECT_TRUE(result.is_ok()) << input;
    return result.value();
}

[[nodiscard]] bool has(const CASContext& ctx, DomainConditionKind kind, ExprPtr subject) {
    return ctx.last_side_conditions().contains(DomainCondition{kind, subject});
}

// Uno scope non committato annulla le condizioni emesse al suo interno e lascia
// intatte quelle che c'erano prima. I due emittenti sono DISTINTI (soggetti
// diversi): con un solo emittente un rollback che azzerasse tutto passerebbe
// comunque il test — la trappola gia' vista in A31 fase 2.
TEST(A53ConditionRollbackTest, UncommittedScopeDropsOnlyItsOwnConditions) {
    CASContext ctx;
    ExprPtr a = parse_expr("a", ctx);
    ExprPtr b = parse_expr("b", ctx);

    ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::NonZero, a).is_ok());
    {
        CASContext::SideConditionRollback guard(ctx);
        ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::Positive, b).is_ok());
        EXPECT_TRUE(has(ctx, DomainConditionKind::Positive, b));
    }
    EXPECT_TRUE(has(ctx, DomainConditionKind::NonZero, a)) << "condizione precedente persa";
    EXPECT_FALSE(has(ctx, DomainConditionKind::Positive, b)) << "condizione del ramo morto rimasta";
    EXPECT_EQ(ctx.last_side_conditions().size(), 1U);
}

TEST(A53ConditionRollbackTest, CommittedScopeKeepsItsConditions) {
    CASContext ctx;
    ExprPtr a = parse_expr("a", ctx);
    ExprPtr b = parse_expr("b", ctx);

    ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::NonZero, a).is_ok());
    {
        CASContext::SideConditionRollback guard(ctx);
        ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::Positive, b).is_ok());
        guard.commit();
    }
    EXPECT_TRUE(has(ctx, DomainConditionKind::NonZero, a));
    EXPECT_TRUE(has(ctx, DomainConditionKind::Positive, b));
    EXPECT_EQ(ctx.last_side_conditions().size(), 2U);
}

// Il mark e' una COPIA, non una dimensione: `add` rimuove la condizione piu'
// debole quando ne entra una che la subsume (Positive subsume NonZero sullo
// stesso soggetto, spec §3.4). Ripristinare per troncamento lascerebbe qui un
// set di dimensione 1 ma con la condizione SBAGLIATA dentro.
TEST(A53ConditionRollbackTest, RestoresConditionSubsumedInsideTheScope) {
    CASContext ctx;
    ExprPtr a = parse_expr("a", ctx);

    ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::NonZero, a).is_ok());
    {
        CASContext::SideConditionRollback guard(ctx);
        ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::Positive, a).is_ok());
        EXPECT_FALSE(has(ctx, DomainConditionKind::NonZero, a)) << "subsunta, come da §3.4";
    }
    EXPECT_TRUE(has(ctx, DomainConditionKind::NonZero, a)) << "la condizione subsunta non e' tornata";
    EXPECT_FALSE(has(ctx, DomainConditionKind::Positive, a));
    EXPECT_EQ(ctx.last_side_conditions().size(), 1U);
}

// Gli scope si annidano: il rollback esterno riporta allo stato del proprio
// mark anche se quello interno aveva committato.
TEST(A53ConditionRollbackTest, NestedScopesRollBackToTheirOwnMark) {
    CASContext ctx;
    ExprPtr a = parse_expr("a", ctx);
    ExprPtr b = parse_expr("b", ctx);
    ExprPtr c = parse_expr("c", ctx);

    ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::NonZero, a).is_ok());
    {
        CASContext::SideConditionRollback outer(ctx);
        ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::NonZero, b).is_ok());
        {
            CASContext::SideConditionRollback inner(ctx);
            ASSERT_TRUE(ctx.emit_side_condition(DomainConditionKind::Positive, c).is_ok());
            inner.commit();
        }
        EXPECT_TRUE(has(ctx, DomainConditionKind::Positive, c));
    }
    EXPECT_EQ(ctx.last_side_conditions().size(), 1U);
    EXPECT_TRUE(has(ctx, DomainConditionKind::NonZero, a));
}

// Il rollback tocca l'accumulatore, NON la cache di simplify: la' le condizioni
// restano attribuite alla loro entry e vanno ri-emesse a ogni hit (spec §4.2).
// Senza questa separazione un ramo morto "avvelenerebbe" la cache per tutti i
// chiamanti successivi, che e' il difetto opposto e altrettanto silenzioso.
TEST(A53ConditionRollbackTest, RollbackDoesNotPoisonTheSimplifyCache) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("log(x*y)", ctx);

    SideConditionSet from_dead_branch;
    {
        CASContext::SideConditionRollback guard(ctx);
        auto first = ctx.simplify(expr);
        ASSERT_TRUE(first.is_ok());
        from_dead_branch = ctx.last_side_conditions();
        ASSERT_FALSE(from_dead_branch.empty()) << "la regola R3 deve emettere condizioni";
    }
    EXPECT_TRUE(ctx.last_side_conditions().empty());

    auto second = ctx.simplify(expr);  // cache hit
    ASSERT_TRUE(second.is_ok());
    EXPECT_EQ(ctx.last_side_conditions().size(), from_dead_branch.size())
        << "la cache ha perso le condizioni della sua entry";
}

}  // namespace
}  // namespace cas::symbolic
