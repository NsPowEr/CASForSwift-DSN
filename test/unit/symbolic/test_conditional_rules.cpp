// A31 fase 2 — conditional domain rules behind
// CASContext::set_conditional_domain_rules (Domain_Conditions_Propagation.md
// §10, test plan §10.6). Default flag OFF: simplify() output must stay
// bit-per-bit identical to fase 1 (§10.7.2); flag ON: each B.1 rule rewrites
// AND registers the exact conditions of §10.3.

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

[[nodiscard]] bool has_condition(const SideConditionSet& set, DomainConditionKind kind,
                                 ExprPtr subject) {
    return set.contains(DomainCondition{kind, subject});
}

// Structural equality of the two inputs' canonical forms under `ctx`.
[[nodiscard]] bool simplify_to_equal(ExprPtr a, ExprPtr b, CASContext& ctx) {
    auto ra = ctx.simplify(a);
    auto rb = ctx.simplify(b);
    EXPECT_TRUE(ra.is_ok());
    EXPECT_TRUE(rb.is_ok());
    if (ra.is_error() || rb.is_error()) return false;
    return ExprEqual{}(ra.value(), rb.value());
}

// ── §10.6.1 Flag OFF: output invariato, set vuoto ──────────────────────────

TEST(ConditionalRulesTest, FlagOff_AllB1InputsStayStructural) {
    const char* inputs[] = {"exp(log(x))", "log(x^n)",   "log(x*y)",
                            "log(x/y)",    "abs(x^2)",   "0^x",
                            "sqrt(x^2)"};
    for (const char* input : inputs) {
        CASContext ctx;  // conditional_domain_rules() defaults to false.
        ASSERT_FALSE(ctx.conditional_domain_rules());
        ExprPtr expr = parse_expr(input, ctx);
        auto result = ctx.simplify(expr);
        ASSERT_TRUE(result.is_ok()) << input;
        // No B.1 rewrite: the result keeps the refused top-level shape
        // (FuncCall for exp/log/abs/sqrt, Pow for 0^x) instead of collapsing
        // to the conditional form.
        const bool kept_funcall = expr_is<FuncCall>(result.value());
        const bool kept_pow = expr_is<Binary>(result.value());
        EXPECT_TRUE(kept_funcall || kept_pow) << input;
        EXPECT_TRUE(ctx.last_side_conditions().empty()) << input;
    }
}

// ── §10.6.2 Flag ON: riscrittura + set esatto per regola ───────────────────

TEST(ConditionalRulesTest, R1_ExpLog_RewritesAndEmitsNonZero) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("exp(log(x))", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(ExprEqual{}(result.value(), xexpr));
    EXPECT_EQ(ctx.last_side_conditions().size(), 1U);
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::NonZero, xexpr));
}

TEST(ConditionalRulesTest, R1_EPowLn_RewritesAndEmitsNonZero) {
    // The Binary Pow(E, ln x) site in simplify_arithmetic_power.cpp — the
    // second §10.3.R1 site, distinct from FuncCall(Exp).
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr xexpr = parse_expr("x", ctx);
    ExprPtr ln_x = ctx.arena().make<FuncCall>(BuiltinOp::Ln,
        std::vector<ExprPtr>{xexpr});
    ExprPtr pow = ctx.arena().make<Binary>(BinaryOp::Pow,
        ctx.arena().make<Constant>(MathConstant::E), ln_x);

    auto result = ctx.simplify(pow);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(ExprEqual{}(result.value(), xexpr));
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::NonZero, xexpr));
}

TEST(ConditionalRulesTest, R2a_LogPower_RewritesAndEmitsPositiveAndReal) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("log(x^n)", ctx);
    ExprPtr expected = parse_expr("n * ln(x)", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);
    ExprPtr nexpr = parse_expr("n", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    // Snapshot BEFORE simplify_to_equal: its own top-level simplify() calls
    // reset the accumulator (§4.1 / §10.1).
    const SideConditionSet conds = ctx.last_side_conditions();
    EXPECT_TRUE(simplify_to_equal(result.value(), expected, ctx));
    EXPECT_TRUE(has_condition(conds,
                              DomainConditionKind::Positive, xexpr));
    EXPECT_TRUE(has_condition(conds,
                              DomainConditionKind::Real, nexpr));
}

TEST(ConditionalRulesTest, R2b_LogProduct_RewritesAndEmitsPositivePerFactor) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("log(x*y)", ctx);
    ExprPtr expected = parse_expr("ln(x) + ln(y)", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);
    ExprPtr yexpr = parse_expr("y", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    // Snapshot BEFORE simplify_to_equal: its own top-level simplify() calls
    // reset the accumulator (§4.1 / §10.1).
    const SideConditionSet conds = ctx.last_side_conditions();
    EXPECT_TRUE(simplify_to_equal(result.value(), expected, ctx));
    EXPECT_TRUE(has_condition(conds,
                              DomainConditionKind::Positive, xexpr));
    EXPECT_TRUE(has_condition(conds,
                              DomainConditionKind::Positive, yexpr));
}

TEST(ConditionalRulesTest, R2c_LogQuotient_RewritesAndEmitsPositiveBothSides) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("log(x/y)", ctx);
    ExprPtr expected = parse_expr("ln(x) - ln(y)", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);
    ExprPtr yexpr = parse_expr("y", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    // Snapshot BEFORE simplify_to_equal: its own top-level simplify() calls
    // reset the accumulator (§4.1 / §10.1).
    const SideConditionSet conds = ctx.last_side_conditions();
    EXPECT_TRUE(simplify_to_equal(result.value(), expected, ctx));
    EXPECT_TRUE(has_condition(conds,
                              DomainConditionKind::Positive, xexpr));
    EXPECT_TRUE(has_condition(conds,
                              DomainConditionKind::Positive, yexpr));
}

TEST(ConditionalRulesTest, R3_AbsEvenPower_RewritesAndEmitsReal) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("abs(x^2)", ctx);
    ExprPtr expected = parse_expr("x^2", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    // Snapshot BEFORE simplify_to_equal: its own top-level simplify() calls
    // reset the accumulator (§4.1 / §10.1).
    const SideConditionSet conds = ctx.last_side_conditions();
    EXPECT_TRUE(simplify_to_equal(result.value(), expected, ctx));
    EXPECT_TRUE(has_condition(conds,
                              DomainConditionKind::Real, xexpr));
}

TEST(ConditionalRulesTest, R4_ZeroPowSymbolic_RewritesAndEmitsPositive) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("0^x", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    const auto* lit = expr_cast<IntegerLit>(result.value());
    ASSERT_NE(lit, nullptr);
    EXPECT_TRUE(lit->value.is_zero());
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::Positive, xexpr));
}

TEST(ConditionalRulesTest, R4_ZeroPowNegativeLiteral_StillRefused) {
    // The flag must NOT swallow 0^(-1): literal exponents keep the exact
    // pre-existing semantics (division by zero / structural refusal).
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("0^(-1)", ctx);
    auto result = ctx.simplify(expr);
    if (result.is_ok()) {
        const auto* lit = expr_cast<IntegerLit>(result.value());
        EXPECT_TRUE(lit == nullptr || !lit->value.is_zero());
    }
    // Either a structured error or a non-zero structural form is acceptable;
    // silently returning 0 is not.
}

TEST(ConditionalRulesTest, R5_SqrtSquare_RewritesToAbsAndEmitsReal) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    ExprPtr expr = parse_expr("sqrt(x^2)", ctx);
    ExprPtr expected = parse_expr("abs(x)", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    // Snapshot BEFORE simplify_to_equal: its own top-level simplify() calls
    // reset the accumulator (§4.1 / §10.1).
    const SideConditionSet conds = ctx.last_side_conditions();
    EXPECT_TRUE(simplify_to_equal(result.value(), expected, ctx));
    EXPECT_TRUE(has_condition(conds,
                              DomainConditionKind::Real, xexpr));
}

// ── §10.6.3 Controprova §3.3: condizioni provate non emesse ────────────────

TEST(ConditionalRulesTest, ProvenPositive_R1EmitsNothing) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    Symbol x{"x"};
    ctx.assumptions().assume_positive(x);
    ExprPtr expr = parse_expr("exp(log(x))", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    EXPECT_TRUE(ExprEqual{}(result.value(), xexpr));
    EXPECT_TRUE(ctx.last_side_conditions().empty());
}

TEST(ConditionalRulesTest, ProvenPositiveBase_R2aEmitsOnlyRealExponent) {
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    Symbol x{"x"};
    ctx.assumptions().assume_positive(x);
    ExprPtr expr = parse_expr("log(x^n)", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);
    ExprPtr nexpr = parse_expr("n", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    EXPECT_FALSE(has_condition(ctx.last_side_conditions(),
                               DomainConditionKind::Positive, xexpr));
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::Real, nexpr));
}

// ── §10.6.4 Guardia anti-contraddizione ────────────────────────────────────

TEST(ConditionalRulesTest, ProvenNegativeBase_R2aRefuses) {
    // Base provably negative: emitting Positive(-y) would be a vacuous
    // (empty-domain) condition — the refusal must be preserved.
    CASContext ctx;
    ctx.set_conditional_domain_rules(true);
    Symbol y{"y"};
    ctx.assumptions().assume_positive(y);
    ExprPtr expr = parse_expr("log((-y)^n)", ctx);

    auto result = ctx.simplify(expr);
    ASSERT_TRUE(result.is_ok());
    // No Positive(-y) may ever be registered.
    ExprPtr neg_y = parse_expr("-y", ctx);
    EXPECT_FALSE(has_condition(ctx.last_side_conditions(),
                               DomainConditionKind::Positive, neg_y));
}

// ── §10.6.5 Fix §10.1: cache-hit top-level azzera l'accumulatore ───────────

TEST(ConditionalRulesTest, CacheHitAfterOtherTopLevel_NoCrossPollution) {
    CASContext ctx;
    ExprPtr a_over_a = parse_expr("a / a", ctx);
    ExprPtr b_over_b = parse_expr("b / b", ctx);
    ExprPtr aexpr = parse_expr("a", ctx);
    ExprPtr bexpr = parse_expr("b", ctx);

    auto r1 = ctx.simplify(a_over_a);  // miss: {NonZero a}
    ASSERT_TRUE(r1.is_ok());
    auto r2 = ctx.simplify(b_over_b);  // miss: {NonZero b}
    ASSERT_TRUE(r2.is_ok());
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::NonZero, bexpr));

    // Cache HIT for a/a: the accumulator must contain EXACTLY this call's
    // conditions — pre-fix it merged into b/b's leftover set and reported
    // {NonZero a, NonZero b} (spec §10.1).
    auto r3 = ctx.simplify(a_over_a);
    ASSERT_TRUE(r3.is_ok());
    EXPECT_EQ(ctx.last_side_conditions().size(), 1U);
    EXPECT_TRUE(has_condition(ctx.last_side_conditions(),
                              DomainConditionKind::NonZero, aexpr));
    EXPECT_FALSE(has_condition(ctx.last_side_conditions(),
                               DomainConditionKind::NonZero, bexpr));
}

// ── §10.6.6 simplify_tracked ───────────────────────────────────────────────

TEST(ConditionalRulesTest, SimplifyTracked_MissAndHitReportSameConditions) {
    CASContext ctx;
    ExprPtr expr = parse_expr("x / x", ctx);
    ExprPtr xexpr = parse_expr("x", ctx);

    auto first = ctx.simplify_tracked(expr);  // miss
    ASSERT_TRUE(first.is_ok());
    const auto* one = expr_cast<IntegerLit>(first.value().expr);
    ASSERT_NE(one, nullptr);
    EXPECT_EQ(one->value, BigInt(1));
    EXPECT_TRUE(has_condition(first.value().conditions,
                              DomainConditionKind::NonZero, xexpr));

    auto second = ctx.simplify_tracked(expr);  // cache hit
    ASSERT_TRUE(second.is_ok());
    EXPECT_TRUE(ExprEqual{}(second.value().expr, first.value().expr));
    EXPECT_EQ(second.value().conditions.size(), first.value().conditions.size());
    EXPECT_TRUE(has_condition(second.value().conditions,
                              DomainConditionKind::NonZero, xexpr));
}

// ── §10.6.7 Certificato numerico (equivalenza matematica, non toString) ────

TEST(ConditionalRulesTest, NumericCertificate_R3R5_MultiPoint) {
    // For each real sample point p, input(p) and output(p) must simplify to
    // the same canonical literal — the emitted condition Real(x) holds at
    // every p, so equality is required by the §3.1 contract.
    const char* points[] = {"-3", "1/2", "7"};
    struct Case { const char* input; const char* output; };
    const Case cases[] = {{"abs(x^2)", "x^2"}, {"sqrt(x^2)", "abs(x)"}};
    for (const auto& c : cases) {
        for (const char* p : points) {
            CASContext ctx;
            ctx.set_conditional_domain_rules(true);
            Symbol x{"x"};
            ExprPtr input = parse_expr(c.input, ctx);
            ExprPtr output = parse_expr(c.output, ctx);
            ExprPtr point = parse_expr(p, ctx);
            auto in_sub = ctx.substitute(input, x, point);
            auto out_sub = ctx.substitute(output, x, point);
            ASSERT_TRUE(in_sub.is_ok() && out_sub.is_ok()) << c.input << " @ " << p;
            EXPECT_TRUE(simplify_to_equal(in_sub.value(), out_sub.value(), ctx))
                << c.input << " @ " << p;
        }
    }
}

}  // namespace
}  // namespace cas::symbolic
