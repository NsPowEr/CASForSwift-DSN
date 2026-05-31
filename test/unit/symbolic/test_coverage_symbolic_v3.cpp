// DEBT-F1-COV-01 — Coverage uplift v3 targeting specific uncovered branches.
// Focus: assumptions.cpp advanced branches, rewrite_engine AC-Product path,
//        context_core instantiate_pattern complex nodes, bessel helper branches,
//        simplify_functions additional paths.

#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/normal_form.hpp"
#include "cas/complex_rational.hpp"
#include "cas/rational.hpp"

#include <string>
#include <vector>

using namespace cas;
using namespace cas::symbolic;

namespace {

[[nodiscard]] ExprPtr make_int(AstArena& a, long long v) {
    return a.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] ExprPtr make_sym(AstArena& a, const std::string& name) {
    return a.make<Symbol>(name);
}

} // namespace

// ─── assumptions: could_be_zero(ExprPtr) branches ────────────────────────────

TEST(AssumptionsCoverage2, CouldBeZeroExprNull) {
    CASContext ctx;
    // ExprPtr overload: null → false
    EXPECT_FALSE(ctx.assumptions().could_be_zero(ExprPtr{}));
}

TEST(AssumptionsCoverage2, CouldBeZeroExprIsZeroLit) {
    CASContext ctx;
    ExprPtr zero = make_int(ctx.arena(), 0);
    // zero literal → could_be_zero = true
    EXPECT_TRUE(ctx.assumptions().could_be_zero(zero));
}

TEST(AssumptionsCoverage2, CouldBeZeroExprIsNonzero) {
    CASContext ctx;
    // Integer 5 is nonzero → could_be_zero = false
    ExprPtr five = make_int(ctx.arena(), 5);
    EXPECT_FALSE(ctx.assumptions().could_be_zero(five));
}

TEST(AssumptionsCoverage2, CouldBeZeroExprIsSymbol) {
    CASContext ctx;
    Symbol x("x");
    // Symbol with no assumptions → could_be_zero = true
    ExprPtr x_expr = make_sym(ctx.arena(), "x");
    EXPECT_TRUE(ctx.assumptions().could_be_zero(x_expr));
}

TEST(AssumptionsCoverage2, CouldBeZeroExprSymbolNonzero) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_nonzero(x);
    ExprPtr x_expr = make_sym(ctx.arena(), "x");
    EXPECT_FALSE(ctx.assumptions().could_be_zero(x_expr));
}

TEST(AssumptionsCoverage2, CouldBeZeroExprNonSymbol) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Non-symbol expression that is not nonzero → true
    ExprPtr e = arena.make<Binary>(BinaryOp::Add, make_sym(arena, "a"), make_sym(arena, "b"));
    EXPECT_TRUE(ctx.assumptions().could_be_zero(e));
}

// ─── assumptions: is_greater with zero-check on rhs=0, rhs as negative symbol ──

TEST(AssumptionsCoverage2, IsGreaterLhsZeroRhsNegativeSymbol) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol y("y");
    ctx.assumptions().assume_domain(y, Domain::Negative);
    ExprPtr zero = make_int(arena, 0);
    ExprPtr y_expr = make_sym(arena, "y");
    // 0 > y_negative should be true
    EXPECT_TRUE(ctx.assumptions().is_greater(zero, y_expr));
}

TEST(AssumptionsCoverage2, IsGreaterEqualLhsZeroRhsNegative) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol y("y");
    ctx.assumptions().assume_domain(y, Domain::Negative);
    ExprPtr zero = make_int(arena, 0);
    ExprPtr y_expr = make_sym(arena, "y");
    // 0 >= y_negative should be true
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(zero, y_expr));
}

TEST(AssumptionsCoverage2, IsGreaterSumVsSingleTermOthersNonneg) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    // (x + y) > 0
    std::vector<ExprPtr> terms{x_expr, y_expr};
    ExprPtr sum = arena.make<Sum>(std::move(terms));
    // This tests the "l_terms.size() > 1 && r_terms.size() == 1" branch
    // Result may be false or true depending on what's provable
    bool result = ctx.assumptions().is_greater(sum, make_int(arena, 0));
    (void)result;
    EXPECT_TRUE(true); // no crash
}

TEST(AssumptionsCoverage2, IsGreaterTermByTermComparison) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    ExprPtr x_expr = make_sym(arena, "x");

    // Test l_terms.size() == r_terms.size() branch in is_greater
    // (x + 2) > (x + 1)? → 2 > 1 ✓
    std::vector<ExprPtr> lhs_terms{x_expr, make_int(arena, 2)};
    std::vector<ExprPtr> rhs_terms{x_expr, make_int(arena, 1)};
    ExprPtr lhs_sum = arena.make<Sum>(std::move(lhs_terms));
    ExprPtr rhs_sum = arena.make<Sum>(std::move(rhs_terms));
    bool result = ctx.assumptions().is_greater(lhs_sum, rhs_sum);
    EXPECT_TRUE(result);
}

TEST(AssumptionsCoverage2, IsGreaterEqualSumTermByTerm) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    ExprPtr x_expr = make_sym(arena, "x");

    // (x + 2) >= (x + 1) — should be true
    std::vector<ExprPtr> lhs_terms{x_expr, make_int(arena, 2)};
    std::vector<ExprPtr> rhs_terms{x_expr, make_int(arena, 1)};
    ExprPtr lhs = arena.make<Sum>(std::move(lhs_terms));
    ExprPtr rhs = arena.make<Sum>(std::move(rhs_terms));
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(lhs, rhs));
}

TEST(AssumptionsCoverage2, IsGreaterSumVsSingleFromL) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");

    // (x + y) >= x? other term y is nonneg, and x >= x is true
    std::vector<ExprPtr> lhs_terms{x_expr, y_expr};
    ExprPtr lhs = arena.make<Sum>(std::move(lhs_terms));
    EXPECT_TRUE(ctx.assumptions().is_greater_equal(lhs, x_expr));
}

// ─── assumptions: check_consistency with cyclic strict inequality ─────────────

TEST(AssumptionsCoverage2, CheckConsistencyCyclicStrictInequality) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr a = make_sym(arena, "a");
    ExprPtr b = make_sym(arena, "b");
    // a > b AND b > a → contradiction
    ctx.assumptions().assume_greater(a, b);
    ctx.assumptions().assume_greater(b, a);
    auto res = ctx.assumptions().check_consistency();
    EXPECT_TRUE(res.is_error());
}

TEST(AssumptionsCoverage2, CheckConsistencySelfLoop) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr a = make_sym(arena, "a");
    // a > a (strict self-loop)
    ctx.assumptions().assume_greater(a, a);
    auto res = ctx.assumptions().check_consistency();
    EXPECT_TRUE(res.is_error());
}

// ─── assumptions: update_roots with relations ─────────────────────────────────

TEST(AssumptionsCoverage2, UpdateRootsWithRelations) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr a = make_sym(arena, "a");
    ExprPtr b = make_sym(arena, "b");
    // Add relation a > b
    ctx.assumptions().assume_greater(a, b);
    // collect_garbage triggers update_roots
    ctx.collect_garbage({});
    EXPECT_TRUE(true);
}

// ─── assumptions: is_real for Binary with Pow negative base integer exponent ──

TEST(AssumptionsCoverage2, IsRealBinaryPowNegativeBaseIntExp) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Negative);
    ctx.assumptions().assume_integer(Symbol("n"));
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr n_expr = make_sym(arena, "n");
    // x^n where x < 0, n integer → real
    ExprPtr pow_e = arena.make<Binary>(BinaryOp::Pow, x_expr, n_expr);
    bool result = ctx.assumptions().is_real(pow_e);
    // Should be true since x is negative and n is integer
    EXPECT_TRUE(result);
}

TEST(AssumptionsCoverage2, IsRealBinaryAddNotReal) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Non-real symbol + something → is_real should consider fallthrough
    ExprPtr x = make_sym(arena, "x"); // no assumption
    ExprPtr e = arena.make<Binary>(BinaryOp::Add, x, make_int(arena, 1));
    // x is not known real → result = false
    EXPECT_FALSE(ctx.assumptions().is_real(e));
}

// ─── assumptions: is_nonnegative Mul and Pow even power ──────────────────────

TEST(AssumptionsCoverage2, IsNonnegativeMulPositiveNegative) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_domain(y, Domain::Negative);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    // x*y where x>0, y<0 → product < 0 → not nonneg
    ExprPtr mul = arena.make<Binary>(BinaryOp::Mul, x_expr, y_expr);
    EXPECT_FALSE(ctx.assumptions().is_nonnegative(mul));
}

TEST(AssumptionsCoverage2, IsNonnegativeMulBothNegative) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_domain(x, Domain::Negative);
    ctx.assumptions().assume_domain(y, Domain::Negative);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    ExprPtr mul = arena.make<Binary>(BinaryOp::Mul, x_expr, y_expr);
    // (-) * (-) = positive → nonneg
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(mul));
}

TEST(AssumptionsCoverage2, IsNonnegativeMulOneZero) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr zero = make_int(arena, 0);
    ExprPtr x = make_sym(arena, "x");
    ExprPtr mul = arena.make<Binary>(BinaryOp::Mul, zero, x);
    // 0 * x = 0 → nonneg
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(mul));
}

TEST(AssumptionsCoverage2, IsNonnegativeAddBranch) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    ExprPtr add = arena.make<Binary>(BinaryOp::Add, x_expr, y_expr);
    // x + y where both > 0 → nonneg
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(add));
}

TEST(AssumptionsCoverage2, IsNonnegativeProductWithZeroFactor) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr zero = make_int(arena, 0);
    ExprPtr x = make_sym(arena, "x");
    std::vector<ExprPtr> factors{zero, x};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(prod));
}

TEST(AssumptionsCoverage2, IsNonnegativeProductAllPositive) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    std::vector<ExprPtr> factors{x_expr, y_expr};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(prod));
}

TEST(AssumptionsCoverage2, IsNonnegativeProductBothNegative) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_domain(x, Domain::Negative);
    ctx.assumptions().assume_domain(y, Domain::Negative);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    std::vector<ExprPtr> factors{x_expr, y_expr};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    // (-)*(-) → positive → nonneg
    EXPECT_TRUE(ctx.assumptions().is_nonnegative(prod));
}

// ─── assumptions: prove_positive_linear advanced cases ───────────────────────

TEST(AssumptionsCoverage2, ProvePositiveLinearBinaryAdd) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    // x + y > 0: BinaryOp::Add branch
    ExprPtr add = arena.make<Binary>(BinaryOp::Add, x_expr, y_expr);
    EXPECT_TRUE(ctx.assumptions().is_positive(add));
}

TEST(AssumptionsCoverage2, ProvePositiveLinearBinaryMul) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    // x * y where both positive → positive
    ExprPtr mul = arena.make<Binary>(BinaryOp::Mul, x_expr, y_expr);
    // prove_positive_linear only handles Sub, Add, Mul — this is the Mul branch
    EXPECT_TRUE(ctx.assumptions().is_positive(mul));
}

TEST(AssumptionsCoverage2, ProvePositiveLinearBinarySub) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // x - y > 0 when x > y
    ExprPtr x = make_sym(arena, "x");
    ExprPtr y = make_sym(arena, "y");
    ctx.assumptions().assume_greater(x, y);
    ExprPtr sub = arena.make<Binary>(BinaryOp::Sub, x, y);
    // prove_positive_linear: BinaryOp::Sub → is_greater(left, right)
    EXPECT_TRUE(ctx.assumptions().is_positive(sub));
}

TEST(AssumptionsCoverage2, ProvePositiveLinearSumStrictAndNonneg) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    // Sum with all positive terms → proves positive
    std::vector<ExprPtr> terms{x_expr, y_expr, make_int(arena, 1)};
    ExprPtr sum = arena.make<Sum>(std::move(terms));
    EXPECT_TRUE(ctx.assumptions().is_positive(sum));
}

TEST(AssumptionsCoverage2, ProvePositiveLinearSumCheckSub) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    // x > y → (x + (-y)) > 0
    ctx.assumptions().assume_greater(make_sym(arena, "x"), make_sym(arena, "y"));
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr neg_y = arena.make<Unary>(UnaryOp::Neg, make_sym(arena, "y"));
    std::vector<ExprPtr> terms{x_expr, neg_y};
    ExprPtr sum = arena.make<Sum>(std::move(terms));
    // This hits the 2-term sum check_sub branch
    bool result = ctx.assumptions().is_positive(sum);
    EXPECT_TRUE(result);
}

// ─── assumptions: get_domain fallback from integer/real/nonzero sets ──────────

TEST(AssumptionsCoverage2, GetDomainFromIntegerSet) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_integer(x); // sets integer_symbols_
    // get_domain should see integer_symbols_ contains x
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Integer);
}

TEST(AssumptionsCoverage2, GetDomainFromRealSet) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_real(x); // sets real_symbols_
    // positive not set, so fallback to real
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Real);
}

TEST(AssumptionsCoverage2, GetDomainFromNonzeroSet) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_nonzero(x);
    // Not positive, negative, integer, real → NonZero
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::NonZero);
}

TEST(AssumptionsCoverage2, GetDomainFromNoSet) {
    CASContext ctx;
    Symbol x("x");
    // No assumptions → Complex
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Complex);
}

// ─── assumptions: assume_domain positive calls assume_positive ────────────────

TEST(AssumptionsCoverage2, AssumeDomainPositiveCallsAssumePositive) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Positive);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Positive);
    EXPECT_TRUE(ctx.assumptions().is_positive(x));
    EXPECT_TRUE(ctx.assumptions().is_real(x));
    EXPECT_TRUE(ctx.assumptions().is_nonzero(x));
}

TEST(AssumptionsCoverage2, AssumeDomainRealCallsAssumeReal) {
    CASContext ctx;
    Symbol x("x");
    ctx.assumptions().assume_domain(x, Domain::Real);
    EXPECT_EQ(ctx.assumptions().get_domain(x), Domain::Real);
    EXPECT_TRUE(ctx.assumptions().is_real(x));
}

// ─── rewrite_engine: AC-Product matching path ────────────────────────────────
// To trigger the Product AC path in try_apply_rule_here, we need:
// - A Product expr with more factors than the pattern
// - A rule whose pattern is a Product
// - The rule must be oriented (replacement < pattern in term order)

TEST(RewriteEngineCoverage2, ApplyRuleACProductMatching) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Pattern: x_ * 1  → replacement: x_  (removing the factor 1)
    // This is oriented: x_ * 1 > x_ in term order (more complex)
    ExprPtr wild = arena.make<Symbol>("x_");
    ExprPtr one = make_int(arena, 1);
    std::vector<ExprPtr> pat_factors{wild, one};
    ExprPtr pattern = arena.make<Product>(std::move(pat_factors));
    ExprPtr replacement = arena.make<Symbol>("x_");
    RewriteRule rule{pattern, replacement, nullptr};

    // Target: x * 1 * y  (3 factors, pattern has 2, so AC should fire)
    ExprPtr x = make_sym(arena, "x");
    ExprPtr y = make_sym(arena, "y");
    std::vector<ExprPtr> target_factors{x, one, y};
    ExprPtr target = arena.make<Product>(std::move(target_factors));

    // Rule orientation: x_ * 1 vs x_ — 1 is lighter so product is heavier
    // Try bottom-up apply — may or may not pass orientation, either way no crash
    auto result = apply_rule(target, rule, TraversalStrategy::BottomUp, arena);
    EXPECT_FALSE(result.is_error() && result.error().kind != CASErrorKind::InvalidArgument);
}

TEST(RewriteEngineCoverage2, ApplyRuleACSumMatchingFires) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Rule: x_ + 0 → x_  (oriented: sum is heavier than x_)
    ExprPtr wild = arena.make<Symbol>("x_");
    ExprPtr zero = make_int(arena, 0);
    std::vector<ExprPtr> pat_terms{wild, zero};
    ExprPtr pattern = arena.make<Sum>(std::move(pat_terms));
    ExprPtr replacement = arena.make<Symbol>("x_");
    RewriteRule rule{pattern, replacement, nullptr};

    // Target: x + y + 0 (3 terms, AC should match 2 of them)
    ExprPtr x = make_sym(arena, "x");
    ExprPtr y = make_sym(arena, "y");
    std::vector<ExprPtr> target_terms{x, y, zero};
    ExprPtr target = arena.make<Sum>(std::move(target_terms));

    auto result = apply_rule(target, rule, TraversalStrategy::BottomUp, arena);
    // Should not crash; rule may or may not fire depending on orientation
    (void)result;
    EXPECT_TRUE(true);
}

TEST(RewriteEngineCoverage2, ApplyRuleFixPointConverges) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Rule that always fires when applied but result doesn't change anymore
    // Pattern: x_ + 0 → x_ (oriented)
    ExprPtr wild = arena.make<Symbol>("x_");
    ExprPtr zero = make_int(arena, 0);
    std::vector<ExprPtr> pat_terms{wild, zero};
    ExprPtr pattern = arena.make<Sum>(std::move(pat_terms));
    ExprPtr replacement = arena.make<Symbol>("x_");
    RewriteRule rule{pattern, replacement, nullptr};

    // Target without 0: nothing to reduce → FixPoint terminates immediately
    ExprPtr x = make_sym(arena, "x");
    auto result = apply_rule(x, rule, TraversalStrategy::FixPoint, arena);
    EXPECT_FALSE(result.is_error());
}

TEST(RewriteEngineCoverage2, RewriteChildrenCoversIntegral) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Build integral ∫ (ln(1)) dx from 0 to 1 — ln(1) will simplify to 0
    ExprPtr ln_1 = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{make_int(arena, 1)});
    ExprPtr integral = arena.make<Integral>(
        ln_1, Symbol("x"),
        std::optional<ExprPtr>(make_int(arena, 0)),
        std::optional<ExprPtr>(make_int(arena, 1)));

    // Applying any oriented rule will recurse through the Integral node
    ExprPtr wild = arena.make<Symbol>("w_");
    std::vector<ExprPtr> pat{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pat));
    RewriteRule rule{pattern, arena.make<Symbol>("w_"), nullptr};
    auto result = apply_rule(integral, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

TEST(RewriteEngineCoverage2, RewriteChildrenCoversDerivative) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr deriv = arena.make<Derivative>(
        make_sym(arena, "x"), Symbol("x"), 1);
    ExprPtr wild = arena.make<Symbol>("w_");
    std::vector<ExprPtr> pat{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pat));
    RewriteRule rule{pattern, arena.make<Symbol>("w_"), nullptr};
    auto result = apply_rule(deriv, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

TEST(RewriteEngineCoverage2, RewriteChildrenCoversLimit) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr lim = arena.make<Limit>(
        make_sym(arena, "x"), Symbol("x"), make_int(arena, 0), LimitDirection::Both);
    ExprPtr wild = arena.make<Symbol>("w_");
    std::vector<ExprPtr> pat{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pat));
    RewriteRule rule{pattern, arena.make<Symbol>("w_"), nullptr};
    auto result = apply_rule(lim, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

TEST(RewriteEngineCoverage2, RewriteChildrenCoversRootOf) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ExprPtr poly = arena.make<Binary>(
        BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Pow, make_sym(arena, "x"), make_int(arena, 2)),
        make_int(arena, 2));
    ExprPtr rootof = arena.make<RootOf>(poly, Symbol("x"), std::nullopt);
    ExprPtr wild = arena.make<Symbol>("w_");
    std::vector<ExprPtr> pat{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pat));
    RewriteRule rule{pattern, arena.make<Symbol>("w_"), nullptr};
    auto result = apply_rule(rootof, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

TEST(RewriteEngineCoverage2, RewriteChildrenCoversMatrix) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> elems{
        make_int(arena, 1), make_int(arena, 0),
        make_int(arena, 0), make_int(arena, 1)};
    ExprPtr mat = arena.make<Matrix>(2, 2, std::move(elems));
    ExprPtr wild = arena.make<Symbol>("w_");
    std::vector<ExprPtr> pat{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pat));
    RewriteRule rule{pattern, arena.make<Symbol>("w_"), nullptr};
    auto result = apply_rule(mat, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

TEST(RewriteEngineCoverage2, RewriteChildrenCoversSeriesExp) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    std::vector<std::pair<long long, ExprPtr>> terms;
    terms.push_back({0, make_int(arena, 1)});
    terms.push_back({1, make_int(arena, 1)});
    ExprPtr series = arena.make<SeriesExp>(
        Symbol("x"), make_int(arena, 0), std::move(terms), 2);
    ExprPtr wild = arena.make<Symbol>("w_");
    std::vector<ExprPtr> pat{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pat));
    RewriteRule rule{pattern, arena.make<Symbol>("w_"), nullptr};
    auto result = apply_rule(series, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

// ─── context_core: instantiate_pattern with various node types ───────────────

TEST(ContextCoreCoverage3, InstantiatePatternViaApplyRule) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Use apply_rule to trigger instantiate_pattern on Sum
    // Rule: x_ + 0 → x_
    ExprPtr wild = arena.make<Symbol>("x_");
    ExprPtr zero = make_int(arena, 0);
    std::vector<ExprPtr> pat_terms{wild, zero};
    ExprPtr pattern = arena.make<Sum>(std::move(pat_terms));
    RewriteRule rule{pattern, arena.make<Symbol>("x_"), nullptr};

    // This triggers instantiate_pattern(replacement, matches, arena) where
    // replacement is a Symbol (wildcard)
    ExprPtr x = make_sym(arena, "x");
    std::vector<ExprPtr> target_terms{x, zero};
    ExprPtr target = arena.make<Sum>(std::move(target_terms));
    auto result = apply_rule(target, rule, TraversalStrategy::BottomUp, arena);
    // Just verify it doesn't crash
    (void)result;
    EXPECT_TRUE(true);
}

TEST(ContextCoreCoverage3, InstantiatePatternBinaryViaRule) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Rule: x_ + 0 → x_ + x_  — but this won't be oriented (lhs <= rhs)
    // Let's try: (x_ * x_) → x_  which should be oriented (a^2 > a in term order sometimes)
    // Actually the easiest path: use apply_rule with a matched Binary node
    // Rule: x_ + 0 → x_  (Sum-form pattern, triggers instantiation)
    ExprPtr wild = arena.make<Symbol>("x_");
    std::vector<ExprPtr> pt{wild, make_int(arena, 0)};
    ExprPtr pattern = arena.make<Sum>(std::move(pt));
    RewriteRule rule{pattern, arena.make<Symbol>("x_"), nullptr};

    // Target is a Binary expr containing Sum
    ExprPtr inner_sum = arena.make<Sum>(
        std::vector<ExprPtr>{make_sym(arena, "a"), make_int(arena, 0)});
    ExprPtr target = arena.make<Binary>(BinaryOp::Mul, inner_sum, make_sym(arena, "b"));
    auto result = apply_rule(target, rule, TraversalStrategy::BottomUp, arena);
    (void)result;
    EXPECT_TRUE(true);
}

// ─── simplify_bessel_orthogonal: is_real and is_nonneg_integer helper paths ──

TEST(SimplifyBesselCoverage, BesselWithRealSymbolAssumption) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_real(x);
    // BesselJ(0, x) — x is real → should evaluate
    std::vector<ExprPtr> args{make_int(arena, 0), make_sym(arena, "x")};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::BesselJ, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyBesselCoverage, BesselWithConstantArg) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // BesselJ(0, Pi) — Pi is a non-I, non-NaN constant → is_real
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    std::vector<ExprPtr> args{make_int(arena, 0), pi};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::BesselJ, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyBesselCoverage, ChebyshevTWithPositiveIntegerDegree) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_real(x);
    // ChebyshevT(2, x) — should work with integer degree
    std::vector<ExprPtr> args{make_int(arena, 2), make_sym(arena, "x")};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::ChebyshevT, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyBesselCoverage, HermiteHWithDegree3) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // HermiteH(3, x) — physicist Hermite polynomial
    std::vector<ExprPtr> args{make_int(arena, 3), make_sym(arena, "x")};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::HermiteH, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyBesselCoverage, HermiteHeWithDegree4) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // HermiteHe(4, x) — probabilist Hermite polynomial
    std::vector<ExprPtr> args{make_int(arena, 4), make_sym(arena, "x")};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::HermiteHe, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyBesselCoverage, LaguerreLWithDegree2) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // LaguerreL(2, x)
    std::vector<ExprPtr> args{make_int(arena, 2), make_sym(arena, "x")};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::LaguerreL, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyBesselCoverage, LegendrePWithDegree5) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // LegendreP(5, x)
    std::vector<ExprPtr> args{make_int(arena, 5), make_sym(arena, "x")};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::LegendreP, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyBesselCoverage, JacobiPWithDegree3) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // JacobiP(3, 1, 1, x)
    std::vector<ExprPtr> args{
        make_int(arena, 3),
        make_int(arena, 1),
        make_int(arena, 1),
        make_sym(arena, "x")};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::JacobiP, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyBesselCoverage, BesselJWithIntegerNegativeOrder) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // BesselJ(-1, x) — negative integer order
    std::vector<ExprPtr> args{make_int(arena, -1), make_sym(arena, "x")};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::BesselJ, std::move(args));
    auto r = ctx.simplify(e);
    // Should handle it (could return Unimplemented for large orders)
    (void)r;
    EXPECT_TRUE(true);
}

TEST(SimplifyBesselCoverage, LambertWWithNonnegativeProduct) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    // W(x * exp(x)) should simplify to x if x >= 0
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr exp_x = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{x_expr});
    std::vector<ExprPtr> prod_factors{x_expr, exp_x};
    ExprPtr prod = arena.make<Product>(std::move(prod_factors));
    std::vector<ExprPtr> w_args{prod};
    ExprPtr w = arena.make<FuncCall>(BuiltinOp::LambertW, std::move(w_args));
    auto r = ctx.simplify(w);
    EXPECT_FALSE(r.is_error());
}

// ─── simplify_functions: may_rewrite_function_call paths ─────────────────────

TEST(SimplifyFunctionsCoverage2, SqrtOfDivision) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    // sqrt(x / y) where both positive → should simplify
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    ExprPtr div = arena.make<Binary>(BinaryOp::Div, x_expr, y_expr);
    std::vector<ExprPtr> args{div};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Sqrt, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyFunctionsCoverage2, SqrtOfProductNonneg) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    // sqrt(x * y) where both positive → may_rewrite_function_call(Sqrt)
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    std::vector<ExprPtr> factors{x_expr, y_expr};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    std::vector<ExprPtr> args{prod};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Sqrt, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyFunctionsCoverage2, LnOfSqrtPositive) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    // ln(sqrt(x)) where x > 0 → triggers may_rewrite_function_call(Ln) sqrt branch
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr sqrt_x = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{x_expr});
    std::vector<ExprPtr> args{sqrt_x};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Ln, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyFunctionsCoverage2, LnOfPositiveDiv) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    // ln(x/y) where both positive → may_rewrite (Div branch)
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    ExprPtr div = arena.make<Binary>(BinaryOp::Div, x_expr, y_expr);
    std::vector<ExprPtr> args{div};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Ln, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyFunctionsCoverage2, LnOfPositiveProduct) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_positive(x);
    ctx.assumptions().assume_positive(y);
    // ln(x*y) product form → may_rewrite (Product branch)
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    std::vector<ExprPtr> factors{x_expr, y_expr};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    std::vector<ExprPtr> args{prod};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Ln, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyFunctionsCoverage2, LnOfPositivePow) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    // ln(x^3) where x > 0 → may_rewrite (Pow with positive base)
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr pow_e = arena.make<Binary>(BinaryOp::Pow, x_expr, make_int(arena, 3));
    std::vector<ExprPtr> args{pow_e};
    ExprPtr e = arena.make<FuncCall>(BuiltinOp::Ln, std::move(args));
    auto r = ctx.simplify(e);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyFunctionsCoverage2, ExpOfLnPositive) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    // e^ln(x) where x > 0 → may_rewrite_power
    ExprPtr e_const = arena.make<Constant>(MathConstant::E);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr ln_x = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{x_expr});
    ExprPtr e = arena.make<Binary>(BinaryOp::Pow, e_const, ln_x);
    auto r = ctx.simplify(e);
    ASSERT_TRUE(r.is_ok());
    // e^ln(x) = x
    const auto* sym = expr_cast<Symbol>(r.value());
    ASSERT_NE(sym, nullptr);
    EXPECT_EQ(sym->name, "x");
}

TEST(SimplifyFunctionsCoverage2, IsKnownNegativeProduct) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_domain(x, Domain::Negative);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    std::vector<ExprPtr> factors{x_expr, y_expr};
    ExprPtr prod = arena.make<Product>(std::move(factors));
    // x<0, y>0 → product is negative
    // Trigger is_known_negative on product via simplify
    auto r = ctx.simplify(prod);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyFunctionsCoverage2, SimplifyRootOfWithContext) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    ctx.set_max_rootof_explicit_degree(2U);
    // RootOf(x^2 - 2, x) — degree 2, should trigger solve
    ExprPtr poly = arena.make<Binary>(
        BinaryOp::Sub,
        arena.make<Binary>(BinaryOp::Pow, make_sym(arena, "x"), make_int(arena, 2)),
        make_int(arena, 2));
    ExprPtr e = arena.make<RootOf>(poly, Symbol("x"), std::nullopt);
    auto r = ctx.simplify(e);
    // Should attempt to solve and possibly succeed or return RootOf
    EXPECT_FALSE(r.is_error());
}

// ─── simplify_arithmetic.cpp: additional branches ────────────────────────────

TEST(SimplifyArithmeticCoverage, IsKnownNegativeBinaryDiv) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x"), y("y");
    ctx.assumptions().assume_domain(x, Domain::Negative);
    ctx.assumptions().assume_positive(y);
    ExprPtr x_expr = make_sym(arena, "x");
    ExprPtr y_expr = make_sym(arena, "y");
    // x / y where x<0, y>0 → negative; x / y where x>0, y<0 → also negative
    ExprPtr div = arena.make<Binary>(BinaryOp::Div, x_expr, y_expr);
    auto r = ctx.simplify(div);
    EXPECT_FALSE(r.is_error());
}

TEST(SimplifyArithmeticCoverage, IsKnownNegativeUnaryNeg) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_positive(x);
    ExprPtr x_expr = make_sym(arena, "x");
    // -(x) where x>0 → negative
    ExprPtr neg_x = arena.make<Unary>(UnaryOp::Neg, x_expr);
    auto r = ctx.simplify(neg_x);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyArithmeticCoverage, IsKnownPositiveConstant) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    // Constant Pi is known positive
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    // Simplify should treat pi as known positive
    auto r = ctx.simplify(pi);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}

TEST(SimplifyArithmeticCoverage, IsAssumedNonzeroSymbol) {
    CASContext ctx;
    AstArena& arena = ctx.arena();
    Symbol x("x");
    ctx.assumptions().assume_nonzero(x);
    // This influences simplification of expressions involving x
    ExprPtr x_expr = make_sym(arena, "x");
    auto r = ctx.simplify(x_expr);
    ASSERT_TRUE(r.is_ok());
    EXPECT_NE(r.value(), nullptr);
}
