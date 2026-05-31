// F1.5 — Discrimination Net unit tests.
//
// Covers:
//   1. Insert 100 rules with diverse pattern kinds; lookup returns correct candidates.
//   2. Pattern type filter: x_Integer rejects RationalLit, accepts IntegerLit.
//   3. Pattern type filter: x_Symbol rejects IntegerLit, accepts Symbol.
//   4. AC matching: a+b+c matches pattern x_+y_+z_ without exponential backtracking.
//   5. FuncCall sub-bucket: sin rules not returned for cos expressions.
//   6. Universal wildcard: rules with wildcard root returned for all expressions.
//   7. clear() resets the net.

#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "rewrite/discrimination_net.hpp"

using namespace cas;
using namespace cas::symbolic;
using namespace cas::rewrite;

namespace {

// Build minimal rules for testing.  We don't need valid LPO-oriented rules
// here because we are testing the net's indexing and lookup, not apply_rule.
[[nodiscard]] RewriteRule make_rule(ExprPtr pattern, ExprPtr replacement) {
    return RewriteRule{.pattern = pattern, .replacement = replacement, .condition = {}};
}

// Arena shared by test fixtures (not interned, just allocated).
struct NetFixture : public ::testing::Test {
protected:
    AstArena arena;
    AstArena pat_arena;  // separate arena for patterns (mirrors production usage)

    ExprPtr make_int(long long v) {
        return pat_arena.make<IntegerLit>(BigInt(v));
    }
    ExprPtr make_sym(const std::string& name) {
        return pat_arena.make<Symbol>(name);
    }
    ExprPtr make_func(const std::string& name, std::vector<ExprPtr> args) {
        return pat_arena.make<FuncCall>(name, std::move(args));
    }
    ExprPtr make_sum(std::vector<ExprPtr> terms) {
        return pat_arena.make<Sum>(std::move(terms));
    }
};

} // anonymous namespace

// ── Test 1: 100 rules, lookup returns only matching candidates ────────────────
TEST_F(NetFixture, InsertHundredRules_LookupReturnsCandidates) {
    DiscriminationNet net;

    // Insert 50 FuncCall(sin, x_) rules and 50 FuncCall(cos, x_) rules.
    std::vector<RewriteRule> rules_storage;
    rules_storage.reserve(100U);

    ExprPtr wildcard = make_sym("x_");
    for (int i = 0; i < 50; ++i) {
        ExprPtr pat_sin = make_func("sin", {wildcard});
        ExprPtr rep     = make_int(i);
        rules_storage.push_back(make_rule(pat_sin, rep));
    }
    for (int i = 0; i < 50; ++i) {
        ExprPtr pat_cos = make_func("cos", {wildcard});
        ExprPtr rep     = make_int(100 + i);
        rules_storage.push_back(make_rule(pat_cos, rep));
    }

    for (const auto& r : rules_storage) {
        net.insert(r);
    }

    EXPECT_EQ(net.size(), 100U);

    // sin(x) expression → should return exactly the 50 sin rules.
    ExprPtr sin_x = arena.make<FuncCall>("sin",
        std::vector<ExprPtr>{arena.make<Symbol>("x")});
    const auto& sin_cands = net.lookup(sin_x);
    EXPECT_EQ(sin_cands.size(), 50U)
        << "Expected 50 sin candidates, got " << sin_cands.size();

    // cos(y) expression → should return exactly the 50 cos rules.
    ExprPtr cos_y = arena.make<FuncCall>("cos",
        std::vector<ExprPtr>{arena.make<Symbol>("y")});
    const auto& cos_cands = net.lookup(cos_y);
    EXPECT_EQ(cos_cands.size(), 50U)
        << "Expected 50 cos candidates, got " << cos_cands.size();

    // A Symbol expression → should return 0 candidates (no Symbol-root rules).
    ExprPtr sym_z = arena.make<Symbol>("z");
    const auto& sym_cands = net.lookup(sym_z);
    EXPECT_EQ(sym_cands.size(), 0U);
}

// ── Test 2: Pattern type filter — x_Integer rejects RationalLit ──────────────
TEST(DiscriminationNetWildcardType, IntegerFilterRejectsRational) {
    EXPECT_EQ(wildcard_type_from_name("x_Integer"), WildcardType::Integer);

    AstArena arena;
    ExprPtr rat = arena.make<RationalLit>(BigInt(1), BigInt(3));
    EXPECT_FALSE(expr_satisfies_wildcard_type(rat, WildcardType::Integer));
}

TEST(DiscriminationNetWildcardType, IntegerFilterAcceptsIntegerLit) {
    AstArena arena;
    ExprPtr il = arena.make<IntegerLit>(BigInt(42));
    EXPECT_TRUE(expr_satisfies_wildcard_type(il, WildcardType::Integer));
}

// ── Test 3: Pattern type filter — x_Symbol ───────────────────────────────────
TEST(DiscriminationNetWildcardType, SymbolFilterRejectsIntegerLit) {
    AstArena arena;
    ExprPtr il = arena.make<IntegerLit>(BigInt(5));
    EXPECT_FALSE(expr_satisfies_wildcard_type(il, WildcardType::Symbol));
}

TEST(DiscriminationNetWildcardType, SymbolFilterAcceptsSymbol) {
    AstArena arena;
    ExprPtr sym = arena.make<Symbol>("x");
    EXPECT_TRUE(expr_satisfies_wildcard_type(sym, WildcardType::Symbol));
}

// ── Test 4: AC matching — a+b+c vs x_+y_+z_ without exponential backtracking ─
// This test verifies that MatchMap is populated correctly for a 3-term Sum
// match via the existing match_ac_pattern API.  The discrimination net feeds
// Sum patterns into Sum buckets, ensuring only Sum rules are attempted.
TEST_F(NetFixture, AC_SumPatternMatches_ThreeTerms) {
    DiscriminationNet net;

    // Pattern: x_ + y_ + z_
    ExprPtr wx = make_sym("x_");
    ExprPtr wy = make_sym("y_");
    ExprPtr wz = make_sym("z_");
    ExprPtr pat_sum = make_sum({wx, wy, wz});
    ExprPtr rep = make_int(0);
    RewriteRule rule = make_rule(pat_sum, rep);
    net.insert(rule);

    // Expression: a + b + c
    ExprPtr a = arena.make<Symbol>("a");
    ExprPtr b = arena.make<Symbol>("b");
    ExprPtr c = arena.make<Symbol>("c");
    ExprPtr expr_sum = arena.make<Sum>(std::vector<ExprPtr>{a, b, c});

    // Net returns the Sum rule as a candidate.
    const auto& cands = net.lookup(expr_sum);
    ASSERT_EQ(cands.size(), 1U) << "Net should return the Sum pattern as candidate";

    // Now verify full AC match succeeds.
    MatchMap m;
    EXPECT_TRUE(match_ac_pattern(expr_sum, cands[0]->pattern, m))
        << "AC match of (a+b+c) against (x_+y_+z_) should succeed";

    // Verify all three wildcards were bound.
    EXPECT_EQ(m.count("x_") + m.count("y_") + m.count("z_"), 3U)
        << "All three wildcard slots should be bound";
}

// ── Test 5: FuncCall sub-bucket isolation ────────────────────────────────────
TEST_F(NetFixture, FuncCallBucketIsolation_SinNotReturnedForCos) {
    DiscriminationNet net;

    ExprPtr wc = make_sym("x_");
    RewriteRule sin_rule = make_rule(make_func("sin", {wc}), make_int(1));
    RewriteRule cos_rule = make_rule(make_func("cos", {wc}), make_int(2));
    net.insert(sin_rule);
    net.insert(cos_rule);

    ExprPtr cos_expr = arena.make<FuncCall>("cos",
        std::vector<ExprPtr>{arena.make<Symbol>("t")});
    const auto& cands = net.lookup(cos_expr);
    ASSERT_EQ(cands.size(), 1U);
    EXPECT_EQ(cands[0], &cos_rule) << "Should return only the cos rule";
}

// ── Test 6: Universal wildcard — returned for every ExprKind ─────────────────
TEST_F(NetFixture, UniversalWildcard_ReturnedForAllKinds) {
    DiscriminationNet net;

    // A rule with plain wildcard pattern (root = Symbol "x_").
    ExprPtr wc = make_sym("x_");
    RewriteRule universal_rule = make_rule(wc, make_int(99));
    net.insert(universal_rule);

    // Should appear in lookup for IntegerLit, Symbol, FuncCall, Sum.
    ExprPtr il = arena.make<IntegerLit>(BigInt(5));
    ExprPtr sym = arena.make<Symbol>("y");
    ExprPtr fc = arena.make<FuncCall>("exp",
        std::vector<ExprPtr>{arena.make<Symbol>("z")});
    ExprPtr s = arena.make<Sum>(std::vector<ExprPtr>{il, sym});

    EXPECT_EQ(net.lookup(il).size(), 1U)  << "Universal rule not returned for IntegerLit";
    EXPECT_EQ(net.lookup(sym).size(), 1U) << "Universal rule not returned for Symbol";
    EXPECT_EQ(net.lookup(fc).size(), 1U)  << "Universal rule not returned for FuncCall";
    EXPECT_EQ(net.lookup(s).size(), 1U)   << "Universal rule not returned for Sum";
}

// ── Test 7: clear() resets the net ───────────────────────────────────────────
TEST_F(NetFixture, Clear_ResetsNet) {
    DiscriminationNet net;
    ExprPtr wc = make_sym("x_");
    RewriteRule r = make_rule(make_func("ln", {wc}), make_int(0));
    net.insert(r);
    EXPECT_EQ(net.size(), 1U);

    net.clear();
    EXPECT_EQ(net.size(), 0U);

    ExprPtr expr = arena.make<FuncCall>("ln",
        std::vector<ExprPtr>{arena.make<Symbol>("x")});
    EXPECT_EQ(net.lookup(expr).size(), 0U);
}

// ── Test 8: wildcard_type_from_name corner cases ──────────────────────────────
TEST(DiscriminationNetWildcardType, PlainWildcard_IsAny) {
    EXPECT_EQ(wildcard_type_from_name("x_"), WildcardType::Any);
    EXPECT_EQ(wildcard_type_from_name("abc_"), WildcardType::Any);
}

TEST(DiscriminationNetWildcardType, PositiveFilter) {
    EXPECT_EQ(wildcard_type_from_name("n_Positive"), WildcardType::Positive);
    AstArena arena;
    ExprPtr pos = arena.make<IntegerLit>(BigInt(7));
    ExprPtr zero = arena.make<IntegerLit>(BigInt(0));
    ExprPtr neg  = arena.make<IntegerLit>(BigInt(-3));
    EXPECT_TRUE(expr_satisfies_wildcard_type(pos, WildcardType::Positive));
    EXPECT_FALSE(expr_satisfies_wildcard_type(zero, WildcardType::Positive));
    EXPECT_FALSE(expr_satisfies_wildcard_type(neg, WildcardType::Positive));
}
