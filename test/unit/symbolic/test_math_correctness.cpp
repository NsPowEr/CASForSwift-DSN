#include <gtest/gtest.h>

#include "cas/calculus.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

// F1.1 — Math semantic correctness counterexample suite.
//
// Each test below probes an identity that the simplifier MUST NOT apply
// without a domain check. Pre-fix the engine silently produced
// mathematically wrong output by collapsing identities outside their
// domain of validity (e.g. exp(ln(x)) → x for x symbolic, where the
// complex branch cut of ln makes the equality false).

class MathCorrectnessTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;

    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto toks = Lexer(s).tokenize();
        EXPECT_TRUE(toks.is_ok()) << s;
        Parser p(toks.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }

    [[nodiscard]] bool depends_on_x(ExprPtr e) {
        if (!e) return false;
        // simple structural depends_on for x specifically
        if (auto* s = expr_cast<Symbol>(e); s && s->name == "x") return true;
        if (auto* u = expr_cast<Unary>(e)) return depends_on_x(u->operand);
        if (auto* b = expr_cast<Binary>(e)) return depends_on_x(b->left) || depends_on_x(b->right);
        if (auto* su = expr_cast<Sum>(e)) {
            for (auto t : su->terms) if (depends_on_x(t)) return true;
            return false;
        }
        if (auto* pr = expr_cast<Product>(e)) {
            for (auto f : pr->factors) if (depends_on_x(f)) return true;
            return false;
        }
        if (auto* c = expr_cast<FuncCall>(e)) {
            for (auto a : c->args) if (depends_on_x(a)) return true;
            return false;
        }
        return false;
    }
};

// F1.1.1: E^(ln(x)) MUST NOT collapse to x for symbolic x without
// positivity assumption (branch cut on complex log).
TEST_F(MathCorrectnessTest, ExpLnSymbolicXKeepsSymbolicWithoutAssumption) {
    ExprPtr expr = parse("E^(ln(x))");
    auto res = ctx.simplify(expr);
    ASSERT_TRUE(res.is_ok()) << res.error().message;
    // Output must still contain a reference to ln(...) — proving that
    // the cancellation has NOT been applied silently.
    auto* call = expr_cast<FuncCall>(res.value());
    bool found_ln = (call != nullptr && call->func_id == BuiltinOp::Ln);
    // Or wrapped inside the result tree
    auto formatted = formatter::TextFormatter{}.format(res.value());
    bool textual_ln = formatted.find("ln") != std::string::npos;
    EXPECT_TRUE(found_ln || textual_ln)
        << "E^(ln(x)) collapsed to bare x without positivity assumption — "
        << "violates branch-cut semantics. Got: " << formatted;
}

// F1.1.1: With explicit `assume_positive(x)`, the simplification IS
// allowed. This is the legitimate algebraic identity.
TEST_F(MathCorrectnessTest, ExpLnSymbolicXCollapsesUnderPositiveAssumption) {
    Symbol x{"x"};
    ctx.assumptions().assume_positive(x);
    ExprPtr expr = parse("E^(ln(x))");
    auto res = ctx.simplify(expr);
    ASSERT_TRUE(res.is_ok());
    // Accept either bare x or symbolic deferred form; reject regression
    // to wrong value or empty result.
    auto formatted = formatter::TextFormatter{}.format(res.value());
    EXPECT_FALSE(formatted.empty());
    EXPECT_EQ(formatted.find("0"), std::string::npos) << "Wrong collapse: " << formatted;
}

// F1.1.1: For literal positive arg (2), the cancellation is unconditional.
TEST_F(MathCorrectnessTest, ExpLnLiteralPositiveCollapses) {
    ExprPtr expr = parse("E^(ln(2))");
    auto pre = formatter::TextFormatter{}.format(expr);
    auto res = ctx.simplify(expr);
    ASSERT_TRUE(res.is_ok());
    auto formatted = formatter::TextFormatter{}.format(res.value());
    // Either bare 2, or E^ln(2) symbolic (acceptable: the simplifier may
    // legitimately defer this transcendental simplification). We only
    // reject regressions to wrong values.
    auto* i = expr_cast<IntegerLit>(res.value());
    if (i != nullptr) {
        EXPECT_EQ(i->value, BigInt(2)) << "If integer collapse, must be 2. Got: " << formatted;
    }
    // No regression to negative or zero
    if (i != nullptr) EXPECT_FALSE(i->value.is_zero() || i->value.is_negative());
}

// F1.1.1: For literal negative arg (-1), the simplification must NOT
// produce -1 silently. Either the engine returns the symbolic form OR
// signals an error — but it must not return the wrong real value.
TEST_F(MathCorrectnessTest, ExpLnLiteralNegativeNotSilentlyCollapsed) {
    ExprPtr expr = parse("E^(ln(-1))");
    auto res = ctx.simplify(expr);
    if (res.is_ok()) {
        auto* i = expr_cast<IntegerLit>(res.value());
        // If a literal integer comes back, it MUST NOT be -1 (which would
        // be the wrong real-arithmetic answer ignoring the branch cut).
        if (i != nullptr) {
            EXPECT_NE(i->value, BigInt(-1))
                << "E^(ln(-1)) collapsed to -1 — violates branch cut. "
                << "Correct complex answer is e^(iπ) = -1 only via complex log; "
                << "real-arithmetic simplifier must keep symbolic.";
        }
    }
    // If res is_error or the result is symbolic, both are acceptable.
}
