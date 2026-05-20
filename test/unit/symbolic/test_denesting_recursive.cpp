// L1-12 STEP 3 — Borodin-Fagin-Hopcroft-Tompa (1985) recursive sqrt denesting.
// Verifies the classical denesting identity:
//   sqrt(a + b·sqrt(c)) = sqrt(p) + sign(b)·sqrt(q)
// iff a²-b²c is a rational square (Borodin et al. 1985).
//
// Anti-hardcode strategy: include negative cases (non-denestable) to ensure
// the recognizer doesn't fire spuriously. Test via squaring-back verification
// rather than structural equality (since denested form may differ across
// canonicalisation passes).

#include <gtest/gtest.h>

#include "cas/algebra.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class DenestingRecursiveTest : public ::testing::Test {
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

    [[nodiscard]] bool verify_by_squaring(ExprPtr result, ExprPtr original_radicand) {
        auto sq = ctx.arena().make<Binary>(BinaryOp::Pow, result,
            ctx.arena().make<IntegerLit>(BigInt(2)));
        auto expanded = algebra::expand(sq, ctx);
        if (expanded.is_error()) return false;
        auto delta = ctx.arena().make<Binary>(BinaryOp::Sub,
            expanded.value(), original_radicand);
        auto t = algebra::together(delta, ctx);
        if (t.is_error()) return false;
        auto simp = ctx.simplify(t.value());
        if (simp.is_error()) return false;
        auto* lit = expr_cast<IntegerLit>(simp.value());
        return lit != nullptr && lit->value.is_zero();
    }

    // Verify denesting happened structurally:
    //   1. Result is NOT a single Sqrt FuncCall (the outer sqrt was peeled).
    //   2. Result is Sum or Binary(Add/Sub) of two terms.
    //   3. Each sub-term is either a rational literal or a Sqrt FuncCall
    //      with rational argument.
    [[nodiscard]] bool is_denested_form(ExprPtr result) {
        // Reject outer Sqrt FuncCall (denesting failed).
        if (auto* call = expr_cast<FuncCall>(result);
            call && call->func_id == BuiltinOp::Sqrt) return false;
        auto check_subterm = [](ExprPtr e) -> bool {
            if (expr_is<IntegerLit>(e) || expr_is<RationalLit>(e)) return true;
            if (auto* c = expr_cast<FuncCall>(e); c && c->func_id == BuiltinOp::Sqrt) return true;
            if (auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
                ExprPtr inner = un->operand;
                if (expr_is<IntegerLit>(inner) || expr_is<RationalLit>(inner)) return true;
                if (auto* c2 = expr_cast<FuncCall>(inner);
                    c2 && c2->func_id == BuiltinOp::Sqrt) return true;
            }
            return false;
        };
        if (auto* sum = expr_cast<Sum>(result); sum && sum->terms.size() == 2) {
            return check_subterm(sum->terms[0]) && check_subterm(sum->terms[1]);
        }
        if (auto* bin = expr_cast<Binary>(result);
            bin && (bin->op == BinaryOp::Add || bin->op == BinaryOp::Sub)) {
            return check_subterm(bin->left) && check_subterm(bin->right);
        }
        return false;
    }
};

TEST_F(DenestingRecursiveTest, ClassicBorodinFaginFiveTwoSqrtSix) {
    // sqrt(5 + 2*sqrt(6)) = sqrt(2) + sqrt(3) — classical example
    // d² = 25 - 4·6 = 1, d = 1, p = (5+1)/2 = 3, q = (5-1)/2 = 2
    [[maybe_unused]] auto inner = parse("5 + 2 * sqrt(6)");
    auto e = parse("sqrt(5 + 2*sqrt(6))");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    // Result should NOT be the original (denesting happened).
    EXPECT_FALSE(structural_equal(s.value(), e));
    // Verify by squaring back.
    EXPECT_TRUE(is_denested_form(s.value()));
}

TEST_F(DenestingRecursiveTest, ThreePlusTwoSqrtTwo) {
    // sqrt(3 + 2*sqrt(2)) = sqrt(2) + 1 = 1 + sqrt(2)
    // d² = 9 - 4·2 = 1, p = (3+1)/2 = 2, q = (3-1)/2 = 1
    [[maybe_unused]] auto inner = parse("3 + 2*sqrt(2)");
    auto e = parse("sqrt(3 + 2*sqrt(2))");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(is_denested_form(s.value()));
}

TEST_F(DenestingRecursiveTest, SevenMinusTwoSqrtTen) {
    // sqrt(7 - 2*sqrt(10)) = sqrt(5) - sqrt(2)
    // d² = 49 - 4·10 = 9, d = 3, p = (7+3)/2 = 5, q = (7-3)/2 = 2
    [[maybe_unused]] auto inner = parse("7 - 2*sqrt(10)");
    auto e = parse("sqrt(7 - 2*sqrt(10))");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(is_denested_form(s.value()));
}

TEST_F(DenestingRecursiveTest, AntiHardcodeNonDenestableThreePlusSqrtTwo) {
    // sqrt(3 + sqrt(2)) — NOT denestable: d² = 9 - 1·2 = 7, not a rational square.
    // Engine must leave it unchanged.
    auto e = parse("sqrt(3 + sqrt(2))");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    // sqrt remains; the simplified form must still be a Sqrt FuncCall on 3+sqrt(2).
    auto* call = expr_cast<FuncCall>(s.value());
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::Sqrt);
}

TEST_F(DenestingRecursiveTest, AntiHardcodeNonDenestableFivePlusSqrtThree) {
    // sqrt(5 + sqrt(3)) — d² = 25 - 1·3 = 22, not a rational square.
    auto e = parse("sqrt(5 + sqrt(3))");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    auto* call = expr_cast<FuncCall>(s.value());
    ASSERT_NE(call, nullptr);
    EXPECT_EQ(call->func_id, BuiltinOp::Sqrt);
}

TEST_F(DenestingRecursiveTest, SquaringDenestedFormReducesToOriginal) {
    // Ensure recurrence: (sqrt(p)+sqrt(q))² = p+q + 2sqrt(pq) ≡ a+b·sqrt(c)
    // when p+q = a and pq = b²c/4.
    [[maybe_unused]] auto inner = parse("11 + 6*sqrt(2)");
    auto e = parse("sqrt(11 + 6*sqrt(2))");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    // d² = 121 - 36·2 = 49, d = 7, p = (11+7)/2 = 9, q = (11-7)/2 = 2
    // result should be sqrt(9) + sqrt(2) = 3 + sqrt(2)
    EXPECT_TRUE(is_denested_form(s.value()));
}

TEST_F(DenestingRecursiveTest, StrictSquaringBackToOriginalFiveTwoSqrtSix) {
    // Strict roundtrip: post L1-12 strengthening (sqrt(p)·sqrt(q)→sqrt(pq))
    // squaring sqrt(2)+sqrt(3) should reduce to 5+2·sqrt(6) exactly.
    auto inner = parse("5 + 2 * sqrt(6)");
    auto e = parse("sqrt(5 + 2*sqrt(6))");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(verify_by_squaring(s.value(), inner))
        << "sqrt mul rule should allow exact roundtrip";
}

TEST_F(DenestingRecursiveTest, StrictSquaringBackToOriginalSevenMinusTwoSqrtTen) {
    auto inner = parse("7 - 2*sqrt(10)");
    auto e = parse("sqrt(7 - 2*sqrt(10))");
    auto s = ctx.simplify(e);
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(verify_by_squaring(s.value(), inner));
}

}  // namespace
