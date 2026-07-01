// A28 — regression guard for to_field_generators generator-mapping bug.
//
// Bug: to_field_generators() substituted tower extensions in FORWARD order
// (t_0 = exp(x) before t_1 = exp(exp(x))). Substituting the inner generator
// first rewrote exp(x) -> t_0 *inside* the outer pattern exp(exp(x)),
// turning it into exp(t_0), so the t_1 pattern no longer matched. The outer
// transcendental was left as a free symbol, silently misclassified as a
// constant coefficient by the exponential poly-part solver — producing a
// wrong elementary "antiderivative" for a genuinely non-elementary integrand
// (e.g. integrate() returned x^2/2 * exp(exp(x)) for x*exp(exp(x)), which
// does not differentiate back to the integrand). Fixed by substituting
// outermost-first (reverse tower order), matching from_field_generators.
//
// These integrands are not known to be elementary, so the sound outcome is
// Unimplemented. The regression contract is: if integrate() ever does return
// ok(), the candidate MUST verify by differentiation (never silent-wrong).

#include <gtest/gtest.h>

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class RischNestedTranscendentalTest : public ::testing::Test {
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
    [[nodiscard]] bool verify_antider(ExprPtr F, ExprPtr expr) {
        auto D = calculus::diff(F, x, 1U, ctx);
        if (D.is_error()) return false;
        ExprPtr delta = ctx.arena().make<Binary>(BinaryOp::Sub, D.value(), expr);
        auto simp = ctx.simplify(delta);
        if (simp.is_error()) return false;
        auto* lit = expr_cast<IntegerLit>(simp.value());
        return lit != nullptr && lit->value.is_zero();
    }
    void expect_sound_or_unimplemented(const std::string& src) {
        auto e = parse(src);
        auto r = calculus::integrate(e, x, ctx);
        if (r.is_ok()) {
            EXPECT_TRUE(verify_antider(r.value(), e))
                << "silent-wrong: D(candidate) != integrand for " << src;
        }
        // Unimplemented is an acceptable sound outcome (REGOLA 0.2).
    }
};

TEST_F(RischNestedTranscendentalTest, XTimesExpExpXNeverSilentWrong) {
    expect_sound_or_unimplemented("x * exp(exp(x))");
}

TEST_F(RischNestedTranscendentalTest, XTimesExpExpXTimesExpXNeverSilentWrong) {
    expect_sound_or_unimplemented("x * exp(exp(x)) * exp(x)");
}

TEST_F(RischNestedTranscendentalTest, TwoExp2xExpExpXSquaredNeverSilentWrong) {
    expect_sound_or_unimplemented("2 * exp(2*x) * exp(exp(x^2))");
}

// exp(exp(x)) * exp(x) IS elementary (D(exp(exp(x))) = exp(x)*exp(exp(x)))
// and must stay correctly solved after the fix (generator-mapping change
// must not regress the already-working case).
TEST_F(RischNestedTranscendentalTest, ExpExpXTimesExpXStillSolved) {
    auto e = parse("exp(exp(x)) * exp(x)");
    auto r = calculus::integrate(e, x, ctx);
    ASSERT_TRUE(r.is_ok()) << (r.is_error() ? r.error().message : "");
    EXPECT_TRUE(verify_antider(r.value(), e));
}

}  // namespace
