// F3.6 probe — verify if MRV cancellation tower is actually broken.
// Audit claimed compare_growth(2·e^x, e^x) returns 0 because leading
// coefficients are ignored, leading to wrong limit values for
// cancellation cases. This probe codifies the alleged failure mode.

#include <gtest/gtest.h>

#include "cas/calculus.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;
using namespace cas::calculus;

class MrvCancellationProbeTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};
    ExprPtr inf() { return ctx.arena().make<Constant>(MathConstant::Infinity); }
    [[nodiscard]] ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
};

// Trivial: (e^x - e^x) / x = 0 / x = 0 (or undefined intermediate).
// limit should be 0.
TEST_F(MrvCancellationProbeTest, ExpMinusExpIsZero) {
    auto e = parse("(exp(x) - exp(x)) / x");
    auto r = limit(e, x, inf(), LimitDirection::Both, ctx);
    if (r.is_ok()) {
        auto formatted = formatter::TextFormatter{}.format(r.value());
        EXPECT_NE(formatted.find("0"), std::string::npos)
            << "Expected 0. Got: " << formatted;
    } else {
        std::cout << "Engine bailed: " << r.error().message << std::endl;
    }
}

// Cancellation that depends on leading coefficient:
// lim_{x→∞} (e^x + 1) / e^x = 1
TEST_F(MrvCancellationProbeTest, ExpPlusOneOverExp) {
    auto e = parse("(exp(x) + 1) / exp(x)");
    auto r = limit(e, x, inf(), LimitDirection::Both, ctx);
    ASSERT_TRUE(r.is_ok()) << r.error().message;
    auto formatted = formatter::TextFormatter{}.format(r.value());
    EXPECT_EQ(formatted, "1") << "Expected 1. Got: " << formatted;
}

// Real cancellation: lim_{x→∞} (2·e^x - e^x) = ∞ (NOT 0)
// This tests if compare_growth correctly identifies "same growth class
// but different leading coefficient → result depends on coefficient
// difference, not zero".
TEST_F(MrvCancellationProbeTest, TwoExpMinusExpDoesNotCancelToZero) {
    auto e = parse("2*exp(x) - exp(x)");
    auto r = limit(e, x, inf(), LimitDirection::Both, ctx);
    if (r.is_ok()) {
        auto formatted = formatter::TextFormatter{}.format(r.value());
        // Should be ∞, NOT 0.
        EXPECT_EQ(formatted.find("0"), std::string::npos)
            << "2·e^x - e^x = e^x → ∞, must NOT collapse to 0. Got: "
            << formatted;
    }
}
