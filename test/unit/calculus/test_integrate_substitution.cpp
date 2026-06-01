#include <gtest/gtest.h>
#include "cas/calculus.hpp"
#include "cas/symbolic.hpp"
#include "cas/parser.hpp"
#include "cas/lexer.hpp"
#include "cas/formatter.hpp"

using namespace cas;
using namespace cas::calculus;

class SubstitutionTest : public ::testing::Test {
protected:
    symbolic::CASContext ctx;
    Symbol x{"x"};

    ExprPtr parse(const std::string& s) {
        auto t = Lexer(s).tokenize();
        EXPECT_TRUE(t.is_ok()) << s;
        Parser p(t.value(), ctx.arena());
        auto r = p.parse();
        EXPECT_TRUE(r.is_ok()) << s;
        return r.value();
    }
    
    void verify_integral(const std::string& integrand_str) {
        auto integrand = parse(integrand_str);
        auto result = integrate(integrand, x, ctx);
        ASSERT_TRUE(result.is_ok()) << "Integration failed for: " << integrand_str << " - " << result.error().message;
        
        auto derived = diff(result.value(), x, 1U, ctx);
        ASSERT_TRUE(derived.is_ok());
        
        auto s_derived = simplify(derived.value(), ctx);
        auto s_integrand = simplify(integrand, ctx);
        
        ASSERT_TRUE(s_derived.is_ok());
        ASSERT_TRUE(s_integrand.is_ok());
        
        auto eq_res = symbolic::mathematically_equal(s_derived.value(), s_integrand.value(), ctx);
        ASSERT_TRUE(eq_res.is_ok());
        
        EXPECT_TRUE(eq_res.value())
            << "Verification failed for: " << integrand_str 
            << "\nResult: " << cas::formatter::TextFormatter().format(result.value())
            << "\nDiff: " << cas::formatter::TextFormatter().format(s_derived.value())
            << "\nOriginal: " << cas::formatter::TextFormatter().format(s_integrand.value());
    }
};

TEST_F(SubstitutionTest, PowerChainRule) {
    // ∫ 2*x*(x^2 + 1)^3 dx = 1/4 * (x^2 + 1)^4
    verify_integral("2*x*(x^2 + 1)^3");
}

TEST_F(SubstitutionTest, ExponentialTrig) {
    // ∫ cos(x)*exp(sin(x)) dx = exp(sin(x))
    verify_integral("cos(x)*exp(sin(x))");
}

TEST_F(SubstitutionTest, Logarithmic) {
    // ∫ ln(x)/x dx = 1/2 * ln(x)^2
    verify_integral("ln(x)/x");
}

TEST_F(SubstitutionTest, SineOfSquare) {
    // ∫ 2*x*sin(x^2) dx = -cos(x^2)
    verify_integral("2*x*sin(x^2)");
}

TEST_F(SubstitutionTest, TangentLog) {
    // ∫ tan(x) dx = -ln(abs(cos(x)))
    // Note: simplify(sin(x)/cos(x)) is tan(x) or stays sin/cos?
    // Let's use the explicit form that triggers u-sub.
    verify_integral("sin(x)/cos(x)"); 
}

TEST_F(SubstitutionTest, ComplexNested) {
    // ∫ (exp(x)+1)*cos(exp(x)+x) dx = sin(exp(x)+x)
    verify_integral("(exp(x)+1)*cos(exp(x)+x)");
}
