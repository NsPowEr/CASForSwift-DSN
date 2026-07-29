// A34 — the F7.0-A4.2 strict-canonical checker (`is_strictly_canonical`) must
// use the SAME ordering key as the real sorter `merge_symbolic_factors`
// (order Product factors by their BASE, a Pow-with-integer-exponent
// contributing its base). Before the fix the checker compared whole factors,
// raising a false "not strictly canonical" warning on correctly-ordered output
// such as sin(x)·sqrt(pi)^(-1). These tests pin both directions: no false
// alarm on canonical output, and genuine misordering still caught.

#include <gtest/gtest.h>
#include <string>
#include <vector>

#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

namespace cas {
namespace {

[[nodiscard]] ExprPtr parse_expr(const std::string& s, symbolic::CASContext& ctx) {
    auto toks = Lexer(s).tokenize();
    EXPECT_TRUE(toks.is_ok()) << s;
    Parser p(toks.value(), ctx.arena());
    auto r = p.parse();
    EXPECT_TRUE(r.is_ok()) << s;
    return r.value();
}

// Repro from A34: simplify(sin(x)/sqrt(pi)) = sin(x)·sqrt(pi)^(-1). The factor
// sqrt(pi)^(-1) sorts by its base sqrt(pi); the output is correctly ordered and
// must be accepted as strictly canonical (previously a false alarm).
TEST(StrictCanonical, SinOverSqrtPi_NoFalseAlarm) {
    symbolic::CASContext ctx;
    auto s = ctx.simplify(parse_expr("sin(x)/sqrt(pi)", ctx));
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(symbolic::is_strictly_canonical(s.value()))
        << "simplify(sin(x)/sqrt(pi)) wrongly flagged non-canonical";
}

// Sibling case that never regressed (non-integer exponent kept whole by the
// sorter): sin(x)·5^(1/2) — sort key is the whole 5^(1/2) factor.
TEST(StrictCanonical, SinTimesSqrt5_NoFalseAlarm) {
    symbolic::CASContext ctx;
    auto s = ctx.simplify(parse_expr("sin(x)*5^(1/2)", ctx));
    ASSERT_TRUE(s.is_ok());
    EXPECT_TRUE(symbolic::is_strictly_canonical(s.value()));
}

// The checker must NOT become a no-op: a hand-built Product whose Pow factors
// are in the wrong BASE order is still rejected. Canonical base order puts x
// before y; [y^2, x^3] is therefore misordered.
TEST(StrictCanonical, MisorderedPowFactors_StillRejected) {
    symbolic::CASContext ctx;
    AstArena& a = ctx.arena();
    ExprPtr x = a.make<Symbol>("x");
    ExprPtr y = a.make<Symbol>("y");
    ExprPtr y2 = a.make<Binary>(BinaryOp::Pow, y, a.make<IntegerLit>(BigInt(2)));
    ExprPtr x3 = a.make<Binary>(BinaryOp::Pow, x, a.make<IntegerLit>(BigInt(3)));
    // Determine the canonical base order and build the DELIBERATELY reversed one.
    ExprPtr canonical = a.make<Product>(std::vector<ExprPtr>{x3, y2});
    ExprPtr reversed = a.make<Product>(std::vector<ExprPtr>{y2, x3});
    // Exactly one of the two orders is canonical; the other must be rejected.
    const bool c_ok = symbolic::is_strictly_canonical(canonical);
    const bool r_ok = symbolic::is_strictly_canonical(reversed);
    EXPECT_NE(c_ok, r_ok)
        << "checker accepted BOTH factor orders — ordering no longer enforced";
}

}  // namespace
}  // namespace cas
