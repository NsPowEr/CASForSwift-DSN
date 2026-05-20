#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "../../helpers/property_test.hpp"
#include <gtest/gtest.h>
#include <string>

namespace cas {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tok = Lexer(input).tokenize();
    if (tok.is_error()) return fail<ExprPtr>(tok.error());
    return Parser(tok.value(), arena).parse();
}

} // namespace

// L0-02: simplify must be idempotent — simplify(simplify(e)) == simplify(e)
TEST(RandomizedAntiHardcode, SimplifyIdempotency) {
    std::size_t violations = 0;
    test::run_seeded_cases(0xDEAD'BEEF'1234ULL, 50U, [&](test::DeterministicRng& rng, std::size_t /*index*/) {
        symbolic::CASContext ctx;
        const std::string expr_str = test::generate_polynomial_expr(rng, 3);

        auto parsed = parse_expr(expr_str, ctx.arena());
        if (!parsed.is_ok()) return;

        auto s1 = ctx.simplify(parsed.value());
        if (!s1.is_ok()) return;

        auto s2 = ctx.simplify(s1.value());
        if (!s2.is_ok()) return;

        if (s1.value() != s2.value()) {
            ++violations;
        }
    });
    EXPECT_EQ(violations, 0U) << "Idempotency violated on " << violations << " out of 50 random expressions";
}

// L0-08: f * 1 must simplify to f (neutral element — anti-hardcode metamorphic check)
TEST(RandomizedAntiHardcode, MultiplyByOneNeutral) {
    std::size_t violations = 0;
    test::run_seeded_cases(0xCAFE'BABE'5678ULL, 40U, [&](test::DeterministicRng& rng, std::size_t /*index*/) {
        symbolic::CASContext ctx;
        const std::string expr_str = test::generate_polynomial_expr(rng, 2);

        auto parsed = parse_expr(expr_str, ctx.arena());
        if (!parsed.is_ok()) return;

        auto s_f = ctx.simplify(parsed.value());
        if (!s_f.is_ok()) return;

        ExprPtr one = ctx.arena().make<IntegerLit>(BigInt(1));
        ExprPtr f_times_1 = ctx.arena().make<Product>(std::vector<ExprPtr>{parsed.value(), one});
        auto s_f1 = ctx.simplify(f_times_1);
        if (!s_f1.is_ok()) return;

        if (s_f.value() != s_f1.value()) {
            ++violations;
        }
    });
    EXPECT_EQ(violations, 0U) << "Neutral element f*1==f violated on " << violations << " out of 40 cases";
}

// L0-08: differentiation linearity — d(f+g)/dx == df/dx + dg/dx
TEST(RandomizedAntiHardcode, DifferentiationLinearity) {
    std::size_t violations = 0;
    test::run_seeded_cases(0xF00D'1337'ABCDULL, 30U, [&](test::DeterministicRng& rng, std::size_t /*index*/) {
        symbolic::CASContext ctx;
        const std::string f_str = test::generate_polynomial_expr(rng, 2);
        const std::string g_str = test::generate_polynomial_expr(rng, 2);

        auto f = parse_expr(f_str, ctx.arena());
        auto g = parse_expr(g_str, ctx.arena());
        if (!f.is_ok() || !g.is_ok()) return;

        Symbol x{"x"};

        ExprPtr f_plus_g = ctx.arena().make<Sum>(std::vector<ExprPtr>{f.value(), g.value()});
        auto d_sum = calculus::diff(f_plus_g, x, 1U, ctx);
        if (!d_sum.is_ok()) return;

        auto df = calculus::diff(f.value(), x, 1U, ctx);
        auto dg = calculus::diff(g.value(), x, 1U, ctx);
        if (!df.is_ok() || !dg.is_ok()) return;

        ExprPtr df_plus_dg = ctx.arena().make<Sum>(std::vector<ExprPtr>{df.value(), dg.value()});

        auto eq = symbolic::mathematically_equal(d_sum.value(), df_plus_dg, ctx);
        if (!eq.is_ok()) return;

        if (!eq.value()) {
            ++violations;
        }
    });
    EXPECT_EQ(violations, 0U) << "Differentiation linearity violated on " << violations << " out of 30 cases";
}

// L0-08: substitution commutes with simplification — simplify(sub(f,x,0)) == simplify(sub(simplify(f),x,0))
TEST(RandomizedAntiHardcode, SubstitutionCommutesSimplifiy) {
    std::size_t violations = 0;
    test::run_seeded_cases(0xBEEF'CAFE'9999ULL, 25U, [&](test::DeterministicRng& rng, std::size_t /*index*/) {
        symbolic::CASContext ctx;
        const std::string expr_str = test::generate_polynomial_expr(rng, 2);

        auto parsed = parse_expr(expr_str, ctx.arena());
        if (!parsed.is_ok()) return;

        ExprPtr zero = ctx.arena().make<IntegerLit>(BigInt(0));
        Symbol x{"x"};

        // path A: sub then simplify
        auto sub_a = symbolic::substitute(parsed.value(), x, zero, ctx);
        if (!sub_a.is_ok()) return;
        auto s_a = ctx.simplify(sub_a.value());
        if (!s_a.is_ok()) return;

        // path B: simplify then sub
        auto sf = ctx.simplify(parsed.value());
        if (!sf.is_ok()) return;
        auto sub_b = symbolic::substitute(sf.value(), x, zero, ctx);
        if (!sub_b.is_ok()) return;
        auto s_b = ctx.simplify(sub_b.value());
        if (!s_b.is_ok()) return;

        auto eq = symbolic::mathematically_equal(s_a.value(), s_b.value(), ctx);
        if (!eq.is_ok()) return;

        if (!eq.value()) {
            ++violations;
        }
    });
    EXPECT_EQ(violations, 0U) << "Substitution commutativity violated on " << violations << " out of 25 cases";
}

} // namespace cas
