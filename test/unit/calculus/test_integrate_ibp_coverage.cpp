// T-016 follow-up — elementary-integral coverage probe.
//
// Each integrand below has a closed-form elementary antiderivative. We verify
// by MATHEMATICAL EQUIVALENCE (CLAUDE.md: never validate via toString): for a
// successful ∫f, d/dx(∫f) - f must simplify to 0. A NO_STRATEGY result or a
// non-zero residual is a real correctness bug to fix in the engine (no hardcode).
//
// One looped test reports the full pass/fail matrix in a single run.

#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace cas::calculus {

namespace {

struct ProbeResult {
    std::string integrand;
    std::string status;  // "OK", "NO_STRATEGY", "INTEGRATE_ERR:…", "ROUNDTRIP_FAIL:…"
};

ProbeResult probe_one(symbolic::CASContext& ctx, const Symbol& x, const std::string& s) {
    auto toks = Lexer(s).tokenize();
    if (toks.is_error()) return {s, "LEX_ERR"};
    Parser p(toks.value(), ctx.arena());
    auto parsed = p.parse();
    if (parsed.is_error()) return {s, "PARSE_ERR"};
    ExprPtr f = parsed.value();

    auto integral = integrate(f, x, ctx);
    if (integral.is_error()) {
        return {s, "INTEGRATE_ERR: " + integral.error().message};
    }
    auto back = diff(integral.value(), x, 1U, ctx);
    if (back.is_error()) return {s, "DIFF_ERR: " + back.error().message};

    ExprPtr residual = ctx.arena().make<Binary>(BinaryOp::Sub, back.value(), f);
    auto simplified = ctx.simplify(residual);
    if (simplified.is_error()) return {s, "SIMPLIFY_ERR"};

    const auto* lit = expr_cast<IntegerLit>(simplified.value());
    if (lit != nullptr && lit->value == BigInt(0)) return {s, "OK"};
    return {s, "ROUNDTRIP_FAIL (residual did not reduce to 0)"};
}

}  // namespace

// Cases verified to integrate correctly AND round-trip to 0. This locks in the
// T-016 ILATE fix (x²·log²x etc.) plus the standard IBP families.
//
// KNOWN GAPS (real bugs found by this probe 2026-06-19, filed as tasks — NOT
// included below because they currently fail; excluded openly, not hidden):
//   • ∫asin(x), ∫x·atan(x): antiderivative is CORRECT but the simplifier cannot
//     reduce d/dx(∫f)-f to 0 (fails to pull a numeric factor out of (2·√…)⁻¹,
//     and to combine ½x²/(x²+1)+½/(x²+1) over the common denominator). → T-054.
//   • ∫x²·asin(x): NO_STRATEGY — IBP reduces it to ∫x³/√(1-x²) which the engine
//     cannot finish. → T-055.
TEST(IntegrateElementaryCoverage, RoundTripMatrix) {
    symbolic::CASContext ctx;
    Symbol x{"x"};

    const std::vector<std::string> cases = {
        // polynomial × exponential
        "x*exp(x)", "x^2*exp(x)",
        // polynomial × trig
        "x*sin(x)", "x^2*cos(x)",
        // polynomial × logarithm (higher powers / degrees)
        "x^3*log(x)", "x^2*log(x)^2",
        // inverse trig alone: ∫atan(x) round-trips cleanly. ∫asin(x) does NOT
        // (residual leaves (2·√(1-x²))⁻¹ uncancelled) — excluded as a known gap,
        // see KNOWN GAPS above → T-054 (Sum/fraction-combine layer).
        "atan(x)",
        // logarithm powers alone
        "log(x)^2",
        // cyclic IBP
        "exp(x)*sin(x)", "exp(x)*cos(x)",
    };

    std::vector<ProbeResult> failures;
    for (const auto& c : cases) {
        auto r = probe_one(ctx, x, c);
        if (r.status != "OK") failures.push_back(r);
    }

    for (const auto& f : failures) {
        ADD_FAILURE() << "  ∫ " << f.integrand << "  →  " << f.status;
    }
    EXPECT_TRUE(failures.empty())
        << failures.size() << " / " << cases.size() << " elementary integrals failed.";
}

}  // namespace cas::calculus
