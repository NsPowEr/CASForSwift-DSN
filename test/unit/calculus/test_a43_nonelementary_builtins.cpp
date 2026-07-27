// A43 incremento 1+2 — builtin per le antiderivate non elementari.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Nonelementary_Antiderivatives.md
// Ogni derivata verificata numericamente con mpmath a 30 cifre prima di essere
// scritta a codice — banco riproducibile: scripts/a43_special_fn_check.py.
//
// Incrementi 1+2: nome + round-trip + derivata. Incremento 3: le identita' di
// riduzione della spec §4, riconosciute dentro mathematically_equal.

#include <gtest/gtest.h>

#include <string>
#include <utility>
#include <vector>

#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

using namespace cas;

namespace {

class A43NonelementaryBuiltinsTest : public ::testing::Test {
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
};

// The canonical spelling must resolve back through the builtin table: a name
// that builtin_from_name does not recognise produces an opaque FuncCall whose
// derivative rule never fires (trap closed twice already, A37 for Gamma and
// A44 for factorial).
TEST_F(A43NonelementaryBuiltinsTest, CanonicalNamesResolveToBuiltins) {
    const std::vector<std::pair<std::string, BuiltinOp>> cases{
        {"erfi", BuiltinOp::Erfi},
        {"Ei", BuiltinOp::ExpIntegralEi},
        {"Si", BuiltinOp::SinIntegral},
        {"Ci", BuiltinOp::CosIntegral},
        {"Shi", BuiltinOp::SinhIntegral},
        {"Chi", BuiltinOp::CoshIntegral},
        {"li", BuiltinOp::LogIntegral},
        {"dilog", BuiltinOp::Dilog},
    };
    for (const auto& [name, op] : cases) {
        EXPECT_EQ(get_builtin_op(name), op) << name;
        EXPECT_EQ(std::string(builtin_op_name(op)), name)
            << "builtin_op_name must round-trip through get_builtin_op";
        auto e = parse(name + "(x)");
        const auto* fc = expr_cast<FuncCall>(e);
        ASSERT_NE(fc, nullptr) << name;
        EXPECT_EQ(fc->func_id, op) << name;
    }
}

// Aliases emitted by the oracles must map to the same builtin.
TEST_F(A43NonelementaryBuiltinsTest, OracleAliasesMapToSameBuiltin) {
    EXPECT_EQ(get_builtin_op("ei"), BuiltinOp::ExpIntegralEi);
    EXPECT_EQ(get_builtin_op("expintegral_ei"), BuiltinOp::ExpIntegralEi);
    EXPECT_EQ(get_builtin_op("expintegral_si"), BuiltinOp::SinIntegral);
    EXPECT_EQ(get_builtin_op("expintegral_ci"), BuiltinOp::CosIntegral);
    EXPECT_EQ(get_builtin_op("logarithmic_integral"), BuiltinOp::LogIntegral);
    EXPECT_EQ(get_builtin_op("Li2"), BuiltinOp::Dilog);
}

// Spec §3 — every derivative, checked against the closed form.
TEST_F(A43NonelementaryBuiltinsTest, DerivativesMatchSpec) {
    const std::vector<std::pair<std::string, std::string>> cases{
        {"Ei(x)", "exp(x)/x"},
        {"Si(x)", "sin(x)/x"},
        {"Ci(x)", "cos(x)/x"},
        {"Shi(x)", "sinh(x)/x"},
        {"Chi(x)", "cosh(x)/x"},
        {"li(x)", "1/ln(x)"},
        {"dilog(x)", "-ln(1-x)/x"},
        {"erfi(x)", "2/sqrt(pi) * exp(x^2)"},
    };
    for (const auto& [f, expected] : cases) {
        auto D = calculus::diff(parse(f), x, 1U, ctx);
        ASSERT_TRUE(D.is_ok()) << f << ": " << D.error().message;
        auto eq = symbolic::mathematically_equal(D.value(), parse(expected), ctx);
        ASSERT_TRUE(eq.is_ok()) << f;
        EXPECT_TRUE(eq.value()) << f << " -> expected " << expected;
    }
}

// Chain rule must apply through the new builtins.
TEST_F(A43NonelementaryBuiltinsTest, DerivativesRespectChainRule) {
    const std::vector<std::pair<std::string, std::string>> cases{
        {"Ei(2*x)", "exp(2*x)/x"},           // (e^{2x}/(2x))*2
        {"Si(x^2)", "2*sin(x^2)/x"},         // (sin(x^2)/x^2)*2x
        {"erfi(3*x)", "6/sqrt(pi) * exp(9*x^2)"},
    };
    for (const auto& [f, expected] : cases) {
        auto D = calculus::diff(parse(f), x, 1U, ctx);
        ASSERT_TRUE(D.is_ok()) << f << ": " << D.error().message;
        auto eq = symbolic::mathematically_equal(D.value(), parse(expected), ctx);
        ASSERT_TRUE(eq.is_ok()) << f;
        EXPECT_TRUE(eq.value()) << f << " -> expected " << expected;
    }
}

// The new builtins must not perturb the ordering of existing terms: their LPO
// precedence sits above Exp deliberately (spec rationale), so a product must
// still simplify to a stable canonical form and be idempotent.
TEST_F(A43NonelementaryBuiltinsTest, OrderingIsStableAndIdempotent) {
    for (const char* src : {"Ei(x) * exp(x) * x", "erfi(x) + Si(x) + x^2",
                            "li(x) * dilog(x)"}) {
        auto once = ctx.simplify(parse(src));
        ASSERT_TRUE(once.is_ok()) << src;
        auto twice = ctx.simplify(once.value());
        ASSERT_TRUE(twice.is_ok()) << src;
        EXPECT_TRUE(structural_equal(once.value(), twice.value()))
            << src << ": simplify must be idempotent";
    }
}

// Spec §4 — the reduction identities must be recognised by mathematically_equal
// (li → Ei(ln), Shi/Chi → Ei combinations, erfi → erf of an imaginary argument).
// Without a single representative the structural comparison forks exactly as in
// A44 (x! vs Γ(x+1)).
TEST_F(A43NonelementaryBuiltinsTest, ReductionIdentitiesAreRecognised) {
    const std::vector<std::pair<std::string, std::string>> identities{
        {"li(x)", "Ei(ln(x))"},
        {"Shi(x)", "(Ei(x) - Ei(-x))/2"},
        {"Chi(x)", "(Ei(x) + Ei(-x))/2"},
        {"erfi(x)", "-i*erf(i*x)"},
        // The identities must survive inside a larger expression, not only at
        // the root: the normaliser is applied bottom-up.
        {"sqrt(pi)/2 * erfi(x)", "-i*sqrt(pi)*erf(i*x)/2"},
        {"li(x^2) + 1", "Ei(ln(x^2)) + 1"},
    };
    for (const auto& [lhs, rhs] : identities) {
        auto eq = symbolic::mathematically_equal(parse(lhs), parse(rhs), ctx);
        ASSERT_TRUE(eq.is_ok()) << lhs << " vs " << rhs;
        EXPECT_TRUE(eq.value()) << lhs << " != " << rhs;
    }
}

// The normaliser must not turn the comparison into a rubber stamp: pairs that
// are genuinely different must still come back false. Si/Ci are deliberately
// left alone (their identities toward Ei carry a spurious `i` on the reals),
// so a Si/Ci pair must not be collapsed either.
TEST_F(A43NonelementaryBuiltinsTest, ReductionIdentitiesDoNotOverreach) {
    const std::vector<std::pair<std::string, std::string>> non_identities{
        {"li(x)", "Ei(x)"},
        {"Shi(x)", "Chi(x)"},
        {"erfi(x)", "erf(x)"},
        {"Si(x)", "Ci(x)"},
        {"Si(x)", "Ei(x)"},
    };
    for (const auto& [lhs, rhs] : non_identities) {
        auto eq = symbolic::mathematically_equal(parse(lhs), parse(rhs), ctx);
        ASSERT_TRUE(eq.is_ok()) << lhs << " vs " << rhs;
        EXPECT_FALSE(eq.value()) << lhs << " must not compare equal to " << rhs;
    }
}

}  // namespace
