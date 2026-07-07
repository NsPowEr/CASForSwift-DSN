// L2-01: Unit tests for ODE Frobenius series expansion around x_0.
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/lexer.hpp"
#include "cas/ode.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <string>

namespace cas::test {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

bool contains_ln_subexpr(ExprPtr e, const std::string& var) {
    if (!e) return false;
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        if (fc->name == "ln" || fc->func_id == BuiltinOp::Ln) {
            return true;
        }
    }
    switch (e->kind) {
    case ExprKind::Unary:
        return contains_ln_subexpr(expr_ref<Unary>(e).operand, var);
    case ExprKind::Binary: {
        const auto& b = expr_ref<Binary>(e);
        return contains_ln_subexpr(b.left, var) || contains_ln_subexpr(b.right, var);
    }
    case ExprKind::Sum:
        for (auto t : expr_ref<Sum>(e).terms)
            if (contains_ln_subexpr(t, var)) return true;
        return false;
    case ExprKind::Product:
        for (auto f : expr_ref<Product>(e).factors)
            if (contains_ln_subexpr(f, var)) return true;
        return false;
    case ExprKind::FuncCall:
        for (auto a : expr_ref<FuncCall>(e).args)
            if (contains_ln_subexpr(a, var)) return true;
        return false;
    default:
        return false;
    }
}

}  // namespace

class OdeFrobeniusTest : public ::testing::Test {
protected:
    void SetUp() override {
        ctx = std::make_unique<symbolic::CASContext>();
        ctx->set_timeout(std::chrono::seconds(30));
    }

    [[nodiscard]] ExprPtr E(const std::string& src) {
        auto r = parse_expr(src, ctx->arena());
        EXPECT_TRUE(r.is_ok()) << "parse: " << src
                               << " err=" << (r.is_error() ? r.error().message : std::string{});
        return r.is_ok() ? r.value() : ExprPtr{};
    }

    std::unique_ptr<symbolic::CASContext> ctx;
};

// 1. Bessel of order 1 (nu=1) at x_0=0:
//    x^2 y'' + x y' + (x^2 - 1) y = 0
TEST_F(OdeFrobeniusTest, BesselOrder1AtZero) {
    ExprPtr a2 = E("x^2");
    ExprPtr a1 = E("x");
    ExprPtr a0 = E("x^2 - 1");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius(a2, a1, a0, x, E("0"), 4U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});

    ExprPtr y = result.value();
    EXPECT_TRUE(contains_ln_subexpr(y, "x"))
        << "expected log branch for Bessel nu=1 in: " << debug_print(y);
}

// 2. Legendre differential equation at x_0=1:
//    (1 - x^2) y'' - 2x y' + 2 y = 0   (l = 1 => l(l+1) = 2)
//    Singular point at x_0=1 (regular singular point with double indicial root r=0).
TEST_F(OdeFrobeniusTest, LegendreAtX1) {
    ExprPtr a2 = E("1 - x^2");
    ExprPtr a1 = E("-2*x");
    ExprPtr a0 = E("2");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius(a2, a1, a0, x, E("1"), 4U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});

    ExprPtr y = result.value();
    EXPECT_TRUE(contains_ln_subexpr(y, "x"))
        << "expected log branch for Legendre at x_0=1 in: " << debug_print(y);
}

// 3. Shifting test around x_0 != 0 (e.g. Euler shifted around x_0=2):
//    (x-2)^2 y'' - 6 y = 0
TEST_F(OdeFrobeniusTest, ShiftedEulerAtX2) {
    ExprPtr a2 = E("(x - 2)^2");
    ExprPtr a1 = E("0");
    ExprPtr a0 = E("-6");
    Symbol x("x");

    auto result = calculus::solve_ode_frobenius(a2, a1, a0, x, E("2"), 4U, *ctx);
    ASSERT_TRUE(result.is_ok()) << "err=" << (result.is_error() ? result.error().message : std::string{});
}

}  // namespace cas::test
