#include <gtest/gtest.h>
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast_debug.hpp"
#include "cas/formatter.hpp"

namespace cas {
namespace {

class ComplexSimplificationTest : public ::testing::Test {
protected:
    void SetUp() override {
        arena = std::make_unique<AstArena>();
        ctx = std::make_unique<symbolic::CASContext>();
    }

    ExprPtr simplify_expr(ExprPtr e) {
        auto res = symbolic::simplify(e, *ctx);
        if (res.is_error()) return nullptr;
        return res.value();
    }

    std::string to_string(ExprPtr e) {
        if (!e) return "null";
        return formatter::TextFormatter().format(e);
    }

    ExprPtr make_complex(long long re_n, long long re_d, long long im_n, long long im_d) {
        return arena->make<ComplexLit>(BigInt(re_n), BigInt(re_d), BigInt(im_n), BigInt(im_d));
    }

    ExprPtr make_int(long long n) {
        return arena->make<IntegerLit>(BigInt(n));
    }

    ExprPtr make_sum(std::vector<ExprPtr> terms) {
        return arena->make<Sum>(std::move(terms));
    }

    ExprPtr make_prod(std::vector<ExprPtr> factors) {
        return arena->make<Product>(std::move(factors));
    }

    std::unique_ptr<AstArena> arena;
    std::unique_ptr<symbolic::CASContext> ctx;
};

TEST_F(ComplexSimplificationTest, BasicSum) {
    // (1 + I) + (2 + 2I) = 3 + 3I
    auto c1 = make_complex(1, 1, 1, 1);
    auto c2 = make_complex(2, 1, 2, 1);
    auto sum = make_sum({c1, c2});
    auto res = simplify_expr(sum);
    
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(to_string(res), "3 + 3 * I");
}

TEST_F(ComplexSimplificationTest, SumCancellation) {
    // (1 + I) + (2 - I) = 3
    auto c1 = make_complex(1, 1, 1, 1);
    auto c2 = make_complex(2, 1, -1, 1);
    auto sum = make_sum({c1, c2});
    auto res = simplify_expr(sum);
    
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(to_string(res), "3");
}

TEST_F(ComplexSimplificationTest, BasicProduct) {
    // (1 + I) * (1 - I) = 1 - I^2 = 2
    auto c1 = make_complex(1, 1, 1, 1);
    auto c2 = make_complex(1, 1, -1, 1);
    auto prod = make_prod({c1, c2});
    auto res = simplify_expr(prod);
    
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(to_string(res), "2");
}

TEST_F(ComplexSimplificationTest, ProductWithImaginaryUnit) {
    // (1 + I) * I = I + I^2 = -1 + I
    auto c1 = make_complex(1, 1, 1, 1);
    auto i_unit = arena->make<Constant>(MathConstant::I);
    auto prod = make_prod({c1, i_unit});
    auto res = simplify_expr(prod);
    
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(to_string(res), "-1 + I");
}

TEST_F(ComplexSimplificationTest, ComplexDistribution) {
    // (1 + I) * (x + 1) = (1+I)x + (1+I)
    auto c1 = make_complex(1, 1, 1, 1);
    auto x = arena->make<Symbol>("x");
    auto sum_x_1 = make_sum({x, make_int(1)});
    auto prod = make_prod({c1, sum_x_1});
    auto res = simplify_expr(prod);
    
    ASSERT_NE(res, nullptr);
    EXPECT_TRUE(to_string(res).find("1 + I") != std::string::npos);
}

} // namespace
} // namespace cas
