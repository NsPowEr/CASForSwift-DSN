#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <map>

namespace cas::algebra {

namespace {

Result<std::vector<MultivariateTerm>> parse_multivariate_impl(ExprPtr expr, symbolic::CASContext& ctx) {
    if (!expr) {
        return ok(std::vector<MultivariateTerm>{});
    }

    if (const auto* integer = expr_cast<IntegerLit>(expr)) {
        if (integer->value.is_zero()) {
            return ok(std::vector<MultivariateTerm>{});
        }
        return ok(std::vector<MultivariateTerm>{{
            .coefficient = integer->value,
            .factors = {}
        }});
    }

    if (const auto* rational = expr_cast<RationalLit>(expr)) {
        if (rational->denominator == BigInt(1)) {
            if (rational->numerator.is_zero()) {
                return ok(std::vector<MultivariateTerm>{});
            }
            return ok(std::vector<MultivariateTerm>{{
                .coefficient = rational->numerator,
                .factors = {}
            }});
        }
        return fail<std::vector<MultivariateTerm>>(make_error(CASErrorKind::Unimplemented, "MultivariatePolynomial supporta solo coefficienti interi"));
    }

    if (const auto* symbol = expr_cast<Symbol>(expr)) {
        return ok(std::vector<MultivariateTerm>{{
            .coefficient = BigInt(1),
            .factors = {{*symbol, 1U}}
        }});
    }

    if (const auto* unary = expr_cast<Unary>(expr)) {
        if (unary->op == UnaryOp::Neg) {
            auto inner = parse_multivariate_impl(unary->operand, ctx);
            if (inner.is_error()) return inner;
            for (auto& term : inner.value()) {
                term.coefficient = -term.coefficient;
            }
            return inner;
        }
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        switch (binary->op) {
            case BinaryOp::Add: {
                auto left = parse_multivariate_impl(binary->left, ctx);
                if (left.is_error()) return left;
                auto right = parse_multivariate_impl(binary->right, ctx);
                if (right.is_error()) return right;
                left.value().insert(left.value().end(), right.value().begin(), right.value().end());
                return left;
            }
            case BinaryOp::Sub: {
                auto left = parse_multivariate_impl(binary->left, ctx);
                if (left.is_error()) return left;
                auto right = parse_multivariate_impl(binary->right, ctx);
                if (right.is_error()) return right;
                for (auto& term : right.value()) {
                    term.coefficient = -term.coefficient;
                }
                left.value().insert(left.value().end(), right.value().begin(), right.value().end());
                return left;
            }
            case BinaryOp::Mul: {
                auto left = parse_multivariate_impl(binary->left, ctx);
                if (left.is_error()) return left;
                auto right = parse_multivariate_impl(binary->right, ctx);
                if (right.is_error()) return right;
                
                std::vector<MultivariateTerm> result;
                result.reserve(left.value().size() * right.value().size());
                for (const auto& lt : left.value()) {
                    for (const auto& rt : right.value()) {
                        auto factors = lt.factors;
                        factors.insert(factors.end(), rt.factors.begin(), rt.factors.end());
                        result.push_back({
                            .coefficient = lt.coefficient * rt.coefficient,
                            .factors = std::move(factors)
                        });
                    }
                }
                return ok(std::move(result));
            }
            case BinaryOp::Pow: {
                auto base = parse_multivariate_impl(binary->left, ctx);
                if (base.is_error()) return base;
                auto exponent_res = poly_parse_nonnegative_integer_exponent(binary->right);
                if (exponent_res.is_error()) return fail<std::vector<MultivariateTerm>>(exponent_res.error());
                std::size_t exponent = exponent_res.value();
                
                if (exponent == 0) {
                    return ok(std::vector<MultivariateTerm>{{.coefficient = BigInt(1), .factors = {}}});
                }
                
                std::vector<MultivariateTerm> result = base.value();
                for (std::size_t i = 1; i < exponent; ++i) {
                    std::vector<MultivariateTerm> next;
                    for (const auto& t1 : result) {
                        for (const auto& t2 : base.value()) {
                            auto factors = t1.factors;
                            factors.insert(factors.end(), t2.factors.begin(), t2.factors.end());
                            next.push_back({
                                .coefficient = t1.coefficient * t2.coefficient,
                                .factors = std::move(factors)
                            });
                        }
                    }
                    result = std::move(next);
                }
                return ok(std::move(result));
            }
            default:
                break;
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<MultivariateTerm> result;
        for (ExprPtr term : sum->terms) {
            auto inner = parse_multivariate_impl(term, ctx);
            if (inner.is_error()) return inner;
            result.insert(result.end(), inner.value().begin(), inner.value().end());
        }
        return ok(std::move(result));
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        std::vector<MultivariateTerm> result{{.coefficient = BigInt(1), .factors = {}}};
        for (ExprPtr factor : product->factors) {
            auto inner = parse_multivariate_impl(factor, ctx);
            if (inner.is_error()) return inner;
            
            std::vector<MultivariateTerm> next;
            for (const auto& t1 : result) {
                for (const auto& t2 : inner.value()) {
                    auto factors = t1.factors;
                    factors.insert(factors.end(), t2.factors.begin(), t2.factors.end());
                    next.push_back({
                        .coefficient = t1.coefficient * t2.coefficient,
                        .factors = std::move(factors)
                    });
                }
            }
            result = std::move(next);
        }
        return ok(std::move(result));
    }

    return fail<std::vector<MultivariateTerm>>(make_error(CASErrorKind::Unimplemented, "Espressione non supportata in MultivariatePolynomial"));
}

} // namespace

Result<MultivariatePolynomial> parse_multivariate_polynomial(ExprPtr expr, symbolic::CASContext& ctx) {
    auto terms = parse_multivariate_impl(expr, ctx);
    if (terms.is_error()) return fail<MultivariatePolynomial>(terms.error());
    return ok(MultivariatePolynomial(std::move(terms.value())));
}

Result<ExprPtr> multivariate_to_expr(const MultivariatePolynomial& poly, symbolic::CASContext& ctx) {
    if (poly.is_zero()) {
        return ok(make_integer(ctx.arena(), 0));
    }
    
    std::vector<ExprPtr> sum_terms;
    for (const auto& term : poly.terms()) {
        auto term_expr = build_multivariate_monomial_expr(term, ctx);
        if (term_expr.is_error()) return term_expr;
        sum_terms.push_back(term_expr.value());
    }
    
    if (sum_terms.size() == 1) return ok(sum_terms[0]);
    return ok(ctx.arena().make<Sum>(std::move(sum_terms)));
}

} // namespace cas::algebra
