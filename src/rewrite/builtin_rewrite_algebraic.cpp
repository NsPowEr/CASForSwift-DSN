// builtin_rewrite_algebraic.cpp — Rational and algebraic rewrite rules.
#include "builtin_rewrite_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <vector>
#include <functional>

namespace cas::symbolic {

[[nodiscard]] static std::optional<Rational> get_exact_rational(ExprPtr e) {
    if (!e) return std::nullopt;
    if (const auto* il = expr_cast<IntegerLit>(e)) return Rational(il->value);
    if (const auto* rl = expr_cast<RationalLit>(e)) {
        auto val = Rational::make(rl->numerator, rl->denominator);
        if (val.is_ok()) return val.value();
    }
    if (const auto* un = expr_cast<Unary>(e)) {
        if (un->op == UnaryOp::Neg) {
            if (auto inner = get_exact_rational(un->operand)) return -*inner;
        }
    }
    return std::nullopt;
}

[[nodiscard]] static BigInt get_poly_lcm(const algebra::PolyExpr& poly) {
    BigInt common_lcm(1);
    for (ExprPtr coefficient : poly.coefficients()) {
        if (!coefficient) continue;
        if (const auto* il = expr_cast<IntegerLit>(coefficient); il && il->value.is_zero()) continue;
        if (const auto* rl = expr_cast<RationalLit>(coefficient); rl && rl->numerator.is_zero()) continue;
        auto rat = get_exact_rational(coefficient);
        if (!rat.has_value()) continue;
        if (!rat->denominator().is_zero()) {
            BigInt d = rat->denominator();
            common_lcm = (common_lcm * d) / gcd(common_lcm, d);
        }
    }
    return common_lcm;
}

Result<ExprPtr> try_rewrite_algebraic(
    ExprPtr expr,
    AstArena& arena,
    const Assumptions* assumptions,
    CASContext* context) {
    if (const auto* binary = expr_cast<Binary>(expr); binary != nullptr && binary->op == BinaryOp::Div) {
        if (context != nullptr) {
            // GCD Reduction for Rational Functions
            auto find_var = [](ExprPtr e) -> std::optional<Symbol> {
                std::optional<Symbol> found;
                std::function<void(ExprPtr)> collect = [&](ExprPtr node) {
                    if (found.has_value() || !node) return;
                    if (const auto* s = expr_cast<Symbol>(node)) {
                        found = *s;
                        return;
                    }
                    visit_expr(node, [&](const auto& n) {
                        using T = std::decay_t<decltype(n)>;
                        if constexpr (std::is_same_v<T, Unary>) collect(n.operand);
                        else if constexpr (std::is_same_v<T, Binary>) { collect(n.left); collect(n.right); }
                        else if constexpr (std::is_same_v<T, FuncCall>) { for (auto a : n.args) collect(a); }
                        else if constexpr (std::is_same_v<T, Sum>) { for (auto t : n.terms) collect(t); }
                        else if constexpr (std::is_same_v<T, Product>) { for (auto f : n.factors) collect(f); }
                    });
                };
                collect(e);
                return found;
            };

            auto var = find_var(expr);
            if (var.has_value()) {
                auto p1_res = algebra::parse_polynomial(binary->left, *var, *context);
                auto p2_res = algebra::parse_polynomial(binary->right, *var, *context);
                
                if (p1_res.is_ok() && p2_res.is_ok()) {
                    auto i1_res = algebra::poly_to_integer_poly(p1_res.value());
                    auto i2_res = algebra::poly_to_integer_poly(p2_res.value());
                    
                    if (i1_res.is_ok() && i2_res.is_ok()) {
                        auto gcd_poly = (context != nullptr)
                            ? algebra::gcd_integer_poly_dispatch(i1_res.value(), i2_res.value(), *context).gcd
                            : algebra::gcd_integer_poly_with_subresultant(i1_res.value(), i2_res.value()).gcd;
                        
                        if (gcd_poly.degree() > 0 || (gcd_poly.size() == 1 && gcd_poly[0] != BigInt(1))) {
                            auto exact_div = [](const algebra::IntPoly& num, const algebra::IntPoly& den) -> algebra::IntPoly {
                                if (den.is_zero()) return num;
                                if (num.is_zero()) return algebra::IntPoly{};
                                if (num.size() < den.size()) return algebra::IntPoly{};
                                
                                algebra::IntPoly q;
                                q.resize(num.size() - den.size() + 1, BigInt(0));
                                algebra::IntPoly r = num;
                                BigInt lc_den = den.leading_coeff();
                                
                                for (int i = static_cast<int>(num.size() - den.size()); i >= 0; --i) {
                                    std::size_t r_idx = static_cast<std::size_t>(i) + den.size() - 1;
                                    q[i] = r[r_idx] / lc_den;
                                    for (std::size_t j = 0; j < den.size(); ++j) {
                                        r[i + j] -= q[i] * den[j];
                                    }
                                }
                                algebra::normalize_integer_poly(q);
                                return q;
                            };

                            auto n_poly = exact_div(i1_res.value(), gcd_poly);
                            auto d_poly = exact_div(i2_res.value(), gcd_poly);

                            BigInt lcm1 = get_poly_lcm(p1_res.value());
                            BigInt lcm2 = get_poly_lcm(p2_res.value());
                            algebra::multiply_integer_coefficients_by_scalar(n_poly, lcm2);
                            algebra::multiply_integer_coefficients_by_scalar(d_poly, lcm1);
                            
                            auto n_expr = algebra::integer_coefficients_to_expr(n_poly, *var, *context);
                            auto d_expr = algebra::integer_coefficients_to_expr(d_poly, *var, *context);
                            
                            if (n_expr.is_ok() && d_expr.is_ok()) {
                                if (!d_poly.is_zero() && d_poly.degree() == 0 && d_poly[0] == BigInt(1)) return n_expr;
                                return ok(arena.make<Binary>(BinaryOp::Div, n_expr.value(), d_expr.value()));
                            }
                        }
                    }
                }
            }
        }
    }

    if (const auto* call = expr_cast<FuncCall>(expr); call != nullptr &&
        call->func_id == BuiltinOp::Sqrt && call->args.size() == 1U) {
        if (const auto* quotient = expr_cast<Binary>(call->args.front());
            quotient != nullptr && quotient->op == BinaryOp::Div &&
            expr_is_nonnegative_under_assumptions(quotient->left, assumptions) &&
            expr_is_positive_under_assumptions(quotient->right, assumptions)) {
            return ok(arena.make<Binary>(
                BinaryOp::Div,
                arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{quotient->left}),
                arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{quotient->right})));
        }

        if (const auto* product = expr_cast<Product>(call->args.front())) {
            std::vector<ExprPtr> numerators;
            std::vector<ExprPtr> denominators;
            for (ExprPtr factor : product->factors) {
                if (const auto* bin = expr_cast<Binary>(factor); 
                    bin && bin->op == BinaryOp::Pow && 
                    exact_expr_is_negative(bin->right)) {
                    denominators.push_back(bin->left);
                } else {
                    numerators.push_back(factor);
                }
            }

            if (!denominators.empty()) {
                ExprPtr num = numerators.empty() ? arena.make<IntegerLit>(BigInt(1)) : 
                             (numerators.size() == 1 ? numerators[0] : arena.make<Product>(std::move(numerators)));
                ExprPtr den = denominators.size() == 1 ? denominators[0] : arena.make<Product>(std::move(denominators));
                return ok(arena.make<Binary>(
                    BinaryOp::Div,
                    arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{num}),
                    arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{den})));
            }

            std::vector<ExprPtr> terms;
            terms.reserve(product->factors.size());
            for (ExprPtr factor : product->factors) {
                if (!expr_is_nonnegative_under_assumptions(factor, assumptions)) {
                    terms.clear();
                    break;
                }
                terms.push_back(arena.make<FuncCall>("sqrt", std::vector<ExprPtr>{factor}));
            }
            if (!terms.empty()) {
                return ok(arena.make<Product>(std::move(terms)));
            }
        }
    }

    return ok(expr);
}

}  // namespace cas::symbolic
