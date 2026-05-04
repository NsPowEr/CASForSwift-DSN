#include "cas/normal_form.hpp"
#include "cas/algebra.hpp"
#include <algorithm>

namespace cas::symbolic {

Result<std::map<Monomial, Rational>> collect_polynomial_terms(ExprPtr expr, CASContext& ctx) {
    auto expanded = algebra::expand(expr, ctx);
    if (expanded.is_error()) return fail<std::map<Monomial, Rational>>(expanded.error());

    std::map<Monomial, Rational> terms;

    auto add_term = [&](ExprPtr term_expr) -> Result<void> {
        Rational coeff(BigInt(1));
        Monomial m;

        auto add_factor = [&](ExprPtr f) {
            if (auto int_lit = expr_cast<IntegerLit>(f)) {
                coeff *= Rational(int_lit->value);
            } else if (auto rat_lit = expr_cast<RationalLit>(f)) {
                coeff *= Rational(rat_lit->numerator, rat_lit->denominator);
            } else if (auto pow = expr_cast<Binary>(f); pow && pow->op == cas::BinaryOp::Pow) {
                if (auto int_exp = expr_cast<IntegerLit>(pow->right)) {
                    if (!int_exp->value.is_negative() && !int_exp->value.is_zero()) {
                        m.factors.push_back({pow->left, static_cast<unsigned int>(int_exp->value.to_u64())});
                        return;
                    }
                }
                m.factors.push_back({f, 1});
            } else if (auto neg = expr_cast<Unary>(f); neg && neg->op == cas::UnaryOp::Neg) {
                coeff *= Rational(BigInt(-1));
                // Recurse into operand if needed or just add it as factor
                // Simplified: expand() should handle negations as *-1
                m.factors.push_back({neg->operand, 1});
            } else {
                m.factors.push_back({f, 1});
            }
        };

        if (auto prod = expr_cast<Product>(term_expr)) {
            for (auto f : prod->factors) {
                add_factor(f);
            }
        } else {
            add_factor(term_expr);
        }

        std::sort(m.factors.begin(), m.factors.end(), [](const auto& a, const auto& b) {
            return canonical_compare(a.first, b.first) < 0;
        });

        Monomial res_m;
        for (const auto& f : m.factors) {
            if (!res_m.factors.empty() && structural_equal(res_m.factors.back().first, f.first)) {
                res_m.factors.back().second += f.second;
            } else {
                res_m.factors.push_back(f);
            }
        }

        terms[res_m] += coeff;
        return ok();
    };

    if (auto sum = expr_cast<Sum>(expanded.value())) {
        for (auto t : sum->terms) {
            auto res = add_term(t);
            if (res.is_error()) return fail<std::map<Monomial, Rational>>(res.error());
        }
    } else {
        auto res = add_term(expanded.value());
        if (res.is_error()) return fail<std::map<Monomial, Rational>>(res.error());
    }

    return ok(terms);
}

Result<ExprPtr> polynomial_normal_form(ExprPtr expr, CASContext& ctx) {
    auto terms_res = collect_polynomial_terms(expr, ctx);
    if (terms_res.is_error()) return fail<ExprPtr>(terms_res.error());

    auto& terms = terms_res.value();
    std::vector<ExprPtr> sum_operands;

    for (const auto& [m, c] : terms) {
        if (c.numerator().is_zero()) continue;

        std::vector<ExprPtr> prod_operands;
        if (c.denominator() == BigInt(1)) {
            if (c.numerator() != BigInt(1) || m.factors.empty()) {
                prod_operands.push_back(ctx.arena().make<IntegerLit>(c.numerator()));
            }
        } else {
            prod_operands.push_back(ctx.arena().make<RationalLit>(c.numerator(), c.denominator()));
        }

        for (const auto& f : m.factors) {
            if (f.second == 0) continue;
            if (f.second == 1) {
                prod_operands.push_back(f.first);
            } else {
                prod_operands.push_back(ctx.arena().make<Binary>(BinaryOp::Pow, f.first, ctx.arena().make<IntegerLit>(BigInt(static_cast<long long>(f.second)))));
            }
        }

        if (prod_operands.empty()) {
            sum_operands.push_back(ctx.arena().make<IntegerLit>(BigInt(1)));
        } else if (prod_operands.size() == 1) {
            sum_operands.push_back(prod_operands[0]);
        } else {
            sum_operands.push_back(ctx.arena().make<Product>(std::move(prod_operands)));
        }
    }

    if (sum_operands.empty()) {
        return ok(ctx.arena().make<IntegerLit>(BigInt(0)));
    } else if (sum_operands.size() == 1) {
        return ok(sum_operands[0]);
    } else {
        std::sort(sum_operands.begin(), sum_operands.end(), [](const ExprPtr& a, const ExprPtr& b) {
            return canonical_compare(a, b) < 0;
        });
        return ok(ctx.arena().make<Sum>(std::move(sum_operands)));
    }
}

} // namespace cas::symbolic
