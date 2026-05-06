#include "cas/normal_form.hpp"
#include "cas/algebra.hpp"
#include <algorithm>

namespace cas::symbolic {

Result<std::map<Monomial, Rational>> collect_polynomial_terms(ExprPtr expr, CASContext& ctx) {
    auto expanded = algebra::expand(expr, ctx);
    if (expanded.is_error()) return fail<std::map<Monomial, Rational>>(expanded.error());

    std::map<Monomial, Rational> terms;

    std::function<Result<void>(ExprPtr, Rational)> collect;
    collect = [&](ExprPtr e, Rational current_coeff) -> Result<void> {
        if (auto sum = expr_cast<Sum>(e)) {
            for (auto t : sum->terms) {
                auto res = collect(t, current_coeff);
                if (res.is_error()) return res;
            }
            return ok();
        }
        if (auto bin = expr_cast<Binary>(e)) {
            if (bin->op == BinaryOp::Add) {
                auto res = collect(bin->left, current_coeff);
                if (res.is_error()) return res;
                return collect(bin->right, current_coeff);
            }
            if (bin->op == BinaryOp::Sub) {
                auto res = collect(bin->left, current_coeff);
                if (res.is_error()) return res;
                return collect(bin->right, current_coeff * Rational(BigInt(-1)));
            }
        }
        if (auto neg = expr_cast<Unary>(e); neg && neg->op == UnaryOp::Neg) {
            return collect(neg->operand, current_coeff * Rational(BigInt(-1)));
        }

        // Single term (Product or factor)
        Rational coeff = current_coeff;
        Monomial m;

        auto add_factor = [&](ExprPtr f, auto& self) -> void {
            if (auto int_lit = expr_cast<IntegerLit>(f)) {
                coeff *= Rational(int_lit->value);
            } else if (auto rat_lit = expr_cast<RationalLit>(f)) {
                coeff *= Rational(rat_lit->numerator, rat_lit->denominator);
            } else if (auto neg_f = expr_cast<Unary>(f); neg_f && neg_f->op == UnaryOp::Neg) {
                coeff *= Rational(BigInt(-1));
                self(neg_f->operand, self);
            } else if (auto prod_f = expr_cast<Product>(f)) {
                for (auto fact : prod_f->factors) self(fact, self);
            } else if (auto mul_f = expr_cast<Binary>(f); mul_f && mul_f->op == BinaryOp::Mul) {
                self(mul_f->left, self);
                self(mul_f->right, self);
            } else if (auto pow = expr_cast<Binary>(f); pow && pow->op == BinaryOp::Pow) {
                if (auto int_exp = expr_cast<IntegerLit>(pow->right)) {
                    if (!int_exp->value.is_negative() && !int_exp->value.is_zero()) {
                        m.factors.push_back({pow->left, static_cast<unsigned int>(int_exp->value.to_u64())});
                        return;
                    }
                }
                m.factors.push_back({f, 1});
            } else {
                m.factors.push_back({f, 1});
            }
        };

        add_factor(e, add_factor);

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

    auto res = collect(expanded.value(), Rational(BigInt(1)));
    if (res.is_error()) return fail<std::map<Monomial, Rational>>(res.error());

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

Result<ExprPtr> transcendental_normal_form(ExprPtr expr, CASContext& ctx) {
    if (!expr) return ok(expr);

    auto recurse = [&](ExprPtr e) -> Result<ExprPtr> {
        return transcendental_normal_form(e, ctx);
    };

    if (auto* f = expr_cast<FuncCall>(expr)) {
        if (f->func_id == BuiltinOp::Ln && f->args.size() == 1U) {
            ExprPtr arg = f->args[0];
            // ln(1) → 0
            if (auto* i = expr_cast<IntegerLit>(arg); i && i->value == BigInt(1))
                return ok(ctx.arena().make<IntegerLit>(BigInt(0)));
            // ln(exp(x)) → x
            if (auto* inner = expr_cast<FuncCall>(arg); inner && inner->func_id == BuiltinOp::Exp && inner->args.size() == 1U)
                return recurse(inner->args[0]);
            // ln(a * b) [Binary] → ln(a) + ln(b)
            if (auto* bin = expr_cast<Binary>(arg); bin && bin->op == BinaryOp::Mul) {
                auto la = recurse(ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{bin->left}));
                if (la.is_error()) return la;
                auto lb = recurse(ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{bin->right}));
                if (lb.is_error()) return lb;
                return ok(ctx.arena().make<Binary>(BinaryOp::Add, la.value(), lb.value()));
            }
            // ln(a * b * ...) [Product] → ln(a) + ln(b) + ...
            if (auto* prod = expr_cast<Product>(arg); prod && prod->factors.size() >= 2U) {
                std::vector<ExprPtr> ln_terms;
                ln_terms.reserve(prod->factors.size());
                for (auto fac : prod->factors) {
                    auto r = recurse(ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{fac}));
                    if (r.is_error()) return r;
                    ln_terms.push_back(r.value());
                }
                return ok(ctx.arena().make<Sum>(std::move(ln_terms)));
            }
            // ln(a / b) → ln(a) - ln(b)
            if (auto* bin = expr_cast<Binary>(arg); bin && bin->op == BinaryOp::Div) {
                auto la = recurse(ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{bin->left}));
                if (la.is_error()) return la;
                auto lb = recurse(ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{bin->right}));
                if (lb.is_error()) return lb;
                return ok(ctx.arena().make<Binary>(BinaryOp::Sub, la.value(), lb.value()));
            }
            // ln(a^n) → n * ln(a), n integer
            if (auto* bin = expr_cast<Binary>(arg); bin && bin->op == BinaryOp::Pow && expr_is<IntegerLit>(bin->right)) {
                auto la = recurse(ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{bin->left}));
                if (la.is_error()) return la;
                return ok(ctx.arena().make<Binary>(BinaryOp::Mul, bin->right, la.value()));
            }
            // recurse into arg
            auto r = recurse(arg);
            if (r.is_error()) return r;
            if (r.value() == arg) return ok(expr);
            return ok(ctx.arena().make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{r.value()}));
        }
        if (f->func_id == BuiltinOp::Exp && f->args.size() == 1U) {
            ExprPtr arg = f->args[0];
            // exp(0) → 1
            if (auto* i = expr_cast<IntegerLit>(arg); i && i->value == BigInt(0))
                return ok(ctx.arena().make<IntegerLit>(BigInt(1)));
            // exp(ln(x)) → x
            if (auto* inner = expr_cast<FuncCall>(arg); inner && inner->func_id == BuiltinOp::Ln && inner->args.size() == 1U)
                return recurse(inner->args[0]);
            auto r = recurse(arg);
            if (r.is_error()) return r;
            if (r.value() == arg) return ok(expr);
            return ok(ctx.arena().make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{r.value()}));
        }
        // other FuncCall: recurse into args
        bool changed = false;
        std::vector<ExprPtr> new_args;
        new_args.reserve(f->args.size());
        for (auto a : f->args) {
            auto r = recurse(a);
            if (r.is_error()) return r;
            new_args.push_back(r.value());
            if (r.value() != a) changed = true;
        }
        if (!changed) return ok(expr);
        return ok(ctx.arena().make<FuncCall>(f->func_id, std::move(new_args)));
    }

    if (auto* b = expr_cast<Binary>(expr)) {
        auto l = recurse(b->left);
        if (l.is_error()) return l;
        auto r = recurse(b->right);
        if (r.is_error()) return r;
        if (l.value() == b->left && r.value() == b->right) return ok(expr);
        return ok(ctx.arena().make<Binary>(b->op, l.value(), r.value()));
    }

    if (auto* s = expr_cast<Sum>(expr)) {
        bool changed = false;
        std::vector<ExprPtr> new_terms;
        new_terms.reserve(s->terms.size());
        for (auto t : s->terms) {
            auto r = recurse(t);
            if (r.is_error()) return r;
            new_terms.push_back(r.value());
            if (r.value() != t) changed = true;
        }
        if (!changed) return ok(expr);
        return ok(ctx.arena().make<Sum>(std::move(new_terms)));
    }

    if (auto* p = expr_cast<Product>(expr)) {
        bool changed = false;
        std::vector<ExprPtr> new_factors;
        new_factors.reserve(p->factors.size());
        for (auto fac : p->factors) {
            auto r = recurse(fac);
            if (r.is_error()) return r;
            new_factors.push_back(r.value());
            if (r.value() != fac) changed = true;
        }
        if (!changed) return ok(expr);
        return ok(ctx.arena().make<Product>(std::move(new_factors)));
    }

    if (auto* u = expr_cast<Unary>(expr)) {
        auto r = recurse(u->operand);
        if (r.is_error()) return r;
        if (r.value() == u->operand) return ok(expr);
        return ok(ctx.arena().make<Unary>(u->op, r.value()));
    }

    return ok(expr);
}

} // namespace cas::symbolic
