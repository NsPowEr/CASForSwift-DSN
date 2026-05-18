#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

// ── sqrt helpers ─────────────────────────────────────────────────────────────

// Extract perfect-square factor from n: returns {k, m} with n = k²·m, m squarefree.
[[nodiscard]] static std::pair<BigInt, BigInt> extract_square_factor(BigInt n) {
    BigInt k(1);
    BigInt i(2);
    while (i * i <= n) {
        BigInt i2 = i * i;
        while ((n % i2).is_zero()) {
            k = k * i;
            n = n / i2;
        }
        i = i + BigInt(1);
    }
    return {k, n};
}

// sqrt(r) for rational r ≥ 0: extract perfect-square factors.
// Returns k·sqrt(m), k rational, m squarefree.
[[nodiscard]] static Result<ExprPtr> simplify_rational_sqrt(const Rational& r, AstArena& arena) {
    const BigInt& p = r.numerator();
    const BigInt& q = r.denominator();
    if (p.is_zero()) return ok(arena.make<IntegerLit>(BigInt(0)));
    auto [p_out, p_rem] = extract_square_factor(p);
    auto [q_out, q_rem] = extract_square_factor(q);
    BigInt final_radicand = p_rem * q_rem;
    Rational coeff(p_out, q_out * q_rem);

    ExprPtr coeff_expr;
    if (coeff.denominator() == BigInt(1)) {
        coeff_expr = arena.make<IntegerLit>(coeff.numerator());
    } else {
        coeff_expr = arena.make<RationalLit>(coeff.numerator(), coeff.denominator());
    }

    if (final_radicand == BigInt(1)) return ok(coeff_expr);

    ExprPtr sqrt_expr = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{arena.make<IntegerLit>(final_radicand)});
    if (coeff.numerator() == BigInt(1) && coeff.denominator() == BigInt(1))
        return ok(sqrt_expr);
    return ok(arena.make<Binary>(BinaryOp::Mul, coeff_expr, sqrt_expr));
}

[[nodiscard]] static BigInt integer_sqrt(const BigInt& n) {
    if (n.is_zero()) return BigInt(0);
    static const BigInt one(1);
    if (n == one) return one;
    BigInt x = one.shift_left_bits((n.bit_length() + 1) / 2);
    while (true) {
        BigInt y = (x + n / x) / BigInt(2);
        if (y >= x) return x;
        x = std::move(y);
    }
}

// ── Simplifier::simplify_funcall_exp_log_sqrt ─────────────────────────────────

Result<ExprPtr> Simplifier::simplify_funcall_exp_log_sqrt(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before) {

    if (op == BuiltinOp::Exp && args.size() == 1U) {
        if (is_zero_expr(args.front()))
            return traced_result(RuleId::SimplifyExpZero, target_before, make_integer(arena_, BigInt(1)));
        if (is_one_expr(args.front()))
            return traced_result(RuleId::SimplifyExpOne, target_before, make_constant(arena_, MathConstant::E));
        if (const auto* ln_call = expr_cast<FuncCall>(args.front());
            ln_call && ln_call->func_id == BuiltinOp::Ln && ln_call->args.size() == 1U)
            return ok(ln_call->args[0]);
        if (is_constant_expr(args.front(), MathConstant::Infinity))
            return traced_result(RuleId::Unknown, target_before, make_constant(arena_, MathConstant::Infinity));
        if (expr_is<Unary>(args.front())
            && expr_ref<Unary>(args.front()).op == UnaryOp::Neg
            && is_constant_expr(expr_ref<Unary>(args.front()).operand, MathConstant::Infinity))
            return traced_result(RuleId::Unknown, target_before, make_integer(arena_, BigInt(0)));
        if (const auto* sum = expr_cast<Sum>(args.front())) {
            std::vector<ExprPtr> factors;
            for (ExprPtr term : sum->terms)
                factors.push_back(arena_.make<FuncCall>(BuiltinOp::Exp,
                    std::vector<ExprPtr>{term}));
            auto rewritten = simplify_product_factors(factors, arena_.make<Product>(factors));
            if (rewritten.is_ok()) {
                append_trace(RuleId::SimplifyExpSum, target_before, rewritten.value());
                return rewritten;
            }
        }
    }

    if (op == BuiltinOp::Ln && args.size() == 1U) {
        if (is_zero_expr(args.front()))
            return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "ln(0) is undefined"));
        if (is_one_expr(args.front()))
            return traced_result(RuleId::SimplifyLnOne, target_before, make_integer(arena_, BigInt(0)));
        if (is_constant_expr(args.front(), MathConstant::E))
            return traced_result(RuleId::SimplifyLnE, target_before, make_integer(arena_, BigInt(1)));
        if (is_constant_expr(args.front(), MathConstant::Infinity))
            return traced_result(RuleId::Unknown, target_before, make_constant(arena_, MathConstant::Infinity));

        // ln(e^x) -> x
        if (const auto* power = expr_cast<Binary>(args.front());
            power != nullptr && power->op == BinaryOp::Pow
            && is_constant_expr(power->left, MathConstant::E))
            return traced_result(RuleId::SimplifyLnExp, target_before, power->right);
        // ln(exp(x)) -> x
        if (const auto* exp_call = expr_cast<FuncCall>(args.front());
            exp_call && exp_call->func_id == BuiltinOp::Exp && exp_call->args.size() == 1U)
            return traced_result(RuleId::SimplifyLnExp, target_before, exp_call->args[0]);

        // ln(a*b) -> ln(a) + ln(b) for a,b > 0
        if (const auto* prod = expr_cast<Product>(args.front())) {
            bool all_pos = true;
            for (auto f : prod->factors) if (!is_known_positive(f)) { all_pos = false; break; }
            if (all_pos) {
                std::vector<ExprPtr> ln_factors;
                for (auto f : prod->factors) {
                    auto res = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln,
                        std::vector<ExprPtr>{f}));
                    if (res.is_error()) return res;
                    ln_factors.push_back(res.value());
                }
                return simplify_expr(arena_.make<Sum>(std::move(ln_factors)));
            }
        }
        // ln(sqrt(x)) = (1/2)*ln(x) — identità esatta
        if (const auto* sqrt_call = expr_cast<FuncCall>(args.front());
            sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt
            && sqrt_call->args.size() == 1U) {
            ExprPtr half = make_rational(arena_, Rational(BigInt(1), BigInt(2)));
            ExprPtr ln_inner = arena_.make<FuncCall>(BuiltinOp::Ln, sqrt_call->args);
            return simplify_expr(arena_.make<Binary>(BinaryOp::Mul, half, ln_inner));
        }
        // Branch cut principal value: ln(-x) = ln(x) + I·π for x > 0
        if (const auto* neg = expr_cast<Unary>(args.front());
            neg != nullptr && neg->op == UnaryOp::Neg) {
            ExprPtr inner = neg->operand;
            bool inner_pos = is_known_positive(inner) || is_constant_expr(inner, MathConstant::E);
            if (!inner_pos) {
                LiteralRational rat;
                auto ex = try_get_exact_rational(inner, rat);
                if (ex.is_ok() && ex.value()
                    && !rat.value.numerator().is_negative()
                    && !rat.value.numerator().is_zero())
                    inner_pos = true;
            }
            if (inner_pos) {
                ExprPtr ln_inner;
                if (is_one_expr(inner)) {
                    ln_inner = make_integer(arena_, BigInt(0));
                } else {
                    auto r = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln,
                        std::vector<ExprPtr>{inner}));
                    if (r.is_error()) return r;
                    ln_inner = r.value();
                }
                ExprPtr i_pi = arena_.make<Binary>(BinaryOp::Mul,
                    make_constant(arena_, MathConstant::I),
                    make_constant(arena_, MathConstant::Pi));
                return simplify_expr(arena_.make<Binary>(BinaryOp::Add, ln_inner, i_pi));
            }
        }
        // Complex principal branch: ln(I) = I·π/2
        if (const auto* c = expr_cast<Constant>(args.front()); c && c->value == MathConstant::I) {
            ExprPtr i_pi = arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Constant>(MathConstant::I),
                arena_.make<Constant>(MathConstant::Pi)});
            return simplify_expr(arena_.make<Binary>(
                BinaryOp::Div, i_pi, make_integer(arena_, BigInt(2))));
        }
        // ln(-1) = I·π  (also covered by ln(-x) above, but explicit for clarity)
        if (const auto* il = expr_cast<IntegerLit>(args.front()); il && il->value == BigInt(-1)) {
            return simplify_expr(arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Constant>(MathConstant::I),
                arena_.make<Constant>(MathConstant::Pi)}));
        }
        // L2-08: ln(a + b·I) = ln|z| + I·arg(z)  (principal branch).
        // Triggered only when the argument is a Sum containing an imaginary term b·I.
        // abs(z) and arg(z) dispatch back to simplify_funcall_complex which already
        // handles the quadrant logic — no duplication.
        if (const auto* sum_z = expr_cast<Sum>(args.front())) {
            auto is_i_unit = [](ExprPtr e) -> bool {
                const auto* c = expr_cast<Constant>(e);
                return c != nullptr && c->value == MathConstant::I;
            };
            auto extract_imag_coeff = [&](ExprPtr term) -> ExprPtr {
                if (is_i_unit(term)) return make_integer(arena_, BigInt(1));
                if (const auto* prod = expr_cast<Product>(term)) {
                    bool found_i = false;
                    std::vector<ExprPtr> others;
                    for (ExprPtr f : prod->factors) {
                        if (!found_i && is_i_unit(f)) { found_i = true; continue; }
                        others.push_back(f);
                    }
                    if (found_i) {
                        if (others.empty()) return make_integer(arena_, BigInt(1));
                        if (others.size() == 1U) return others[0];
                        return arena_.make<Product>(std::move(others));
                    }
                }
                return nullptr;
            };
            ExprPtr real_part = nullptr, imag_part = nullptr;
            for (ExprPtr term : sum_z->terms) {
                ExprPtr b = extract_imag_coeff(term);
                if (b && !imag_part) imag_part = b;
                else if (!b && !real_part) real_part = term;
            }
            if (real_part && imag_part) {
                // ln(a + b·I) = ln(abs(a+b·I)) + I·arg(a+b·I)
                ExprPtr z = args.front();
                ExprPtr abs_z = arena_.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{z});
                ExprPtr arg_z = arena_.make<FuncCall>(BuiltinOp::Arg, std::vector<ExprPtr>{z});
                ExprPtr ln_abs = arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{abs_z});
                ExprPtr i_arg = arena_.make<Product>(std::vector<ExprPtr>{
                    arena_.make<Constant>(MathConstant::I), arg_z});
                ExprPtr result = arena_.make<Sum>(std::vector<ExprPtr>{ln_abs, i_arg});
                return simplify_expr(result);
            }
        }
    }

    if (op == BuiltinOp::Sqrt && args.size() == 1U) {
        // Denesting sqrt(a + b*sqrt(c))
        if (const auto* sum = expr_cast<Sum>(args.front()); sum && sum->terms.size() == 2) {
            LiteralRational rat_a, rat_b, rat_c;
            ExprPtr a_ptr = nullptr, b_ptr = nullptr, c_ptr = nullptr;
            for (auto term : sum->terms) {
                if (auto ex = try_get_exact_rational(term, rat_a); ex.is_ok() && ex.value()) {
                    a_ptr = term;
                } else if (const auto* prod = expr_cast<Product>(term)) {
                    Rational b_coeff(1);
                    ExprPtr c_val = nullptr;
                    bool found_sqrt = false;
                    for (ExprPtr f : prod->factors) {
                        LiteralRational lr;
                        if (auto ex = try_get_exact_rational(f, lr); ex.is_ok() && ex.value()) {
                            b_coeff *= lr.value;
                        } else if (const auto* sqrt_c = expr_cast<FuncCall>(f);
                            sqrt_c && sqrt_c->func_id == BuiltinOp::Sqrt && !found_sqrt) {
                            c_val = sqrt_c->args[0];
                            found_sqrt = true;
                        } else {
                            found_sqrt = false;
                            break;
                        }
                    }
                    if (found_sqrt && c_val) {
                        LiteralRational lr_c;
                        if (auto ex_c = try_get_exact_rational(c_val, lr_c); ex_c.is_ok() && ex_c.value()) {
                            rat_b.value = b_coeff;
                            rat_c.value = lr_c.value;
                            b_ptr = make_rational(arena_, b_coeff);
                            c_ptr = c_val;
                        }
                    }
                }
            }
            if (a_ptr && b_ptr && c_ptr) {
                Rational a = rat_a.value;
                Rational b = rat_b.value;
                Rational c = rat_c.value;
                Rational discriminant = a*a - b*b*c;
                if (discriminant >= Rational(0)) {
                    BigInt d_num = discriminant.numerator();
                    BigInt d_den = discriminant.denominator();
                    BigInt s_num = integer_sqrt(d_num);
                    BigInt s_den = integer_sqrt(d_den);
                    if (s_num * s_num == d_num && s_den * s_den == d_den) {
                        Rational s(s_num, s_den);
                        Rational x = (a + s) / Rational(2);
                        Rational y = (a - s) / Rational(2);
                        auto sqrt_x = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Sqrt,
                            std::vector<ExprPtr>{make_rational(arena_, x)}));
                        auto sqrt_y = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Sqrt,
                            std::vector<ExprPtr>{make_rational(arena_, y)}));
                        if (sqrt_x.is_ok() && sqrt_y.is_ok()) {
                            ExprPtr res;
                            if (b >= Rational(0))
                                res = arena_.make<Sum>(std::vector<ExprPtr>{sqrt_x.value(), sqrt_y.value()});
                            else
                                res = arena_.make<Binary>(BinaryOp::Sub, sqrt_x.value(), sqrt_y.value());
                            return simplify_expr(res);
                        }
                    }
                }
            }
        }

        LiteralRational rat;
        auto exact = try_get_exact_rational(args.front(), rat);
        if (exact.is_error()) return fail<ExprPtr>(exact.error());
        if (exact.is_ok() && exact.value()) {
            if (rat.value.numerator().is_zero())
                return traced_result(RuleId::Unknown, target_before, make_integer(arena_, BigInt(0)));
            if (rat.value == Rational(BigInt(1)))
                return traced_result(RuleId::Unknown, target_before, make_integer(arena_, BigInt(1)));
            if (rat.value.numerator().is_negative()) {
                auto pos_rat = -rat.value;
                auto sqrt_pos = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Sqrt,
                    std::vector<ExprPtr>{make_rational(arena_, pos_rat)}));
                if (sqrt_pos.is_ok()) {
                    auto product = simplify_expr(arena_.make<Binary>(
                        BinaryOp::Mul,
                        arena_.make<Constant>(MathConstant::I),
                        sqrt_pos.value()));
                    if (product.is_error()) return product;
                    return traced_result(RuleId::Unknown, target_before, product.value());
                }
            }
            auto num_sqrt = integer_sqrt(rat.value.numerator());
            auto den_sqrt = integer_sqrt(rat.value.denominator());
            if (num_sqrt * num_sqrt == rat.value.numerator()
                && den_sqrt * den_sqrt == rat.value.denominator())
                return traced_result(RuleId::Unknown, target_before,
                    make_rational(arena_, Rational(num_sqrt, den_sqrt)));
            auto denested = simplify_rational_sqrt(rat.value, arena_);
            if (denested.is_ok())
                return traced_result(RuleId::Unknown, target_before, denested.value());
        }

        if (is_known_negative(args.front())) {
            auto negated_arg = simplify_expr(arena_.make<Unary>(UnaryOp::Neg, args.front()));
            if (negated_arg.is_ok()) {
                auto sqrt_pos = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Sqrt,
                    std::vector<ExprPtr>{negated_arg.value()}));
                if (sqrt_pos.is_ok()) {
                    ExprPtr res = arena_.make<Binary>(BinaryOp::Mul,
                        arena_.make<Constant>(MathConstant::I), sqrt_pos.value());
                    return traced_result(RuleId::Unknown, target_before, res);
                }
            }
        }
        // sqrt(x^2) = x (if x >= 0) or abs(x)
        if (const auto* power = expr_cast<Binary>(args.front());
            power != nullptr && power->op == BinaryOp::Pow) {
            if (auto exp = try_get_integer_exponent(power->right);
                exp.has_value() && *exp == BigInt(2)) {
                if (is_known_nonnegative(power->left)) {
                    append_assumption(target_before);
                    return traced_result(RuleId::SimplifySqrtSquare,
                        target_before, power->left);
                }
                return traced_result(RuleId::SimplifySqrtSquare, target_before,
                    arena_.make<FuncCall>(BuiltinOp::Abs,
                        std::vector<ExprPtr>{power->left}));
            }
        }
        // sqrt(sqrt(x)) -> x^(1/4)
        if (const auto* inner = expr_cast<FuncCall>(args.front());
            inner && inner->func_id == BuiltinOp::Sqrt) {
            return simplify_expr(arena_.make<Binary>(BinaryOp::Pow,
                inner->args[0],
                make_rational(arena_, Rational(BigInt(1), BigInt(4)))));
        }
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

} // namespace cas::symbolic::detail
