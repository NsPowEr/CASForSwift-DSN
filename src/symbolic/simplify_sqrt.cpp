// W9.3 split: Simplifier::simplify_funcall_sqrt — the Sqrt branch extracted from
// the former monolithic simplify_exp_log.cpp (854 LOC). Rational/denesting
// helpers it relies on live in simplify_sqrt_helpers.cpp (declared in
// simplify_sqrt_internal.hpp).

#include "simplify_sqrt_internal.hpp"
#include "simplify_branch_cut.hpp"

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_funcall_sqrt(
    ExprPtr original, std::vector<ExprPtr> args, ExprPtr target_before) {

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
        const std::size_t trial_bound = context_
            ? context_->simplify_sqrt_trial_division_bound()
            : 10000U;
        auto denested = simplify_rational_sqrt(rat.value, arena_, trial_bound);
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
    // sqrt(x^2) with branch-cut awareness (F4.K / Branch_Cut_Propagation.md):
    //   x known nonneg     → x          (always safe)
    //   x known real       → abs(x)     (real branch)
    //   complex generic    : if ctx.strict_branch_cuts() → keep structural
    //                        else                         → abs(x) (legacy)
    // Reference: Kahan 1987; Corless-Davenport-Jeffrey 2000.
    if (const auto* power = expr_cast<Binary>(args.front());
        power != nullptr && power->op == BinaryOp::Pow) {
        if (auto exp = try_get_integer_exponent(power->right);
            exp.has_value() && *exp == BigInt(2)) {
            if (is_known_nonnegative(power->left)) {
                append_assumption(target_before);
                return traced_result(RuleId::SimplifySqrtSquare,
                    target_before, power->left);
            }
            const bool x_is_real = (context_ != nullptr)
                && context_->assumptions().is_real(power->left);
            if (x_is_real) {
                return traced_result(RuleId::SimplifySqrtSquare, target_before,
                    arena_.make<FuncCall>(BuiltinOp::Abs,
                        std::vector<ExprPtr>{power->left}));
            }
            // x not provably real.
            const bool strict = (context_ != nullptr) && context_->strict_branch_cuts();
            if (strict) {
                // F8.0-6.2 / Task 20 BC-1 (Branch_Cut_Propagation.md §2 rule 1):
                //   sqrt(z²) = z · (-1)^K(2·ln(z))
                // Emit the explicit unwinding correction so the identity
                // stays algebraically exact in the complex plane.
                ExprPtr correction = symbolic::branch_cut::make_sqrt_of_square_correction(
                    power->left, arena_);
                ExprPtr corrected = arena_.make<Binary>(BinaryOp::Mul,
                    power->left, correction);
                return traced_result(RuleId::SimplifySqrtSquare, target_before, corrected);
            }
            // Legacy default for complex generic x: emit abs(x).
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
    // HC-KV-04: sqrt(Pow(a, b)) with rational b, a non-literal — legacy mode.
    //   sqrt(a^b) → a^(b/2)
    // The b=2 case is handled above with branch-cut awareness.  This rule
    // covers b != 2 for symbolic bases (e.g. Kovacic Case 2 sqrt(x^{-1}) →
    // x^{-1/2}).  Integer-literal bases are routed through the existing
    // integer-sqrt and prime-power factorization paths so that downstream
    // tests expecting canonical `2^{k}` forms keep working.
    if (const auto* power = expr_cast<Binary>(args.front());
        power != nullptr && power->op == BinaryOp::Pow) {
        const bool base_is_literal = expr_is<IntegerLit>(power->left)
            || expr_is<RationalLit>(power->left);
        LiteralRational exp_rat;
        auto ex_res = try_get_exact_rational(power->right, exp_rat);
        const bool strict = (context_ != nullptr)
            && context_->strict_branch_cuts();
        if (!strict && !base_is_literal
            && ex_res.is_ok() && ex_res.value()
            && exp_rat.value != Rational(BigInt(2))) {
            Rational half_b = exp_rat.value / Rational(BigInt(2));
            if (half_b.numerator().is_zero())
                return traced_result(RuleId::Unknown, target_before,
                    make_integer(arena_, BigInt(1)));
            return simplify_expr(arena_.make<Binary>(BinaryOp::Pow,
                power->left, make_rational(arena_, half_b)));
        }
    }
    // HC-KV-04: sqrt(N / D) → sqrt(N) / sqrt(D) — legacy mode only.
    if (const auto* div = expr_cast<Binary>(args.front());
        div != nullptr && div->op == BinaryOp::Div) {
        const bool strict = (context_ != nullptr)
            && context_->strict_branch_cuts();
        if (!strict) {
            auto sqrt_n = simplify_expr(arena_.make<FuncCall>(
                BuiltinOp::Sqrt, std::vector<ExprPtr>{div->left}));
            auto sqrt_d = simplify_expr(arena_.make<FuncCall>(
                BuiltinOp::Sqrt, std::vector<ExprPtr>{div->right}));
            if (sqrt_n.is_ok() && sqrt_d.is_ok()) {
                return simplify_expr(arena_.make<Binary>(BinaryOp::Div,
                    sqrt_n.value(), sqrt_d.value()));
            }
        }
    }
    // HC-KV-04: sqrt(c · Pow(x, n)) with c ≥ 0 rational, x symbolic →
    //   sqrt(c) · Pow(x, n/2).
    // Restricted to the Kovacic Case 2 / ω-construction shape: a single
    // Pow factor whose base is non-literal.  This avoids disturbing
    // existing tests that depend on the canonical `2^k` form when the
    // radicand is purely rational-number arithmetic.
    if (const auto* prod = expr_cast<Product>(args.front()); prod) {
        const bool strict = (context_ != nullptr)
            && context_->strict_branch_cuts();
        if (!strict) {
            std::vector<ExprPtr> nonrational;
            Rational scalar(BigInt(1));
            bool found_scalar = false;
            bool scalar_negative = false;
            for (ExprPtr f : prod->factors) {
                LiteralRational lr;
                if (auto ex = try_get_exact_rational(f, lr);
                    ex.is_ok() && ex.value()) {
                    if (lr.value.numerator().is_negative())
                        scalar_negative = true;
                    scalar = scalar * lr.value;
                    found_scalar = true;
                } else {
                    nonrational.push_back(f);
                }
            }
            // Require: single Pow(non-literal-base, anything) factor.
            bool shape_ok = (nonrational.size() == 1U);
            if (shape_ok) {
                auto* pw = expr_cast<Binary>(nonrational[0]);
                shape_ok = pw && pw->op == BinaryOp::Pow
                    && !(expr_is<IntegerLit>(pw->left)
                         || expr_is<RationalLit>(pw->left));
            }
            if (found_scalar && shape_ok && !scalar_negative
                && !scalar.numerator().is_zero()) {
                auto sqrt_c = simplify_expr(arena_.make<FuncCall>(
                    BuiltinOp::Sqrt,
                    std::vector<ExprPtr>{make_rational(arena_, scalar)}));
                auto sqrt_rest = simplify_expr(arena_.make<FuncCall>(
                    BuiltinOp::Sqrt,
                    std::vector<ExprPtr>{nonrational[0]}));
                if (sqrt_c.is_ok() && sqrt_rest.is_ok()) {
                    return simplify_expr(arena_.make<Binary>(BinaryOp::Mul,
                        sqrt_c.value(), sqrt_rest.value()));
                }
            }
        }
    }
    // Borodin-Fagin-Hopcroft-Tompa denesting:
    //   sqrt(a + b·sqrt(c)) → sqrt(p) ± sqrt(q)
    // when a²-b²c is a rational square.
    if (auto denested = try_denest_borodin_fagin(args.front(), arena_)) {
        // Recurse simplify on result to cascade any inner reductions.
        auto recursed = simplify_expr(*denested);
        if (recursed.is_ok())
            return traced_result(RuleId::Unknown, target_before, recursed.value());
        return traced_result(RuleId::Unknown, target_before, *denested);
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(BuiltinOp::Sqrt, std::move(args)));
}

}  // namespace cas::symbolic::detail
