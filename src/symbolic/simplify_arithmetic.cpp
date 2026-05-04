#include "simplify_impl.hpp"
#include "cas/linalg/Matrix.hpp"
#include <algorithm>

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const IntegerLit&) {
    return ok(original);
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const RationalLit& node) {
    auto value = Rational::make(node.numerator, node.denominator);
    if (value.is_error()) return fail<ExprPtr>(value.error());
    if (value.value().numerator() == node.numerator && value.value().denominator() == node.denominator) return ok(original);
    return make_rational_result(arena_, std::move(value.value()));
}

Result<ExprPtr> Simplifier::simplify_node(const DecimalLit& node) {
    auto rational = decimal_to_rational(node);
    if (rational.denominator() == BigInt(1)) {
        return ok(arena_.make<IntegerLit>(rational.numerator()));
    }
    return make_rational_result(arena_, std::move(rational));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const DecimalLit& node) {
    auto rational = decimal_to_rational(node);
    if (rational.denominator() == BigInt(1)) {
        if (rational.numerator() == BigInt(0)) return ok(original); // For completeness, might want to just make new node
        return ok(arena_.make<IntegerLit>(rational.numerator()));
    }
    return make_rational_result(arena_, std::move(rational));
}

Result<ExprPtr> Simplifier::simplify_node(const ExprNode&) {
    return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot simplify null expression node"));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Symbol&) {
    return ok(original);
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Constant&) {
    return ok(original);
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Unary& node) {
    Result<ExprPtr> operand = fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
    if (trace_enabled_) {
        ScopedFrame frame(*this, [this, op = node.op](ExprPtr current_operand) {
            return arena_.make<Unary>(op, current_operand);
        });
        operand = simplify_expr(node.operand);
    } else {
        operand = simplify_expr(node.operand);
    }
    if (operand.is_error()) return operand;

    const ExprPtr target_before = (operand.value() == node.operand) ? original : arena_.make<Unary>(node.op, operand.value());

    if (node.op == UnaryOp::Neg) {
        LiteralRational rational;
        auto exact = try_get_exact_rational(operand.value(), rational);
        if (exact.is_error()) return fail<ExprPtr>(exact.error());
        if (exact.value()) return make_rational_result(arena_, -rational.value);
        
        if (const auto* nested = expr_cast<Unary>(operand.value()); nested != nullptr && nested->op == UnaryOp::Neg) {
            return traced_result(RuleId::SimplifyCollectLikeTerms, target_before, nested->operand);
        }

        if (const auto* sum = expr_cast<Sum>(operand.value())) {
            std::vector<ExprPtr> negated_terms;
            negated_terms.reserve(sum->terms.size());
            for (ExprPtr term : sum->terms) {
                auto negated = simplify_expr(arena_.make<Unary>(UnaryOp::Neg, term));
                if (negated.is_error()) return negated;
                negated_terms.push_back(negated.value());
            }
            return simplify_sum_terms(negated_terms, target_before, true);
        }

        if (const auto* product = expr_cast<Product>(operand.value())) {
            std::vector<ExprPtr> factors = product->factors;
            factors.push_back(make_integer(arena_, BigInt(-1)));
            return simplify_product_factors(factors, target_before, false);
        }
    }

    if (operand.value() == node.operand) return ok(original);
    return ok(arena_.make<Unary>(node.op, operand.value()));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Binary& node) {
    if (!trace_enabled_ && (node.op == BinaryOp::Add || node.op == BinaryOp::Sub)) {
        return simplify_additive_chain_fast(original);
    }

    Result<ExprPtr> lhs = fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
    if (trace_enabled_) {
        ScopedFrame left_frame(*this, [this, op = node.op, right = node.right](ExprPtr left) {
            return arena_.make<Binary>(op, left, right);
        });
        lhs = simplify_expr(node.left);
    } else {
        lhs = simplify_expr(node.left);
    }
    if (lhs.is_error()) return lhs;

    Result<ExprPtr> rhs = fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
    if (trace_enabled_) {
        ScopedFrame right_frame(*this, [this, op = node.op, left = lhs.value()](ExprPtr right) {
            return arena_.make<Binary>(op, left, right);
        });
        rhs = simplify_expr(node.right);
    } else {
        rhs = simplify_expr(node.right);
    }
    if (rhs.is_error()) return rhs;

    const ExprPtr target_before = (lhs.value() == node.left && rhs.value() == node.right) ? original : 
        (trace_enabled_ ? arena_.make<Binary>(node.op, lhs.value(), rhs.value()) : ExprPtr{});

    if (rewrite_provider_ != nullptr) {
        ExprPtr rewrite_target = target_before ? target_before : arena_.make<Binary>(node.op, lhs.value(), rhs.value());
        auto rewritten = rewrite_provider_->try_rewrite(rewrite_target, arena_, assumptions_, context_);
        if (rewritten.is_error()) return rewritten;
        if (rewritten.value() != rewrite_target) {
            append_trace(RuleId::RewriteProviderApplied, rewrite_target, rewritten.value());
            return simplify_expr(rewritten.value());
        }
    }

    switch (node.op) {
    case BinaryOp::Add: return simplify_sum_terms({lhs.value(), rhs.value()}, target_before, true);
    case BinaryOp::Sub: {
        auto negated_rhs = simplify_expr(arena_.make<Unary>(UnaryOp::Neg, rhs.value()));
        if (negated_rhs.is_error()) return negated_rhs;
        return simplify_sum_terms({lhs.value(), negated_rhs.value()}, target_before, true);
    }
    case BinaryOp::Mul: return simplify_product_factors({lhs.value(), rhs.value()}, target_before, true);
    case BinaryOp::Div: {
        LiteralRational left_rat, right_rat;
        auto l_exact = try_get_exact_rational(lhs.value(), left_rat);
        if (l_exact.is_error()) return fail<ExprPtr>(l_exact.error());
        auto r_exact = try_get_exact_rational(rhs.value(), right_rat);
        if (r_exact.is_error()) return fail<ExprPtr>(r_exact.error());
        if (l_exact.value() && r_exact.value()) {
            auto quot = checked_divide(left_rat.value, right_rat.value);
            if (quot.is_error()) return fail<ExprPtr>(quot.error());
            return make_rational_result(arena_, std::move(quot.value()));
        }
        return simplify_product_factors({lhs.value(), arena_.make<Binary>(BinaryOp::Pow, rhs.value(), make_integer(arena_, BigInt(-1)))}, target_before, false);
    }
    case BinaryOp::Pow: return simplify_power(lhs.value(), rhs.value(), target_before);
    case BinaryOp::Mod: {
        LiteralRational l_rat, r_rat;
        auto l_exact = try_get_exact_rational(lhs.value(), l_rat);
        if (l_exact.is_error()) return fail<ExprPtr>(l_exact.error());
        auto r_exact = try_get_exact_rational(rhs.value(), r_rat);
        if (r_exact.is_error()) return fail<ExprPtr>(r_exact.error());
        if (l_exact.value() && r_exact.value() && l_rat.value.is_integer() && r_rat.value.is_integer()) {
            auto rem = checked_mod(l_rat.value.numerator(), r_rat.value.numerator());
            if (rem.is_error()) return fail<ExprPtr>(rem.error());
            return ok(make_integer(arena_, std::move(rem.value())));
        }
        if (lhs.value() == node.left && rhs.value() == node.right) return ok(original);
        return ok(arena_.make<Binary>(BinaryOp::Mod, lhs.value(), rhs.value()));
    }
    case BinaryOp::Equal:
        if (lhs.value() == node.left && rhs.value() == node.right) return ok(original);
        return ok(arena_.make<Binary>(BinaryOp::Equal, lhs.value(), rhs.value()));
    }
    return fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unsupported binary operator"));
}

void Simplifier::collect_additive_operands(ExprPtr expr, bool negate, std::vector<std::pair<ExprPtr, bool>>& operands) {
    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Add) {
            collect_additive_operands(binary->left, negate, operands);
            collect_additive_operands(binary->right, negate, operands);
            return;
        }
        if (binary->op == BinaryOp::Sub) {
            collect_additive_operands(binary->left, negate, operands);
            collect_additive_operands(binary->right, !negate, operands);
            return;
        }
    }
    operands.push_back({expr, negate});
}

Result<ExprPtr> Simplifier::simplify_additive_chain_fast(ExprPtr original) {
    std::vector<std::pair<ExprPtr, bool>> raw_operands;
    collect_additive_operands(original, false, raw_operands);
    std::vector<ExprPtr> simplified_terms;
    simplified_terms.reserve(raw_operands.size());
    for (const auto& [operand, negate] : raw_operands) {
        auto simplified = simplify_expr(operand);
        if (simplified.is_error()) return simplified;
        if (!negate) {
            simplified_terms.push_back(simplified.value());
        } else {
            auto negated = simplify_expr(arena_.make<Unary>(UnaryOp::Neg, simplified.value()));
            if (negated.is_error()) return negated;
            simplified_terms.push_back(negated.value());
        }
    }
    return simplify_sum_terms(simplified_terms, original, true);
}

// ── Trig power linearization (Chebyshev / DeMoivre) ──────────────────────
[[nodiscard]] static BigInt trig_binomial(unsigned int n, unsigned int k) {
    if (k > n) return BigInt(0);
    if (k == 0U || k == n) return BigInt(1);
    if (k > n - k) k = n - k;
    BigInt result(1);
    for (unsigned int i = 0U; i < k; ++i) {
        result = result * BigInt(static_cast<long long>(n - i));
        result = result / BigInt(static_cast<long long>(i + 1U));
    }
    return result;
}

[[nodiscard]] static BigInt trig_pow4(unsigned int m) {
    BigInt result(1);
    const BigInt four(4);
    for (unsigned int i = 0U; i < m; ++i) result = result * four;
    return result;
}

// Returns linearized form of sin^n or cos^n using multiple-angle identities.
// Odd n=2m+1:  trig^n = (1/4^m) * sum_{j=0}^{m} c_j * trig((n-2j)*arg)
// Even n=2m:   trig^n = (1/4^m) * [C(n,m) + 2*sum_{j=0}^{m-1} c_j * cos((n-2j)*arg)]
[[nodiscard]] static ExprPtr try_trig_linearize(AstArena& arena, BuiltinOp func, ExprPtr arg, unsigned int n) {
    const bool is_sin = (func == BuiltinOp::Sin);
    std::vector<ExprPtr> terms;

    if (n % 2U == 1U) {
        // odd
        const unsigned int m = (n - 1U) / 2U;
        const BigInt denom = trig_pow4(m);
        for (unsigned int j = 0U; j <= m; ++j) {
            const unsigned int freq = n - 2U * j;
            const BigInt c = trig_binomial(n, j);
            const bool negate_coeff = is_sin && ((m - j) % 2U == 1U);
            const BigInt num = negate_coeff ? -c : c;
            ExprPtr angle_arg = (freq == 1U)
                ? arg
                : arena.make<Binary>(BinaryOp::Mul,
                    arena.make<IntegerLit>(BigInt(static_cast<long long>(freq))), arg);
            ExprPtr trig_node = arena.make<FuncCall>(func, std::vector<ExprPtr>{angle_arg});
            terms.push_back(arena.make<Binary>(BinaryOp::Mul,
                arena.make<RationalLit>(num, denom), trig_node));
        }
    } else {
        // even
        const unsigned int m = n / 2U;
        const BigInt denom = trig_pow4(m);
        terms.push_back(arena.make<RationalLit>(trig_binomial(n, m), denom));
        for (unsigned int j = 0U; j < m; ++j) {
            const unsigned int freq = n - 2U * j;
            const BigInt c = trig_binomial(n, j);
            const bool negate_coeff = is_sin && ((m - j) % 2U == 1U);
            const BigInt num = negate_coeff ? BigInt(-2) * c : BigInt(2) * c;
            ExprPtr angle_arg = (freq == 1U)
                ? arg
                : arena.make<Binary>(BinaryOp::Mul,
                    arena.make<IntegerLit>(BigInt(static_cast<long long>(freq))), arg);
            ExprPtr trig_node = arena.make<FuncCall>(BuiltinOp::Cos, std::vector<ExprPtr>{angle_arg});
            terms.push_back(arena.make<Binary>(BinaryOp::Mul,
                arena.make<RationalLit>(num, denom), trig_node));
        }
    }

    if (terms.empty()) return ExprPtr{};
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}
// ──────────────────────────────────────────────────────────────────────────

Result<ExprPtr> Simplifier::simplify_power(ExprPtr base, ExprPtr exponent, ExprPtr target_before) {
    if (!target_before && trace_enabled_) {
        target_before = arena_.make<Binary>(BinaryOp::Pow, base, exponent);
    }

    if (rewrite_provider_ != nullptr && may_rewrite_power(base, exponent)) {
        ExprPtr rewrite_target = target_before ? target_before : arena_.make<Binary>(BinaryOp::Pow, base, exponent);
        auto rewritten = rewrite_provider_->try_rewrite(rewrite_target, arena_, assumptions_, context_);
        if (rewritten.is_error()) return rewritten;
        if (rewritten.value() != rewrite_target) {
            append_trace(RuleId::RewriteProviderApplied, rewrite_target, rewritten.value());
            return simplify_expr(rewritten.value());
        }
    }

    LiteralRational base_rat, exp_rat;
    auto b_exact = try_get_exact_rational(base, base_rat);
    if (b_exact.is_error()) return fail<ExprPtr>(b_exact.error());
    auto e_exact = try_get_exact_rational(exponent, exp_rat);
    if (e_exact.is_error()) return fail<ExprPtr>(e_exact.error());

    if (b_exact.value() && e_exact.value()) {
        if (exp_rat.value.is_integer()) {
            const BigInt power = exp_rat.value.numerator();
            if (!power.is_negative()) {
                return make_rational_result(arena_, pow_rational_nonnegative(base_rat.value, power));
            }
            Rational result = pow_rational_nonnegative(base_rat.value, -power);
            auto inverse = checked_divide(Rational(BigInt(1)), result);
            if (inverse.is_error()) return fail<ExprPtr>(inverse.error());
            return make_rational_result(arena_, std::move(inverse.value()));
        }
    }

    if (is_constant_expr(base, MathConstant::Infinity) && e_exact.value() && exp_rat.value.is_integer()) {
        const BigInt power = exp_rat.value.numerator();
        if (power.is_negative()) {
            return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(0)));
        }
        if (!power.is_zero()) {
            return traced_result(RuleId::SimplifyPowerOne, target_before, base);
        }
    }

    if (is_zero_expr(exponent)) return traced_result(RuleId::SimplifyPowerZero, target_before, make_integer(arena_, BigInt(1)));
    if (is_one_expr(exponent)) return traced_result(RuleId::SimplifyPowerOne, target_before, base);
    if (is_one_expr(base)) return ok(base);

    // x^(1/2) → sqrt(x)  (canonical form)
    {
        LiteralRational exp_half;
        auto e_exact = try_get_exact_rational(exponent, exp_half);
        if (e_exact.is_ok() && e_exact.value() &&
            exp_half.value.numerator() == BigInt(1) && exp_half.value.denominator() == BigInt(2)) {
            ExprPtr sqrt_expr = arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{base});
            return simplify_expr(sqrt_expr);
        }
    }
    if (is_zero_expr(base)) {
        LiteralRational exp_rat_check;
        auto exp_check = try_get_exact_rational(exponent, exp_rat_check);
        if (exp_check.is_ok() && exp_check.value() && exp_rat_check.value.is_integer() && !exp_rat_check.value.numerator().is_negative() && !exp_rat_check.value.numerator().is_zero()) {
            return traced_result(RuleId::SimplifyZeroPowerPositive, target_before, make_integer(arena_, BigInt(0)));
        }
    }

    if (is_constant_expr(base, MathConstant::E)) {
        // E^(a + b + ...) = E^a * E^b * ...  (linearizzazione esponenziale)
        if (const auto* sum_exp = expr_cast<Sum>(exponent); sum_exp != nullptr && !sum_exp->terms.empty()) {
            std::vector<ExprPtr> factors;
            factors.reserve(sum_exp->terms.size());
            for (ExprPtr term : sum_exp->terms) {
                factors.push_back(arena_.make<Binary>(BinaryOp::Pow, base, term));
            }
            return simplify_expr(arena_.make<Product>(std::move(factors)));
        }

        // E^(ln(x)) = x: valido per x non noto negativo/zero
        const auto* call = expr_cast<FuncCall>(exponent);
        if (call != nullptr && call->func_id == BuiltinOp::Ln && call->args.size() == 1U) {
            ExprPtr arg = call->args.front();
            LiteralRational arg_rat;
            auto arg_exact = try_get_exact_rational(arg, arg_rat);
            bool is_known_nonpositive = arg_exact.is_ok() && arg_exact.value() &&
                (arg_rat.value.numerator().is_zero() || arg_rat.value.numerator().is_negative());
            if (!is_known_nonpositive) {
                return traced_result(RuleId::SimplifyExpLnPositive, target_before, arg);
            }
        }
    }

    // I^n simplification
    if (const auto* constant = expr_cast<Constant>(base); constant != nullptr && constant->value == MathConstant::I) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            BigInt rem = n % BigInt(4);
            if (rem.is_negative()) rem += BigInt(4);
            const long long r = rem.to_u64();
            if (r == 0) return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(1)));
            if (r == 1) return traced_result(RuleId::SimplifyPowerOne, target_before, base);
            if (r == 2) return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(-1)));
            if (r == 3) return traced_result(RuleId::SimplifyPowerOne, target_before, arena_.make<Unary>(UnaryOp::Neg, base));
        }
    }
    if (const auto* unary = expr_cast<Unary>(base);
        unary != nullptr && unary->op == UnaryOp::Neg && is_constant_expr(unary->operand, MathConstant::I)) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            BigInt rem = n % BigInt(4);
            if (rem.is_negative()) rem += BigInt(4);
            const long long r = rem.to_u64();
            if (r == 0) return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(1)));
            if (r == 1) return traced_result(RuleId::SimplifyPowerOne, target_before, base);
            if (r == 2) return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(-1)));
            if (r == 3) return traced_result(RuleId::SimplifyPowerOne, target_before, unary->operand);
        }
    }

    // sqrt(A)^n simplification
    if (const auto* sqrt_call = expr_cast<FuncCall>(base);
        sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            if (!n.is_negative() && !n.is_zero()) {
                ExprPtr arg = sqrt_call->args.front();
                BigInt k = n / BigInt(2);
                BigInt rem = n % BigInt(2);
                
                if (k.is_zero()) {
                    // sqrt(A)^1 = sqrt(A), handled by default return
                } else if (rem.is_zero()) {
                    // (sqrt(A))^(2k) = A^k
                    return simplify_power(arg, make_integer(arena_, k));
                } else {
                    // (sqrt(A))^(2k+1) = A^k * sqrt(A)
                    auto ak = simplify_power(arg, make_integer(arena_, k));
                    if (ak.is_error()) return ak;
                    return simplify_product_factors({ak.value(), base});
                }
            }
        }
    }

    if (const auto* product = expr_cast<Product>(base); product != nullptr) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            if (!n.is_zero() && !n.is_negative() && n <= BigInt(20)) {
                std::vector<ExprPtr> factors;
                for (ExprPtr factor : product->factors) {
                    auto p = simplify_power(factor, exponent);
                    if (p.is_error()) return p;
                    factors.push_back(p.value());
                }
                return simplify_product_factors(factors);
            }
        }
    }

    if (const auto* outer = expr_cast<Binary>(base); outer != nullptr && outer->op == BinaryOp::Pow) {
        bool safe = false;
        auto c_int = try_get_integer_exponent(exponent);
        if (c_int.has_value()) {
            safe = true;
        } else if (is_known_positive(outer->left)) {
            safe = true;
        }

        if (safe) {
            auto new_exp = simplify_expr(arena_.make<Binary>(BinaryOp::Mul, outer->right, exponent));
            if (new_exp.is_ok()) {
                auto rewritten = simplify_power(outer->left, new_exp.value());
                if (rewritten.is_ok()) {
                    append_trace(RuleId::SimplifyFlattenNestedPowers, target_before, rewritten.value());
                    return rewritten;
                }
            }
        } else {
            // Fallback for literal rationals (already covers some cases)
            LiteralRational l_exp, r_exp;
            auto l_exact = try_get_exact_rational(outer->right, l_exp);
            auto r_exact = try_get_exact_rational(exponent, r_exp);
            if (l_exact.is_ok() && l_exact.value() && r_exact.is_ok() && r_exact.value()) {
                // (a^b)^c where b, c are rational. Still needs a > 0 or specific conditions,
                // but CAS often flattens these if they are literal rationals unless it's a very strict mode.
                // For now, we keep the previous behavior for literal rationals but we could be more precise.
                auto rewritten = simplify_power(outer->left, make_rational(arena_, l_exp.value * r_exp.value));
                if (rewritten.is_ok()) {
                    append_trace(RuleId::SimplifyFlattenNestedPowers, target_before, rewritten.value());
                    return rewritten;
                }
            }
        }
    }

    // sqrt(A)^2 = A for any real A (principal branch: sqrt(-k) = i*sqrt(k), (i*sqrt(k))^2 = -k = A).
    // Only apply when A is a known real constant (nonneg or known-negative rational).
    if (const auto* sqrt_call = expr_cast<FuncCall>(base);
        sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U) {
        LiteralRational exp_r;
        auto e_ex = try_get_exact_rational(exponent, exp_r);
        if (e_ex.is_ok() && e_ex.value() && exp_r.value.is_integer() && exp_r.value.numerator() == BigInt(2)) {
            // (sqrt(x))^2 = x is always true for principal complex branch.
            return traced_result(RuleId::SimplifyFlattenNestedPowers, target_before, sqrt_call->args.front());
        }
    }

    // Trig power linearization: sin(u)^n / cos(u)^n for integer n >= 2
    if (const auto* trig_call = expr_cast<FuncCall>(base);
        trig_call != nullptr &&
        (trig_call->func_id == BuiltinOp::Sin || trig_call->func_id == BuiltinOp::Cos) &&
        trig_call->args.size() == 1U) {
        if (auto maybe_n = try_get_integer_exponent(exponent);
            maybe_n.has_value() && !maybe_n->is_negative() &&
            *maybe_n >= BigInt(2) && *maybe_n <= BigInt(20)) {
            const auto n_val = static_cast<unsigned int>(maybe_n->to_u64());
            if (auto lin = try_trig_linearize(arena_, trig_call->func_id, trig_call->args.front(), n_val)) {
                return simplify_expr(lin);
            }
        }
    }

    // GAP #6: exp(a)^b -> exp(a*b)
    ExprPtr exp_base = base;
    bool neg_sign = false;
    if (const auto* neg = expr_cast<Unary>(base); neg != nullptr && neg->op == UnaryOp::Neg) {
        exp_base = neg->operand;
        neg_sign = true;
    }

    if (const auto* exp_call = expr_cast<FuncCall>(exp_base);
        exp_call != nullptr && exp_call->func_id == BuiltinOp::Exp && exp_call->args.size() == 1U) {
        auto new_arg = simplify_expr(arena_.make<Binary>(BinaryOp::Mul, exp_call->args.front(), exponent));
        if (new_arg.is_ok()) {
            ExprPtr res = arena_.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{new_arg.value()});
            if (neg_sign) {
                if (auto exp_val = try_get_integer_exponent(exponent); exp_val.has_value()) {
                    if ((*exp_val % BigInt(2)) != BigInt(0)) res = arena_.make<Unary>(UnaryOp::Neg, res);
                    // if even, remains positive
                } else {
                    // complex case, skip for now
                    return ok(arena_.make<Binary>(BinaryOp::Pow, base, exponent));
                }
            }
            return simplify_expr(res);
        }
    }

    if (expr_is<Binary>(target_before)) {
        const auto& before = expr_ref<Binary>(target_before);
        if (before.op == BinaryOp::Pow && before.left == base && before.right == exponent) return ok(target_before);
    }
    return ok(arena_.make<Binary>(BinaryOp::Pow, base, exponent));
}

bool Simplifier::monomial_keys_equal(const MonomialKey& lhs, const MonomialKey& rhs) {
    if (lhs.factors.size() != rhs.factors.size()) return false;
    for (std::size_t i = 0; i < lhs.factors.size(); ++i) {
        if (!structural_equal(lhs.factors[i].first, rhs.factors[i].first) || lhs.factors[i].second != rhs.factors[i].second) return false;
    }
    return true;
}

Result<std::optional<MonomialTerm>> Simplifier::extract_monomial(ExprPtr expr) {
    LiteralRational rational;
    auto exact = try_get_exact_rational(expr, rational);
    if (exact.is_error()) return fail<std::optional<MonomialTerm>>(exact.error());
    if (exact.value()) return ok(std::optional<MonomialTerm>(MonomialTerm{.coefficient = rational.value, .key = MonomialKey{}}));

    if (expr_is<Symbol>(expr)) return ok(std::optional<MonomialTerm>(MonomialTerm{.coefficient = Rational(BigInt(1)), .key = MonomialKey{{{expr, BigInt(1)}}}}));

    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) {
        auto inner = extract_monomial(unary->operand);
        if (inner.is_ok() && inner.value().has_value()) {
            inner.value()->coefficient = -inner.value()->coefficient;
            return inner;
        }
    }

    if (const auto* binary = expr_cast<Binary>(expr)) {
        if (binary->op == BinaryOp::Div) {
            auto lhs_res = extract_monomial(binary->left);
            auto rhs_res = extract_monomial(binary->right);
            if (lhs_res.is_ok() && lhs_res.value().has_value() && rhs_res.is_ok() && rhs_res.value().has_value()) {
                auto lhs = *lhs_res.value();
                auto rhs = *rhs_res.value();
                auto coeff = checked_divide(lhs.coefficient, rhs.coefficient);
                if (coeff.is_error()) return fail<std::optional<MonomialTerm>>(coeff.error());
                MonomialKey key = lhs.key;
                for (const auto& rhs_factor : rhs.key.factors) {
                    bool found = false;
                    for (auto& lhs_factor : key.factors) {
                        if (structural_equal(rhs_factor.first, lhs_factor.first)) {
                            lhs_factor.second -= rhs_factor.second;
                            found = true;
                            break;
                        }
                    }
                    if (!found) key.factors.push_back({rhs_factor.first, -rhs_factor.second});
                }
                key.factors.erase(std::remove_if(key.factors.begin(), key.factors.end(), [](const auto& p) { return p.second.is_zero(); }), key.factors.end());
                return ok(std::optional<MonomialTerm>(MonomialTerm{.coefficient = coeff.value(), .key = std::move(key)}));
            }
        }
        if (binary->op == BinaryOp::Pow) {
            if (auto exponent = try_get_integer_exponent(binary->right); exponent.has_value()) {
                return ok(std::optional<MonomialTerm>(MonomialTerm{.coefficient = Rational(BigInt(1)), .key = MonomialKey{{{binary->left, *exponent}}}}));
            }
        }
    }

    if (const auto* product = expr_cast<Product>(expr)) {
        Rational coefficient(BigInt(1));
        std::vector<std::pair<ExprPtr, BigInt>> factors;
        for (ExprPtr factor : product->factors) {
            auto factor_exact = try_get_exact_rational(factor, rational);
            if (factor_exact.is_ok() && factor_exact.value()) {
                coefficient *= rational.value;
                continue;
            }
            if (expr_is<Symbol>(factor)) {
                factors.push_back({factor, BigInt(1)});
                continue;
            }
            if (const auto* power = expr_cast<Binary>(factor); power != nullptr && power->op == BinaryOp::Pow) {
                if (auto exponent = try_get_integer_exponent(power->right); exponent.has_value() && !exponent->is_negative()) {
                    factors.push_back({power->left, *exponent});
                    continue;
                }
            }
            factors.push_back({factor, BigInt(1)});
        }
        merge_symbolic_factors(factors);
        return ok(std::optional<MonomialTerm>(MonomialTerm{.coefficient = coefficient, .key = MonomialKey{std::move(factors)}}));
    }

    if (!expr_is<Sum>(expr)) return ok(std::optional<MonomialTerm>(MonomialTerm{.coefficient = Rational(BigInt(1)), .key = MonomialKey{{{expr, BigInt(1)}}}}));
    return ok(std::optional<MonomialTerm>{});
}

ExprPtr Simplifier::build_monomial(const MonomialKey& key, const Rational& coefficient) {
    std::vector<ExprPtr> factors;
    bool is_negative_one = (coefficient == Rational(BigInt(-1)));
    
    if (!is_negative_one && (!(coefficient == Rational(BigInt(1))) || key.factors.empty())) {
        factors.push_back(make_rational(arena_, coefficient));
    }
    
    for (const auto& [base, exponent] : key.factors) {
        if (exponent == BigInt(1)) {
            factors.push_back(base);
        } else {
            factors.push_back(arena_.make<Binary>(BinaryOp::Pow, base, make_integer(arena_, exponent)));
        }
    }
    
    if (is_negative_one) {
        if (factors.empty()) return make_integer(arena_, BigInt(-1));
        if (factors.size() == 1U) return arena_.make<Unary>(UnaryOp::Neg, factors.front());
        return arena_.make<Unary>(UnaryOp::Neg, arena_.make<Product>(std::move(factors)));
    }

    if (factors.size() == 1U) return factors.front();
    return arena_.make<Product>(std::move(factors));
}

void Simplifier::merge_symbolic_factors(std::vector<std::pair<ExprPtr, BigInt>>& factors) {
    std::sort(factors.begin(), factors.end(), [](const auto& lhs, const auto& rhs) { return canonical_compare(lhs.first, rhs.first) < 0; });
    std::vector<std::pair<ExprPtr, BigInt>> merged;
    for (const auto& factor : factors) {
        // structural_equal -> O(1) pointer comparison thanks to interning
        if (!merged.empty() && merged.back().first == factor.first) {
            merged.back().second += factor.second;
        } else {
            merged.push_back(factor);
        }
    }
    factors = std::move(merged);
}

} // namespace cas::symbolic::detail
