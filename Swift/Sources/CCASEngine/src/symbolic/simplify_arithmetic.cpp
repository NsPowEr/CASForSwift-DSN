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

Result<ExprPtr> Simplifier::simplify_node(const DecimalLit&) {
    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Decimal literals are preserved for lossless parsing but are not supported by the symbolic core"));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr, const DecimalLit& node) {
    return simplify_node(node);
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
        auto rewritten = rewrite_provider_->try_rewrite(rewrite_target, arena_, assumptions_);
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

Result<ExprPtr> Simplifier::simplify_power(ExprPtr base, ExprPtr exponent, ExprPtr target_before) {
    if (!target_before && trace_enabled_) {
        target_before = arena_.make<Binary>(BinaryOp::Pow, base, exponent);
    }

    if (rewrite_provider_ != nullptr && may_rewrite_power(base, exponent)) {
        ExprPtr rewrite_target = target_before ? target_before : arena_.make<Binary>(BinaryOp::Pow, base, exponent);
        auto rewritten = rewrite_provider_->try_rewrite(rewrite_target, arena_, assumptions_);
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

    if (is_zero_expr(exponent)) return traced_result(RuleId::SimplifyPowerZero, target_before, make_integer(arena_, BigInt(1)));
    if (is_one_expr(exponent)) return traced_result(RuleId::SimplifyPowerOne, target_before, base);
    if (is_one_expr(base)) return ok(base);
    if (is_zero_expr(base)) {
        LiteralRational exp_rat_check;
        auto exp_check = try_get_exact_rational(exponent, exp_rat_check);
        if (exp_check.is_ok() && exp_check.value() && exp_rat_check.value.is_integer() && !exp_rat_check.value.numerator().is_negative() && !exp_rat_check.value.numerator().is_zero()) {
            return traced_result(RuleId::SimplifyZeroPowerPositive, target_before, make_integer(arena_, BigInt(0)));
        }
    }

    if (is_constant_expr(base, MathConstant::E)) {
        const auto* call = expr_cast<FuncCall>(exponent);
        if (call != nullptr && call->name == "ln" && call->args.size() == 1U && is_known_positive(call->args.front())) {
            return traced_result(RuleId::SimplifyExpLnPositive, target_before, call->args.front());
        }
    }

    if (const auto* outer = expr_cast<Binary>(base); outer != nullptr && outer->op == BinaryOp::Pow) {
        LiteralRational l_exp, r_exp;
        auto l_exact = try_get_exact_rational(outer->right, l_exp);
        auto r_exact = try_get_exact_rational(exponent, r_exp);
        if (l_exact.is_ok() && l_exact.value() && r_exact.is_ok() && r_exact.value()) {
            auto rewritten = simplify_power(outer->left, make_rational(arena_, l_exp.value * r_exp.value));
            if (rewritten.is_error()) return rewritten;
            append_trace(RuleId::SimplifyFlattenNestedPowers, target_before, rewritten.value());
            return rewritten;
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
        if (!merged.empty() && structural_equal(merged.back().first, factor.first)) {
            merged.back().second += factor.second;
        } else {
            merged.push_back(factor);
        }
    }
    factors = std::move(merged);
}

} // namespace cas::symbolic::detail
