#include "simplify_impl.hpp"
#include <algorithm>
#include <map>

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Sum& node) {
    return simplify_sum_terms(node.terms, original);
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Product& node) {
    return simplify_product_factors(node.factors, original);
}

Result<ExprPtr> Simplifier::simplify_sum_terms(const std::vector<ExprPtr>& terms, ExprPtr target_before, bool inputs_are_simplified) {
    if (!target_before && trace_enabled_) target_before = make_sum_target(terms);
    std::vector<ExprPtr> flat_terms;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        auto simplify_current = [&]() -> Result<ExprPtr> { return inputs_are_simplified ? ok(terms[i]) : simplify_expr(terms[i]); };
        Result<ExprPtr> simplified = fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
        if (trace_enabled_) {
            std::vector<ExprPtr> current_terms = terms;
            for (std::size_t j = 0; j < flat_terms.size() && j < current_terms.size(); ++j) current_terms[j] = flat_terms[j];
            ScopedFrame frame(*this, [this, current_terms = std::move(current_terms), i](ExprPtr value) mutable {
                current_terms[i] = value;
                return make_sum_target(current_terms);
            });
            simplified = simplify_current();
        } else {
            simplified = simplify_current();
        }
        if (simplified.is_error()) return simplified;
        if (const auto* nested = expr_cast<Sum>(simplified.value())) {
            flat_terms.insert(flat_terms.end(), nested->terms.begin(), nested->terms.end());
        } else {
            flat_terms.push_back(simplified.value());
        }
    }

    while (rewrite_provider_ != nullptr && flat_terms.size() >= 2U && may_rewrite_sum_terms(flat_terms)) {
        ExprPtr rewrite_target = make_sum_target(flat_terms);
        auto rewritten = rewrite_provider_->try_rewrite(rewrite_target, arena_, assumptions_, context_);
        if (rewritten.is_error()) return rewritten;
        if (rewritten.value() == rewrite_target) break;
        append_trace(RuleId::RewriteProviderApplied, rewrite_target, rewritten.value());
        if (const auto* rewritten_sum = expr_cast<Sum>(rewritten.value())) {
            flat_terms.assign(rewritten_sum->terms.begin(), rewritten_sum->terms.end());
        } else {
            return simplify_expr(rewritten.value());
        }
    }

    Rational constant(BigInt(0));
    std::map<MonomialKey, Rational> collected;
    for (ExprPtr term : flat_terms) {
        auto timeout = check_timeout();
        if (timeout.is_error()) return fail<ExprPtr>(timeout.error());
        LiteralRational rational;
        auto exact = try_get_exact_rational(term, rational);
        if (exact.is_ok() && exact.value()) {
            constant += rational.value;
            continue;
        }
        auto monomial_res = extract_monomial(term);
        if (monomial_res.is_error()) return fail<ExprPtr>(monomial_res.error());
        if (!monomial_res.value().has_value()) {
            collected[MonomialKey{.factors = {{term, BigInt(1)}}}] += Rational(BigInt(1));
            continue;
        }
        auto& monomial = *monomial_res.value();
        collected[monomial.key] += monomial.coefficient;
    }

    std::vector<ExprPtr> normalized;
    for (const auto& [key, coeff] : collected) {
        if (!coeff.numerator().is_zero()) normalized.push_back(build_monomial(key, coeff));
    }
    if (!constant.numerator().is_zero()) normalized.push_back(make_rational(arena_, constant));

    if (normalized.empty()) return traced_result(RuleId::SimplifyCollectLikeTerms, target_before, make_integer(arena_, BigInt(0)));
    std::sort(normalized.begin(), normalized.end(), [](ExprPtr lhs, ExprPtr rhs) {
        int l_deg = polynomial_degree(lhs), r_deg = polynomial_degree(rhs);
        return l_deg != r_deg ? l_deg > r_deg : canonical_compare(lhs, rhs) < 0;
    });
    if (normalized.size() == 1U) {
        if (std::any_of(terms.begin(), terms.end(), [](ExprPtr expr) { return is_zero_expr(expr); })) return traced_result(RuleId::SimplifyAddZero, target_before, normalized.front());
        return traced_result(RuleId::SimplifyCollectLikeTerms, target_before, normalized.front());
    }
    if (expr_is<Sum>(target_before) && expr_ptr_sequence_identical(normalized, expr_ref<Sum>(target_before).terms)) return ok(target_before);
    return ok(arena_.make<Sum>(std::move(normalized)));
}

Result<ExprPtr> Simplifier::simplify_product_factors(const std::vector<ExprPtr>& factors, ExprPtr target_before, bool inputs_are_simplified) {
    if (!target_before && trace_enabled_) target_before = make_product_target(factors);
    std::vector<ExprPtr> initial_factors;
    bool has_zero = false;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        auto simplify_current = [&]() -> Result<ExprPtr> { return inputs_are_simplified ? ok(factors[i]) : simplify_expr(factors[i]); };
        Result<ExprPtr> simplified = fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
        if (trace_enabled_) {
            std::vector<ExprPtr> current_factors = factors;
            for (std::size_t j = 0; j < initial_factors.size() && j < current_factors.size(); ++j) current_factors[j] = initial_factors[j];
            ScopedFrame frame(*this, [this, current_factors = std::move(current_factors), i](ExprPtr value) mutable {
                current_factors[i] = value;
                return make_product_target(current_factors);
            });
            simplified = simplify_current();
        } else {
            simplified = simplify_current();
        }
        if (simplified.is_error()) return simplified;
        if (is_zero_expr(simplified.value())) {
            has_zero = true;
        }
        initial_factors.push_back(simplified.value());
    }

    if (has_zero) {
        for (ExprPtr factor : initial_factors) {
            if (const auto* c = expr_cast<Constant>(factor)) {
                if (c->value == MathConstant::Infinity) {
                    return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "0 * Infinity is undefined"));
                }
            }
        }
        return traced_result(RuleId::SimplifyMultiplyByZero, target_before, make_integer(arena_, BigInt(0)));
    }

    std::vector<ExprPtr> flat_factors;
    std::function<Result<void>(ExprPtr, bool)> flatten = [&](ExprPtr f, bool invert) -> Result<void> {
        if (const auto* prod = expr_cast<Product>(f)) {
            for (ExprPtr factor : prod->factors) {
                auto res = flatten(factor, invert);
                if (res.is_error()) return res;
            }
        } else if (const auto* div = expr_cast<Binary>(f); div != nullptr && div->op == BinaryOp::Div) {
            if (auto res_l = flatten(div->left, invert); res_l.is_error()) return res_l;
            if (auto res_r = flatten(div->right, !invert); res_r.is_error()) return res_r;
        } else if (invert) {
            auto inv = simplify_power(f, make_integer(arena_, BigInt(-1)));
            if (inv.is_error()) return fail<void>(inv.error());
            flat_factors.push_back(inv.value());
        } else {
            flat_factors.push_back(f);
        }
        return ok();
    };
    for (ExprPtr f : initial_factors) if (auto res = flatten(f, false); res.is_error()) return fail<ExprPtr>(res.error());

    /*
    for (std::size_t i = 0; i < flat_factors.size(); ++i) {
        const auto* sum = expr_cast<Sum>(flat_factors[i]);
        if (sum == nullptr) continue;
        std::vector<ExprPtr> dist_terms;
        for (ExprPtr term : sum->terms) {
            std::vector<ExprPtr> dist_factors;
            for (std::size_t j = 0; j < flat_factors.size(); ++j) dist_factors.push_back(j == i ? term : flat_factors[j]);
            auto dist = simplify_product_factors(dist_factors, ExprPtr{}, true);
            if (dist.is_error()) return dist;
            dist_terms.push_back(dist.value());
        }
        auto rewritten = simplify_sum_terms(dist_terms, make_sum_target(dist_terms));
        if (rewritten.is_error()) return rewritten;
        append_trace(RuleId::SimplifyDistributeProductOverSum, target_before, rewritten.value());
        return rewritten;
    }
    */

    Rational coefficient(BigInt(1));
    std::vector<std::pair<ExprPtr, BigInt>> symbolic;
    BigInt i_count(0);
    for (ExprPtr f : flat_factors) {
        if (auto timeout = check_timeout(); timeout.is_error()) return fail<ExprPtr>(timeout.error());
        LiteralRational rat;
        auto exact = try_get_exact_rational(f, rat);
        if (exact.is_ok() && exact.value()) { coefficient *= rat.value; continue; }
        if (const auto* unary = expr_cast<Unary>(f); unary != nullptr && unary->op == UnaryOp::Neg) { coefficient *= Rational(BigInt(-1)); f = unary->operand; }
        
        if (is_constant_expr(f, MathConstant::I)) {
            i_count += BigInt(1);
            continue;
        }

        if (const auto* binary = expr_cast<Binary>(f); binary != nullptr && binary->op == BinaryOp::Pow) {
            if (is_constant_expr(binary->left, MathConstant::I)) {
                if (auto exponent = try_get_integer_exponent(binary->right); exponent.has_value()) {
                    i_count += *exponent;
                    continue;
                }
            }
            if (auto exponent = try_get_integer_exponent(binary->right); exponent.has_value()) { symbolic.push_back({binary->left, *exponent}); continue; }
        }
        symbolic.push_back({f, BigInt(1)});
    }
    
    if (!i_count.is_zero()) {
        BigInt rem = i_count % BigInt(4);
        if (rem.is_negative()) rem += BigInt(4);
        long long r = rem.to_u64();
        if (r == 1) symbolic.push_back({arena_.make<Constant>(MathConstant::I), BigInt(1)});
        else if (r == 2) coefficient *= Rational(BigInt(-1));
        else if (r == 3) {
            coefficient *= Rational(BigInt(-1));
            symbolic.push_back({arena_.make<Constant>(MathConstant::I), BigInt(1)});
        }
    }

    if (coefficient.numerator().is_zero()) return traced_result(RuleId::SimplifyMultiplyByZero, target_before, make_integer(arena_, BigInt(0)));
    merge_symbolic_factors(symbolic);

    // GAP #4: Branch Cuts / Sqrt Product Merging
    // Collect sqrt calls to merge them safely: sqrt(A)*sqrt(B) -> sqrt(A*B) if A>=0 or B>=0.
    std::vector<ExprPtr> sqrt_args;
    std::vector<std::pair<ExprPtr, BigInt>> other_symbolic;
    for (auto& factor : symbolic) {
        if (factor.second == BigInt(1)) {
            if (const auto* call = expr_cast<FuncCall>(factor.first);
                call != nullptr && call->func_id == BuiltinOp::Sqrt && call->args.size() == 1U) {
                sqrt_args.push_back(call->args.front());
                continue;
            }
        }
        other_symbolic.push_back(std::move(factor));
    }

    if (sqrt_args.size() >= 2U) {
        // Safe to merge if at most one is not known non-negative.
        std::size_t known_nonnegative_count = 0;
        for (const auto& arg : sqrt_args) {
            if (is_known_nonnegative(arg)) {
                known_nonnegative_count++;
            }
        }

        if (known_nonnegative_count >= sqrt_args.size() - 1) {
            auto inner_res = simplify_product_factors(sqrt_args, ExprPtr{}, true);
            if (inner_res.is_ok()) {
                ExprPtr merged_sqrt = arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{inner_res.value()});
                other_symbolic.push_back({merged_sqrt, BigInt(1)});
            }
        } else {
            // Merge not safe: put sqrt factors back
            for (ExprPtr arg : sqrt_args) {
                other_symbolic.push_back({arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{arg}), BigInt(1)});
            }
        }
    } else if (!sqrt_args.empty()) {
        // Single sqrt: put it back
        other_symbolic.push_back({arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{sqrt_args.front()}), BigInt(1)});
    }
    // Always restore symbolic from other_symbolic (move invalidated original entries)
    symbolic = std::move(other_symbolic);
    merge_symbolic_factors(symbolic);

    std::vector<ExprPtr> normalized;
    bool is_neg = (coefficient == Rational(BigInt(-1)));
    if (!is_neg && (!(coefficient == Rational(BigInt(1))) || symbolic.empty())) {
        normalized.push_back(make_rational(arena_, coefficient));
    }
    for (const auto& [base, exp] : symbolic) {
        if (exp.is_zero()) continue;
        normalized.push_back(exp == BigInt(1) ? base : arena_.make<Binary>(BinaryOp::Pow, base, make_integer(arena_, exp)));
    }
    normalized.erase(std::remove_if(normalized.begin(), normalized.end(), [](ExprPtr e) { return is_one_expr(e); }), normalized.end());
    
    if (normalized.empty()) {
        return traced_result(RuleId::SimplifyMultiplyByOne, target_before, make_rational(arena_, coefficient));
    }

    ExprPtr result;
    if (normalized.size() == 1U) {
        result = normalized.front();
    } else {
        result = arena_.make<Product>(std::move(normalized));
    }

    if (is_neg) {
        return traced_result(RuleId::SimplifyMultiplyByOne, target_before, arena_.make<Unary>(UnaryOp::Neg, result));
    }

    return ok(result);
}

} // namespace cas::symbolic::detail
