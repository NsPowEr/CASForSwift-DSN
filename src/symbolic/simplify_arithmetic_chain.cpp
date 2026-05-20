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
    bool has_infinity = false;
    bool has_neg_infinity = false;

    for (ExprPtr term : flat_terms) {
        auto timeout = check_timeout();
        if (timeout.is_error()) return fail<ExprPtr>(timeout.error());

        if (is_constant_expr(term, MathConstant::Infinity)) { has_infinity = true; continue; }
        if (const auto* u = expr_cast<Unary>(term); u != nullptr && u->op == UnaryOp::Neg) {
            if (is_constant_expr(u->operand, MathConstant::Infinity)) {
                has_neg_infinity = true; continue;
            }
        }

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

    if (has_infinity && has_neg_infinity) return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "Infinity - Infinity is undefined"));
    if (has_infinity) return ok(arena_.make<Constant>(MathConstant::Infinity));
    if (has_neg_infinity) return ok(arena_.make<Unary>(UnaryOp::Neg, arena_.make<Constant>(MathConstant::Infinity)));

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
            if (is_constant_expr(factor, MathConstant::Infinity)) return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "0 * Infinity is undefined"));
            if (const auto* u = expr_cast<Unary>(factor); u != nullptr && u->op == UnaryOp::Neg) {
                if (is_constant_expr(u->operand, MathConstant::Infinity)) return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "0 * Infinity is undefined"));
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

    Rational coefficient(BigInt(1));
    std::vector<std::pair<ExprPtr, BigInt>> symbolic;
    BigInt i_count(0);
    bool has_infinity = false;
    int infinity_sign = 1;

    for (ExprPtr f : flat_factors) {
        if (auto timeout = check_timeout(); timeout.is_error()) return fail<ExprPtr>(timeout.error());
        
        int sign = 1;
        ExprPtr core = f;
        if (const auto* u = expr_cast<Unary>(f); u != nullptr && u->op == UnaryOp::Neg) {
            sign = -1;
            core = u->operand;
        }

        if (is_constant_expr(core, MathConstant::Infinity)) {
            has_infinity = true;
            infinity_sign *= sign;
            continue;
        }

        LiteralRational rat;
        auto exact = try_get_exact_rational(f, rat);
        if (exact.is_ok() && exact.value()) { 
            coefficient *= rat.value; 
            continue; 
        }
        
        if (const auto* unary = expr_cast<Unary>(f); unary != nullptr && unary->op == UnaryOp::Neg) { 
            coefficient *= Rational(BigInt(-1)); 
            f = unary->operand; 
        }
        
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
            if (auto exponent = try_get_integer_exponent(binary->right); exponent.has_value()) { 
                symbolic.push_back({binary->left, *exponent}); 
                continue; 
            }
        }
        symbolic.push_back({f, BigInt(1)});
    }

    if (has_infinity) {
        if (coefficient.numerator().is_negative()) infinity_sign *= -1;
        // Check if i_count affects sign (I^2 = -1)
        BigInt rem = i_count % BigInt(4);
        if (rem == BigInt(2) || rem == BigInt(3)) infinity_sign *= -1;
        
        ExprPtr inf = arena_.make<Constant>(MathConstant::Infinity);
        if (infinity_sign < 0) return ok(arena_.make<Unary>(UnaryOp::Neg, inf));
        return ok(inf);
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

    // L3-04 Gamma reflection identity (and 0-sum variant):
    //   Γ(z) · Γ(1 − z) = π / sin(π z)
    //   Γ(z) · Γ(−z)    = −π / (z · sin(π z))     (derived: Γ(z+1)=z·Γ(z) +
    //                                              reflection on (z+1, −z))
    // Scan symbolic factors (exponent 1) for pairs of Gamma calls whose
    // argument sum is an integer m ∈ {0, 1} and rewrite via the
    // corresponding closed form.  Pure algorithmic identity, no lookup.
    {
        auto gamma_arg = [](const std::pair<ExprPtr, BigInt>& p) -> ExprPtr {
            if (p.second != BigInt(1)) return nullptr;
            const auto* fc = expr_cast<FuncCall>(p.first);
            if (fc == nullptr || fc->func_id != BuiltinOp::Gamma || fc->args.size() != 1U) {
                return nullptr;
            }
            return fc->args[0];
        };
        bool reflected_any = false;
        bool keep_scanning = true;
        while (keep_scanning) {
            keep_scanning = false;
            std::vector<std::size_t> gamma_idx;
            for (std::size_t i = 0; i < symbolic.size(); ++i) {
                if (gamma_arg(symbolic[i]) != nullptr) gamma_idx.push_back(i);
            }
            for (std::size_t a = 0; a < gamma_idx.size() && !keep_scanning; ++a) {
                for (std::size_t b = a + 1; b < gamma_idx.size() && !keep_scanning; ++b) {
                    ExprPtr za = gamma_arg(symbolic[gamma_idx[a]]);
                    ExprPtr zb = gamma_arg(symbolic[gamma_idx[b]]);
                    if (za == nullptr || zb == nullptr) continue;
                    ExprPtr sum_zz = arena_.make<Sum>(std::vector<ExprPtr>{za, zb});
                    auto sum_simp = simplify_expr(sum_zz);
                    if (sum_simp.is_error()) continue;
                    const auto* m_lit = expr_cast<IntegerLit>(sum_simp.value());
                    if (m_lit == nullptr) continue;
                    // Only m ∈ {0, 1} have a closed form via reflection
                    // alone — higher |m| would need a finite descending
                    // chain via the functional equation, leave inert.
                    const bool m_is_one = (m_lit->value == BigInt(1));
                    const bool m_is_zero = m_lit->value.is_zero();
                    if (!m_is_one && !m_is_zero) continue;
                    ExprPtr pi_const = arena_.make<Constant>(MathConstant::Pi);
                    ExprPtr pi_z = arena_.make<Product>(std::vector<ExprPtr>{pi_const, za});
                    auto pi_z_simp = simplify_expr(pi_z);
                    if (pi_z_simp.is_error()) continue;
                    ExprPtr sin_call = arena_.make<FuncCall>(
                        BuiltinOp::Sin, std::vector<ExprPtr>{pi_z_simp.value()});
                    auto sin_simp = simplify_expr(sin_call);
                    if (sin_simp.is_error()) continue;
                    std::size_t ia = gamma_idx[a];
                    std::size_t ib = gamma_idx[b];
                    if (ia > ib) std::swap(ia, ib);
                    symbolic.erase(symbolic.begin() + ib);
                    symbolic.erase(symbolic.begin() + ia);
                    symbolic.push_back({pi_const, BigInt(1)});
                    symbolic.push_back({sin_simp.value(), BigInt(-1)});
                    if (m_is_zero) {
                        // Extra factor: -1 / z_a.
                        coefficient *= Rational(BigInt(-1));
                        symbolic.push_back({za, BigInt(-1)});
                    }
                    reflected_any = true;
                    keep_scanning = true;
                }
            }
        }
        if (reflected_any) merge_symbolic_factors(symbolic);
    }

    // L2-07: Double-angle compaction — sin(x)·cos(x) → (1/2)·sin(2x).
    // Fires at most once per product pass; sin(x)^n·cos(x)^m with n,m≥1 not handled here
    // (would require exponent matching beyond BigInt(1)).  Clean form: C·sin·cos = (C/2)·sin(2x).
    {
        std::optional<std::pair<std::size_t, std::size_t>> sc_pair;
        ExprPtr sc_arg = nullptr;
        for (std::size_t ii = 0; ii < symbolic.size() && !sc_pair; ++ii) {
            if (symbolic[ii].second != BigInt(1)) continue;
            const auto* fi = expr_cast<FuncCall>(symbolic[ii].first);
            if (!fi || fi->func_id != BuiltinOp::Sin || fi->args.size() != 1U) continue;
            for (std::size_t jj = ii + 1; jj < symbolic.size(); ++jj) {
                if (symbolic[jj].second != BigInt(1)) continue;
                const auto* fj = expr_cast<FuncCall>(symbolic[jj].first);
                if (!fj || fj->func_id != BuiltinOp::Cos || fj->args.size() != 1U) continue;
                if (structural_equal(fi->args[0], fj->args[0])) {
                    sc_pair = {ii, jj};
                    sc_arg = fi->args[0];
                    break;
                }
            }
        }
        if (sc_pair) {
            const auto [si, ci] = *sc_pair;
            coefficient *= Rational(BigInt(1), BigInt(2));
            ExprPtr two_arg = arena_.make<Binary>(BinaryOp::Mul,
                make_integer(arena_, BigInt(2)), sc_arg);
            ExprPtr sin2x = arena_.make<FuncCall>(BuiltinOp::Sin,
                std::vector<ExprPtr>{two_arg});
            // Remove sin(x) and cos(x) from symbolic (erase larger index first)
            symbolic.erase(symbolic.begin() + ci);
            symbolic.erase(symbolic.begin() + si);
            auto sin2x_s = simplify_expr(sin2x);
            symbolic.push_back({sin2x_s.is_ok() ? sin2x_s.value() : sin2x, BigInt(1)});
        }
    }

    // L1-12 strengthening: merge sqrt(rat_pos) · sqrt(rat_pos) → sqrt(prod).
    // Strict to non-negative rational arguments to avoid branch-cut
    // violations (sqrt(-1)·sqrt(-1) ≠ sqrt(1)).
    {
        auto get_rat_pos = [](ExprPtr e) -> std::optional<Rational> {
            const auto* fc = expr_cast<FuncCall>(e);
            if (!fc || fc->func_id != BuiltinOp::Sqrt || fc->args.size() != 1U)
                return std::nullopt;
            if (const auto* il = expr_cast<IntegerLit>(fc->args[0])) {
                if (il->value.is_negative()) return std::nullopt;
                return Rational(il->value, BigInt(1));
            }
            if (const auto* rl = expr_cast<RationalLit>(fc->args[0])) {
                if (rl->numerator.is_negative()) return std::nullopt;
                return Rational(rl->numerator, rl->denominator);
            }
            return std::nullopt;
        };
        bool merged_any = true;
        while (merged_any) {
            merged_any = false;
            for (std::size_t i = 0; i < symbolic.size() && !merged_any; ++i) {
                if (symbolic[i].second != BigInt(1)) continue;
                auto ra = get_rat_pos(symbolic[i].first);
                if (!ra) continue;
                for (std::size_t j = i + 1; j < symbolic.size(); ++j) {
                    if (symbolic[j].second != BigInt(1)) continue;
                    auto rb = get_rat_pos(symbolic[j].first);
                    if (!rb) continue;
                    Rational prod = (*ra) * (*rb);
                    ExprPtr arg = (prod.denominator() == BigInt(1))
                        ? static_cast<ExprPtr>(arena_.make<IntegerLit>(prod.numerator()))
                        : static_cast<ExprPtr>(arena_.make<RationalLit>(
                            prod.numerator(), prod.denominator()));
                    ExprPtr new_sqrt_raw = arena_.make<FuncCall>(BuiltinOp::Sqrt,
                        std::vector<ExprPtr>{arg});
                    auto new_sqrt = simplify_expr(new_sqrt_raw);
                    ExprPtr replacement = new_sqrt.is_ok() ? new_sqrt.value() : new_sqrt_raw;
                    symbolic.erase(symbolic.begin() + j);
                    symbolic.erase(symbolic.begin() + i);
                    LiteralRational rep_rat;
                    auto rep_check = try_get_exact_rational(replacement, rep_rat);
                    if (rep_check.is_ok() && rep_check.value()) {
                        coefficient *= rep_rat.value;
                    } else {
                        symbolic.push_back({replacement, BigInt(1)});
                    }
                    merged_any = true;
                    break;
                }
            }
        }
    }

    // Distribute over Sum
    ExprPtr sum_factor = nullptr;
    std::size_t sum_idx = 0;
    bool found_sum = false;
    for (std::size_t i = 0; i < symbolic.size(); ++i) {
        if (symbolic[i].second == BigInt(1) && expr_is<Sum>(symbolic[i].first)) {
            sum_factor = symbolic[i].first;
            sum_idx = i;
            found_sum = true;
            break;
        }
    }

    if (found_sum || (!(coefficient == Rational(BigInt(1))) && std::any_of(symbolic.begin(), symbolic.end(), [](const auto& p) { return expr_is<Sum>(p.first) && p.second == BigInt(1); }))) {
        // Re-find if we only found it via coefficient check
        if (!found_sum) {
            for (std::size_t i = 0; i < symbolic.size(); ++i) {
                if (symbolic[i].second == BigInt(1) && expr_is<Sum>(symbolic[i].first)) {
                    sum_factor = symbolic[i].first;
                    sum_idx = i;
                    found_sum = true;
                    break;
                }
            }
        }

        if (found_sum) {
            const auto* sum = expr_cast<Sum>(sum_factor);
            std::vector<std::pair<ExprPtr, BigInt>> other_symbolic = symbolic;
            other_symbolic.erase(other_symbolic.begin() + sum_idx);
            
            std::vector<ExprPtr> distributed_terms;
            for (ExprPtr term : sum->terms) {
                std::vector<ExprPtr> factors_for_term;
                if (!(coefficient == Rational(BigInt(1)))) {
                    factors_for_term.push_back(make_rational(arena_, coefficient));
                }
                factors_for_term.push_back(term);
                for (const auto& [base, exp] : other_symbolic) {
                    factors_for_term.push_back(exp == BigInt(1) ? base : arena_.make<Binary>(BinaryOp::Pow, base, make_integer(arena_, exp)));
                }
                
                auto prod = simplify_product_factors(factors_for_term, ExprPtr{}, false);
                if (prod.is_error()) return prod;
                distributed_terms.push_back(prod.value());
            }
            return simplify_sum_terms(distributed_terms, target_before, true);
        }
    }

    std::vector<ExprPtr> normalized;

    bool is_neg = (coefficient == Rational(BigInt(-1)));
    if (!is_neg && (!(coefficient == Rational(BigInt(1))) || symbolic.empty())) normalized.push_back(make_rational(arena_, coefficient));
    for (const auto& [base, exp] : symbolic) {
        if (exp.is_zero()) continue;
        normalized.push_back(exp == BigInt(1) ? base : arena_.make<Binary>(BinaryOp::Pow, base, make_integer(arena_, exp)));
    }
    normalized.erase(std::remove_if(normalized.begin(), normalized.end(), [](ExprPtr e) { return is_one_expr(e); }), normalized.end());
    
    if (normalized.empty()) return traced_result(RuleId::SimplifyMultiplyByOne, target_before, make_rational(arena_, coefficient));

    ExprPtr result;
    if (normalized.size() == 1U) result = normalized.front();
    else result = arena_.make<Product>(std::move(normalized));

    if (is_neg) return traced_result(RuleId::SimplifyMultiplyByOne, target_before, arena_.make<Unary>(UnaryOp::Neg, result));
    return ok(result);
}

} // namespace cas::symbolic::detail
