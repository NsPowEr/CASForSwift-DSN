#include "simplify_impl.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/algebra.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <algorithm>

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const IntegerLit&) {
    return ok(original);
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const ComplexLit&) {
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

        LiteralComplex complex_val;
        auto exact_c = try_get_exact_complex(operand.value(), complex_val);
        if (exact_c.is_error()) return fail<ExprPtr>(exact_c.error());
        if (exact_c.value()) {
            return traced_result(RuleId::SimplifyCollectLikeTerms, target_before,
                                 make_complex(arena_, -complex_val.value));
        }

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
    case BinaryOp::Less:
    case BinaryOp::Greater:
    case BinaryOp::LessEqual:
    case BinaryOp::GreaterEqual:
        if (lhs.value() == node.left && rhs.value() == node.right) return ok(original);
        return ok(arena_.make<Binary>(node.op, lhs.value(), rhs.value()));
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

// `Simplifier::simplify_power` (+ Chebyshev/DeMoivre trig power linearization)
// estratto in `simplify_arithmetic_power.cpp` per HC-F8-MONOLITH-WAIVER
// tier-1 split (2026-06-11). Vedi `ANTI_MONOLITHIC_REPORT.md`.


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
        // Canonical marker for the imaginary unit in monomial keys: a hash-
        // consed Constant::I.  Both `ComplexLit(0, n/d)` and `Constant::I`
        // map to the same marker, so monomials that differ only by the sign
        // or magnitude of a purely-imaginary factor collapse correctly in
        // the Sum collector.
        const ExprPtr i_marker = arena_.make<Constant>(MathConstant::I);
        // num^e for integer e (e may be negative): reciprocal when e < 0.
        auto pow_rat_int = [](const Rational& b, const BigInt& e) -> Rational {
            const bool neg = e.is_negative();
            Rational p = pow_rational_nonnegative(b, neg ? -e : e);
            return neg ? Rational(p.denominator(), p.numerator()) : p;
        };
        for (ExprPtr factor : product->factors) {
            // T-054: pull rational numeric factors out of a (c·rest)^e base so
            // (2·√R)^-1 and (√R)^-1 produce the SAME monomial key — the numeric
            // part c^e folds into the coefficient. (c·e)^n = c^n·e^n is exact for
            // integer n. This is SAFE at the Sum/like-term layer (it only reshapes
            // coefficients); doing it in simplify_power instead breaks Gruntz mrv,
            // which needs the single un-split Pow node (see simplify_arithmetic_power).
            if (const auto* pw = expr_cast<Binary>(factor);
                pw != nullptr && pw->op == BinaryOp::Pow) {
                if (auto e = try_get_integer_exponent(pw->right); e.has_value()) {
                    if (const auto* pbase = expr_cast<Product>(pw->left)) {
                        Rational num(BigInt(1));
                        std::vector<ExprPtr> rest;
                        rest.reserve(pbase->factors.size());
                        for (ExprPtr pf : pbase->factors) {
                            LiteralRational lr;
                            auto ex = try_get_exact_rational(pf, lr);
                            if (ex.is_ok() && ex.value() && !lr.value.numerator().is_zero())
                                num *= lr.value;
                            else
                                rest.push_back(pf);
                        }
                        // Restrict to SURD bases: every residual factor must be a
                        // sqrt(...). Pulling a numeric factor out of a polynomial /
                        // exp-tower base (e.g. (c·wᵏ)^n) re-keys mrv intermediates and
                        // breaks Gruntz (AcidTest Test1 → ComplexRational div by zero).
                        // The integration artifacts that need this normalization
                        // (∫xᵏ/√(…), ∫asin, …) always have radical denominators.
                        bool all_sqrt = true;
                        for (ExprPtr rf : rest) {
                            const auto* fc = expr_cast<FuncCall>(rf);
                            if (fc == nullptr || fc->func_id != BuiltinOp::Sqrt) {
                                all_sqrt = false;
                                break;
                            }
                        }
                        if (!(num == Rational(BigInt(1))) && !rest.empty() && all_sqrt) {
                            coefficient *= pow_rat_int(num, *e);
                            ExprPtr rest_base = rest.size() == 1U
                                ? rest[0]
                                : arena_.make<Product>(std::move(rest));
                            factor = arena_.make<Binary>(BinaryOp::Pow, rest_base,
                                make_integer(arena_, *e));
                        }
                    }
                }
            }
            auto factor_exact = try_get_exact_rational(factor, rational);
            if (factor_exact.is_ok() && factor_exact.value()) {
                coefficient *= rational.value;
                continue;
            }
            // Purely-imaginary ComplexLit: split into rational coefficient
            // and a single shared `I` marker so monomial keys are stable
            // under sign/magnitude variations.
            if (const auto* cl = expr_cast<ComplexLit>(factor);
                cl != nullptr && cl->re_num.is_zero() && !cl->im_num.is_zero()) {
                auto rat = Rational::make(cl->im_num, cl->im_den);
                if (rat.is_ok()) {
                    coefficient *= rat.value();
                    bool found = false;
                    for (auto& f : factors) {
                        if (f.first == i_marker) { f.second += BigInt(1); found = true; break; }
                    }
                    if (!found) factors.push_back({i_marker, BigInt(1)});
                    continue;
                }
            }
            // Real-only ComplexLit folds into coefficient.
            if (const auto* cl = expr_cast<ComplexLit>(factor);
                cl != nullptr && cl->im_num.is_zero() && !cl->re_num.is_zero()) {
                auto rat = Rational::make(cl->re_num, cl->re_den);
                if (rat.is_ok()) {
                    coefficient *= rat.value();
                    continue;
                }
            }
            // Constant::I shares the marker with imaginary ComplexLit so
            // mixed-form products (e.g. `Pi * I` vs `Pi * ComplexLit(0,1)`)
            // collapse to the same key.
            if (const auto* c = expr_cast<Constant>(factor);
                c != nullptr && c->value == MathConstant::I) {
                bool found = false;
                for (auto& f : factors) {
                    if (f.first == i_marker) { f.second += BigInt(1); found = true; break; }
                }
                if (!found) factors.push_back({i_marker, BigInt(1)});
                continue;
            }
            if (expr_is<Symbol>(factor)) {
                factors.push_back({factor, BigInt(1)});
                continue;
            }
            if (const auto* power = expr_cast<Binary>(factor); power != nullptr && power->op == BinaryOp::Pow) {
                if (auto exponent = try_get_integer_exponent(power->right); exponent.has_value()) {
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
