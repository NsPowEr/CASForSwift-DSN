#include "simplify_arithmetic_chain_impl.hpp"
#include <algorithm>
#include <map>

// simplify_node(Sum) + simplify_sum_terms implementation.
// Product handling lives in simplify_arithmetic_chain.cpp.

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Sum& node) {
    return simplify_sum_terms(node.terms, original);
}

Result<ExprPtr> Simplifier::simplify_sum_terms(
    const std::vector<ExprPtr>& terms,
    ExprPtr target_before,
    bool inputs_are_simplified)
{
    if (!target_before && trace_enabled_) target_before = make_sum_target(terms);

    // Step 1: flatten and recursively simplify each term.
    std::vector<ExprPtr> flat_terms;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        auto simplify_current = [&]() -> Result<ExprPtr> {
            return inputs_are_simplified ? ok(terms[i]) : simplify_expr(terms[i]);
        };
        Result<ExprPtr> simplified = fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
        if (trace_enabled_) {
            std::vector<ExprPtr> current_terms = terms;
            for (std::size_t j = 0; j < flat_terms.size() && j < current_terms.size(); ++j)
                current_terms[j] = flat_terms[j];
            ScopedFrame frame(*this,
                [this, current_terms = std::move(current_terms), i](ExprPtr value) mutable {
                    current_terms[i] = value;
                    return make_sum_target(current_terms);
                });
            simplified = simplify_current();
        } else {
            simplified = simplify_current();
        }
        if (simplified.is_error()) return simplified;
        if (const auto* nested = expr_cast<Sum>(simplified.value())) {
            flat_terms.insert(flat_terms.end(),
                nested->terms.begin(), nested->terms.end());
        } else {
            flat_terms.push_back(simplified.value());
        }
    }

    // Step 2: optional rewrite-provider pass.
    while (rewrite_provider_ != nullptr
        && flat_terms.size() >= 2U
        && may_rewrite_sum_terms(flat_terms))
    {
        ExprPtr rewrite_target = make_sum_target(flat_terms);
        auto rewritten = rewrite_provider_->try_rewrite(
            rewrite_target, arena_, assumptions_, context_);
        if (rewritten.is_error()) return rewritten;
        if (rewritten.value() == rewrite_target) break;
        append_trace(RuleId::RewriteProviderApplied, rewrite_target, rewritten.value());
        if (const auto* rewritten_sum = expr_cast<Sum>(rewritten.value())) {
            flat_terms.assign(
                rewritten_sum->terms.begin(), rewritten_sum->terms.end());
        } else {
            return simplify_expr(rewritten.value());
        }
    }

    // Step 3: L3-08 Quantity addition — group by SI dimension, sum values.
    // F6.6-T1: enforce SI dimensional analysis on the Sum.  A term mixing
    // two distinct SI dimensions (e.g. `1·m + 1·s`) violates the physical
    // law that only commensurable quantities can be summed.  Once detected
    // we surface CASErrorKind::Undefined with the offending dimension pair
    // so callers can react instead of silently producing a half-collapsed
    // Quantity expression.
    {
        std::map<SIDimensions, std::vector<ExprPtr>> by_dim;
        std::vector<ExprPtr> non_qty;
        for (ExprPtr t : flat_terms) {
            if (auto* q = expr_cast<Quantity>(t)) {
                by_dim[q->dimensions].push_back(q->value);
            } else {
                non_qty.push_back(t);
            }
        }
        if (by_dim.size() > 1) {
            return fail<ExprPtr>(make_error(CASErrorKind::Undefined,
                "Sum mixes incompatible SI dimensions — addition is "
                "undefined across distinct quantity dimensions (F6.6)"));
        }
        if (!by_dim.empty() && by_dim.begin()->second.size() > 1)
        {
            std::vector<ExprPtr> rebuilt = std::move(non_qty);
            for (auto& [dim, vals] : by_dim) {
                ExprPtr summed;
                if (vals.size() == 1U) {
                    summed = vals[0];
                } else {
                    auto s = simplify_sum_terms(vals, ExprPtr{}, false);
                    summed = s.is_ok() ? s.value() : arena_.make<Sum>(std::move(vals));
                }
                rebuilt.push_back(arena_.make<Quantity>(summed, dim));
            }
            flat_terms = std::move(rebuilt);
        }
    }
// Step 4: collect numeric-coefficient like-terms into MonomialKey map.
ComplexRational constant = ComplexRational::zero();
std::map<MonomialKey, Rational> collected;
bool has_infinity = false;
bool has_neg_infinity = false;

for (ExprPtr term : flat_terms) {
    auto timeout = check_timeout();
    if (timeout.is_error()) return fail<ExprPtr>(timeout.error());

    if (is_constant_expr(term, MathConstant::Infinity)) {
        has_infinity = true; continue;
    }
    if (const auto* u = expr_cast<Unary>(term);
        u != nullptr && u->op == UnaryOp::Neg)
    {
        if (is_constant_expr(u->operand, MathConstant::Infinity)) {
            has_neg_infinity = true; continue;
        }
    }

    LiteralComplex complex;
    auto exact = try_get_exact_complex(term, complex);
    if (exact.is_ok() && exact.value()) {
        constant = constant + complex.value;
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

    if (has_infinity && has_neg_infinity)
        return fail<ExprPtr>(make_error(CASErrorKind::Undefined,
            "Infinity - Infinity is undefined"));
    if (has_infinity)  return ok(arena_.make<Constant>(MathConstant::Infinity));
    if (has_neg_infinity)
        return ok(arena_.make<Unary>(UnaryOp::Neg,
            arena_.make<Constant>(MathConstant::Infinity)));

    std::vector<ExprPtr> normalized;
    for (const auto& [key, coeff] : collected) {
        if (!coeff.numerator().is_zero())
            normalized.push_back(build_monomial(key, coeff));
    }
    if (!constant.is_zero() || (collected.empty() && normalized.empty())) {
        normalized.insert(normalized.begin(), make_complex(arena_, constant));
    }

    // Step 5: F1.4 symbolic like-term pass (a·x + b·x → (a+b)·x).
    // Numeric like-terms already merged above; this handles symbolic coefficients.
    // The residual-overlap guard inside try_merge_symbolic_like_terms prevents
    // merging polynomial terms that differ only in degree (e.g. x^3 + x).
    //
    // NOTE: We do NOT call simplify_expr on the merged product here — doing so
    // would re-enter simplify_product_factors and trigger the distribute step,
    // which would expand (a+b)*x back to a*x + b*x, creating an infinite
    // factor/distribute cycle.  Idempotency is instead ensured by sorting
    // the coefficient Sum's terms canonically inside try_merge_symbolic_like_terms.
    {
        bool merged = true;
        int pass_limit = static_cast<int>(normalized.size()) + 1;
        while (merged && pass_limit-- > 0) {
            merged = try_merge_symbolic_like_terms(normalized, arena_);
        }
    }

    if (normalized.empty())
        return traced_result(RuleId::SimplifyCollectLikeTerms,
            target_before, make_integer(arena_, BigInt(0)));

    std::sort(normalized.begin(), normalized.end(), [](ExprPtr lhs, ExprPtr rhs) {
        int l_deg = polynomial_degree(lhs), r_deg = polynomial_degree(rhs);
        return l_deg != r_deg ? l_deg > r_deg : canonical_compare(lhs, rhs) < 0;
    });

    if (normalized.size() == 1U) {
        if (std::any_of(terms.begin(), terms.end(),
                [](ExprPtr expr) { return is_zero_expr(expr); }))
            return traced_result(RuleId::SimplifyAddZero,
                target_before, normalized.front());
        return traced_result(RuleId::SimplifyCollectLikeTerms,
            target_before, normalized.front());
    }

    if (expr_is<Sum>(target_before)
        && expr_ptr_sequence_identical(
            normalized, expr_ref<Sum>(target_before).terms))
        return ok(target_before);

    return ok(arena_.make<Sum>(std::move(normalized)));
}

} // namespace cas::symbolic::detail
