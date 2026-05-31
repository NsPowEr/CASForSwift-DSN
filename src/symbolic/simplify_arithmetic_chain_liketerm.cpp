#include "simplify_arithmetic_chain_impl.hpp"
#include <algorithm>

// F1.4: symbolic like-term collection helpers.
//
// Collects a·x + b·x → (a+b)·x for SYMBOLIC coefficients a, b (not just
// rational literals).
//
// Algorithm (coefficient / monomial decomposition):
//   Given a Sum {t₁, t₂, …, tₙ}, for each pair (ti, tj) compute the
//   intersection of their symbolic factor sets.  If the intersection is
//   non-empty AND no residual overlaps the shared monomial (canonicality
//   guard), replace with (qi + qj)·shared_monomial.
//
// The residual-overlap guard is critical for idempotence: without it,
//   x^3 + x → (x^2+1)*x which re-expands to x^3+x on the next simplify.
// The guard rejects any merge where a residual base appears in the shared
// monomial, preventing non-canonical factored forms for polynomial sums.

namespace cas::symbolic::detail {

// Internal struct (definition-local; the public API is the three functions
// declared in simplify_arithmetic_chain_impl.hpp).
struct SymbolicMonomial {
    std::vector<std::pair<ExprPtr, BigInt>> factors;
    ExprPtr coefficient;  // nullptr means coefficient is 1
};

bool decompose_term(
    ExprPtr expr,
    Rational& coeff_out,
    std::vector<std::pair<ExprPtr, BigInt>>& factors_out)
{
    coeff_out = Rational(BigInt(1));
    factors_out.clear();

    // Pure rational scalar
    {
        LiteralRational lr;
        auto ex = try_get_exact_rational(expr, lr);
        if (ex.is_ok() && ex.value()) {
            coeff_out = lr.value;
            return true;
        }
    }

    auto push_factor = [&](ExprPtr f) {
        if (const auto* bin = expr_cast<Binary>(f);
            bin && bin->op == BinaryOp::Pow) {
            if (auto exp = try_get_integer_exponent(bin->right); exp.has_value()) {
                factors_out.push_back({bin->left, *exp});
                return;
            }
        }
        factors_out.push_back({f, BigInt(1)});
    };

    // Single non-numeric factor
    if (expr_is<Symbol>(expr) || expr_is<FuncCall>(expr) || expr_is<Constant>(expr)) {
        push_factor(expr);
        return true;
    }

    if (const auto* u = expr_cast<Unary>(expr); u && u->op == UnaryOp::Neg) {
        coeff_out = Rational(BigInt(-1));
        expr = u->operand;
        LiteralRational lr;
        auto ex = try_get_exact_rational(expr, lr);
        if (ex.is_ok() && ex.value()) { coeff_out *= lr.value; return true; }
        push_factor(expr);
        return true;
    }

    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) {
            LiteralRational lr;
            auto ex = try_get_exact_rational(f, lr);
            if (ex.is_ok() && ex.value()) {
                coeff_out *= lr.value;
            } else {
                push_factor(f);
            }
        }
        std::sort(factors_out.begin(), factors_out.end(),
            [](const auto& a, const auto& b) {
                int c = canonical_compare(a.first, b.first);
                return c != 0 ? c < 0 : a.second < b.second;
            });
        return true;
    }

    if (const auto* bin = expr_cast<Binary>(expr);
        bin && bin->op == BinaryOp::Pow) {
        if (auto exp = try_get_integer_exponent(bin->right); exp.has_value()) {
            factors_out.push_back({bin->left, *exp});
            return true;
        }
    }

    // Cannot decompose further — leave as single opaque factor with exponent 1.
    factors_out.push_back({expr, BigInt(1)});
    return true;
}

ExprPtr build_coeff_monomial(
    ExprPtr coeff_expr,
    Rational numeric_coeff,
    const std::vector<std::pair<ExprPtr, BigInt>>& mono_factors,
    AstArena& arena)
{
    std::vector<ExprPtr> parts;
    bool neg = numeric_coeff.numerator().is_negative();
    Rational abs_coeff = neg ? -numeric_coeff : numeric_coeff;

    if (coeff_expr) {
        if (!(abs_coeff == Rational(BigInt(1)))) {
            parts.push_back(make_rational(arena, abs_coeff));
        }
        parts.push_back(coeff_expr);
    } else {
        if (!(abs_coeff == Rational(BigInt(1))) || mono_factors.empty()) {
            parts.push_back(make_rational(arena, abs_coeff));
        }
    }

    for (const auto& [base, exp] : mono_factors) {
        if (exp == BigInt(1)) {
            parts.push_back(base);
        } else {
            parts.push_back(arena.make<Binary>(BinaryOp::Pow, base,
                make_integer(arena, exp)));
        }
    }

    ExprPtr result;
    if (parts.empty()) {
        result = make_integer(arena, BigInt(1));
    } else if (parts.size() == 1U) {
        result = parts[0];
    } else {
        result = arena.make<Product>(std::move(parts));
    }
    return neg ? arena.make<Unary>(UnaryOp::Neg, result) : result;
}

bool try_merge_symbolic_like_terms(
    std::vector<ExprPtr>& terms,
    AstArena& arena)
{
    const std::size_t n = terms.size();
    for (std::size_t i = 0; i < n; ++i) {
        Rational ci;
        std::vector<std::pair<ExprPtr, BigInt>> fi;
        if (!decompose_term(terms[i], ci, fi)) continue;
        if (fi.empty()) continue;  // pure scalar, skip

        for (std::size_t j = i + 1; j < n; ++j) {
            Rational cj;
            std::vector<std::pair<ExprPtr, BigInt>> fj;
            if (!decompose_term(terms[j], cj, fj)) continue;
            if (fj.empty()) continue;

            // Compute intersection of factor sets (structural equality for
            // interned nodes is pointer comparison — O(1)).
            std::vector<std::pair<ExprPtr, BigInt>> shared;
            for (const auto& [base_i, exp_i] : fi) {
                for (const auto& [base_j, exp_j] : fj) {
                    if (structural_equal(base_i, base_j)) {
                        BigInt min_exp = exp_i < exp_j ? exp_i : exp_j;
                        if (!min_exp.is_zero() && !min_exp.is_negative()) {
                            shared.push_back({base_i, min_exp});
                        }
                        break;
                    }
                }
            }
            if (shared.empty()) continue;

            // Polynomial-context guard: if ALL shared bases are plain Symbols
            // or powers of Symbols, skip.  In that context the monomial
            // collector in Step 4 handles like-term collection; running F1.4
            // on top would produce factored forms (e.g. (b+z)*x) that the
            // distribute step (Step 8) immediately expands back, causing
            // non-confluence with the distributive property tests.
            // F1.4 is only useful when at least one shared base is a
            // transcendental / radical expression (FuncCall / non-Symbol Power).
            {
                bool any_nonpoly_shared = false;
                for (const auto& [sb, _se] : shared) {
                    if (!expr_is<Symbol>(sb)) {
                        any_nonpoly_shared = true;
                        break;
                    }
                }
                if (!any_nonpoly_shared) continue;
            }

            // Compute residual factors after removing the shared part.
            auto subtract_shared = [&](const std::vector<std::pair<ExprPtr, BigInt>>& src)
                -> std::vector<std::pair<ExprPtr, BigInt>>
            {
                std::vector<std::pair<ExprPtr, BigInt>> res;
                for (const auto& [base, exp] : src) {
                    BigInt remaining = exp;
                    for (const auto& [sb, se] : shared) {
                        if (structural_equal(base, sb)) {
                            remaining = remaining - se;
                            break;
                        }
                    }
                    if (!remaining.is_zero()) res.push_back({base, remaining});
                }
                return res;
            };

            auto ri = subtract_shared(fi);
            auto rj = subtract_shared(fj);

            // Canonicality guard: reject merge if any residual base appears in
            // the shared monomial.  Prevents x^3 + x → (x^2+1)*x (would break
            // idempotence and the polynomial coefficient-extraction invariant).
            // E.g. x^3+x: shared={x:1}, ri={x:2}. x is in both → reject.
            auto residual_overlaps_shared =
                [&](const std::vector<std::pair<ExprPtr, BigInt>>& r) -> bool {
                for (const auto& [rb, _r] : r) {
                    for (const auto& [sb, _s] : shared) {
                        if (structural_equal(rb, sb)) return true;
                    }
                }
                return false;
            };
            if (residual_overlaps_shared(ri) || residual_overlaps_shared(rj)) continue;

            // Build quotients qi = ci*∏ri, qj = cj*∏rj.
            auto build_quotient = [&](Rational c,
                const std::vector<std::pair<ExprPtr, BigInt>>& r) -> ExprPtr
            {
                return build_coeff_monomial(nullptr, c, r, arena);
            };

            ExprPtr qi = build_quotient(ci, ri);
            ExprPtr qj = build_quotient(cj, rj);

            // Form (qi + qj) * shared_monomial.
            // Sort the coefficient sum canonically so the result is in the
            // same normal form on repeated simplify passes (idempotency).
            std::vector<ExprPtr> coeff_terms{qi, qj};
            std::sort(coeff_terms.begin(), coeff_terms.end(),
                [](ExprPtr a, ExprPtr b) {
                    int d = polynomial_degree(a) - polynomial_degree(b);
                    return d != 0 ? d > 0 : canonical_compare(a, b) < 0;
                });
            ExprPtr coeff_sum = arena.make<Sum>(std::move(coeff_terms));

            // Build the list of factor expressions for the merged product.
            // Each entry is (base_for_sort, full_factor_expr).
            // base_for_sort is what merge_symbolic_factors uses as sort key:
            // for a plain expr it is the expr itself; for base^exp it is base.
            // Sorting here mirrors merge_symbolic_factors so that the Product
            // node created below already has canonical factor order.  Without
            // this sort pass-1 produces Product{coeff_sum, base^exp} while
            // pass-2 (via simplify_product_factors → merge_symbolic_factors)
            // produces Product{base^exp, coeff_sum}, violating idempotency.
            struct FactorEntry {
                ExprPtr sort_key;   // base used by canonical_compare
                ExprPtr factor_expr;
            };
            std::vector<FactorEntry> factor_entries;
            factor_entries.push_back({coeff_sum, coeff_sum});
            for (const auto& [base, exp] : shared) {
                ExprPtr full = (exp == BigInt(1))
                    ? base
                    : arena.make<Binary>(BinaryOp::Pow, base, make_integer(arena, exp));
                factor_entries.push_back({base, full});
            }
            // Stable sort preserves relative order among equal-key entries,
            // matching the stable deterministic output of merge_symbolic_factors.
            std::stable_sort(factor_entries.begin(), factor_entries.end(),
                [](const FactorEntry& lhs, const FactorEntry& rhs) {
                    return canonical_compare(lhs.sort_key, rhs.sort_key) < 0;
                });
            std::vector<ExprPtr> merged_parts;
            merged_parts.reserve(factor_entries.size());
            for (const auto& fe : factor_entries)
                merged_parts.push_back(fe.factor_expr);
            ExprPtr merged = merged_parts.size() == 1U
                ? merged_parts[0]
                : arena.make<Product>(std::move(merged_parts));

            terms[i] = merged;
            terms.erase(terms.begin() + static_cast<std::ptrdiff_t>(j));
            return true;
        }
    }
    return false;
}

} // namespace cas::symbolic::detail
