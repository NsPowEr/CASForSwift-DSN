#include "polynomial_internal.hpp"

#include <optional>
#include <vector>

namespace cas::algebra {
namespace {

[[nodiscard]] BigInt mod_positive(const BigInt& value, const BigInt& modulus) {
    BigInt result = value % modulus;
    if (result.is_negative()) {
        result += modulus.abs();
    }
    return result;
}

[[nodiscard]] IntPoly multiply_mod(const IntPoly& lhs, const IntPoly& rhs, const BigInt& modulus) {
    if (lhs.is_zero() || rhs.is_zero()) {
        return IntPoly{};
    }

    IntPoly result;
    result.resize(lhs.size() + rhs.size() - 1U, BigInt(0));
    for (std::size_t i = 0U; i < lhs.size(); ++i) {
        if (lhs[i].is_zero()) {
            continue;
        }
        for (std::size_t j = 0U; j < rhs.size(); ++j) {
            result[i + j] = mod_positive(result[i + j] + lhs[i] * rhs[j], modulus);
        }
    }
    normalize_integer_poly(result);
    return result;
}

[[nodiscard]] IntPoly center_modular_coefficients(IntPoly poly, const BigInt& modulus) {
    for (BigInt& coefficient : poly.coefficients()) {
        coefficient = mod_positive(coefficient, modulus);
        if (coefficient * BigInt(2) > modulus) {
            coefficient -= modulus;
        }
    }
    return primitive_integer_poly(std::move(poly));
}

[[nodiscard]] bool divides_exactly(const IntPoly& dividend, const IntPoly& divisor) {
    if (divisor.is_zero() || divisor.degree() == 0U || divisor.degree() >= dividend.degree()) {
        return false;
    }

    auto remainder = pseudo_remainder_integer_poly(dividend, divisor);
    normalize_integer_poly(remainder);
    return remainder.is_zero();
}

// Landau-Mignotte factor bound (Knuth TAOCP vol2 §4.6.2 Theorem F):
// if h | f in Z[x] with deg(h) ≤ d_h and ||f||_inf = c, then
//
//     ||h||_inf ≤ 2^d_h · c
//
// This is the strict envelope used to prune candidate subsets whose
// lifted product cannot be a Z-factor of `f`. Using BigInt arithmetic
// throughout to avoid overflow on high-degree inputs.
[[nodiscard]] BigInt mignotte_factor_bound(const IntPoly& f, std::size_t max_factor_degree) {
    BigInt max_abs(0);
    for (std::size_t k = 0; k < f.size(); ++k) {
        const BigInt a = f[k].abs();
        if (a > max_abs) max_abs = a;
    }
    BigInt bound = max_abs;
    for (std::size_t k = 0; k < max_factor_degree; ++k) {
        bound *= BigInt(2);
    }
    // +1 to absorb the strict inequality.
    return bound + BigInt(1);
}

[[nodiscard]] bool exceeds_factor_bound(const IntPoly& candidate, const BigInt& bound) {
    for (std::size_t k = 0; k < candidate.size(); ++k) {
        if (candidate[k].abs() > bound) return true;
    }
    return false;
}

// Zassenhaus subset test on ALREADY-LIFTED factors: the candidate is
//     primitive(center(lc(f) · ∏_{i∈S} g̃_i  mod p^a)),
// where the g̃_i are the Hensel factors mod p^a (∏ g̃_i ≡ f mod p^a is the
// invariant of hensel_lift_multi). No per-subset Hensel lifting — the
// historic implementation re-lifted a (left, right) pair for EVERY subset,
// which (a) cost one full quadratic lift per node of the search tree and
// (b) after dropping large factors from the pool violated
// left·right ≡ f (mod p), silently producing garbage candidates — real
// factors were missed and composites declared irreducible (the
// "recombination misses a factor" silent-wrong found via A6/Φ₇-R₃).
[[nodiscard]] std::optional<IntPoly> check_subset(
    const IntPoly& f,
    const std::vector<IntPoly>& lifted,
    const std::vector<bool>& selected,
    const BigInt& modulus,
    std::size_t max_degree,
    const BigInt& mignotte_bound) {
    IntPoly product(std::vector<BigInt>{mod_positive(f.leading_coeff(), modulus)});
    for (std::size_t index = 0U; index < lifted.size(); ++index) {
        if (selected[index]) {
            product = multiply_mod(product, lifted[index], modulus);
        }
    }
    if (product.is_zero() || product.degree() == 0U ||
        product.degree() > max_degree) {
        return std::nullopt;
    }

    IntPoly candidate = center_modular_coefficients(std::move(product), modulus);
    if (candidate.degree() > max_degree) return std::nullopt;
    // Landau-Mignotte pruning: any factor h of f satisfies
    //     ||h||_inf ≤ 2^deg(h) · ||f||_inf.
    // If the centered candidate exceeds this bound, the subset cannot be a
    // Z-factor of f. Cheap O(deg) test that sidesteps the expensive
    // pseudo-division in divides_exactly.
    if (exceeds_factor_bound(candidate, mignotte_bound)) {
        return std::nullopt;
    }
    if (divides_exactly(f, candidate)) {
        return candidate;
    }

    return std::nullopt;
}

// Subset enumeration with Mignotte-based pruning (Lecerf 2007 §3 +
// Landau-Mignotte). The legacy 32768-subset cap is replaced by the
// natural enumeration bound 2^r where r = |usable_factors|, and the
// Mignotte coefficient bound check inside check_subset prunes the
// exponential tree to polynomial-time in practice on inputs that are
// not pathological Swinnerton-Dyer-style products. (A full van Hoeij
// knapsack-lattice replacement remains as a follow-up for the
// pathological case; see Step 5b in /Users/davidesaba/.claude/plans/.)
//
// Termination: bounded structurally by 2^r enumeration depth + the
// fact that selected_degree > max_degree prunes whole sub-trees.
[[nodiscard]] std::optional<IntPoly> recombine_from(
    const IntPoly& f,
    const std::vector<IntPoly>& lifted,
    const BigInt& modulus,
    std::size_t max_degree,
    std::size_t start,
    std::vector<bool>& selected,
    std::size_t selected_degree,
    const BigInt& mignotte_bound) {
    for (std::size_t index = start; index < lifted.size(); ++index) {
        const std::size_t next_degree = selected_degree + lifted[index].degree();
        if (next_degree > max_degree) {
            continue;
        }

        selected[index] = true;
        auto candidate = check_subset(
            f, lifted, selected, modulus, max_degree, mignotte_bound);
        if (candidate.has_value()) {
            return candidate;
        }

        auto nested = recombine_from(
            f, lifted, modulus, max_degree,
            index + 1U, selected, next_degree, mignotte_bound);
        if (nested.has_value()) {
            return nested;
        }
        selected[index] = false;
    }
    return std::nullopt;
}

} // namespace

std::optional<IntPoly> find_factor_by_hensel_recombination(
    const IntPoly& f,
    const std::vector<IntPoly>& lifted_factors,
    const BigInt& modulus,
    std::size_t max_degree) {
    if (f.is_zero() || modulus <= BigInt(1) || lifted_factors.empty()) {
        return std::nullopt;
    }

    // Every nonconstant lifted factor stays in the pool: a factor of degree
    // > max_degree is never *selected* (recombine_from prunes it via
    // next_degree > max_degree), but keeping it preserves the product
    // invariant that makes complements of found factors exact.
    std::vector<IntPoly> usable_factors;
    usable_factors.reserve(lifted_factors.size());
    bool any_selectable = false;
    for (IntPoly factor : lifted_factors) {
        normalize_integer_poly(factor);
        if (!factor.is_zero() && factor.degree() > 0U) {
            if (factor.degree() <= max_degree) any_selectable = true;
            usable_factors.push_back(std::move(factor));
        }
    }
    if (!any_selectable) {
        return std::nullopt;
    }

    std::vector<bool> selected(usable_factors.size(), false);
    const BigInt mignotte_bound = mignotte_factor_bound(f, max_degree);
    return recombine_from(
        f,
        usable_factors,
        modulus,
        max_degree,
        0U,
        selected,
        0U,
        mignotte_bound);
}

} // namespace cas::algebra
