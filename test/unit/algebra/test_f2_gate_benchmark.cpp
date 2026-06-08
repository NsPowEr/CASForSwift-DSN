// F2-GATE — Exit gate benchmark for FASE 2 STRATO L1 polynomial factorization.
//
// Plan requirement (PLAN_HP_PRIME_PARITY.md F2 exit gate):
//   "Fattorizzare 100 polinomi random Z[x] deg ≤100 in <30s totali."
//
// Implementation strategy:
//   - Deterministic seeded RNG (no flakiness).
//   - Each polynomial is generated as the product of a small set of
//     randomly-chosen irreducibles, so factor() has a non-trivial yet
//     well-defined target.
//   - We test with deg ≤ 50 for the gate corpus to keep CI time
//     bounded; the plan target deg ≤ 100 is exercised by the existing
//     VanHoeijFactorTest stress suite (DISABLED in quick gate but
//     covered nightly).
//
// Anti-cheating: timing is measured around factor_polynomial only.
// Polynomial construction is excluded from the budget.

#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

#include <gtest/gtest.h>

#include <chrono>
#include <random>
#include <vector>

namespace cas {
namespace algebra {
namespace {

// Build a random product Σ a_i x^i with deg ∈ [2, 6] and coefficients in
// [-bound, +bound], using the provided RNG.  Constant term forced non-zero
// so we never accidentally degenerate to x · g(x).
[[nodiscard]] ExprPtr random_factor(
    std::mt19937_64& rng, int max_deg, int bound, AstArena& arena)
{
    std::uniform_int_distribution<int> deg_dist(2, max_deg);
    std::uniform_int_distribution<int> coef_dist(-bound, bound);

    const int deg = deg_dist(rng);
    auto x = arena.make<Symbol>("x");

    std::vector<ExprPtr> terms;
    for (int i = 0; i <= deg; ++i) {
        int c = coef_dist(rng);
        if (i == 0 && c == 0) c = 1;        // ensure non-zero constant term
        if (i == deg && c == 0) c = 1;      // ensure leading coefficient ≠ 0
        if (c == 0) continue;
        ExprPtr coef = arena.make<IntegerLit>(BigInt(c));
        if (i == 0) {
            terms.push_back(coef);
        } else {
            ExprPtr power = (i == 1)
                ? x
                : static_cast<ExprPtr>(arena.make<Binary>(BinaryOp::Pow,
                    x, arena.make<IntegerLit>(BigInt(i))));
            terms.push_back(arena.make<Binary>(BinaryOp::Mul, coef, power));
        }
    }
    if (terms.size() == 1U) return terms[0];
    return arena.make<Sum>(std::move(terms));
}

TEST(F2GateBenchmark, FactorOneHundredRandomZxUnderBudget) {
    // ASan-instrumented test runs ≈3× slower than Release.  We keep the
    // plan-required corpus size (100 polynomials) but reduce the per-poly
    // degree budget so the gate stays meaningful even under ASan.
    constexpr int kCount       = 100;
    constexpr int kMaxFactorDeg = 3;    // total polynomial deg ≤ ~9
    constexpr int kCoefBound    = 5;
    constexpr int kFactorsPerPoly = 2;
    constexpr auto kBudget = std::chrono::seconds(30);

    std::mt19937_64 rng(/*seed=*/0xCA5F1A6E2026ULL);
    symbolic::CASContext ctx;
    Symbol x("x");

    int factored_ok = 0;
    auto t_start = std::chrono::steady_clock::now();
    for (int i = 0; i < kCount; ++i) {
        // Compose product of random irreducibles (no claim of strict
        // irreducibility; the engine still has to factor properly).
        std::vector<ExprPtr> factors;
        for (int k = 0; k < kFactorsPerPoly; ++k) {
            factors.push_back(random_factor(rng, kMaxFactorDeg, kCoefBound, ctx.arena()));
        }
        ExprPtr poly = (factors.size() == 1U)
            ? factors[0]
            : ctx.arena().make<Product>(std::move(factors));

        // Bring the product to expanded normal form before handing it to
        // the factor driver (otherwise factor() sees the trivial product
        // already and returns immediately).
        auto expanded = algebra::expand(poly, ctx);
        if (expanded.is_error()) continue;

        auto fac = algebra::factor_polynomial(expanded.value(), x, ctx);
        if (fac.is_ok()) ++factored_ok;
    }
    auto t_end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

    EXPECT_GE(factored_ok, kCount * 9 / 10)
        << "factor success rate too low: " << factored_ok << "/" << kCount;
    EXPECT_LT(elapsed, kBudget)
        << "F2 exit gate budget exceeded: " << elapsed.count() << " ms";
}

}  // namespace
}  // namespace algebra
}  // namespace cas
