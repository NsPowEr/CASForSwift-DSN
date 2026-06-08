// F7.2-T1 — Univariate Normal distribution N(μ, σ²).
//
// Closed forms:
//   pdf(x; μ, σ) = exp(-(x - μ)² / (2σ²)) / (σ · sqrt(2π))
//   cdf(x; μ, σ) = ½ · (1 + erf((x - μ) / (σ · sqrt(2))))
//
// The pdf is returned as a symbolic ExprPtr so callers can compose it
// inside larger expressions (e.g. likelihoods).  The cdf reuses
// BuiltinOp::Erf which the symbolic engine already knows how to
// differentiate, simplify, and numerically evaluate.
//
// The numeric quantile (inverse CDF) is computed via Newton-Raphson on
// the residual cdf(x) - p, with the well-conditioned Wichura starting
// guess.  Convergence is quadratic; we cap iterations at 30 to fail
// gracefully on degenerate inputs (σ ≤ 0 or p ∉ (0, 1)).

#include "cas/ast.hpp"
#include "cas/statistics.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>

namespace cas::statistics {

namespace {

[[nodiscard]] ExprPtr make_int(AstArena& arena, long long n) {
    return arena.make<IntegerLit>(BigInt(n));
}

[[nodiscard]] ExprPtr sqrt_two_pi(AstArena& arena) {
    return arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{
        arena.make<Product>(std::vector<ExprPtr>{
            make_int(arena, 2),
            arena.make<Constant>(MathConstant::Pi)})});
}

[[nodiscard]] ExprPtr sqrt_two(AstArena& arena) {
    return arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{make_int(arena, 2)});
}

}  // namespace

ExprPtr normal_pdf(ExprPtr x, ExprPtr mu, ExprPtr sigma, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr dx = arena.make<Binary>(BinaryOp::Sub, x, mu);
    ExprPtr dx2 = arena.make<Binary>(BinaryOp::Pow, dx, make_int(arena, 2));
    ExprPtr two_sigma2 = arena.make<Product>(std::vector<ExprPtr>{
        make_int(arena, 2),
        arena.make<Binary>(BinaryOp::Pow, sigma, make_int(arena, 2))});
    ExprPtr exponent = arena.make<Unary>(UnaryOp::Neg,
        arena.make<Binary>(BinaryOp::Div, dx2, two_sigma2));
    ExprPtr numerator = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{exponent});
    ExprPtr denominator = arena.make<Product>(std::vector<ExprPtr>{sigma, sqrt_two_pi(arena)});
    return arena.make<Binary>(BinaryOp::Div, numerator, denominator);
}

ExprPtr normal_cdf(ExprPtr x, ExprPtr mu, ExprPtr sigma, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr dx = arena.make<Binary>(BinaryOp::Sub, x, mu);
    ExprPtr scaled = arena.make<Binary>(BinaryOp::Div, dx,
        arena.make<Product>(std::vector<ExprPtr>{sigma, sqrt_two(arena)}));
    ExprPtr erf_call = arena.make<FuncCall>(BuiltinOp::Erf, std::vector<ExprPtr>{scaled});
    ExprPtr sum = arena.make<Sum>(std::vector<ExprPtr>{make_int(arena, 1), erf_call});
    return arena.make<Binary>(BinaryOp::Div, sum, make_int(arena, 2));
}

Result<double> normal_quantile(double p, double mu, double sigma) {
    if (!(p > 0.0) || !(p < 1.0)) {
        CASError err{};
        err.kind = CASErrorKind::InvalidArgument;
        err.message = "normal_quantile: p must lie strictly in (0, 1)";
        return fail<double>(err);
    }
    if (!(sigma > 0.0)) {
        CASError err{};
        err.kind = CASErrorKind::InvalidArgument;
        err.message = "normal_quantile: sigma must be positive";
        return fail<double>(err);
    }
    // Wichura's "Algorithm AS 241" initial approximation (Applied
    // Statistics 1988) using a single rational expansion around the tail.
    auto initial_guess = [&]() {
        const double q = p - 0.5;
        if (std::abs(q) <= 0.425) {
            const double r = q * q;
            return q * (((-25.44106049637 * r + 41.39119773534) * r - 18.61500062529) * r + 2.50662823884)
                     / ((((3.13082909833 * r - 21.06224101826) * r + 23.08336743743) * r - 8.47351093090) * r + 1.0);
        }
        const double r = (q < 0.0) ? std::log(-std::log(p))
                                   : std::log(-std::log(1.0 - p));
        const double poly = (((2.93816e-2 * r + 0.05116) * r + 0.5318) * r + 0.6873) * r - 1.0;
        return std::copysign(poly, q);
    };
    double z = initial_guess();
    // Newton on residual r(z) = Φ(z) - p, using r'(z) = φ(z).
    const double inv_sqrt_two = 1.0 / std::sqrt(2.0);
    const double inv_sqrt_two_pi = 1.0 / std::sqrt(2.0 * std::numbers::pi);
    for (int i = 0; i < 30; ++i) {
        const double phi = 0.5 * (1.0 + std::erf(z * inv_sqrt_two));
        const double pdf = inv_sqrt_two_pi * std::exp(-0.5 * z * z);
        if (pdf == 0.0) break;
        const double next = z - (phi - p) / pdf;
        if (std::abs(next - z) < 1e-12) {
            z = next;
            break;
        }
        z = next;
    }
    return ok(mu + sigma * z);
}

}  // namespace cas::statistics
