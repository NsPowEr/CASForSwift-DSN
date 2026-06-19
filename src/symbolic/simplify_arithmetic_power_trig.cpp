// Chebyshev / DeMoivre trig power linearization — extracted from
// simplify_arithmetic_power.cpp (2026-06-19, anti-monolith split T-054 follow-up).
// Self-contained section flagged in the original file header; no behavioural
// change. Dependencies via simplify_impl.hpp.

#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

// Returns the linearized form of sin^n / cos^n using multiple-angle identities,
// or std::nullopt when the rule does not apply (caller falls through).
//
//   Even n=2m:   trig^n = (1/4^m)*[C(n,m) + 2*Σ_{j=0}^{m-1} s_j*C(n,j)*cos((n-2j)*arg)]
//   Odd  n=2m+1: trig^n = (1/4^m)*Σ_{j=0}^{m} s_j*C(n,j)*trig((n-2j)*arg)
//     sin: s_j = (-1)^(m-j),  cos: s_j = 1
//
// Limit: ctx.max_trig_power_reduction (default 32) — returns nullopt if exceeded.
std::optional<Result<ExprPtr>> Simplifier::try_linearize_trig_power(ExprPtr base, const BigInt& n) {
    if (!(n >= BigInt(2) && !n.is_negative())) return std::nullopt;

    const auto* func = expr_cast<FuncCall>(base);
    if (func == nullptr) return std::nullopt;
    const bool is_sin = (func->func_id == BuiltinOp::Sin);
    const bool is_cos = (func->func_id == BuiltinOp::Cos);
    if (!((is_sin || is_cos) && func->args.size() == 1U)) return std::nullopt;

    const long long max_n = context_ ? static_cast<long long>(context_->max_trig_power_reduction()) : 32LL;
    if (n > BigInt(max_n)) return std::nullopt;

    const long long n_ll = static_cast<long long>(n.to_u64());
    ExprPtr arg = func->args[0];
    // Pascal's triangle for C(n, j), j = 0..n
    std::vector<BigInt> binom(static_cast<std::size_t>(n_ll + 1), BigInt(0));
    binom[0] = BigInt(1);
    for (long long i = 1; i <= n_ll; ++i)
        for (long long j = i; j >= 1; --j)
            binom[static_cast<std::size_t>(j)] = binom[static_cast<std::size_t>(j)] + binom[static_cast<std::size_t>(j - 1)];

    const long long m = n_ll / 2;
    BigInt denom(1);
    for (long long i = 0; i < m; ++i) denom = denom * BigInt(4);

    std::vector<ExprPtr> terms;
    if (n_ll % 2 == 0) {
        // Constant: C(n,m)/4^m
        terms.push_back(make_rational(arena_, Rational(binom[static_cast<std::size_t>(m)], denom)));
        for (long long j = 0; j < m; ++j) {
            BigInt c2 = binom[static_cast<std::size_t>(j)] * BigInt(2);
            Rational coeff(c2, denom);
            if (is_sin && ((m - j) % 2 == 1)) coeff = -coeff;
            const long long k = n_ll - 2 * j;
            ExprPtr ka = arena_.make<Binary>(BinaryOp::Mul, make_integer(arena_, BigInt(k)), arg);
            ExprPtr cos_ka = arena_.make<FuncCall>(BuiltinOp::Cos, std::vector<ExprPtr>{ka});
            terms.push_back(arena_.make<Binary>(BinaryOp::Mul, make_rational(arena_, coeff), cos_ka));
        }
    } else {
        const BuiltinOp trig_op = is_sin ? BuiltinOp::Sin : BuiltinOp::Cos;
        for (long long j = 0; j <= m; ++j) {
            Rational coeff(binom[static_cast<std::size_t>(j)], denom);
            if (is_sin && ((m - j) % 2 == 1)) coeff = -coeff;
            const long long k = n_ll - 2 * j;
            ExprPtr ka = arena_.make<Binary>(BinaryOp::Mul, make_integer(arena_, BigInt(k)), arg);
            ExprPtr trig_ka = arena_.make<FuncCall>(trig_op, std::vector<ExprPtr>{ka});
            terms.push_back(arena_.make<Binary>(BinaryOp::Mul, make_rational(arena_, coeff), trig_ka));
        }
    }
    return simplify_expr(arena_.make<Sum>(std::move(terms)));
}

}  // namespace cas::symbolic::detail
