// W9.3 split: rational-sqrt / denesting helpers extracted from the former
// monolithic simplify_exp_log.cpp (854 LOC). Definitions only; declarations and
// rationale live in simplify_sqrt_internal.hpp. Consumed by simplify_sqrt.cpp.

#include "simplify_sqrt_internal.hpp"

namespace cas::symbolic::detail {

// ── sqrt helpers ─────────────────────────────────────────────────────────────

// Extract perfect-square factor from n: returns {k, m} with n = k²·m.
// Phase 1: trial-divide squares i² for i ∈ [2, trial_bound]; this is the
// "small squarefull" extraction step.
// Phase 2: check whether the residue is itself a perfect square via
// integer_sqrt (Newton). If so, fold it entirely into k.
//
// Without a bound, the unbounded O(sqrt(n)) loop hung on big rational
// radicands (norm² for QR Householder on 8×8 random Q rationals).
// Reference: HC-F4-QR-SYMBOLIC-TIMEOUT.
[[nodiscard]] std::pair<BigInt, BigInt> extract_square_factor(BigInt n, std::size_t trial_bound) {
    // Regola 1: no int64_t/double arithmetic in symbolic core; loop counter is
    // BigInt.  trial_bound is a CASContext-configurable budget (boundary
    // conversion to BigInt is permitted; arithmetic on the value is BigInt).
    BigInt k(1);
    BigInt bound(static_cast<std::int64_t>(trial_bound));
    BigInt i(2);
    while (i <= bound) {
        BigInt i2 = i * i;
        if (i2 > n) break;
        while ((n % i2).is_zero()) {
            k = k * i;
            n = n / i2;
        }
        i = i + BigInt(1);
    }
    // Fallback: if residue is itself a perfect square, absorb it.
    BigInt s = integer_sqrt(n);
    if (!(s * s == n)) return {k, n};
    k = k * s;
    return {k, BigInt(1)};
}

// sqrt(r) for rational r ≥ 0: extract perfect-square factors.
// Returns k·sqrt(m), k rational, m squarefree.
[[nodiscard]] Result<ExprPtr> simplify_rational_sqrt(const Rational& r, AstArena& arena, std::size_t trial_bound) {
    const BigInt& p = r.numerator();
    const BigInt& q = r.denominator();
    if (p.is_zero()) return ok(arena.make<IntegerLit>(BigInt(0)));
    auto [p_out, p_rem] = extract_square_factor(p, trial_bound);
    auto [q_out, q_rem] = extract_square_factor(q, trial_bound);
    BigInt final_radicand = p_rem * q_rem;
    Rational coeff(p_out, q_out * q_rem);

    ExprPtr coeff_expr;
    if (coeff.denominator() == BigInt(1)) {
        coeff_expr = arena.make<IntegerLit>(coeff.numerator());
    } else {
        coeff_expr = arena.make<RationalLit>(coeff.numerator(), coeff.denominator());
    }

    if (final_radicand == BigInt(1)) return ok(coeff_expr);

    ExprPtr sqrt_expr = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{arena.make<IntegerLit>(final_radicand)});
    if (coeff.numerator() == BigInt(1) && coeff.denominator() == BigInt(1))
        return ok(sqrt_expr);
    return ok(arena.make<Binary>(BinaryOp::Mul, coeff_expr, sqrt_expr));
}

// Try to extract a rational sqrt: if r = (p/q)² returns p/q, else nullopt.
[[nodiscard]] std::optional<Rational> try_rational_sqrt(const Rational& r) {
    if (r.numerator().is_negative()) return std::nullopt;
    if (r.numerator().is_zero()) return Rational(BigInt(0));
    auto isqrt = [](const BigInt& n) -> std::optional<BigInt> {
        if (n.is_zero()) return BigInt(0);
        // Newton-Raphson integer sqrt
        BigInt x = n;
        BigInt y = (x + BigInt(1)) / BigInt(2);
        while (y < x) {
            x = y;
            y = (x + n / x) / BigInt(2);
        }
        if (x * x == n) return x;
        return std::nullopt;
    };
    auto num_sqrt = isqrt(r.numerator());
    auto den_sqrt = isqrt(r.denominator());
    if (!num_sqrt || !den_sqrt) return std::nullopt;
    return Rational(*num_sqrt, *den_sqrt);
}

// Borodin-Fagin-Hopcroft-Tompa (1985) denesting:
//   sqrt(a + b·sqrt(c)) = sqrt(p) + sign(b)·sqrt(q)
// iff a² - b²·c is a rational square. Then p = (a+d)/2, q = (a-d)/2
// with d = sqrt(a² - b²c).
//
// Detects argument shapes:
//   Sum([rat_a, Product([rat_b, sqrt(rat_c)])]) — generic a + b·sqrt(c)
//   Sum([rat_a, sqrt(rat_c)])                   — b = 1
//   Sum([rat_a, Unary(Neg, Product([rat_b, sqrt(rat_c)]))]) — negative b
//   Sum([rat_a, Unary(Neg, sqrt(rat_c))])       — b = -1
//
// Returns the denested form on match; nullopt on no match or non-denestable.
[[nodiscard]] std::optional<ExprPtr> try_denest_borodin_fagin(
    ExprPtr radicand, AstArena& arena) {
    const auto* sum = expr_cast<Sum>(radicand);
    if (!sum || sum->terms.size() != 2) return std::nullopt;

    auto extract_rational = [](ExprPtr e) -> std::optional<Rational> {
        if (auto* il = expr_cast<IntegerLit>(e))
            return Rational(il->value, BigInt(1));
        if (auto* rl = expr_cast<RationalLit>(e))
            return Rational(rl->numerator, rl->denominator);
        return std::nullopt;
    };

    // Extract b·sqrt(c) from a term: returns {b, c} or nullopt.
    auto extract_b_sqrt_c =
        [&](ExprPtr term) -> std::optional<std::pair<Rational, Rational>> {
        bool negate = false;
        if (auto* un = expr_cast<Unary>(term); un && un->op == UnaryOp::Neg) {
            negate = true;
            term = un->operand;
        }
        Rational b(BigInt(1), BigInt(1));
        ExprPtr sqrt_node = nullptr;
        if (auto* call = expr_cast<FuncCall>(term);
            call && call->func_id == BuiltinOp::Sqrt && call->args.size() == 1) {
            sqrt_node = term;
        } else if (auto* prod = expr_cast<Product>(term)) {
            std::vector<ExprPtr> coeff_factors;
            for (ExprPtr f : prod->factors) {
                if (auto* call = expr_cast<FuncCall>(f);
                    call && call->func_id == BuiltinOp::Sqrt && call->args.size() == 1) {
                    if (sqrt_node) return std::nullopt;
                    sqrt_node = f;
                } else {
                    coeff_factors.push_back(f);
                }
            }
            if (!sqrt_node) return std::nullopt;
            if (coeff_factors.size() == 1) {
                auto r = extract_rational(coeff_factors[0]);
                if (!r) return std::nullopt;
                b = *r;
            } else if (!coeff_factors.empty()) {
                return std::nullopt;
            }
        } else if (auto* bin = expr_cast<Binary>(term);
                   bin && bin->op == BinaryOp::Mul) {
            // Binary Mul: rat * sqrt(c)
            ExprPtr lhs = bin->left, rhs = bin->right;
            auto rl = extract_rational(lhs);
            auto rr = extract_rational(rhs);
            if (rl && !rr) {
                b = *rl;
                sqrt_node = rhs;
            } else if (rr && !rl) {
                b = *rr;
                sqrt_node = lhs;
            } else {
                return std::nullopt;
            }
            auto* call = expr_cast<FuncCall>(sqrt_node);
            if (!call || call->func_id != BuiltinOp::Sqrt || call->args.size() != 1)
                return std::nullopt;
        } else {
            return std::nullopt;
        }
        auto* call = expr_cast<FuncCall>(sqrt_node);
        if (!call) return std::nullopt;
        auto c = extract_rational(call->args[0]);
        if (!c) return std::nullopt;
        if (c->numerator().is_negative()) return std::nullopt;  // complex outside scope
        if (negate) b = -b;
        return std::make_pair(b, *c);
    };

    // Try both orderings: term[0]=a, term[1]=b·sqrt(c); and swapped.
    for (int swap = 0; swap < 2; ++swap) {
        ExprPtr a_term = sum->terms[swap];
        ExprPtr bsc_term = sum->terms[1 - swap];
        auto a = extract_rational(a_term);
        if (!a) continue;
        auto bsc = extract_b_sqrt_c(bsc_term);
        if (!bsc) continue;
        auto [b, c] = *bsc;
        // Discriminant: d² = a² - b²·c
        Rational disc_sq = (*a) * (*a) - b * b * c;
        if (disc_sq.numerator().is_negative()) continue;
        auto d = try_rational_sqrt(disc_sq);
        if (!d) continue;
        Rational p = ((*a) + *d) / Rational(BigInt(2), BigInt(1));
        Rational q = ((*a) - *d) / Rational(BigInt(2), BigInt(1));
        if (p.numerator().is_negative()) continue;
        if (q.numerator().is_negative()) continue;
        // Build sqrt(p) and sqrt(q) (could simplify if perfect square)
        ExprPtr sqrt_p, sqrt_q;
        if (auto pr = try_rational_sqrt(p)) {
            sqrt_p = arena.make<RationalLit>(pr->numerator(), pr->denominator());
        } else {
            sqrt_p = arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{arena.make<RationalLit>(p.numerator(), p.denominator())});
        }
        if (auto qr = try_rational_sqrt(q)) {
            sqrt_q = arena.make<RationalLit>(qr->numerator(), qr->denominator());
        } else {
            sqrt_q = arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{arena.make<RationalLit>(q.numerator(), q.denominator())});
        }
        // sign(b) determines + or - on sqrt(q)
        if (b.numerator().is_negative()) {
            return arena.make<Binary>(BinaryOp::Sub, sqrt_p, sqrt_q);
        }
        return arena.make<Binary>(BinaryOp::Add, sqrt_p, sqrt_q);
    }
    return std::nullopt;
}

[[nodiscard]] BigInt integer_sqrt(const BigInt& n) {
    if (n.is_zero()) return BigInt(0);
    static const BigInt one(1);
    if (n == one) return one;
    BigInt x = one.shift_left_bits((n.bit_length() + 1) / 2);
    while (true) {
        BigInt y = (x + n / x) / BigInt(2);
        if (y >= x) return x;
        x = std::move(y);
    }
}

}  // namespace cas::symbolic::detail
