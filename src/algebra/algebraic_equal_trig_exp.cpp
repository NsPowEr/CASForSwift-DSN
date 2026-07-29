// F7.5 (B.2) — exponential (Laurent-in-e^{ix}) trig-RATIONAL equivalence
// helper for `mathematically_equal`. Complements `weierstrass_zero_diff`.
//
// For a trig RATIONAL (sums/products/integer-powers — positive OR negative — of
// sin/cos of rational multiples of a single variable, including trig in the
// denominator) the substitution z = e^{ix} maps the expression to a *rational*
// Laurent function N(z)/D(z), where N, D are Laurent polynomials
//   Σ_k c_k z^k     (k ∈ Q,  c_k ∈ Q(i))
// built via the exact identities
//   cos(αx) = (z^α + z^{-α})/2,   sin(αx) = (z^α − z^{-α})/(2i) = −(i/2)z^α + (i/2)z^{-α}.
// A Laurent polynomial vanishes on an arc of the unit circle (z = e^{ix},
// x over any interval) only if every coefficient is zero. Hence
//   N(z)/D(z) ≡ 0  ⟺  N ≡ 0   (with D ≢ 0).
// So a trig-rational difference is identically zero iff its combined numerator
// (after clearing denominators by exact cross-multiplication of Laurent maps)
// has an empty coefficient map.
//
// Why this complements the Weierstrass path: t = tan(x/2) turns the difference
// into a *rational* function of t whose numerator must then be collapsed to
// zero — which can exceed the canonicaliser's budget on larger linearisable
// identities (e.g. sin³x·cos²x vs its Fourier expansion, or the derivative of a
// `ln|tan(x/2)|`-style antiderivative which carries negative half-angle powers).
// The exponential form is a *direct* canonical form (coefficient collection on
// a single cleared numerator), so it succeeds where the t-form stalls.
//
// SOUND (REGOLA ZERO / A14): the conversion uses only exact, unconditional
// identities. Any sub-term that is not a single-variable trig rational — a
// non-integer/symbolic power (e.g. sqrt(sin)), a non-linear or offset or
// multi-variable argument, a non-rational coefficient, a non-trig function
// (abs, ln, …) — makes the whole conversion bail and the helper returns
// `false`. A degenerate 0/0 (empty denominator) also bails. It can therefore
// only ever PROVE equality (empty numerator), never assert a false equality.

#include "cas/algebra.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/complex_rational.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"

#include <map>
#include <optional>
#include <set>
#include <string>

namespace cas::symbolic {

namespace {

using LaurentMap = std::map<Rational, ComplexRational>;

// A trig rational in z = e^{ix}: num / den, both Laurent polynomials.
struct LaurentRat {
    LaurentMap num;
    LaurentMap den;  // invariant: never empty (≡0 denom) by construction
};

void add_term(LaurentMap& m, const Rational& freq, const ComplexRational& c) {
    if (c.is_zero()) return;
    auto it = m.find(freq);
    if (it == m.end()) {
        m.emplace(freq, c);
        return;
    }
    ComplexRational s = it->second + c;
    if (s.is_zero()) m.erase(it);
    else it->second = s;
}

LaurentMap convolve(const LaurentMap& a, const LaurentMap& b) {
    LaurentMap r;
    for (const auto& [fa, ca] : a)
        for (const auto& [fb, cb] : b)
            add_term(r, fa + fb, ca * cb);
    return r;
}

LaurentMap one_map() {
    LaurentMap m;
    add_term(m, Rational(BigInt(0)), ComplexRational::one());
    return m;
}

LaurentRat lr_const(const ComplexRational& c) {
    LaurentMap n;
    add_term(n, Rational(BigInt(0)), c);
    return LaurentRat{std::move(n), one_map()};
}

LaurentRat lr_mul(const LaurentRat& a, const LaurentRat& b) {
    return LaurentRat{convolve(a.num, b.num), convolve(a.den, b.den)};
}

LaurentRat lr_div(const LaurentRat& a, const LaurentRat& b) {
    return LaurentRat{convolve(a.num, b.den), convolve(a.den, b.num)};
}

LaurentRat lr_neg(const LaurentRat& a) {
    LaurentMap n;
    for (const auto& [f, c] : a.num) add_term(n, f, -c);
    return LaurentRat{std::move(n), a.den};
}

// a + sign·b  (sign = +1 or −1), combined over a common denominator.
LaurentRat lr_add(const LaurentRat& a, const LaurentRat& b, bool subtract) {
    LaurentMap n = convolve(a.num, b.den);
    LaurentMap cross = convolve(b.num, a.den);
    for (const auto& [f, c] : cross) add_term(n, f, subtract ? -c : c);
    return LaurentRat{std::move(n), convolve(a.den, b.den)};
}

// a^n for integer exponent (negative ⇒ reciprocal). Repeated multiply.
LaurentRat lr_pow(const LaurentRat& a, bool negative, std::uint64_t e) {
    LaurentRat base = negative ? LaurentRat{a.den, a.num} : a;
    LaurentRat r = lr_const(ComplexRational::one());
    for (std::uint64_t i = 0; i < e; ++i) r = lr_mul(r, base);
    return r;
}

// arg must be exactly `c·var` (rational c, no constant offset, single variable).
[[nodiscard]] std::optional<Rational> linear_freq(ExprPtr arg, const std::string& var) {
    if (const auto* s = expr_cast<Symbol>(arg))
        return s->name == var ? std::optional<Rational>(Rational(BigInt(1))) : std::nullopt;
    auto as_rational = [](ExprPtr k) -> std::optional<Rational> {
        if (const auto* il = expr_cast<IntegerLit>(k)) return Rational(il->value);
        if (const auto* rl = expr_cast<RationalLit>(k)) return Rational(rl->numerator, rl->denominator);
        return std::nullopt;
    };
    auto try_pair = [&](ExprPtr k, ExprPtr v) -> std::optional<Rational> {
        const auto* sv = expr_cast<Symbol>(v);
        if (!sv || sv->name != var) return std::nullopt;
        return as_rational(k);
    };
    if (const auto* b = expr_cast<Binary>(arg); b != nullptr && b->op == BinaryOp::Mul) {
        if (auto c = try_pair(b->left, b->right)) return c;
        if (auto c = try_pair(b->right, b->left)) return c;
    }
    if (const auto* p = expr_cast<Product>(arg); p != nullptr && p->factors.size() == 2U) {
        if (auto c = try_pair(p->factors[0], p->factors[1])) return c;
        if (auto c = try_pair(p->factors[1], p->factors[0])) return c;
    }
    return std::nullopt;
}

// Convert e to a rational Laurent function N(z)/D(z) in z = e^{ix}; nullopt if
// e is not a single-variable trig rational (see file header for the bail set).
[[nodiscard]] std::optional<LaurentRat> to_laurent_rat(ExprPtr e, const std::string& var) {
    const Rational zero_freq(BigInt(0));
    if (const auto* il = expr_cast<IntegerLit>(e))
        return lr_const(ComplexRational(Rational(il->value)));
    if (const auto* rl = expr_cast<RationalLit>(e))
        return lr_const(ComplexRational(Rational(rl->numerator, rl->denominator)));
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        if ((fc->func_id == BuiltinOp::Sin || fc->func_id == BuiltinOp::Cos) && fc->args.size() == 1U) {
            auto freq = linear_freq(fc->args[0], var);
            if (!freq) return std::nullopt;
            LaurentMap m;
            const Rational half(BigInt(1), BigInt(2));
            if (fc->func_id == BuiltinOp::Cos) {
                add_term(m, *freq, ComplexRational(half));
                add_term(m, -*freq, ComplexRational(half));
            } else {
                add_term(m, *freq, ComplexRational(zero_freq, -half));
                add_term(m, -*freq, ComplexRational(zero_freq, half));
            }
            return LaurentRat{std::move(m), one_map()};
        }
        return std::nullopt;  // non-trig function (abs, ln, arctan, …) → bail
    }
    if (const auto* u = expr_cast<Unary>(e); u != nullptr && u->op == UnaryOp::Neg) {
        auto inner = to_laurent_rat(u->operand, var);
        if (!inner) return std::nullopt;
        return lr_neg(*inner);
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        LaurentRat acc = lr_const(ComplexRational(zero_freq));
        for (ExprPtr t : sum->terms) {
            auto tm = to_laurent_rat(t, var);
            if (!tm) return std::nullopt;
            acc = lr_add(acc, *tm, /*subtract=*/false);
        }
        return acc;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        LaurentRat acc = lr_const(ComplexRational::one());
        for (ExprPtr f : prod->factors) {
            auto fm = to_laurent_rat(f, var);
            if (!fm) return std::nullopt;
            acc = lr_mul(acc, *fm);
        }
        return acc;
    }
    if (const auto* bin = expr_cast<Binary>(e)) {
        if (bin->op == BinaryOp::Mul) {
            auto l = to_laurent_rat(bin->left, var); if (!l) return std::nullopt;
            auto r = to_laurent_rat(bin->right, var); if (!r) return std::nullopt;
            return lr_mul(*l, *r);
        }
        if (bin->op == BinaryOp::Div) {
            auto l = to_laurent_rat(bin->left, var); if (!l) return std::nullopt;
            auto r = to_laurent_rat(bin->right, var); if (!r) return std::nullopt;
            if (r->num.empty()) return std::nullopt;  // division by ≡0 → bail
            return lr_div(*l, *r);
        }
        if (bin->op == BinaryOp::Add || bin->op == BinaryOp::Sub) {
            auto l = to_laurent_rat(bin->left, var); if (!l) return std::nullopt;
            auto r = to_laurent_rat(bin->right, var); if (!r) return std::nullopt;
            return lr_add(*l, *r, /*subtract=*/bin->op == BinaryOp::Sub);
        }
        if (bin->op == BinaryOp::Pow) {
            const auto* il = expr_cast<IntegerLit>(bin->right);
            if (il == nullptr) return std::nullopt;  // non-integer power (sqrt, symbolic) → bail
            auto base = to_laurent_rat(bin->left, var); if (!base) return std::nullopt;
            const bool neg = il->value.is_negative();
            if (neg && base->num.empty()) return std::nullopt;  // 0^{-n} → bail
            return lr_pow(*base, neg, il->value.abs().to_u64());
        }
        return std::nullopt;
    }
    return std::nullopt;
}

void collect_symbols(ExprPtr e, std::set<std::string>& out) {
    if (!e) return;
    if (const auto* s = expr_cast<Symbol>(e)) { out.insert(s->name); return; }
    if (const auto* u = expr_cast<Unary>(e)) { collect_symbols(u->operand, out); return; }
    if (const auto* b = expr_cast<Binary>(e)) { collect_symbols(b->left, out); collect_symbols(b->right, out); return; }
    if (const auto* sum = expr_cast<Sum>(e)) { for (ExprPtr t : sum->terms) collect_symbols(t, out); return; }
    if (const auto* p = expr_cast<Product>(e)) { for (ExprPtr f : p->factors) collect_symbols(f, out); return; }
    if (const auto* fc = expr_cast<FuncCall>(e)) { for (ExprPtr a : fc->args) collect_symbols(a, out); return; }
}

}  // namespace

// True iff `diff_expr` is identically zero as a single-variable trig rational.
// Returns false (no claim) on anything outside that class — see file header.
[[nodiscard]] bool trig_exponential_zero_diff(ExprPtr diff_expr, CASContext& ctx) {
    (void)ctx;
    std::set<std::string> syms;
    collect_symbols(diff_expr, syms);
    if (syms.size() != 1U) return false;  // need exactly one free variable
    auto r = to_laurent_rat(diff_expr, *syms.begin());
    if (!r) return false;
    if (r->den.empty()) return false;  // degenerate 0/0 → no claim
    return r->num.empty();  // numerator cancelled ⇒ N/D ≡ 0 (D ≢ 0) ⇒ diff ≡ 0
}

}  // namespace cas::symbolic
