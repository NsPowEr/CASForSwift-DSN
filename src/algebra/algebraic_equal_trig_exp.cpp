// F7.5 (B.2) — exponential (Laurent-in-e^{ix}) trig-polynomial equivalence
// helper for `mathematically_equal`. Complements `weierstrass_zero_diff`.
//
// For a trig POLYNOMIAL (sums/products/non-negative powers of sin/cos of
// rational multiples of a single variable, no trig in any denominator) the
// substitution z = e^{ix} maps the expression to a Laurent polynomial
//   Σ_k c_k z^k     (k ∈ Q,  c_k ∈ Q(i))
// via the exact identities
//   cos(αx) = (z^α + z^{-α})/2,   sin(αx) = (z^α − z^{-α})/(2i) = −(i/2)z^α + (i/2)z^{-α}.
// Two trig polynomials are equal iff their Laurent coefficient maps coincide,
// so a difference is identically zero iff every coefficient is zero.
//
// Why this complements the Weierstrass path: t = tan(x/2) turns the difference
// into a *rational* function of t whose numerator must then be collapsed to
// zero — which can exceed the canonicaliser's budget on larger linearisable
// identities (e.g. sin³x·cos²x vs its Fourier expansion). The exponential form
// is a *direct* canonical form (coefficient collection only, no rational
// numerator), so it succeeds on the trig-polynomial cases where the t-form
// stalls.
//
// SOUND (REGOLA ZERO / A14): the conversion uses only exact, unconditional
// identities. Any sub-term that is not a single-variable trig polynomial — trig
// in a denominator (negative power / div by non-constant), a non-linear or
// offset or multi-variable argument, a non-rational coefficient — makes the
// whole conversion bail and the helper returns `false`. It therefore can only
// ever PROVE equality (empty residual map), never assert a false equality.

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

// Convert e to a Laurent polynomial in z = e^{ix}; nullopt if e is not a
// single-variable trig polynomial.
[[nodiscard]] std::optional<LaurentMap> to_laurent(ExprPtr e, const std::string& var) {
    const Rational zero_freq(BigInt(0));
    if (const auto* il = expr_cast<IntegerLit>(e)) {
        LaurentMap m; add_term(m, zero_freq, ComplexRational(Rational(il->value))); return m;
    }
    if (const auto* rl = expr_cast<RationalLit>(e)) {
        LaurentMap m; add_term(m, zero_freq, ComplexRational(Rational(rl->numerator, rl->denominator))); return m;
    }
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
            return m;
        }
        return std::nullopt;
    }
    if (const auto* u = expr_cast<Unary>(e); u != nullptr && u->op == UnaryOp::Neg) {
        auto inner = to_laurent(u->operand, var);
        if (!inner) return std::nullopt;
        LaurentMap m;
        for (const auto& [f, c] : *inner) add_term(m, f, -c);
        return m;
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        LaurentMap m;
        for (ExprPtr t : sum->terms) {
            auto tm = to_laurent(t, var);
            if (!tm) return std::nullopt;
            for (const auto& [f, c] : *tm) add_term(m, f, c);
        }
        return m;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        LaurentMap m; add_term(m, zero_freq, ComplexRational::one());
        for (ExprPtr f : prod->factors) {
            auto fm = to_laurent(f, var);
            if (!fm) return std::nullopt;
            m = convolve(m, *fm);
        }
        return m;
    }
    if (const auto* bin = expr_cast<Binary>(e)) {
        if (bin->op == BinaryOp::Mul) {
            auto l = to_laurent(bin->left, var); if (!l) return std::nullopt;
            auto r = to_laurent(bin->right, var); if (!r) return std::nullopt;
            return convolve(*l, *r);
        }
        if (bin->op == BinaryOp::Sub) {
            auto l = to_laurent(bin->left, var); if (!l) return std::nullopt;
            auto r = to_laurent(bin->right, var); if (!r) return std::nullopt;
            LaurentMap m = *l;
            for (const auto& [f, c] : *r) add_term(m, f, -c);
            return m;
        }
        if (bin->op == BinaryOp::Add) {
            auto l = to_laurent(bin->left, var); if (!l) return std::nullopt;
            auto r = to_laurent(bin->right, var); if (!r) return std::nullopt;
            LaurentMap m = *l;
            for (const auto& [f, c] : *r) add_term(m, f, c);
            return m;
        }
        if (bin->op == BinaryOp::Pow) {
            const auto* il = expr_cast<IntegerLit>(bin->right);
            if (il == nullptr || il->value.is_negative()) return std::nullopt;  // trig in denom → bail
            auto base = to_laurent(bin->left, var); if (!base) return std::nullopt;
            LaurentMap m; add_term(m, zero_freq, ComplexRational::one());
            const std::uint64_t n = il->value.to_u64();
            for (std::uint64_t i = 0; i < n; ++i) m = convolve(m, *base);
            return m;
        }
        if (bin->op == BinaryOp::Div) {
            // Only division by a pure non-zero constant (single freq-0 term) is
            // a trig polynomial; trig in the denominator bails.
            auto r = to_laurent(bin->right, var); if (!r) return std::nullopt;
            if (r->size() != 1U) return std::nullopt;
            auto it = r->begin();
            if (it->first != zero_freq) return std::nullopt;
            const ComplexRational denom = it->second;
            auto l = to_laurent(bin->left, var); if (!l) return std::nullopt;
            LaurentMap m;
            for (const auto& [f, c] : *l) {
                auto q = c.divide(denom);
                if (q.is_error()) return std::nullopt;
                add_term(m, f, q.value());
            }
            return m;
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

// True iff `diff_expr` is identically zero as a single-variable trig polynomial.
// Returns false (no claim) on anything outside that class — see file header.
[[nodiscard]] bool trig_exponential_zero_diff(ExprPtr diff_expr, CASContext& ctx) {
    (void)ctx;
    std::set<std::string> syms;
    collect_symbols(diff_expr, syms);
    if (syms.size() != 1U) return false;  // need exactly one free variable
    auto m = to_laurent(diff_expr, *syms.begin());
    if (!m) return false;
    return m->empty();  // all Laurent coefficients cancelled ⇒ diff ≡ 0
}

}  // namespace cas::symbolic
