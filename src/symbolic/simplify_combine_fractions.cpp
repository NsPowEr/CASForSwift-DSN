// T-054b: combine sum terms over a common denominator and cancel.
//
// Non-recursive by design: detection and cancellation run on local univariate
// Rational-polynomial arithmetic, never calling simplify_expr. (An earlier
// speculative-simplify version was correct but re-simplified per candidate on the
// default path and timed out diff / hung Taylor — see HANDOFF / T-054b notes.)
//
// Scope (conservative — never produces a wrong or larger form; on anything outside
// it the sum is left untouched, which is always a valid simplification):
//   • single variable v;
//   • every term, after removing the common base g, is a polynomial in v;
//   • g is either a univariate polynomial Q(v) (any denominator power), or a
//     surd √(P(v)) with the denominator power |m| = 1 (the form integration of
//     x^k/√(…) and ∫asin/atan produce).
// Commits only when the denominator fully cancels, so non-cancelling fraction
// sums (e.g. 1/(x+1)+1/(x+2)) are deliberately left alone.

#include "simplify_arithmetic_chain_impl.hpp"

#include <algorithm>
#include <optional>

namespace cas::symbolic::detail {
namespace {

using Poly = std::vector<Rational>;  // ascending coefficients; never empty

void poly_trim(Poly& p) {
    while (p.size() > 1U && p.back().numerator().is_zero()) p.pop_back();
    if (p.empty()) p.push_back(Rational(BigInt(0)));
}
bool poly_is_zero(const Poly& p) {
    for (const auto& c : p) if (!c.numerator().is_zero()) return false;
    return true;
}
Poly poly_add(const Poly& a, const Poly& b) {
    Poly r(std::max(a.size(), b.size()), Rational(BigInt(0)));
    for (std::size_t i = 0; i < a.size(); ++i) r[i] = r[i] + a[i];
    for (std::size_t i = 0; i < b.size(); ++i) r[i] = r[i] + b[i];
    poly_trim(r);
    return r;
}
Poly poly_mul(const Poly& a, const Poly& b) {
    if (poly_is_zero(a) || poly_is_zero(b)) return {Rational(BigInt(0))};
    Poly r(a.size() + b.size() - 1U, Rational(BigInt(0)));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            r[i + j] = r[i + j] + a[i] * b[j];
    poly_trim(r);
    return r;
}
Poly poly_pow(const Poly& b, std::size_t n) {
    Poly r{Rational(BigInt(1))};
    for (std::size_t i = 0; i < n; ++i) r = poly_mul(r, b);
    return r;
}
// a = q·b + rem with deg(rem) < deg(b); b ≠ 0. Returns false if b is zero.
bool poly_divmod(const Poly& a, const Poly& b, Poly& q, Poly& rem) {
    Poly d = b; poly_trim(d);
    if (poly_is_zero(d)) return false;
    Poly r = a; poly_trim(r);
    const std::size_t db = d.size() - 1U;
    const Rational lead = d.back();
    q.assign((r.size() >= d.size()) ? (r.size() - db) : 1U, Rational(BigInt(0)));
    while (!poly_is_zero(r) && r.size() >= d.size()) {
        const std::size_t shift = (r.size() - 1U) - db;
        const Rational f = r.back() / lead;
        q[shift] = f;
        for (std::size_t i = 0; i < d.size(); ++i)
            r[shift + i] = r[shift + i] - f * d[i];
        poly_trim(r);
    }
    poly_trim(q);
    rem = r;
    return true;
}

// First Symbol encountered in a depth-first walk, or nullptr.
ExprPtr first_symbol(ExprPtr e) {
    if (e == nullptr) return nullptr;
    if (expr_is<Symbol>(e)) return e;
    if (const auto* u = expr_cast<Unary>(e)) return first_symbol(u->operand);
    if (const auto* b = expr_cast<Binary>(e)) {
        if (ExprPtr s = first_symbol(b->left)) return s;
        return first_symbol(b->right);
    }
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms) if (ExprPtr r = first_symbol(t)) return r;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr f : p->factors) if (ExprPtr r = first_symbol(f)) return r;
    }
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        for (ExprPtr a : fc->args) if (ExprPtr r = first_symbol(a)) return r;
    }
    return nullptr;
}

// expr → univariate polynomial in v (Rational coeffs), or nullopt if it is not a
// polynomial in v (other symbol, negative/fractional power, surd, transcendental).
// Fully structural so coeff·(a+b), products of polynomials, and powers all reduce.
std::optional<Poly> expr_to_poly(ExprPtr e, ExprPtr v) {
    {
        LiteralRational lr;
        auto r = try_get_exact_rational(e, lr);
        if (r.is_error()) return std::nullopt;
        if (r.value()) return Poly{lr.value};
    }
    if (structural_equal(e, v)) return Poly{Rational(BigInt(0)), Rational(BigInt(1))};
    if (const auto* u = expr_cast<Unary>(e); u && u->op == UnaryOp::Neg) {
        auto p = expr_to_poly(u->operand, v);
        if (!p) return std::nullopt;
        for (auto& c : *p) c = c * Rational(BigInt(-1));
        return p;
    }
    if (const auto* s = expr_cast<Sum>(e)) {
        Poly acc{Rational(BigInt(0))};
        for (ExprPtr t : s->terms) {
            auto p = expr_to_poly(t, v);
            if (!p) return std::nullopt;
            acc = poly_add(acc, *p);
        }
        return acc;
    }
    if (const auto* pr = expr_cast<Product>(e)) {
        Poly acc{Rational(BigInt(1))};
        for (ExprPtr f : pr->factors) {
            auto p = expr_to_poly(f, v);
            if (!p) return std::nullopt;
            acc = poly_mul(acc, *p);
        }
        return acc;
    }
    if (const auto* b = expr_cast<Binary>(e); b && b->op == BinaryOp::Pow) {
        auto n = try_get_integer_exponent(b->right);
        if (!n.has_value() || n->is_negative()) return std::nullopt;
        auto base = expr_to_poly(b->left, v);
        if (!base) return std::nullopt;
        return poly_pow(*base, static_cast<std::size_t>(n->to_u64()));
    }
    return std::nullopt;
}

}  // namespace

ExprPtr Simplifier::poly_to_expr_for_combine(const std::vector<Rational>& p, ExprPtr v) {
    std::vector<ExprPtr> terms;
    for (std::size_t k = 0; k < p.size(); ++k) {
        if (p[k].numerator().is_zero()) continue;
        ExprPtr coeff = make_rational(arena_, p[k]);
        if (k == 0) { terms.push_back(coeff); continue; }
        ExprPtr mono = (k == 1)
            ? v
            : arena_.make<Binary>(BinaryOp::Pow, v,
                  make_integer(arena_, BigInt(static_cast<long long>(k))));
        if (p[k] == Rational(BigInt(1))) terms.push_back(mono);
        else terms.push_back(arena_.make<Product>(std::vector<ExprPtr>{coeff, mono}));
    }
    if (terms.empty()) return make_integer(arena_, BigInt(0));
    if (terms.size() == 1U) return terms[0];
    return arena_.make<Sum>(std::move(terms));
}

bool Simplifier::try_combine_common_denominator(std::vector<ExprPtr>& terms) {
    if (terms.size() < 2U) return false;

    // Decompose every term once into (coeff, factors). Bail on any term we cannot
    // read as coeff·∏ base^exp (decompose_term never fails on valid AST).
    struct Dec { Rational coeff; std::vector<std::pair<ExprPtr, BigInt>> factors; };
    std::vector<Dec> dec(terms.size());
    for (std::size_t i = 0; i < terms.size(); ++i)
        if (!decompose_term(terms[i], dec[i].coeff, dec[i].factors)) return false;

    // Candidate denominator bases: any base appearing with a negative exponent.
    std::vector<ExprPtr> candidates;
    for (const auto& d : dec)
        for (const auto& [base, exp] : d.factors)
            if (exp.is_negative()) {
                bool seen = false;
                for (ExprPtr c : candidates) if (structural_equal(c, base)) { seen = true; break; }
                if (!seen) candidates.push_back(base);
            }

    for (ExprPtr g : candidates) {
        const bool is_surd = [&] {
            const auto* fc = expr_cast<FuncCall>(g);
            return fc != nullptr && fc->func_id == BuiltinOp::Sqrt && fc->args.size() == 1U;
        }();
        ExprPtr radicand = is_surd ? expr_cast<FuncCall>(g)->args[0] : g;

        ExprPtr v = first_symbol(radicand);
        if (v == nullptr) continue;
        auto base_poly = expr_to_poly(is_surd ? radicand : g, v);   // P or Q
        if (!base_poly.has_value()) continue;

        // Per term: net exponent e of g, and the polynomial rest = term / g^e
        // (rebuilt as an expression and converted — rest may be any v-polynomial,
        // e.g. (3/16·x²+3/32) for ∫x³·asin). nullopt if rest is not a v-polynomial.
        auto term_rest = [&](const Dec& d, BigInt& e_out) -> std::optional<Poly> {
            e_out = BigInt(0);
            std::vector<ExprPtr> factors{make_rational(arena_, d.coeff)};
            for (const auto& [base, exp] : d.factors) {
                if (structural_equal(base, g)) { e_out += exp; continue; }
                if (exp.is_negative()) return std::nullopt;   // a second denominator
                factors.push_back(exp == BigInt(1)
                    ? base
                    : arena_.make<Binary>(BinaryOp::Pow, base, make_integer(arena_, exp)));
            }
            ExprPtr rest_expr = factors.size() == 1U
                ? factors[0]
                : arena_.make<Product>(std::move(factors));
            return expr_to_poly(rest_expr, v);
        };

        std::vector<std::pair<BigInt, Poly>> parts;
        parts.reserve(dec.size());
        BigInt m(0);
        bool shape_ok = true;
        for (const auto& d : dec) {
            BigInt e;
            auto rp = term_rest(d, e);
            if (!rp) { shape_ok = false; break; }
            if (e < m) m = e;
            parts.emplace_back(e, std::move(*rp));
        }
        if (!shape_ok || !m.is_negative()) continue;

        if (!is_surd) {
            // Polynomial denominator Q^{|m|}. N = Σ restᵢ · Q^{eᵢ − m}.
            const std::size_t mm = static_cast<std::size_t>((-m).to_u64());
            Poly N{Rational(BigInt(0))};
            for (const auto& [e, r] : parts)
                N = poly_add(N, poly_mul(r, poly_pow(*base_poly,
                        static_cast<std::size_t>((e - m).to_u64()))));
            Poly q, rem;
            if (!poly_divmod(N, poly_pow(*base_poly, mm), q, rem)) continue;
            if (!poly_is_zero(rem)) continue;             // Q does not fully cancel
            std::vector<ExprPtr> out{poly_to_expr_for_combine(q, v)};
            terms.swap(out);
            return true;
        }

        // Surd denominator √P with |m| = 1. Split by sqrt parity:
        //   eᵢ = −1: rest/√P     eᵢ = +1: rest·√P = rest·P/√P     eᵢ = 0: rest (poly)
        //   anything else (|e| ≥ 2 / even ≠ 0): out of scope.
        if (m != BigInt(-1)) continue;
        Poly num_over_sqrt{Rational(BigInt(0))};   // collected numerator of (…)/√P
        Poly poly_part{Rational(BigInt(0))};       // the √-free remainder S
        bool ok_shape = true;
        for (const auto& [e, r] : parts) {
            if (e == BigInt(-1))      num_over_sqrt = poly_add(num_over_sqrt, r);
            else if (e == BigInt(1))  num_over_sqrt = poly_add(num_over_sqrt, poly_mul(r, *base_poly));
            else if (e == BigInt(0))  poly_part = poly_add(poly_part, r);
            else { ok_shape = false; break; }
        }
        if (!ok_shape) continue;

        // (num_over_sqrt)/√P + poly_part. The √ cancels iff num_over_sqrt is 0
        // (→ poly_part) or divisible by P (→ quotient·√P + poly_part).
        if (poly_is_zero(num_over_sqrt)) {
            std::vector<ExprPtr> out{poly_to_expr_for_combine(poly_part, v)};
            terms.swap(out);
            return true;
        }
        Poly q, rem;
        if (!poly_divmod(num_over_sqrt, *base_poly, q, rem)) continue;
        if (!poly_is_zero(rem)) continue;                 // √ denominator remains
        std::vector<ExprPtr> out;
        ExprPtr qsqrt = arena_.make<Product>(std::vector<ExprPtr>{
            poly_to_expr_for_combine(q, v), g});
        out.push_back(qsqrt);
        if (!poly_is_zero(poly_part)) out.push_back(poly_to_expr_for_combine(poly_part, v));
        terms.swap(out);
        return true;
    }
    return false;
}

}  // namespace cas::symbolic::detail
