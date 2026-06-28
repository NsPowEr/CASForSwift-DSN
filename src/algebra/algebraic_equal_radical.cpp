// F7.5 (B.2 algebraic) — single square-root extension equivalence helper for
// `mathematically_equal`. The radical analogue of `weierstrass_zero_diff` /
// `trig_exponential_zero_diff`.
//
// Antiderivatives of an integrand containing one radical s = √p(x) frequently
// come out in different but equal shapes: e.g. ∫√(x²−1) dx, where the CAS emits
//   ½·x·√(x²−1) − ½·acosh(x)         (acosh(x) = ln(x+√(x²−1)))
// and Maxima emits
//   ½·x·√(x²−1) − ½·log(2x+2√(x²−1)).
// These differ only by a constant, so their DERIVATIVES are equal; but the
// derivative of the Maxima form carries √(x²−1) in a denominator, which neither
// `simplify` (no conjugate rationalisation) nor `polynomial_normal_form` (treats
// √(x²−1) as an opaque atom, never applying √p·√p = p) reduces to the bare
// √(x²−1) the CAS side produces.
//
// Method: every sub-expression of a single-radical difference is an element of
// the field Q(x)(s) with s² = p(x). Represent it as
//   (A + B·s) / (C + D·s),     A,B,C,D ∈ Q[x] (polynomials, no radical)
// closed under +,−,·,/ and integer powers, using s² → p on every product. The
// difference is identically zero iff its numerator A + B·s ≡ 0.
//
// SOUND (REGOLA ZERO / A14): the zero test requires BOTH A ≡ 0 AND B ≡ 0 as
// polynomials (via `polynomial_normal_form`). A ≡ 0 ∧ B ≡ 0 ⟹ A + B·s ≡ 0
// *unconditionally* — independent of whether s is irrational over Q(x) — so it
// can only ever PROVE equality, never assert a false one. Anything outside the
// single-radical class (two distinct radicands, a radical nested in a non-sqrt
// function, a non-integer/symbolic power, a DecimalLit) bails to `false` and the
// caller continues its fallback chain.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/normal_form.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"

#include <optional>
#include <vector>

namespace cas::symbolic {

namespace {

// Is `e` a square root √q? Returns the radicand q (else nullptr). Handles both
// FuncCall sqrt(q) and the Pow(q, 1/2) spelling.
[[nodiscard]] ExprPtr as_sqrt(ExprPtr e) {
    if (const auto* fc = expr_cast<FuncCall>(e);
        fc != nullptr && fc->func_id == BuiltinOp::Sqrt && fc->args.size() == 1U) {
        return fc->args[0];
    }
    if (const auto* b = expr_cast<Binary>(e); b != nullptr && b->op == BinaryOp::Pow) {
        if (const auto* rl = expr_cast<RationalLit>(b->right);
            rl != nullptr && rl->numerator == BigInt(1) && rl->denominator == BigInt(2)) {
            return b->left;
        }
    }
    return nullptr;
}

bool contains_sqrt(ExprPtr e) {
    if (!e) return false;
    if (as_sqrt(e) != nullptr) return true;
    if (const auto* u = expr_cast<Unary>(e)) return contains_sqrt(u->operand);
    if (const auto* b = expr_cast<Binary>(e)) return contains_sqrt(b->left) || contains_sqrt(b->right);
    if (const auto* s = expr_cast<Sum>(e)) {
        for (ExprPtr t : s->terms) if (contains_sqrt(t)) return true;
        return false;
    }
    if (const auto* p = expr_cast<Product>(e)) {
        for (ExprPtr f : p->factors) if (contains_sqrt(f)) return true;
        return false;
    }
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        for (ExprPtr a : fc->args) if (contains_sqrt(a)) return true;
        return false;
    }
    return false;
}

// Collect distinct radicands (by structural identity). Returns false if any
// radical is nested inside another radical (not a single flat extension).
void collect_radicands(ExprPtr e, std::vector<ExprPtr>& out) {
    if (!e) return;
    if (ExprPtr q = as_sqrt(e); q != nullptr) {
        bool seen = false;
        for (ExprPtr r : out) if (structural_equal(q, r)) { seen = true; break; }
        if (!seen) out.push_back(q);
        collect_radicands(q, out);  // detect nested radicals (→ >1 radicand)
        return;
    }
    if (const auto* u = expr_cast<Unary>(e)) { collect_radicands(u->operand, out); return; }
    if (const auto* b = expr_cast<Binary>(e)) { collect_radicands(b->left, out); collect_radicands(b->right, out); return; }
    if (const auto* s = expr_cast<Sum>(e)) { for (ExprPtr t : s->terms) collect_radicands(t, out); return; }
    if (const auto* p = expr_cast<Product>(e)) { for (ExprPtr f : p->factors) collect_radicands(f, out); return; }
    if (const auto* fc = expr_cast<FuncCall>(e)) { for (ExprPtr a : fc->args) collect_radicands(a, out); return; }
}

// Element of Q(x)(s): num/den with each component a + b·s (a,b ∈ Q[x]).
struct R2 { ExprPtr a; ExprPtr b; };       // a + b·s
struct RatR { R2 num; R2 den; };

class RadicalField {
public:
    RadicalField(AstArena& arena, ExprPtr p) : arena_(arena), p_(p) {}

    ExprPtr zero() { return arena_.make<IntegerLit>(BigInt(0)); }
    ExprPtr one()  { return arena_.make<IntegerLit>(BigInt(1)); }

    ExprPtr add(ExprPtr x, ExprPtr y) { return arena_.make<Binary>(BinaryOp::Add, x, y); }
    ExprPtr sub(ExprPtr x, ExprPtr y) { return arena_.make<Binary>(BinaryOp::Sub, x, y); }
    ExprPtr mul(ExprPtr x, ExprPtr y) { return arena_.make<Binary>(BinaryOp::Mul, x, y); }

    // (a1+b1 s)(a2+b2 s) = (a1 a2 + b1 b2 p) + (a1 b2 + a2 b1) s
    R2 r2_mul(const R2& x, const R2& y) {
        return R2{ add(mul(x.a, y.a), mul(mul(x.b, y.b), p_)),
                   add(mul(x.a, y.b), mul(x.b, y.a)) };
    }
    R2 r2_add(const R2& x, const R2& y) { return R2{ add(x.a, y.a), add(x.b, y.b) }; }
    R2 r2_neg(const R2& x) { return R2{ sub(zero(), x.a), sub(zero(), x.b) }; }

    RatR lit(ExprPtr scalar) { return RatR{ R2{scalar, zero()}, R2{one(), zero()} }; }
    RatR radical() { return RatR{ R2{zero(), one()}, R2{one(), zero()} }; }  // s itself

    RatR rmul(const RatR& f, const RatR& g) {
        return RatR{ r2_mul(f.num, g.num), r2_mul(f.den, g.den) };
    }
    RatR rdiv(const RatR& f, const RatR& g) {
        return RatR{ r2_mul(f.num, g.den), r2_mul(f.den, g.num) };
    }
    RatR radd(const RatR& f, const RatR& g, bool subtract) {
        R2 cross = r2_mul(subtract ? r2_neg(g.num) : g.num, f.den);
        return RatR{ r2_add(r2_mul(f.num, g.den), cross), r2_mul(f.den, g.den) };
    }
    RatR rneg(const RatR& f) { return RatR{ r2_neg(f.num), f.den }; }

private:
    AstArena& arena_;
    ExprPtr p_;
};

// Recursively map `e` into Q(x)(s). nullopt on anything outside the class.
[[nodiscard]] std::optional<RatR> to_ratr(ExprPtr e, ExprPtr p, RadicalField& F) {
    if (ExprPtr q = as_sqrt(e); q != nullptr) {
        if (!structural_equal(q, p)) return std::nullopt;  // different radicand
        return F.radical();
    }
    switch (e->kind) {
        case ExprKind::IntegerLit:
        case ExprKind::RationalLit:
        case ExprKind::Constant:
        case ExprKind::ComplexLit:
        case ExprKind::Symbol:
            return F.lit(e);  // radical-free scalar leaf (polynomial coefficient)
        case ExprKind::DecimalLit:
            return std::nullopt;  // symbolic/numeric boundary (CLAUDE.md §5)
        case ExprKind::Unary: {
            const auto* u = expr_cast<Unary>(e);
            if (u->op != UnaryOp::Neg) return std::nullopt;
            auto inner = to_ratr(u->operand, p, F);
            if (!inner) return std::nullopt;
            return F.rneg(*inner);
        }
        case ExprKind::Sum: {
            const auto* s = expr_cast<Sum>(e);
            RatR acc = F.lit(F.zero());
            for (ExprPtr t : s->terms) {
                auto tm = to_ratr(t, p, F);
                if (!tm) return std::nullopt;
                acc = F.radd(acc, *tm, /*subtract=*/false);
            }
            return acc;
        }
        case ExprKind::Product: {
            const auto* pr = expr_cast<Product>(e);
            RatR acc = F.lit(F.one());
            for (ExprPtr f : pr->factors) {
                auto fm = to_ratr(f, p, F);
                if (!fm) return std::nullopt;
                acc = F.rmul(acc, *fm);
            }
            return acc;
        }
        case ExprKind::Binary: {
            const auto* b = expr_cast<Binary>(e);
            if (b->op == BinaryOp::Add || b->op == BinaryOp::Sub) {
                auto l = to_ratr(b->left, p, F); if (!l) return std::nullopt;
                auto r = to_ratr(b->right, p, F); if (!r) return std::nullopt;
                return F.radd(*l, *r, /*subtract=*/b->op == BinaryOp::Sub);
            }
            if (b->op == BinaryOp::Mul) {
                auto l = to_ratr(b->left, p, F); if (!l) return std::nullopt;
                auto r = to_ratr(b->right, p, F); if (!r) return std::nullopt;
                return F.rmul(*l, *r);
            }
            if (b->op == BinaryOp::Div) {
                auto l = to_ratr(b->left, p, F); if (!l) return std::nullopt;
                auto r = to_ratr(b->right, p, F); if (!r) return std::nullopt;
                return F.rdiv(*l, *r);
            }
            if (b->op == BinaryOp::Pow) {
                const auto* il = expr_cast<IntegerLit>(b->right);
                if (il == nullptr) return std::nullopt;  // non-integer power (other than the √ handled above)
                auto base = to_ratr(b->left, p, F); if (!base) return std::nullopt;
                const bool neg = il->value.is_negative();
                RatR unit = F.lit(F.one());
                RatR acc = unit;
                const std::uint64_t n = il->value.abs().to_u64();
                for (std::uint64_t i = 0; i < n; ++i) acc = F.rmul(acc, *base);
                return neg ? F.rdiv(unit, acc) : acc;
            }
            return std::nullopt;
        }
        case ExprKind::FuncCall: {
            // A non-sqrt function is an opaque scalar leaf ONLY if it contains no
            // radical; a radical nested in a transcendental cannot be linearised.
            if (contains_sqrt(e)) return std::nullopt;
            return F.lit(e);
        }
        default:
            return std::nullopt;
    }
}

}  // namespace

// True iff `diff_expr` is identically zero as a rational function over a single
// square-root extension Q(x)(√p). Returns false (no claim) otherwise.
[[nodiscard]] bool radical_zero_diff(ExprPtr diff_expr, CASContext& ctx) {
    std::vector<ExprPtr> radicands;
    collect_radicands(diff_expr, radicands);
    if (radicands.size() != 1U) return false;  // need exactly one radical extension
    ExprPtr p = radicands.front();

    RadicalField F(ctx.arena(), p);
    auto r = to_ratr(diff_expr, p, F);
    if (!r) return false;

    auto is_zero_poly = [&](ExprPtr poly) {
        auto nf = polynomial_normal_form(poly, ctx);
        if (nf.is_error()) return false;
        const auto* il = expr_cast<IntegerLit>(nf.value());
        return il != nullptr && il->value.is_zero();
    };
    // diff = (A + B·s)/(C + D·s) ≡ 0  ⟺  A ≡ 0 ∧ B ≡ 0  (numerator vanishes).
    // Guard against a degenerate 0/0: require the denominator ≢ 0.
    if (is_zero_poly(r->den.a) && is_zero_poly(r->den.b)) return false;
    return is_zero_poly(r->num.a) && is_zero_poly(r->num.b);
}

}  // namespace cas::symbolic
