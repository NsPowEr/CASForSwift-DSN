// A7 step 5 — Meijer G fallback integrator (Meijer_G_Slater.md §8-§9).
// Runs ONLY after Risch and Weierstrass have failed (wired as the last
// resort in Integrator::integrate — spec §8 "wiring POST-Risch").
//
// Handled family (general in the parameters, structural in the shape):
//   int  K * x^mu * f(x) dx     with K var-free, mu a rational literal,
// where to_meijerg(f) = P * G(c * x^r) with P var-free-or-monomial,
// c var-free and r a positive rational literal. Chain (all spec-verified):
//   u-substitution t = x^r  ->  (1/r) int t^{(mu+1)/r - 1} G(c t) dt,
//   §6.2 power shift absorbs t^nu (constant c^{-nu}),
//   §6.6 antiderivative  int G(ct) dt = t * G^{m,n+1}_{p+1,q+1}(ct|0,a;b,-1),
//   §6.3 parameter cancellation, then expand_meijerg_nodes: the result is
// elementary when the inverse table / Slater + the engine's pFq identities
// can fold it (e.g. int e^{-x^2} dx -> (sqrt(pi)/2) erf(x)); otherwise it
// legitimately STAYS a Meijer G / pFq closed form (§9.4 — a first-class
// antiderivative, not an error). Anything outside the shape family fails
// with a structured Unimplemented, never a wrong result.

#include "integrate_engine.hpp"

#include "cas/error_helpers.hpp"
#include "cas/meijerg.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

using symbolic::CASContext;
using integrate_detail::depends_on;
using integrate_detail::is_same_symbol;
using integrate_detail::make_integer;

struct VarMonomial {
    ExprPtr coeff;    // var-free (nullptr = 1)
    Rational power;   // exponent of var
};

// expr == c * var^r with c var-free and r a rational literal (r may be 0 =>
// pure var-free expression). Conservative: unknown shapes return nullopt.
[[nodiscard]] std::optional<VarMonomial> match_var_monomial(
    ExprPtr expr, const Symbol& var, AstArena& arena) {
    if (expr == nullptr) return std::nullopt;
    if (!depends_on(expr, var)) return VarMonomial{expr, Rational(BigInt(0))};
    if (is_same_symbol(expr, var)) {
        return VarMonomial{nullptr, Rational(BigInt(1))};
    }
    if (const auto* un = expr_cast<Unary>(expr); un != nullptr && un->op == UnaryOp::Neg) {
        auto inner = match_var_monomial(un->operand, var, arena);
        if (!inner.has_value()) return std::nullopt;
        ExprPtr minus_one = make_integer(arena, -1);
        inner->coeff = inner->coeff == nullptr
            ? minus_one
            : arena.make<Product>(std::vector<ExprPtr>{minus_one, inner->coeff});
        return inner;
    }
    if (const auto* pw = expr_cast<Binary>(expr);
        pw != nullptr && pw->op == BinaryOp::Pow && is_same_symbol(pw->left, var)) {
        if (const auto* il = expr_cast<IntegerLit>(pw->right)) {
            return VarMonomial{nullptr, Rational(il->value)};
        }
        if (const auto* rl = expr_cast<RationalLit>(pw->right)) {
            return VarMonomial{nullptr, Rational(rl->numerator, rl->denominator)};
        }
        return std::nullopt;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        Rational power(BigInt(0));
        std::vector<ExprPtr> coeff_factors;
        for (ExprPtr f : prod->factors) {
            auto part = match_var_monomial(f, var, arena);
            if (!part.has_value()) return std::nullopt;
            power = power + part->power;
            if (part->coeff != nullptr) coeff_factors.push_back(part->coeff);
        }
        ExprPtr coeff = nullptr;
        if (coeff_factors.size() == 1U) coeff = coeff_factors.front();
        else if (!coeff_factors.empty())
            coeff = arena.make<Product>(std::move(coeff_factors));
        return VarMonomial{coeff, power};
    }
    return std::nullopt;
}

[[nodiscard]] ExprPtr rational_lit(AstArena& arena, const Rational& r) {
    if (r.is_integer()) return arena.make<IntegerLit>(r.numerator());
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

[[nodiscard]] std::optional<Rational> rational_of(ExprPtr e) {
    if (const auto* il = expr_cast<IntegerLit>(e)) return Rational(il->value);
    if (const auto* rl = expr_cast<RationalLit>(e))
        return Rational(rl->numerator, rl->denominator);
    return std::nullopt;
}

// True when the expression contains a node of the pFq/G/Bessel world — the
// ops whose §5/§3.1 to_meijerg entries are SPECIAL-function sources (not the
// elementary rows). Extend together with the table (e.g. incomplete gamma
// when §5.9 lands). Governs the §9.4 result policy below: a Meijer G
// antiderivative is a first-class RESULT only when the integrand already
// lives in this world; for purely elementary integrands an unfolded G means
// "inverse table cannot fold this yet" -> structured Unimplemented (the
// pre-A7 behaviour, never a readability downgrade vs the elementary form
// that may exist).
[[nodiscard]] bool contains_special_fn_node(ExprPtr e) {
    if (e == nullptr) return false;
    if (const auto* call = expr_cast<FuncCall>(e)) {
        switch (call->func_id) {
            case BuiltinOp::MeijerG:
            case BuiltinOp::BesselJ:
            case BuiltinOp::Hypergeometric0F1:
            case BuiltinOp::Hypergeometric1F1:
            case BuiltinOp::Hypergeometric2F1:
                return true;
            default: break;
        }
        for (ExprPtr a : call->args)
            if (contains_special_fn_node(a)) return true;
        return false;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr f : prod->factors)
            if (contains_special_fn_node(f)) return true;
        return false;
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        for (ExprPtr t : sum->terms)
            if (contains_special_fn_node(t)) return true;
        return false;
    }
    if (const auto* bin = expr_cast<Binary>(e))
        return contains_special_fn_node(bin->left)
            || contains_special_fn_node(bin->right);
    if (const auto* un = expr_cast<Unary>(e))
        return contains_special_fn_node(un->operand);
    return false;
}

[[nodiscard]] bool contains_meijerg_node(ExprPtr e) {
    if (e == nullptr) return false;
    if (const auto* call = expr_cast<FuncCall>(e)) {
        if (call->func_id == BuiltinOp::MeijerG) return true;
        for (ExprPtr a : call->args)
            if (contains_meijerg_node(a)) return true;
        return false;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr f : prod->factors)
            if (contains_meijerg_node(f)) return true;
        return false;
    }
    if (const auto* sum = expr_cast<Sum>(e)) {
        for (ExprPtr t : sum->terms)
            if (contains_meijerg_node(t)) return true;
        return false;
    }
    if (const auto* bin = expr_cast<Binary>(e))
        return contains_meijerg_node(bin->left) || contains_meijerg_node(bin->right);
    if (const auto* un = expr_cast<Unary>(e))
        return contains_meijerg_node(un->operand);
    return false;
}

// Post-pass over the expanded antiderivative: fold (c*x^r)^s -> c^s * x^{r*s}
// wherever the base is a var-monomial and s is a rational literal (Sqrt = s
// 1/2). This is the explicit implementation of the Mellin half-line
// convention (Meijer_G_Slater.md §4, §9): with x > 0 the fold is EXACT on
// principal branches for any var-free c != 0 (arg(x^r) = 0), and it is what
// turns the raw expansion's sqrt(x^2)-type artifacts into the intended
// x-powers. Unconditionally exact only for integer s AND integer r; every
// other fold is recorded by the caller as a Positive(var) side condition
// (A31) — never applied silently.
struct MonomialPowerFold {
    AstArena& arena;
    const Symbol& var;
    bool conditional = false;

    [[nodiscard]] ExprPtr fold_pow(ExprPtr base, const Rational& s, ExprPtr orig) {
        auto mono = match_var_monomial(base, var, arena);
        if (!mono.has_value() || mono->power.numerator().is_zero()) return orig;
        const Rational rs = mono->power * s;
        const bool exact = s.is_integer() && mono->power.is_integer();
        if (!exact) conditional = true;
        std::vector<ExprPtr> factors;
        if (mono->coeff != nullptr) {
            factors.push_back(arena.make<Binary>(BinaryOp::Pow, mono->coeff,
                                                 rational_lit(arena, s)));
        }
        if (rs == Rational(BigInt(1))) {
            factors.push_back(arena.make<Symbol>(var));
        } else if (!rs.numerator().is_zero()) {
            factors.push_back(arena.make<Binary>(BinaryOp::Pow,
                arena.make<Symbol>(var), rational_lit(arena, rs)));
        }
        if (factors.empty()) return make_integer(arena, 1);
        if (factors.size() == 1U) return factors.front();
        return arena.make<Product>(std::move(factors));
    }

    // Structural sharing: children are rebuilt only when a fold changed them.
    [[nodiscard]] ExprPtr walk(ExprPtr e) {
        if (e == nullptr) return e;
        if (const auto* bin = expr_cast<Binary>(e)) {
            if (bin->op == BinaryOp::Pow) {
                if (auto s = rational_of(bin->right); s.has_value()) {
                    ExprPtr base = walk(bin->left);
                    ExprPtr folded = fold_pow(base,
                        *s, base == bin->left ? e : arena.make<Binary>(
                            BinaryOp::Pow, base, bin->right));
                    return folded;
                }
            }
            ExprPtr l = walk(bin->left);
            ExprPtr r = walk(bin->right);
            if (l == bin->left && r == bin->right) return e;
            return arena.make<Binary>(bin->op, l, r);
        }
        if (const auto* call = expr_cast<FuncCall>(e)) {
            if (call->func_id == BuiltinOp::Sqrt && call->args.size() == 1U) {
                ExprPtr inner = walk(call->args.front());
                ExprPtr orig = inner == call->args.front()
                    ? e
                    : arena.make<FuncCall>(BuiltinOp::Sqrt,
                                           std::vector<ExprPtr>{inner});
                return fold_pow(inner, Rational(BigInt(1), BigInt(2)), orig);
            }
            std::vector<ExprPtr> args;
            bool changed = false;
            args.reserve(call->args.size());
            for (ExprPtr a : call->args) {
                ExprPtr w = walk(a);
                changed = changed || (w != a);
                args.push_back(w);
            }
            if (!changed) return e;
            return arena.make<FuncCall>(call->func_id, std::move(args));
        }
        if (const auto* prod = expr_cast<Product>(e)) {
            std::vector<ExprPtr> fs;
            bool changed = false;
            fs.reserve(prod->factors.size());
            for (ExprPtr f : prod->factors) {
                ExprPtr w = walk(f);
                changed = changed || (w != f);
                fs.push_back(w);
            }
            if (!changed) return e;
            return arena.make<Product>(std::move(fs));
        }
        if (const auto* sum = expr_cast<Sum>(e)) {
            std::vector<ExprPtr> ts;
            bool changed = false;
            ts.reserve(sum->terms.size());
            for (ExprPtr t : sum->terms) {
                ExprPtr w = walk(t);
                changed = changed || (w != t);
                ts.push_back(w);
            }
            if (!changed) return e;
            return arena.make<Sum>(std::move(ts));
        }
        if (const auto* un = expr_cast<Unary>(e)) {
            ExprPtr w = walk(un->operand);
            if (w == un->operand) return e;
            return arena.make<Unary>(un->op, w);
        }
        return e;
    }
};

}  // namespace

Result<ExprPtr> integrate_meijerg_fallback(
    ExprPtr expr, const Symbol& var, CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto unsupported = [&](const char* shape) {
        return make_unimplemented<ExprPtr>(
            "calculus", "integrate_meijerg_fallback", shape,
            cas::error::reason_codes::GENERIC,
            "Integrand outside the K*x^mu*f(c*x^r) Meijer-G family "
            "(Meijer_G_Slater.md §9)",
            "A7");
    };

    // 1. Split the integrand: var-free coefficient, x^mu monomial, one core.
    std::vector<ExprPtr> coeff_factors;
    Rational mu(BigInt(0));
    std::vector<ExprPtr> core;
    auto classify = [&](ExprPtr f) {
        auto mono = match_var_monomial(f, var, arena);
        if (mono.has_value()) {
            mu = mu + mono->power;
            if (mono->coeff != nullptr) coeff_factors.push_back(mono->coeff);
        } else {
            core.push_back(f);
        }
    };
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) classify(f);
    } else {
        classify(expr);
    }
    if (core.size() != 1U) {
        return unsupported(core.empty()
            ? "pure monomial (already Risch territory)"
            : "more than one non-monomial factor");
    }

    // 2. Convert the core to prefactor * G.
    auto g_form = symbolic::to_meijerg(ctx, core.front());
    if (g_form.is_error()) return g_form;
    const FuncCall* g_node = nullptr;
    auto absorb = [&](ExprPtr f) -> bool {
        if (const auto* call = expr_cast<FuncCall>(f);
            call != nullptr && call->func_id == BuiltinOp::MeijerG
            && g_node == nullptr) {
            g_node = call;
            return true;
        }
        auto mono = match_var_monomial(f, var, arena);
        if (!mono.has_value()) return false;
        mu = mu + mono->power;
        if (mono->coeff != nullptr) coeff_factors.push_back(mono->coeff);
        return true;
    };
    if (const auto* prod = expr_cast<Product>(g_form.value())) {
        for (ExprPtr f : prod->factors) {
            if (!absorb(f)) return unsupported("non-monomial G prefactor");
        }
    } else if (!absorb(g_form.value())) {
        return unsupported("to_meijerg result shape");
    }
    if (g_node == nullptr) return unsupported("no MeijerG node produced");

    // 3. The G argument must be c * x^r, r a positive rational literal.
    auto view = symbolic::view_meijerg(*g_node);
    if (view.is_error()) return fail<ExprPtr>(view.error());
    auto z_mono = match_var_monomial(view.value().z, var, arena);
    if (!z_mono.has_value() || !z_mono->power.numerator().is_positive()) {
        return unsupported("G argument is not c * x^r with r > 0");
    }
    const Rational r = z_mono->power;
    ExprPtr c = z_mono->coeff;  // nullptr = 1

    // 4. nu = (mu+1)/r - 1; absorb t^nu via §6.2 (constant factor c^{-nu}).
    const Rational nu = (mu + Rational(BigInt(1))) / r - Rational(BigInt(1));
    const FuncCall* work = g_node;
    ExprPtr shifted_holder;
    if (!nu.numerator().is_zero()) {
        auto shifted = symbolic::meijerg_power_shift(
            ctx, *work, rational_lit(arena, nu));
        if (shifted.is_error()) return shifted;
        shifted_holder = shifted.value();
        work = expr_cast<FuncCall>(shifted_holder);
        if (work == nullptr) return unsupported("power shift result shape");
    }

    // 5. §6.6 antiderivative: gives (c x^r) * G'(c x^r); the 1/c and the
    //    substitution bookkeeping fold into the constant below.
    auto antider = symbolic::meijerg_antiderivative(ctx, *work);
    if (antider.is_error()) return antider;

    // 6. Constant: K * (1/r) * c^{-nu} * (1/c) = K * (1/r) * c^{-(nu+1)}.
    std::vector<ExprPtr> result_factors = std::move(coeff_factors);
    const Rational inv_r = Rational(BigInt(1)) / r;
    if (!(inv_r == Rational(BigInt(1)))) {
        result_factors.push_back(rational_lit(arena, inv_r));
    }
    if (c != nullptr) {
        const Rational neg_nu_minus_one =
            (Rational(BigInt(0)) - nu) - Rational(BigInt(1));
        if (!neg_nu_minus_one.numerator().is_zero()) {
            result_factors.push_back(arena.make<Binary>(BinaryOp::Pow, c,
                rational_lit(arena, neg_nu_minus_one)));
        }
    }
    result_factors.push_back(antider.value());
    ExprPtr assembled = result_factors.size() == 1U
        ? result_factors.front()
        : arena.make<Product>(std::move(result_factors));

    // 7. §6.3 cleanup inside, then try to fold back to elementary/pFq form.
    auto expanded = symbolic::expand_meijerg_nodes(ctx, assembled);
    if (expanded.is_error()) return expanded;

    // 8. Mellin half-line back-substitution (§4): fold (c x^r)^s artifacts
    //    into plain x-powers. Non-trivially-exact folds are declared via a
    //    Positive(var) side condition (A31) — emitted AFTER the final
    //    top-level simplify, which clears the accumulator on entry.
    MonomialPowerFold fold{arena, var};
    ExprPtr folded = fold.walk(expanded.value());
    auto simplified = ctx.simplify(folded);
    if (simplified.is_error()) return simplified;

    // §9.4 result policy (see contains_special_fn_node): an antiderivative
    // still carrying a Meijer G node is returned as-is only for integrands
    // of the special world; on purely elementary integrands it would be a
    // readability downgrade (an elementary form may exist that the inverse
    // table cannot reach yet) -> structured refusal, golden stays a SKIP.
    if (contains_meijerg_node(simplified.value())
        && !contains_special_fn_node(expr)) {
        return make_unimplemented<ExprPtr>(
            "calculus", "integrate_meijerg_fallback",
            "unfolded Meijer G on purely elementary integrand",
            cas::error::reason_codes::GENERIC,
            "Extend from_meijerg inverse table / Slater folds for this "
            "G shape (Meijer_G_Slater.md §5 inverse, A7 residuo)",
            "A7");
    }
    if (fold.conditional) {
        auto cond = ctx.emit_side_condition(
            symbolic::DomainConditionKind::Positive,
            arena.make<Symbol>(var));
        if (cond.is_error()) return fail<ExprPtr>(cond.error());
    }
    return simplified;
}

}  // namespace cas::calculus
