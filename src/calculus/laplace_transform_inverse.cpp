// CAS-L3-07 — Inverse Laplace transform L⁻¹{F(s)}(t).
// Split out of laplace_transform.cpp (T-050 anti-monolith).

#include "calculus_internal.hpp"
#include "laplace_transform_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <vector>

namespace cas::calculus {

using laplace_detail::make_int;
using laplace_detail::is_sym;
using laplace_detail::factorial_expr;
using laplace_detail::depends_on_t;

namespace {

// Detect s² + a² where a² rational positive square or rational positive.
// Returns a (as expr) on match, else nullptr.
[[nodiscard]] ExprPtr match_s_squared_plus_const(ExprPtr expr, const Symbol& s,
                                                  AstArena& arena) {
    const auto* sum = expr_cast<Sum>(expr);
    if (!sum || sum->terms.size() != 2) return nullptr;
    ExprPtr s_sq_part = nullptr;
    ExprPtr const_part = nullptr;
    for (ExprPtr term : sum->terms) {
        if (const auto* bin = expr_cast<Binary>(term);
            bin && bin->op == BinaryOp::Pow) {
            if (const auto* sym = expr_cast<Symbol>(bin->left);
                sym && sym->name == s.name) {
                if (const auto* il = expr_cast<IntegerLit>(bin->right);
                    il && il->value == BigInt(2)) {
                    if (s_sq_part) return nullptr;
                    s_sq_part = term;
                    continue;
                }
            }
        }
        const_part = term;
    }
    if (!s_sq_part || !const_part) return nullptr;
    // a² extracted; need sqrt(a²) = a as expr. For integer/rational c>0:
    if (const auto* il = expr_cast<IntegerLit>(const_part)) {
        if (il->value.is_negative() || il->value.is_zero()) return nullptr;
        // Try integer sqrt.
        BigInt c = il->value;
        BigInt x = c;
        BigInt y = (x + BigInt(1)) / BigInt(2);
        while (y < x) { x = y; y = (x + c / x) / BigInt(2); }
        if (x * x == c) return arena.make<IntegerLit>(std::move(x));
        // Else a = sqrt(c) symbolic.
        return arena.make<FuncCall>(BuiltinOp::Sqrt,
            std::vector<ExprPtr>{const_part});
    }
    return nullptr;
}

}  // namespace

// Inverse Laplace L⁻¹{F(s)}(t) — elementary pattern table.
// Covers:
//   1/s → 1
//   1/s^(n+1) → t^n / n!
//   1/(s-a) → exp(a·t)
//   a/(s²+a²) → sin(a·t)
//   s/(s²+a²) → cos(a·t)
//   Linearity over Sum
[[nodiscard]] Result<ExprPtr> inverse_laplace_transform(
    ExprPtr expr, const Symbol& s, const Symbol& t,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr t_expr = arena.make<Symbol>(t);

    // Linearity over Sum.
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        for (ExprPtr term : sum->terms) {
            auto inv = inverse_laplace_transform(term, s, t, ctx);
            if (inv.is_error()) return inv;
            terms.push_back(inv.value());
        }
        return ctx.simplify(arena.make<Sum>(std::move(terms)));
    }

    // Unary(Neg, X) → -L⁻¹{X}
    if (const auto* un = expr_cast<Unary>(expr); un && un->op == UnaryOp::Neg) {
        auto inner = inverse_laplace_transform(un->operand, s, t, ctx);
        if (inner.is_error()) return inner;
        return ctx.simplify(arena.make<Unary>(UnaryOp::Neg, inner.value()));
    }

    // Constant scalar factor: c · F(s) → c · L⁻¹{F(s)}
    if (const auto* prod = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> consts, vars;
        for (ExprPtr f : prod->factors) {
            if (depends_on_t(f, s)) vars.push_back(f);
            else consts.push_back(f);
        }
        if (!consts.empty() && !vars.empty()) {
            ExprPtr f_only = vars.size() == 1 ? vars[0]
                : static_cast<ExprPtr>(arena.make<Product>(std::move(vars)));
            auto inv = inverse_laplace_transform(f_only, s, t, ctx);
            if (inv.is_error()) return inv;
            consts.push_back(inv.value());
            return ctx.simplify(arena.make<Product>(std::move(consts)));
        }
    }

    // 1 / Q(s) patterns. Detect both Binary(Div, num, den) AND
    // Binary(Pow, den, -1) (== 1/den).
    ExprPtr num_extracted = nullptr;
    ExprPtr den_extracted = nullptr;
    if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Div) {
        num_extracted = bin->left;
        den_extracted = bin->right;
    } else if (const auto* bin = expr_cast<Binary>(expr); bin && bin->op == BinaryOp::Pow) {
        if (auto* exp_lit = expr_cast<IntegerLit>(bin->right);
            exp_lit && exp_lit->value == BigInt(-1)) {
            num_extracted = arena.make<IntegerLit>(BigInt(1));
            den_extracted = bin->left;
        }
    } else if (const auto* prod = expr_cast<Product>(expr)) {
        // Detect Product with Pow(_, -1) factor as denominator.
        std::vector<ExprPtr> num_factors;
        ExprPtr den = nullptr;
        bool ok2 = true;
        for (auto f : prod->factors) {
            if (auto* pb = expr_cast<Binary>(f); pb && pb->op == BinaryOp::Pow) {
                if (auto* el = expr_cast<IntegerLit>(pb->right);
                    el && el->value == BigInt(-1)) {
                    if (den) { ok2 = false; break; }
                    den = pb->left;
                    continue;
                }
            }
            num_factors.push_back(f);
        }
        if (ok2 && den) {
            num_extracted = num_factors.empty()
                ? static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(1)))
                : (num_factors.size() == 1 ? num_factors[0]
                   : static_cast<ExprPtr>(arena.make<Product>(std::move(num_factors))));
            den_extracted = den;
        }
    }
    if (num_extracted && den_extracted) {
        ExprPtr num = num_extracted;
        ExprPtr den = den_extracted;

        // num must be constant in s for most patterns; or s for cos pattern.
        bool num_is_s = false;
        if (const auto* sym = expr_cast<Symbol>(num); sym && sym->name == s.name) {
            num_is_s = true;
        }

        // 1/s → 1
        if (is_sym(den, s) && !depends_on_t(num, s)) {
            return ctx.simplify(num);
        }

        // 1/s^(n+1) → num · t^n / n!
        if (const auto* den_pow = expr_cast<Binary>(den);
            den_pow && den_pow->op == BinaryOp::Pow) {
            if (is_sym(den_pow->left, s)) {
                if (const auto* exp_lit = expr_cast<IntegerLit>(den_pow->right);
                    exp_lit && exp_lit->value > BigInt(0)
                    && !depends_on_t(num, s)) {
                    const long long n_plus_1 = static_cast<long long>(exp_lit->value.to_u64());
                    const long long n = n_plus_1 - 1;
                    ExprPtr t_pow_n = (n == 0)
                        ? static_cast<ExprPtr>(make_int(ctx, 1))
                        : static_cast<ExprPtr>(arena.make<Binary>(BinaryOp::Pow,
                              t_expr, make_int(ctx, n)));
                    ExprPtr fact = factorial_expr(ctx, n);
                    ExprPtr ratio = arena.make<Binary>(BinaryOp::Div, t_pow_n, fact);
                    ExprPtr scaled = arena.make<Product>(std::vector<ExprPtr>{num, ratio});
                    return ctx.simplify(scaled);
                }
            }
        }

        // 1/(s-a) → exp(a·t).  Detect den == Binary(Sub, s, a) OR
        // den == Sum([s + c]) where c constant (then a = -c).
        auto try_match_s_minus_a = [&](ExprPtr d) -> ExprPtr {
            if (const auto* db = expr_cast<Binary>(d);
                db && db->op == BinaryOp::Sub
                && is_sym(db->left, s) && !depends_on_t(db->right, s)) {
                return db->right;
            }
            if (const auto* sm = expr_cast<Sum>(d); sm && sm->terms.size() == 2) {
                ExprPtr s_term = nullptr;
                ExprPtr c_term = nullptr;
                for (auto term : sm->terms) {
                    if (is_sym(term, s)) {
                        if (s_term) return nullptr;
                        s_term = term;
                    } else if (!depends_on_t(term, s)) {
                        if (c_term) return nullptr;
                        c_term = term;
                    } else {
                        return nullptr;
                    }
                }
                if (!s_term || !c_term) return nullptr;
                // den = s + c → a = -c
                return arena.make<Unary>(UnaryOp::Neg, c_term);
            }
            return nullptr;
        };
        if (!depends_on_t(num, s)) {
            ExprPtr a_extract = try_match_s_minus_a(den);
            if (a_extract) {
                ExprPtr at = arena.make<Product>(std::vector<ExprPtr>{a_extract, t_expr});
                ExprPtr exp_at = arena.make<FuncCall>(BuiltinOp::Exp,
                    std::vector<ExprPtr>{at});
                ExprPtr scaled = arena.make<Product>(std::vector<ExprPtr>{num, exp_at});
                return ctx.simplify(scaled);
            }
        }

        // a/(s²+a²) → sin(a·t)   or   s/(s²+a²) → cos(a·t)
        ExprPtr a_extracted = match_s_squared_plus_const(den, s, arena);
        if (a_extracted) {
            ExprPtr a = a_extracted;
            ExprPtr at = arena.make<Product>(std::vector<ExprPtr>{a, t_expr});
            // cos pattern: num == s
            if (num_is_s) {
                return ctx.simplify(arena.make<FuncCall>(BuiltinOp::Cos,
                    std::vector<ExprPtr>{at}));
            }
            // sin pattern: num == a (or constant proportional to a)
            // Verify num structurally equals a (simple case).
            if (structural_equal(num, a)) {
                return ctx.simplify(arena.make<FuncCall>(BuiltinOp::Sin,
                    std::vector<ExprPtr>{at}));
            }
            // num = c (constant), → (c/a) · sin(a·t)
            if (!depends_on_t(num, s)) {
                ExprPtr c_over_a = arena.make<Binary>(BinaryOp::Div, num, a);
                ExprPtr sin_at = arena.make<FuncCall>(BuiltinOp::Sin,
                    std::vector<ExprPtr>{at});
                return ctx.simplify(arena.make<Product>(
                    std::vector<ExprPtr>{c_over_a, sin_at}));
            }
        }
    }

    // PFD fallback: if F(s) is rational in s, decompose via partial
    // fractions and apply pattern table per term.
    if (num_extracted && den_extracted) {
        // Build expr = num/den as Binary(Div) for partial_fractions.
        ExprPtr rational_expr = arena.make<Binary>(BinaryOp::Div,
            num_extracted, den_extracted);
        auto pfd_res = algebra::partial_fractions(rational_expr, s, ctx);
        if (pfd_res.is_ok()) {
            const auto& terms = pfd_res.value();
            if (terms.size() >= 2U) {
                std::vector<ExprPtr> inverse_terms;
                bool all_ok = true;
                for (auto term : terms) {
                    auto inv = inverse_laplace_transform(term, s, t, ctx);
                    if (inv.is_error()) { all_ok = false; break; }
                    inverse_terms.push_back(inv.value());
                }
                if (all_ok) {
                    return ctx.simplify(arena.make<Sum>(std::move(inverse_terms)));
                }
            }
        }
    }

    // Fallback finale: Bronstein residue formula (F5.8 / Task #17).
    //   L⁻¹{F(s)}(t) = Σ_k Res_{s=p_k}[F(s)·e^(s·t)].
    //   Più generale del pattern table elementare; copre F razionale con
    //   polos qualsiasi (anche multipli).  Restituisce Unimplemented solo
    //   se anche il residue path fallisce.
    {
        auto via_residue = inverse_laplace_residue_q(expr, s, t, ctx);
        if (via_residue.is_ok()) return via_residue;
    }

    // F0.8-MIGRATED
    return make_unimplemented<ExprPtr>(
        "calculus", "inverse_laplace_transform",
        "expression not in elementary inverse Laplace table nor residue-decomposable",
        cas::error::reason_codes::LAPLACE_UNKNOWN_FORM,
        "Extend the inverse Laplace table or implement Bromwich integral / PFD bridge",
        "F0.8");
}

}  // namespace cas::calculus
