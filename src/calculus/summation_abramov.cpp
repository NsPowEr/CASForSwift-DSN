// F5.7 — Abramov rational summation helpers.
//
// Sub-block 0 (Gosper):  definite sum via Newton-Leibniz finite calculus.
// Sub-block 1 (Polygamma): A/(linear(k))^m atoms via ψ^(m-1) antidifference.
// Sub-block 2 (Abramov-Full): partial-fraction decomposition + per-atom routing.
// B6-bis (Quadratic-RootOf): Q-irreducible quadratic atoms (B₁k+B₀)/Q(k) via
//   C·ψ(k−α) + D·ψ(k−β) where α,β = RootOf(Q̃,t,0/1) and C,D are residues.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/builtin_functions.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "summation_internal.hpp"
#include "../symbolic/summation_gosper.hpp"
#include <optional>
#include <vector>

namespace cas::calculus {

// ── static helpers ────────────────────────────────────────────────────────────

[[nodiscard]] static bool extract_rational_const(ExprPtr e, Rational& out) {
    if (const auto* il = expr_cast<IntegerLit>(e)) {
        out = Rational(il->value);
        return true;
    }
    if (const auto* rl = expr_cast<RationalLit>(e)) {
        out = Rational(rl->numerator, rl->denominator);
        return true;
    }
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (extract_rational_const(un->operand, out)) {
            out = -out;
            return true;
        }
    }
    return false;
}

[[nodiscard]] static bool parse_affine_in_k(
    ExprPtr e, const Symbol& k, ExprPtr& c0_out, ExprPtr& c1_out,
    symbolic::CASContext& ctx) {
    auto coeffs_res = algebra::univariate_coefficients(e, k, ctx);
    if (coeffs_res.is_error()) return false;
    const auto& cs = coeffs_res.value();
    if (cs.empty() || cs.size() > 2U) return false;
    c0_out = cs[0];
    c1_out = cs.size() == 2U ? cs[1] : ctx.arena().make<IntegerLit>(BigInt(0));
    auto c1_simp = ctx.simplify(c1_out);
    if (c1_simp.is_error()) return false;
    Rational c1_rat;
    if (!extract_rational_const(c1_simp.value(), c1_rat)) return false;
    if (c1_rat.numerator().is_zero()) return false;
    c1_out = c1_simp.value();
    return true;
}

// Build the antidifference of A/(k+a)^m as A·(−1)^(m−1)/(m−1)!·ψ^(m−1)(k+a).
// All arguments are ExprPtrs — symbolic or algebraic coefficients pass through.
[[nodiscard]] static ExprPtr polygamma_antidiff(
    ExprPtr A, ExprPtr k_plus_a, unsigned int m, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    BigInt sign = (m % 2U == 1U) ? BigInt(1) : BigInt(-1);
    BigInt fact(1);
    for (unsigned int i = 2U; i < m; ++i)
        fact *= BigInt(static_cast<long long>(i));
    ExprPtr scale_expr = rational_expr(arena, Rational(sign, fact));
    ExprPtr coefficient = arena.make<Binary>(BinaryOp::Mul, A, scale_expr);

    ExprPtr ant;
    if (m == 1U) {
        std::vector<ExprPtr> args{k_plus_a};
        ant = arena.make<FuncCall>(BuiltinOp::Digamma, std::move(args));
    } else {
        std::vector<ExprPtr> args{
            arena.make<IntegerLit>(BigInt(static_cast<long long>(m - 1U))),
            k_plus_a,
        };
        ant = arena.make<FuncCall>(BuiltinOp::Polygamma, std::move(args));
    }
    return arena.make<Binary>(BinaryOp::Mul, coefficient, ant);
}

// Fast structural check: does `expr` contain a Div node or Pow with negative
// exponent?  If not, the polygamma/Abramov path is skipped.
[[nodiscard]] static bool has_rational_dependency(
    ExprPtr expr, const std::string& k_name) {
    (void)k_name;
    if (!expr) return false;
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Div) return true;
        if (bin->op == BinaryOp::Pow) {
            if (const auto* exp_lit = expr_cast<IntegerLit>(bin->right);
                exp_lit && exp_lit->value.is_negative()) {
                return true;
            }
        }
        return has_rational_dependency(bin->left, k_name) ||
               has_rational_dependency(bin->right, k_name);
    }
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr t : sum->terms) {
            if (has_rational_dependency(t, k_name)) return true;
        }
        return false;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr t : prod->factors) {
            if (has_rational_dependency(t, k_name)) return true;
        }
        return false;
    }
    if (const auto* un = expr_cast<Unary>(expr)) {
        return has_rational_dependency(un->operand, k_name);
    }
    return false;
}

// Try to interpret term as A/(linear(k))^m and return the polygamma antidiff.
[[nodiscard]] static std::optional<ExprPtr> try_polygamma_antidiff(
    ExprPtr term, const Symbol& k, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    if (!has_rational_dependency(term, k.name)) return std::nullopt;

    auto together_res = algebra::together(term, ctx);
    if (together_res.is_error()) return std::nullopt;
    auto parts = algebra::apart_num_den(together_res.value(), ctx);
    if (parts.is_error()) return std::nullopt;

    auto N_simp = ctx.simplify(parts.value().numerator);
    auto D_simp = ctx.simplify(parts.value().denominator);
    if (N_simp.is_error() || D_simp.is_error()) return std::nullopt;
    ExprPtr N = N_simp.value();
    ExprPtr D = D_simp.value();

    auto N_coeffs_res = algebra::univariate_coefficients(N, k, ctx);
    if (N_coeffs_res.is_error()) return std::nullopt;
    if (N_coeffs_res.value().size() != 1U) return std::nullopt;
    ExprPtr A_const = N_coeffs_res.value()[0];

    ExprPtr base = D;
    unsigned int m = 1U;
    if (const auto* bin = expr_cast<Binary>(D); bin && bin->op == BinaryOp::Pow) {
        const auto* exp_lit = expr_cast<IntegerLit>(bin->right);
        if (!exp_lit) return std::nullopt;
        if (exp_lit->value <= BigInt(0)) return std::nullopt;
        if (exp_lit->value.bit_length() > 16U) return std::nullopt;
        m = static_cast<unsigned int>(exp_lit->value.to_u64());
        base = bin->left;
    }

    ExprPtr c0, c1;
    if (!parse_affine_in_k(base, k, c0, c1, ctx)) return std::nullopt;

    ExprPtr a = arena.make<Binary>(BinaryOp::Div, c0, c1);
    auto a_simp = ctx.simplify(a);
    if (a_simp.is_error()) return std::nullopt;
    a = a_simp.value();

    ExprPtr k_sym = arena.make<Symbol>(k);
    ExprPtr k_plus_a = arena.make<Binary>(BinaryOp::Add, k_sym, a);
    auto kpa_simp = ctx.simplify(k_plus_a);
    if (kpa_simp.is_error()) return std::nullopt;
    k_plus_a = kpa_simp.value();

    ExprPtr c1_m = (m == 1U)
        ? c1
        : arena.make<Binary>(BinaryOp::Pow, c1,
            arena.make<IntegerLit>(BigInt(static_cast<long long>(m))));
    ExprPtr A_eff_raw = arena.make<Binary>(BinaryOp::Div, A_const, c1_m);
    auto A_eff = ctx.simplify(A_eff_raw);
    if (A_eff.is_error()) return std::nullopt;

    return polygamma_antidiff(A_eff.value(), k_plus_a, m, ctx);
}

// B6-bis: RootOf-aware digamma antidifference for Q-irreducible quadratic atoms.
//
// Given atom = (B₁·k + B₀) / Q(k)  with Q = d₂·k² + d₁·k + d₀ irreducible over Q
// (m = 1 only; m > 1 returns nullopt — future extension for higher-order polygamma
// at algebraic shifts):
//
//   Let Q̃(t) = t² + p·t + q  (monic, p = d₁/d₂, q = d₀/d₂),
//       α = RootOf(Q̃, t, 0),  β = RootOf(Q̃, t, 1).
//
//   Partial fractions over Q(α):
//     (A₁k + A₀) / Q(k) = C/(k−α) + D/(k−β)
//   where  C = (A₁α + A₀)/(α−β),   D = (A₁β + A₀)/(β−α),   A_i = B_i/d₂.
//
//   Antidifference:  S(k) = C·ψ(k−α) + D·ψ(k−β).
//
// Irreducibility is guaranteed by the caller: atoms come from
// partial_fractions/factor_over_integers which returns Q-irreducible factors.
[[nodiscard]] static std::optional<ExprPtr> try_quadratic_atom_antidiff(
    ExprPtr term, const Symbol& k, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    if (!has_rational_dependency(term, k.name)) return std::nullopt;

    // Normalise to single fraction N/D.
    auto together_res = algebra::together(term, ctx);
    if (together_res.is_error()) return std::nullopt;
    auto parts = algebra::apart_num_den(together_res.value(), ctx);
    if (parts.is_error()) return std::nullopt;

    auto N_s = ctx.simplify(parts.value().numerator);
    auto D_s = ctx.simplify(parts.value().denominator);
    if (N_s.is_error() || D_s.is_error()) return std::nullopt;
    ExprPtr N = N_s.value();
    ExprPtr D = D_s.value();

    // Only handle m = 1 (plain quadratic denominator, not a power).
    if (const auto* bin = expr_cast<Binary>(D); bin && bin->op == BinaryOp::Pow) {
        return std::nullopt;
    }

    // Denominator must be exactly degree 2 in k with rational coefficients.
    auto D_coeffs_res = algebra::univariate_coefficients(D, k, ctx);
    if (D_coeffs_res.is_error() || D_coeffs_res.value().size() != 3U)
        return std::nullopt;
    const auto& dc = D_coeffs_res.value();  // [d0, d1, d2] ascending degree

    auto d2_s = ctx.simplify(dc[2]);
    auto d1_s = ctx.simplify(dc[1]);
    auto d0_s = ctx.simplify(dc[0]);
    if (d2_s.is_error() || d1_s.is_error() || d0_s.is_error()) return std::nullopt;

    Rational d2_rat, d1_rat, d0_rat;
    if (!extract_rational_const(d2_s.value(), d2_rat) || d2_rat.numerator().is_zero())
        return std::nullopt;
    if (!extract_rational_const(d1_s.value(), d1_rat)) return std::nullopt;
    if (!extract_rational_const(d0_s.value(), d0_rat)) return std::nullopt;

    // Numerator must be at most linear in k.
    auto N_coeffs_res = algebra::univariate_coefficients(N, k, ctx);
    if (N_coeffs_res.is_error() || N_coeffs_res.value().empty() ||
        N_coeffs_res.value().size() > 2U)
        return std::nullopt;
    const auto& nc = N_coeffs_res.value();

    // Monic form: p = d₁/d₂,  q = d₀/d₂.
    Rational p_rat = d1_rat / d2_rat;
    Rational q_rat = d0_rat / d2_rat;

    // Build Q̃(t) = t² + p·t + q in a fresh variable t (bound inside RootOf).
    Symbol t_sym = ctx.make_fresh_symbol("t");
    ExprPtr t_e   = arena.make<Symbol>(t_sym);
    ExprPtr t_sq  = arena.make<Binary>(BinaryOp::Pow, t_e,
                        arena.make<IntegerLit>(BigInt(2)));
    ExprPtr p_e   = rational_expr(arena, p_rat);
    ExprPtr q_e   = rational_expr(arena, q_rat);
    ExprPtr p_t   = arena.make<Binary>(BinaryOp::Mul, p_e, t_e);
    ExprPtr q_tilde_raw = arena.make<Sum>(
        std::vector<ExprPtr>{t_sq, p_t, q_e});
    auto q_tilde_s = ctx.simplify(q_tilde_raw);
    if (q_tilde_s.is_error()) return std::nullopt;
    ExprPtr q_tilde = q_tilde_s.value();

    // α = first root, β = second root of Q̃.
    ExprPtr alpha = arena.make<RootOf>(q_tilde, t_sym,
        std::optional<std::size_t>{0U});
    ExprPtr beta  = arena.make<RootOf>(q_tilde, t_sym,
        std::optional<std::size_t>{1U});

    // Effective numerator coefficients: A_i = nc[i] / d₂.
    ExprPtr d2_e  = d2_s.value();
    ExprPtr A0_raw = arena.make<Binary>(BinaryOp::Div, nc[0], d2_e);
    auto A0_s = ctx.simplify(A0_raw);
    if (A0_s.is_error()) return std::nullopt;
    ExprPtr A0 = A0_s.value();

    ExprPtr A1 = arena.make<IntegerLit>(BigInt(0));
    if (nc.size() == 2U) {
        ExprPtr A1_raw = arena.make<Binary>(BinaryOp::Div, nc[1], d2_e);
        auto A1_simp = ctx.simplify(A1_raw);
        if (A1_simp.is_error()) return std::nullopt;
        A1 = A1_simp.value();
    }

    // Residue formula:  C = (A₁α + A₀)/(α−β),  D = (A₁β + A₀)/(β−α).
    ExprPtr alpha_minus_beta = arena.make<Binary>(BinaryOp::Sub, alpha, beta);
    ExprPtr beta_minus_alpha = arena.make<Unary>(UnaryOp::Neg, alpha_minus_beta);

    ExprPtr A1_alpha     = arena.make<Binary>(BinaryOp::Mul, A1, alpha);
    ExprPtr num_at_alpha = arena.make<Binary>(BinaryOp::Add, A1_alpha, A0);
    ExprPtr C = arena.make<Binary>(BinaryOp::Div, num_at_alpha, alpha_minus_beta);

    ExprPtr A1_beta     = arena.make<Binary>(BinaryOp::Mul, A1, beta);
    ExprPtr num_at_beta = arena.make<Binary>(BinaryOp::Add, A1_beta, A0);
    ExprPtr D_coeff = arena.make<Binary>(BinaryOp::Div, num_at_beta, beta_minus_alpha);

    // Antidifference: S(k) = C·ψ(k−α) + D·ψ(k−β).
    ExprPtr k_sym_e      = arena.make<Symbol>(k);
    ExprPtr k_minus_alpha = arena.make<Binary>(BinaryOp::Sub, k_sym_e, alpha);
    ExprPtr k_minus_beta  = arena.make<Binary>(BinaryOp::Sub, k_sym_e, beta);

    ExprPtr S1 = polygamma_antidiff(C,       k_minus_alpha, 1U, ctx);
    ExprPtr S2 = polygamma_antidiff(D_coeff, k_minus_beta,  1U, ctx);
    return arena.make<Binary>(BinaryOp::Add, S1, S2);
}

// ── public functions (declared in summation_internal.hpp) ─────────────────────

Result<ExprPtr> try_gosper_definite(
    ExprPtr term, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto antidiff = symbolic::gosper_sum(term, var, ctx);
    if (antidiff.is_error()) return fail<ExprPtr>(antidiff.error());
    if (!antidiff.value().has_value()) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Gosper: term is not Gosper-summable",
        });
    }
    AstArena& arena = ctx.arena();
    ExprPtr S = antidiff.value().value();
    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));
    auto S_upper = ctx.substitute(S, var, upper_plus_one);
    if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
    auto S_lower = ctx.substitute(S, var, lower);
    if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
    ExprPtr diff = arena.make<Binary>(BinaryOp::Sub,
        S_upper.value(), S_lower.value());
    return ctx.simplify(diff);
}

Result<ExprPtr> try_polygamma_definite(
    ExprPtr term, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    auto antidiff_opt = try_polygamma_antidiff(term, var, ctx);
    if (!antidiff_opt.has_value()) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Polygamma summation: term is not of the form "
                       "A/(linear(k))^m",
        });
    }
    AstArena& arena = ctx.arena();
    ExprPtr S = antidiff_opt.value();
    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));
    auto S_upper = ctx.substitute(S, var, upper_plus_one);
    if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
    auto S_lower = ctx.substitute(S, var, lower);
    if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
    ExprPtr diff = arena.make<Binary>(BinaryOp::Sub,
        S_upper.value(), S_lower.value());
    return ctx.simplify(diff);
}

Result<ExprPtr> try_abramov_definite(
    ExprPtr term, const Symbol& var,
    ExprPtr lower, ExprPtr upper,
    symbolic::CASContext& ctx) {
    if (!has_rational_dependency(term, var.name)) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Abramov: term has no rational dependency on the "
                       "summation variable",
        });
    }

    auto together_res = algebra::together(term, ctx);
    if (together_res.is_error()) return fail<ExprPtr>(together_res.error());
    auto parts = algebra::apart_num_den(together_res.value(), ctx);
    if (parts.is_error()) return fail<ExprPtr>(parts.error());

    auto deg_N_res = algebra::polynomial_degree(parts.value().numerator, var, ctx);
    auto deg_D_res = algebra::polynomial_degree(parts.value().denominator, var, ctx);
    if (deg_N_res.is_error() || deg_D_res.is_error()) {
        return fail<ExprPtr>(CASError{
            .kind    = CASErrorKind::Unimplemented,
            .message = "Abramov: rational term must be polynomial-over-"
                       "polynomial in the summation variable",
        });
    }
    const std::size_t deg_N = deg_N_res.value();
    const std::size_t deg_D = deg_D_res.value();
    AstArena& arena = ctx.arena();
    ExprPtr N = parts.value().numerator;
    ExprPtr D = parts.value().denominator;

    ExprPtr polynomial_part  = arena.make<IntegerLit>(BigInt(0));
    ExprPtr proper_remainder = term;

    if (deg_N >= deg_D) {
        auto divmod = algebra::polynomial_divmod(N, D, var, ctx);
        if (divmod.is_error()) return fail<ExprPtr>(divmod.error());
        polynomial_part  = divmod.value().quotient;
        proper_remainder = arena.make<Binary>(BinaryOp::Div,
            divmod.value().remainder, D);
        auto simp_r = ctx.simplify(proper_remainder);
        if (simp_r.is_error()) return fail<ExprPtr>(simp_r.error());
        proper_remainder = simp_r.value();
    }

    auto fractions = algebra::partial_fractions(proper_remainder, var, ctx);
    if (fractions.is_error()) return fail<ExprPtr>(fractions.error());

    ExprPtr upper_plus_one = arena.make<Binary>(BinaryOp::Add, upper,
        arena.make<IntegerLit>(BigInt(1)));

    std::vector<ExprPtr> total_terms;

    bool poly_part_is_zero = false;
    if (const auto* il = expr_cast<IntegerLit>(polynomial_part))
        poly_part_is_zero = il->value.is_zero();
    if (!poly_part_is_zero) {
        auto poly_def = try_gosper_definite(polynomial_part, var, lower, upper, ctx);
        if (poly_def.is_ok()) {
            total_terms.push_back(poly_def.value());
        } else {
            return poly_def;
        }
    }

    for (ExprPtr atom : fractions.value()) {
        auto atom_simp = ctx.simplify(atom);
        if (atom_simp.is_error()) return fail<ExprPtr>(atom_simp.error());

        // Q-linear atoms: A/(c₁k+c₀)^m → polygamma antidifference.
        auto antidiff_opt = try_polygamma_antidiff(atom_simp.value(), var, ctx);

        // Q-irreducible quadratic atoms: (B₁k+B₀)/Q(k), m=1 → RootOf digamma.
        if (!antidiff_opt.has_value())
            antidiff_opt = try_quadratic_atom_antidiff(atom_simp.value(), var, ctx);

        if (!antidiff_opt.has_value()) {
            return fail<ExprPtr>(CASError{
                .kind    = CASErrorKind::Unimplemented,
                .message = "Abramov: partial-fraction atom not of supported form. "
                           "Handled: Q-linear A/(ck+d)^m and Q-irreducible "
                           "quadratic (B₁k+B₀)/Q(k) with m=1. "
                           "Unhandled: Q-irreducible deg≥3 or quadratic^m (m>1).",
            });
        }

        ExprPtr S = antidiff_opt.value();
        auto S_upper = ctx.substitute(S, var, upper_plus_one);
        if (S_upper.is_error()) return fail<ExprPtr>(S_upper.error());
        auto S_lower = ctx.substitute(S, var, lower);
        if (S_lower.is_error()) return fail<ExprPtr>(S_lower.error());
        total_terms.push_back(
            arena.make<Binary>(BinaryOp::Sub,
                S_upper.value(), S_lower.value()));
    }

    if (total_terms.empty())
        return ok(arena.make<IntegerLit>(BigInt(0)));

    ExprPtr combined = (total_terms.size() == 1U)
        ? total_terms.front()
        : arena.make<Sum>(std::move(total_terms));
    return ctx.simplify(combined);
}

} // namespace cas::calculus
