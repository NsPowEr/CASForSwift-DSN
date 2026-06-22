#include "summation_abramov_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas::calculus::abramov_detail {

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

std::optional<ExprPtr> try_polygamma_antidiff(
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

[[nodiscard]] static Rational compute_hermite_factor(unsigned int m, unsigned int r) {
    if (r == 0U) return Rational(1);
    
    auto fact = [](unsigned int n) -> BigInt {
        BigInt res(1);
        for (unsigned int i = 2U; i <= n; ++i) {
            res *= BigInt(static_cast<long long>(i));
        }
        return res;
    };
    
    BigInt num = fact(m + r - 2U);
    BigInt den = fact(r) * fact(m - 1U);
    return Rational(num, den);
}

std::optional<ExprPtr> try_quadratic_atom_antidiff(
    ExprPtr term, const Symbol& k, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    if (!has_rational_dependency(term, k.name)) return std::nullopt;

    auto together_res = algebra::together(term, ctx);
    if (together_res.is_error()) return std::nullopt;
    auto parts = algebra::apart_num_den(together_res.value(), ctx);
    if (parts.is_error()) return std::nullopt;

    auto N_s = ctx.simplify(parts.value().numerator);
    auto D_s = ctx.simplify(parts.value().denominator);
    if (N_s.is_error() || D_s.is_error()) return std::nullopt;
    ExprPtr N = N_s.value();
    ExprPtr D = D_s.value();

    unsigned int m = 1U;
    ExprPtr base = D;
    if (const auto* bin = expr_cast<Binary>(D); bin && bin->op == BinaryOp::Pow) {
        const auto* exp_lit = expr_cast<IntegerLit>(bin->right);
        if (!exp_lit || exp_lit->value.is_negative() || exp_lit->value.bit_length() > 64U) return std::nullopt;
        m = static_cast<unsigned int>(exp_lit->value.to_u64());
        base = bin->left;
    }

    auto D_coeffs_res = algebra::univariate_coefficients(base, k, ctx);
    if (D_coeffs_res.is_error() || D_coeffs_res.value().size() != 3U)
        return std::nullopt;
    const auto& dc = D_coeffs_res.value();

    auto d2_s = ctx.simplify(dc[2]);
    auto d1_s = ctx.simplify(dc[1]);
    auto d0_s = ctx.simplify(dc[0]);
    if (d2_s.is_error() || d1_s.is_error() || d0_s.is_error()) return std::nullopt;

    Rational d2_rat, d1_rat, d0_rat;
    if (!extract_rational_const(d2_s.value(), d2_rat) || d2_rat.numerator().is_zero())
        return std::nullopt;
    if (!extract_rational_const(d1_s.value(), d1_rat)) return std::nullopt;
    if (!extract_rational_const(d0_s.value(), d0_rat)) return std::nullopt;

    auto N_coeffs_res = algebra::univariate_coefficients(N, k, ctx);
    if (N_coeffs_res.is_error() || N_coeffs_res.value().empty() ||
        N_coeffs_res.value().size() > 2U)
        return std::nullopt;
    const auto& nc = N_coeffs_res.value();

    Rational p_rat = d1_rat / d2_rat;
    Rational q_rat = d0_rat / d2_rat;

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

    ExprPtr alpha = arena.make<RootOf>(q_tilde, t_sym,
        std::optional<std::size_t>{0U});
    ExprPtr beta  = arena.make<RootOf>(q_tilde, t_sym,
        std::optional<std::size_t>{1U});

    ExprPtr d2_e  = d2_s.value();
    ExprPtr d2_m = (m == 1U)
        ? d2_e
        : arena.make<Binary>(BinaryOp::Pow, d2_e,
            arena.make<IntegerLit>(BigInt(static_cast<long long>(m))));
    auto d2_m_s = ctx.simplify(d2_m);
    if (d2_m_s.is_error()) return std::nullopt;
    ExprPtr d2_m_eff = d2_m_s.value();

    ExprPtr A0_raw = arena.make<Binary>(BinaryOp::Div, nc[0], d2_m_eff);
    auto A0_s = ctx.simplify(A0_raw);
    if (A0_s.is_error()) return std::nullopt;
    ExprPtr A0 = A0_s.value();

    ExprPtr A1 = arena.make<IntegerLit>(BigInt(0));
    if (nc.size() == 2U) {
        ExprPtr A1_raw = arena.make<Binary>(BinaryOp::Div, nc[1], d2_m_eff);
        auto A1_simp = ctx.simplify(A1_raw);
        if (A1_simp.is_error()) return std::nullopt;
        A1 = A1_simp.value();
    }

    std::vector<ExprPtr> S_terms;
    S_terms.reserve(2U * m);

    ExprPtr alpha_minus_beta = arena.make<Binary>(BinaryOp::Sub, alpha, beta);
    ExprPtr beta_minus_alpha = arena.make<Unary>(UnaryOp::Neg, alpha_minus_beta);

    for (unsigned int j = 1U; j <= m; ++j) {
        unsigned int r = m - j;
        ExprPtr C_coeff, D_coeff;
        if (r == 0U) {
            ExprPtr denom_C = (m == 1U)
                ? alpha_minus_beta
                : arena.make<Binary>(BinaryOp::Pow, alpha_minus_beta,
                    arena.make<IntegerLit>(BigInt(static_cast<long long>(m))));
            ExprPtr denom_D = (m == 1U)
                ? beta_minus_alpha
                : arena.make<Binary>(BinaryOp::Pow, beta_minus_alpha,
                    arena.make<IntegerLit>(BigInt(static_cast<long long>(m))));

            ExprPtr A1_alpha = arena.make<Binary>(BinaryOp::Mul, A1, alpha);
            ExprPtr num_C = arena.make<Binary>(BinaryOp::Add, A1_alpha, A0);
            C_coeff = arena.make<Binary>(BinaryOp::Div, num_C, denom_C);

            ExprPtr A1_beta = arena.make<Binary>(BinaryOp::Mul, A1, beta);
            ExprPtr num_D = arena.make<Binary>(BinaryOp::Add, A1_beta, A0);
            D_coeff = arena.make<Binary>(BinaryOp::Div, num_D, denom_D);
        } else {
            Rational factor_rat = compute_hermite_factor(m, r);
            if (r % 2U == 1U) {
                factor_rat = -factor_rat;
            }
            ExprPtr scale = rational_expr(arena, factor_rat);

            ExprPtr denom_C = arena.make<Binary>(BinaryOp::Pow, alpha_minus_beta,
                arena.make<IntegerLit>(BigInt(static_cast<long long>(m + r))));
            ExprPtr denom_D = arena.make<Binary>(BinaryOp::Pow, beta_minus_alpha,
                arena.make<IntegerLit>(BigInt(static_cast<long long>(m + r))));

            ExprPtr A1_alpha = arena.make<Binary>(BinaryOp::Mul, A1, alpha);
            ExprPtr A1_alpha_plus_A0 = arena.make<Binary>(BinaryOp::Add, A1_alpha, A0);
            ExprPtr m_plus_r_minus_1_expr = arena.make<IntegerLit>(BigInt(static_cast<long long>(m + r - 1U)));
            ExprPtr term_C_1 = arena.make<Binary>(BinaryOp::Mul, m_plus_r_minus_1_expr, A1_alpha_plus_A0);
            ExprPtr r_A1 = arena.make<Binary>(BinaryOp::Mul, arena.make<IntegerLit>(BigInt(static_cast<long long>(r))), A1);
            ExprPtr term_C_2 = arena.make<Binary>(BinaryOp::Mul, r_A1, alpha_minus_beta);
            ExprPtr expr_C = arena.make<Binary>(BinaryOp::Sub, term_C_1, term_C_2);
            ExprPtr scale_expr_C = arena.make<Binary>(BinaryOp::Mul, scale, expr_C);
            C_coeff = arena.make<Binary>(BinaryOp::Div, scale_expr_C, denom_C);

            ExprPtr A1_beta = arena.make<Binary>(BinaryOp::Mul, A1, beta);
            ExprPtr A1_beta_plus_A0 = arena.make<Binary>(BinaryOp::Add, A1_beta, A0);
            ExprPtr term_D_1 = arena.make<Binary>(BinaryOp::Mul, m_plus_r_minus_1_expr, A1_beta_plus_A0);
            ExprPtr term_D_2 = arena.make<Binary>(BinaryOp::Mul, r_A1, beta_minus_alpha);
            ExprPtr expr_D = arena.make<Binary>(BinaryOp::Sub, term_D_1, term_D_2);
            ExprPtr scale_expr_D = arena.make<Binary>(BinaryOp::Mul, scale, expr_D);
            D_coeff = arena.make<Binary>(BinaryOp::Div, scale_expr_D, denom_D);
        }

        ExprPtr k_sym_e = arena.make<Symbol>(k);
        ExprPtr k_minus_alpha = arena.make<Binary>(BinaryOp::Sub, k_sym_e, alpha);
        ExprPtr k_minus_beta = arena.make<Binary>(BinaryOp::Sub, k_sym_e, beta);

        ExprPtr S1 = polygamma_antidiff(C_coeff, k_minus_alpha, j, ctx);
        ExprPtr S2 = polygamma_antidiff(D_coeff, k_minus_beta,  j, ctx);
        S_terms.push_back(S1);
        S_terms.push_back(S2);
    }

    ExprPtr combined = (S_terms.size() == 1U)
        ? S_terms.front()
        : arena.make<Sum>(std::move(S_terms));
    return combined;
}

} // namespace cas::calculus::abramov_detail
