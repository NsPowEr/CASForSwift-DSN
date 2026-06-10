// F7.5.E2: split from simplify_bessel_orthogonal.cpp to honor the 500-line
// anti-monolith limit. Hosts canonical recurrences for orthogonal polynomial
// families — Chebyshev T/U, Hermite H/He, Legendre P, Laguerre L, Jacobi P^(α,β)
// — plus the LambertW principal-branch helper.

#include "simplify_impl.hpp"
#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_funcall_orthogonal_polys(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before) {
    (void)target_before;

    // Chebyshev T_n — three-term recurrence, scales to any n
    if (op == BuiltinOp::ChebyshevT && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && !il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_orthogonal_polys",
                    "ChebyshevT with degree > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Cap degree via ctx.max_polynomial_degree or request numeric evaluation",
                    "F1.x");
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            if (n == 1U) return ok(x);
            ExprPtr t_prev = make_integer(arena_, BigInt(1));
            ExprPtr t_curr = x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr two_x_tk = arena_.make<Product>(std::vector<ExprPtr>{
                    make_integer(arena_, BigInt(2)), x, t_curr});
                ExprPtr neg_prev = arena_.make<Unary>(UnaryOp::Neg, t_prev);
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{two_x_tk, neg_prev});
                auto s = simplify_expr(next);
                if (s.is_error()) return s;
                t_prev = t_curr;
                t_curr = s.value();
            }
            return ok(t_curr);
        }
    }

    // Chebyshev U_n — three-term recurrence
    if (op == BuiltinOp::ChebyshevU && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && !il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_orthogonal_polys",
                    "ChebyshevU with degree > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Cap degree via ctx.max_polynomial_degree or request numeric evaluation",
                    "F1.x");
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            ExprPtr two_x = arena_.make<Product>(std::vector<ExprPtr>{
                make_integer(arena_, BigInt(2)), x});
            if (n == 1U) return simplify_expr(two_x);
            ExprPtr u_prev = make_integer(arena_, BigInt(1));
            ExprPtr u_curr = two_x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr two_x_uk = arena_.make<Product>(std::vector<ExprPtr>{
                    make_integer(arena_, BigInt(2)), x, u_curr});
                ExprPtr neg_prev = arena_.make<Unary>(UnaryOp::Neg, u_prev);
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{two_x_uk, neg_prev});
                auto s = simplify_expr(next);
                if (s.is_error()) return s;
                u_prev = u_curr;
                u_curr = s.value();
            }
            return ok(u_curr);
        }
    }

    // Hermite H_n (physicist) — three-term recurrence
    if (op == BuiltinOp::HermiteH && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && !il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_orthogonal_polys",
                    "HermiteH (physicist) with degree > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Cap degree via ctx.max_polynomial_degree or request numeric evaluation",
                    "F1.x");
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            ExprPtr two_x = arena_.make<Product>(std::vector<ExprPtr>{
                make_integer(arena_, BigInt(2)), x});
            if (n == 1U) return simplify_expr(two_x);
            ExprPtr h_prev = make_integer(arena_, BigInt(1));
            ExprPtr h_curr = two_x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr two_x_hk = arena_.make<Product>(std::vector<ExprPtr>{
                    make_integer(arena_, BigInt(2)), x, h_curr});
                ExprPtr two_k = make_integer(arena_, BigInt(static_cast<std::int64_t>(2U * k)));
                ExprPtr neg = arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Product>(std::vector<ExprPtr>{two_k, h_prev}));
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{two_x_hk, neg});
                auto s = simplify_expr(next);
                if (s.is_error()) return s;
                h_prev = h_curr;
                h_curr = s.value();
            }
            return ok(h_curr);
        }
    }

    // Hermite He_n (probabilist) — three-term recurrence
    if (op == BuiltinOp::HermiteHe && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && !il->value.is_negative()) {
            if (il->value.bit_length() > 16)
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_orthogonal_polys",
                    "HermiteHe (probabilist) with degree > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Cap degree via ctx.max_polynomial_degree or request numeric evaluation",
                    "F1.x");
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            if (n == 1U) return ok(x);
            ExprPtr h_prev = make_integer(arena_, BigInt(1));
            ExprPtr h_curr = x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr x_hk = arena_.make<Product>(std::vector<ExprPtr>{x, h_curr});
                ExprPtr k_const = make_integer(arena_, BigInt(static_cast<std::int64_t>(k)));
                ExprPtr neg = arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Product>(std::vector<ExprPtr>{k_const, h_prev}));
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{x_hk, neg});
                auto s = simplify_expr(next);
                if (s.is_error()) return s;
                h_prev = h_curr;
                h_curr = s.value();
            }
            return ok(h_curr);
        }
    }

    // LambertW(z) — principal branch W₀ (moved here to consolidate single-arg
    // special functions away from the Bessel file under anti-monolith).
    //   W(0) = 0,  W(e) = 1,  W(x · exp(x)) = x  per x ≥ 0.
    if (op == BuiltinOp::LambertW && args.size() == 1U) {
        ExprPtr arg = args[0];
        if (is_zero_expr(arg))
            return ok(make_integer(arena_, BigInt(0)));
        if (is_constant_expr(arg, MathConstant::E))
            return ok(make_integer(arena_, BigInt(1)));
        if (const auto* prod = expr_cast<Product>(arg); prod && prod->factors.size() == 2U) {
            ExprPtr a = prod->factors[0];
            ExprPtr b = prod->factors[1];
            auto match = [&](ExprPtr x_cand, ExprPtr exp_cand) -> ExprPtr {
                const auto* fc = expr_cast<FuncCall>(exp_cand);
                if (fc && fc->func_id == BuiltinOp::Exp && fc->args.size() == 1U) {
                    if (structural_equal(fc->args[0], x_cand)
                        && is_known_nonnegative(x_cand)) {
                        return x_cand;
                    }
                }
                return nullptr;
            };
            if (ExprPtr m = match(a, b)) return ok(m);
            if (ExprPtr m = match(b, a)) return ok(m);
        }
    }

    // Legendre P_n — Bonnet recurrence
    if (op == BuiltinOp::LegendreP && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && il->value >= BigInt(0)) {
            if (il->value.bit_length() > 16)
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_orthogonal_polys",
                    "LegendreP with degree > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Cap degree via ctx.max_polynomial_degree or use asymptotic Legendre formula",
                    "F1.x");
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            if (n == 1U) return ok(x);
            ExprPtr p_prev = make_integer(arena_, BigInt(1));
            ExprPtr p_curr = x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                Rational a(BigInt(static_cast<std::int64_t>(2U * k + 1U)),
                           BigInt(static_cast<std::int64_t>(k + 1U)));
                Rational b(BigInt(static_cast<std::int64_t>(k)),
                           BigInt(static_cast<std::int64_t>(k + 1U)));
                ExprPtr term1 = arena_.make<Product>(std::vector<ExprPtr>{
                    make_rational(arena_, a), x, p_curr});
                ExprPtr neg_term2 = arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Product>(std::vector<ExprPtr>{make_rational(arena_, b), p_prev}));
                ExprPtr next = arena_.make<Sum>(std::vector<ExprPtr>{term1, neg_term2});
                auto simp = simplify_expr(next);
                if (simp.is_error()) return simp;
                p_prev = p_curr;
                p_curr = simp.value();
            }
            return ok(p_curr);
        }
    }

    // Laguerre L_n(x) — classical recurrence (Abramowitz-Stegun 22.7.12).
    if (op == BuiltinOp::LaguerreL && args.size() == 2U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && il->value >= BigInt(0)) {
            if (il->value.bit_length() > 16)
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_orthogonal_polys",
                    "LaguerreL with degree > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Cap degree via ctx.max_polynomial_degree or request numeric evaluation",
                    "F1.x");
            const std::uint64_t n = il->value.to_u64();
            ExprPtr x = args[1];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            ExprPtr one_minus_x = arena_.make<Binary>(BinaryOp::Sub,
                make_integer(arena_, BigInt(1)), x);
            if (n == 1U) return ok(one_minus_x);
            ExprPtr l_prev = make_integer(arena_, BigInt(1));
            ExprPtr l_curr = one_minus_x;
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr two_k_plus_1 = make_integer(arena_,
                    BigInt(static_cast<std::int64_t>(2*k + 1)));
                ExprPtr two_k_plus_1_minus_x = arena_.make<Binary>(BinaryOp::Sub,
                    two_k_plus_1, x);
                ExprPtr term1 = arena_.make<Product>(std::vector<ExprPtr>{
                    two_k_plus_1_minus_x, l_curr});
                ExprPtr k_e = make_integer(arena_, BigInt(static_cast<std::int64_t>(k)));
                ExprPtr term2 = arena_.make<Product>(std::vector<ExprPtr>{k_e, l_prev});
                ExprPtr numer = arena_.make<Binary>(BinaryOp::Sub, term1, term2);
                ExprPtr k_plus_1 = make_integer(arena_,
                    BigInt(static_cast<std::int64_t>(k + 1)));
                ExprPtr next = arena_.make<Binary>(BinaryOp::Div, numer, k_plus_1);
                auto simp = simplify_expr(next);
                if (simp.is_error()) return simp;
                l_prev = l_curr;
                l_curr = simp.value();
            }
            return ok(l_curr);
        }
    }

    // Jacobi P_n^{(α,β)}(x) — Bonnet recurrence (Abramowitz-Stegun 22.7.1).
    if (op == BuiltinOp::JacobiP && args.size() == 4U) {
        if (const auto* il = expr_cast<IntegerLit>(args[0]);
            il != nullptr && il->value >= BigInt(0)) {
            if (il->value.bit_length() > 16)
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_orthogonal_polys",
                    "JacobiP P_n^(α,β) with degree > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Cap degree via ctx.max_polynomial_degree or request numeric evaluation",
                    "F1.x");
            const std::uint64_t n = il->value.to_u64();
            ExprPtr alpha = args[1];
            ExprPtr beta = args[2];
            ExprPtr x = args[3];
            if (n == 0U) return ok(make_integer(arena_, BigInt(1)));
            ExprPtr alpha_minus_beta = arena_.make<Binary>(BinaryOp::Sub, alpha, beta);
            ExprPtr alpha_plus_beta_plus_2 = arena_.make<Sum>(std::vector<ExprPtr>{
                alpha, beta, make_integer(arena_, BigInt(2))});
            ExprPtr half = make_rational(arena_, Rational(BigInt(1), BigInt(2)));
            ExprPtr p1_const = arena_.make<Product>(std::vector<ExprPtr>{half, alpha_minus_beta});
            ExprPtr p1_xterm = arena_.make<Product>(std::vector<ExprPtr>{
                half, alpha_plus_beta_plus_2, x});
            ExprPtr p1_raw = arena_.make<Sum>(std::vector<ExprPtr>{p1_const, p1_xterm});
            auto p1_simp = simplify_expr(p1_raw);
            if (p1_simp.is_error()) return p1_simp;
            if (n == 1U) return p1_simp;
            ExprPtr p_prev = make_integer(arena_, BigInt(1));
            ExprPtr p_curr = p1_simp.value();
            for (std::uint64_t k = 1U; k < n; ++k) {
                ExprPtr k_e = make_integer(arena_, BigInt(static_cast<std::int64_t>(k)));
                ExprPtr two = make_integer(arena_, BigInt(2));
                ExprPtr k_plus_1 = make_integer(arena_, BigInt(static_cast<std::int64_t>(k + 1)));
                ExprPtr k_plus_alpha = arena_.make<Sum>(std::vector<ExprPtr>{k_e, alpha});
                ExprPtr k_plus_beta = arena_.make<Sum>(std::vector<ExprPtr>{k_e, beta});
                ExprPtr k_plus_alpha_beta = arena_.make<Sum>(std::vector<ExprPtr>{
                    k_e, alpha, beta});
                ExprPtr k_plus_alpha_beta_plus_1 = arena_.make<Sum>(std::vector<ExprPtr>{
                    k_e, alpha, beta, make_integer(arena_, BigInt(1))});
                ExprPtr two_k_alpha_beta = arena_.make<Product>(std::vector<ExprPtr>{
                    two, k_plus_alpha_beta});
                ExprPtr two_k_alpha_beta_plus_1 = arena_.make<Sum>(std::vector<ExprPtr>{
                    two_k_alpha_beta, make_integer(arena_, BigInt(1))});
                ExprPtr two_k_alpha_beta_plus_2 = arena_.make<Sum>(std::vector<ExprPtr>{
                    two_k_alpha_beta, two});
                ExprPtr a_k = arena_.make<Product>(std::vector<ExprPtr>{
                    two, k_plus_1, k_plus_alpha_beta_plus_1, two_k_alpha_beta});
                ExprPtr alpha_sq = arena_.make<Binary>(BinaryOp::Pow, alpha, two);
                ExprPtr beta_sq = arena_.make<Binary>(BinaryOp::Pow, beta, two);
                ExprPtr alpha_sq_minus_beta_sq = arena_.make<Binary>(BinaryOp::Sub,
                    alpha_sq, beta_sq);
                ExprPtr b_k = arena_.make<Product>(std::vector<ExprPtr>{
                    two_k_alpha_beta_plus_1, alpha_sq_minus_beta_sq});
                ExprPtr c_k = arena_.make<Product>(std::vector<ExprPtr>{
                    two_k_alpha_beta, two_k_alpha_beta_plus_1, two_k_alpha_beta_plus_2});
                ExprPtr d_k = arena_.make<Product>(std::vector<ExprPtr>{
                    two, k_plus_alpha, k_plus_beta, two_k_alpha_beta_plus_2});
                ExprPtr num_term1 = arena_.make<Sum>(std::vector<ExprPtr>{
                    b_k, arena_.make<Product>(std::vector<ExprPtr>{c_k, x})});
                ExprPtr factor1 = arena_.make<Product>(std::vector<ExprPtr>{
                    num_term1, p_curr});
                ExprPtr factor2 = arena_.make<Product>(std::vector<ExprPtr>{d_k, p_prev});
                ExprPtr numer = arena_.make<Binary>(BinaryOp::Sub, factor1, factor2);
                ExprPtr p_next = arena_.make<Binary>(BinaryOp::Div, numer, a_k);
                ExprPtr norm = p_next;
                if (context_ != nullptr) {
                    auto t = algebra::together(p_next, *context_);
                    if (t.is_ok()) norm = t.value();
                }
                auto simp = simplify_expr(norm);
                if (simp.is_error()) return simp;
                p_prev = p_curr;
                p_curr = simp.value();
            }
            return ok(p_curr);
        }
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

}  // namespace cas::symbolic::detail
