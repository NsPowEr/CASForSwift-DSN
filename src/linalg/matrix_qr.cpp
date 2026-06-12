// CAS-F4.1b — QR decomposition via Householder Reflectors (symbolic).
//
// Decomposes A (m×n, m ≥ n, full column rank) into Q (m×n) e R (n×n)
// triangolare superiore tali che A = Q·R e Q^T·Q = I.
//
// Algoritmo (Householder puramente simbolico e razionalizzato):
//   Per risolvere l'instabilità dell'ortogonalizzazione Gram-Schmidt (MGS),
//   viene utilizzato il riflettore di Householder H_k. Per prevenire
//   l'esplosione dell'AST e il timeout del simplifier, le equazioni
//   sono state razionalizzate analiticamente. Il denominatore standard
//   di Householder contiene radici quadrate, ma espandendo (I - 2 v v^T / N_v)
//   con v = x + alpha*e_1, si dimostra che l'aggiornamento per un vettore y
//   può essere scritto come:
//     y_0 = - (x^T y / N_x) * alpha
//     y_i = (y_i - A x_i) - B x_i * alpha    per i > 0
//   dove A e B sono funzioni PURAMENTE RAZIONALI di x e y.
//   Questo garantisce che le radici quadrate appaiano solo nei numeratori.

#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace cas::linalg {

Result<QRDecomposition> qr_decompose(const MatrixExpr& matrix,
                                       symbolic::CASContext& ctx) {
    const std::size_t m = matrix.rows();
    const std::size_t n = matrix.cols();
    if (m == 0U || n == 0U) {
        return fail<QRDecomposition>(CASError{
            CASErrorKind::InvalidArgument, "qr_decompose: empty matrix", std::nullopt});
    }

    MatrixExpr R = matrix; 

    struct Reflector {
        std::vector<ExprPtr> x;
        ExprPtr Nx;
        ExprPtr S;
        ExprPtr alpha;
    };
    std::vector<Reflector> reflectors;
    reflectors.resize(n);

    auto apply_H = [&](std::vector<ExprPtr>& y, const Reflector& ref) -> Result<void> {
        const auto& x = ref.x;
        ExprPtr x0 = x[0];
        ExprPtr y0 = y[0];

        ExprPtr dot_xy = integer(ctx, 0);
        for(std::size_t i=0; i<x.size(); ++i) {
            auto term = mul_expr(ctx, x[i], y[i]);
            if (term.is_error()) return fail<void>(term.error());
            auto sum = add_expr(ctx, dot_xy, term.value());
            if (sum.is_error()) return fail<void>(sum.error());
            dot_xy = sum.value();
        }
        
        auto x0_y0 = mul_expr(ctx, x0, y0);
        if (x0_y0.is_error()) return fail<void>(x0_y0.error());
        auto dot_rest = sub_expr(ctx, dot_xy, x0_y0.value());
        if (dot_rest.is_error()) return fail<void>(dot_rest.error());

        auto A_res = div_expr(ctx, dot_rest.value(), ref.S);
        if (A_res.is_error()) return fail<void>(A_res.error());
        ExprPtr A = A_res.value();

        auto Nx_y0 = mul_expr(ctx, ref.Nx, y0);
        if (Nx_y0.is_error()) return fail<void>(Nx_y0.error());
        auto x0_dot = mul_expr(ctx, x0, dot_xy);
        if (x0_dot.is_error()) return fail<void>(x0_dot.error());
        auto B_num = sub_expr(ctx, Nx_y0.value(), x0_dot.value());
        if (B_num.is_error()) return fail<void>(B_num.error());
        auto Nx_S = mul_expr(ctx, ref.Nx, ref.S);
        if (Nx_S.is_error()) return fail<void>(Nx_S.error());
        auto B_res = div_expr(ctx, B_num.value(), Nx_S.value());
        if (B_res.is_error()) return fail<void>(B_res.error());
        ExprPtr B = B_res.value();

        auto dot_div_Nx = div_expr(ctx, dot_xy, ref.Nx);
        if (dot_div_Nx.is_error()) return fail<void>(dot_div_Nx.error());
        auto neg_dot_div_Nx = negate_expr(ctx, dot_div_Nx.value());
        if (neg_dot_div_Nx.is_error()) return fail<void>(neg_dot_div_Nx.error());
        auto y0_res = mul_expr(ctx, neg_dot_div_Nx.value(), ref.alpha);
        if (y0_res.is_error()) return fail<void>(y0_res.error());
        y[0] = y0_res.value();

        for(std::size_t i=1; i<x.size(); ++i) {
            auto A_xi = mul_expr(ctx, A, x[i]);
            if (A_xi.is_error()) return fail<void>(A_xi.error());
            auto term1 = sub_expr(ctx, y[i], A_xi.value());
            if (term1.is_error()) return fail<void>(term1.error());
            
            auto B_xi = mul_expr(ctx, B, x[i]);
            if (B_xi.is_error()) return fail<void>(B_xi.error());
            auto term2 = mul_expr(ctx, B_xi.value(), ref.alpha);
            if (term2.is_error()) return fail<void>(term2.error());
            
            auto yi_res = sub_expr(ctx, term1.value(), term2.value());
            if (yi_res.is_error()) return fail<void>(yi_res.error());
            y[i] = yi_res.value();
        }
        return ok();
    };

    for (std::size_t k = 0; k < n; ++k) {
        std::vector<ExprPtr> x(m - k);
        for (std::size_t i = k; i < m; ++i) {
            x[i - k] = R(i, k);
        }

        bool already_zero = true;
        for (std::size_t i = 1; i < m - k; ++i) {
            if (!is_zero_expr(x[i])) {
                already_zero = false;
                break;
            }
        }
        if (already_zero) {
            if (is_zero_expr(x[0])) {
                return make_unimplemented<QRDecomposition>(
                    "linalg", "qr_decompose", "matrix column is linearly dependent",
                    error::reason_codes::LINALG_LINEAR_DEPENDENT,
                    "Matrix must have full column rank for this QR implementation",
                    "F4.1b");
            }
            reflectors[k] = Reflector{x, integer(ctx, 0), integer(ctx, 0), integer(ctx, 0)};
            continue;
        }

        auto norm_sq_res = sym_norm_sq(x, ctx);
        if (norm_sq_res.is_error()) return fail<QRDecomposition>(norm_sq_res.error());
        ExprPtr Nx = norm_sq_res.value();

        if (is_zero_expr(Nx)) {
            return make_unimplemented<QRDecomposition>(
                "linalg", "qr_decompose", "matrix column is linearly dependent",
                error::reason_codes::LINALG_LINEAR_DEPENDENT,
                "Matrix must have full column rank for this QR implementation",
                "F4.1b");
        }

        auto x0_sq = mul_expr(ctx, x[0], x[0]);
        if (x0_sq.is_error()) return fail<QRDecomposition>(x0_sq.error());
        auto S_res = sub_expr(ctx, Nx, x0_sq.value());
        if (S_res.is_error()) return fail<QRDecomposition>(S_res.error());
        ExprPtr S = S_res.value();

        ExprPtr alpha = ctx.arena().make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{Nx});
        auto sim_alpha = simplify(ctx, alpha);
        if (sim_alpha.is_ok()) alpha = sim_alpha.value();

        // Se alpha non si semplifica a un numero e contiene simboli, le espressioni generate
        // saranno corrette ma il simplifier attualmente non riesce a dimostrare Q*R == A 
        // (espansione di binomi con Sqrt simboliche). Restituiamo Unimplemented come faceva MGS.
        if (alpha->kind == ExprKind::FuncCall) {
            auto call = expr_cast<FuncCall>(alpha);
            if (call->func_id == BuiltinOp::Sqrt) {
                if (estimate_complexity(Nx) > 2 && !is_one_expr(Nx) && !is_zero_expr(Nx)) {
                    // Controlliamo se ci sono simboli. Un approccio semplice è vedere la complessità.
                    if (Nx->kind != ExprKind::IntegerLit && Nx->kind != ExprKind::RationalLit) {
                        // Per le matrici numeriche come sqrt(5), la complessità è bassa (1 o 2).
                        // Ma per x^2+y^2 la complessità è alta.
                        if (total_degree(Nx) > 0) {
                            return fail<QRDecomposition>(CASError{
                                CASErrorKind::Unimplemented,
                                "Symbolic QR with unresolved symbolic square roots is mathematically computed but rejected to prevent downstream simplification timeouts",
                                std::nullopt
                            });
                        }
                    }
                }
            }
        }
        
        reflectors[k] = Reflector{x, Nx, S, alpha};

        for (std::size_t j = k + 1; j < n; ++j) {
            std::vector<ExprPtr> y(m - k);
            for (std::size_t i = k; i < m; ++i) y[i - k] = R(i, j);
            
            auto res = apply_H(y, reflectors[k]);
            if (res.is_error()) return fail<QRDecomposition>(res.error());

            for (std::size_t i = k; i < m; ++i) R(i, j) = y[i - k];
        }
        
        for (std::size_t i = k + 1; i < m; ++i) {
            R(i, k) = integer(ctx, 0);
        }
        auto neg_alpha = negate_expr(ctx, alpha);
        if (neg_alpha.is_error()) return fail<QRDecomposition>(neg_alpha.error());
        R(k, k) = neg_alpha.value();
    }

    MatrixExpr Q(m, n);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            Q(i, j) = integer(ctx, i == j ? 1 : 0);
        }
    }

    for (int k_int = static_cast<int>(n) - 1; k_int >= 0; --k_int) {
        std::size_t k = static_cast<std::size_t>(k_int);
        if (is_zero_expr(reflectors[k].S)) continue;

        for (std::size_t j = 0; j < n; ++j) {
            std::vector<ExprPtr> y(m - k);
            for (std::size_t i = k; i < m; ++i) y[i - k] = Q(i, j);
            
            auto res = apply_H(y, reflectors[k]);
            if (res.is_error()) return fail<QRDecomposition>(res.error());

            for (std::size_t i = k; i < m; ++i) Q(i, j) = y[i - k];
        }
    }

    MatrixExpr R_ret(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            R_ret(i, j) = R(i, j);
        }
    }

    return ok(QRDecomposition{std::move(Q), std::move(R_ret)});
}

}  // namespace cas::linalg
