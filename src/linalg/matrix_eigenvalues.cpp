#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/algebra.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

// Division-free cofactor expansion — correct for symbolic entries where Bareiss pivots
// may not be decidably nonzero (e.g., entries of the form a_ii - lambda).
[[nodiscard]] Result<ExprPtr> cofactor_det(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n == 0U) return ok(integer(ctx, 1));
    if (n == 1U) return ok(matrix(0U, 0U));

    if (n == 2U) {
        auto ad = mul_expr(ctx, matrix(0U, 0U), matrix(1U, 1U));
        if (ad.is_error()) return ad;
        auto bc = mul_expr(ctx, matrix(0U, 1U), matrix(1U, 0U));
        if (bc.is_error()) return bc;
        return sub_expr(ctx, ad.value(), bc.value());
    }

    // Expand along first row.
    ExprPtr result = integer(ctx, 0);
    for (std::size_t col = 0; col < n; ++col) {
        ExprPtr entry = matrix(0U, col);
        if (is_zero_expr(entry)) continue;

        MatrixExpr minor(n - 1U, n - 1U);
        std::size_t out_col = 0U;
        for (std::size_t c = 0; c < n; ++c) {
            if (c == col) continue;
            for (std::size_t r = 1U; r < n; ++r) {
                minor(r - 1U, out_col) = matrix(r, c);
            }
            ++out_col;
        }

        auto minor_val = cofactor_det(minor, ctx);
        if (minor_val.is_error()) return minor_val;

        auto term = mul_expr(ctx, entry, minor_val.value());
        if (term.is_error()) return term;

        ExprPtr signed_term = term.value();
        if ((col % 2U) != 0U) {
            auto negated = simplify(ctx, ctx.arena().make<Unary>(UnaryOp::Neg, signed_term));
            if (negated.is_error()) return negated;
            signed_term = negated.value();
        }

        auto next = add_expr(ctx, result, signed_term);
        if (next.is_error()) return next;
        result = next.value();
    }
    return ok(result);
}

} // namespace

Result<ExprPtr> characteristic_polynomial(const MatrixExpr& matrix, const Symbol& lambda, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n != matrix.cols()) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "characteristic_polynomial: non-square matrix"));
    }

    if (n <= 4) { // Small matrix, cofactor is fine
        MatrixExpr lambda_i(n, n);
        ExprPtr lambda_expr = ctx.arena().make<Symbol>(lambda.name);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (i == j) {
                    auto diff = sub_expr(ctx, matrix(i, j), lambda_expr);
                    if (diff.is_error()) return diff;
                    lambda_i(i, j) = diff.value();
                } else {
                    lambda_i(i, j) = matrix(i, j);
                }
            }
        }
        return cofactor_det(lambda_i, ctx);
    }

    // Faddeev-Leverrier algorithm (n > 4)
    // p(lambda) = det(lambda*I - A) = lambda^n + c_{n-1}*lambda^(n-1) + ... + c_0
    // M_1 = A, c_{n-1} = -tr(M_1)
    // M_k = A*(M_{k-1} + c_{n-k+1}*I), c_{n-k} = -1/k * tr(M_k)
    
    std::vector<ExprPtr> c(n + 1);
    c[n] = integer(ctx, 1);
    
    MatrixExpr B(n, n);
    B.fill(integer(ctx, 0));
    MatrixExpr M(n, n);
    M.fill(integer(ctx, 0));
    
    for (std::size_t k = 1; k <= n; ++k) {
        // B = M + c_{n-k+1}*I (where M is from previous step)
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (k == 1) {
                    B(i, j) = (i == j) ? integer(ctx, 1) : integer(ctx, 0); // Not used for k=1
                } else {
                    if (i == j) {
                        auto val = add_expr(ctx, M(i, j), c[n - k + 1]);
                        if (val.is_error()) return val;
                        B(i, j) = val.value();
                    } else {
                        B(i, j) = M(i, j);
                    }
                }
            }
        }
        
        // M = A * B (except for k=1, M=A)
        if (k == 1) {
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    M(i, j) = matrix(i, j);
        } else {
            auto prod_res = multiply(matrix, B, ctx);
            if (prod_res.is_error()) return fail<ExprPtr>(prod_res.error());
            M = std::move(prod_res.value());
        }
        
        // c_{n-k} = -1/k * tr(M)
        auto tr_res = trace(M, ctx);
        if (tr_res.is_error()) return tr_res;
        
        auto neg_tr = simplify(ctx, ctx.arena().make<Unary>(UnaryOp::Neg, tr_res.value()));
        if (neg_tr.is_error()) return neg_tr;
        
        auto coeff = div_expr(ctx, neg_tr.value(), integer(ctx, static_cast<long long>(k)));
        if (coeff.is_error()) return coeff;
        
        // M-FIX: Faddeev-Leverrier together
        auto tg = algebra::together(coeff.value(), ctx);
        c[n - k] = tg.is_ok() ? tg.value() : coeff.value();
    }
    
    // Build the polynomial: det(A - lambda*I). 
    // Faddeev-Leverrier gives det(lambda*I - A).
    // det(A - lambda*I) = (-1)^n * det(lambda*I - A)
    
    std::vector<ExprPtr> terms;
    ExprPtr lambda_expr = ctx.arena().make<Symbol>(lambda.name);
    for (std::size_t i = 0; i <= n; ++i) {
        if (is_zero_expr(c[i])) continue;
        
        ExprPtr term;
        if (i == 0) {
            term = c[i];
        } else {
            ExprPtr pow_l;
            if (i == 1) {
                pow_l = lambda_expr;
            } else {
                pow_l = ctx.arena().make<Binary>(BinaryOp::Pow, lambda_expr, integer(ctx, static_cast<long long>(i)));
            }
            auto prod = mul_expr(ctx, c[i], pow_l);
            if (prod.is_error()) return prod;
            term = prod.value();
        }
        terms.push_back(term);
    }
    
    ExprPtr poly = (terms.size() == 1) ? terms[0] : ctx.arena().make<Sum>(std::move(terms));
    auto poly_s = simplify(ctx, poly);
    if (poly_s.is_error()) return poly_s;
    
    if (n % 2 != 0) {
        return simplify(ctx, ctx.arena().make<Unary>(UnaryOp::Neg, poly_s.value()));
    }
    return poly_s;
}

Result<std::vector<ExprPtr>> eigenvalues(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n == 0) return ok(std::vector<ExprPtr>{});

    Symbol lambda{"lambda"};
    auto poly_res = characteristic_polynomial(matrix, lambda, ctx);
    if (poly_res.is_error()) return fail<std::vector<ExprPtr>>(poly_res.error());

    auto roots_res = algebra::solve_polynomial(poly_res.value(), lambda, ctx);
    if (roots_res.is_error()) return fail<std::vector<ExprPtr>>(roots_res.error());

    return roots_res;
}

Result<std::vector<Eigenpair>> eigenvectors(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    auto ev_res = eigenvalues(matrix, ctx);
    if (ev_res.is_error()) return fail<std::vector<Eigenpair>>(ev_res.error());
    auto evs = std::move(ev_res.value());
    
    // Unique-ify eigenvalues (using expr_is_equal would be better, but pointer equality usually enough for simplified)
    // Actually, solve_polynomial returns all roots including multiplicity.
    // For eigenvectors, we want distinct roots.
    std::vector<ExprPtr> distinct_evs;
    for (ExprPtr e : evs) {
        bool found = false;
        for (ExprPtr d : distinct_evs) {
            if (e == d) { found = true; break; }
        }
        if (!found) distinct_evs.push_back(e);
    }
    
    const std::size_t n = matrix.rows();
    std::vector<Eigenpair> result;
    
    for (ExprPtr val : distinct_evs) {
        // Solve (A - val*I)x = 0
        MatrixExpr mat(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (i == j) {
                    auto diff = sub_expr(ctx, matrix(i, j), val);
                    if (diff.is_error()) return fail<std::vector<Eigenpair>>(diff.error());
                    mat(i, j) = diff.value();
                } else {
                    mat(i, j) = matrix(i, j);
                }
            }
        }

        auto ns_res = null_space(mat, ctx);
        if (ns_res.is_error()) return fail<std::vector<Eigenpair>>(ns_res.error());

        for (const auto& vec : ns_res.value()) {
            result.push_back(Eigenpair{.eigenvalue = val, .eigenvector = vec});
        }
    }

    return ok(std::move(result));
}

}  // namespace cas::linalg
