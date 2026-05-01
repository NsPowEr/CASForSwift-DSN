#include "cas/linalg/Matrix.hpp"
#include "cas/algebra.hpp"
#include "cas/ast.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace cas::linalg {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] ExprPtr integer(symbolic::CASContext& ctx, long long value) {
    return ctx.arena().make<IntegerLit>(BigInt(value));
}

[[nodiscard]] Result<ExprPtr> sub_expr(symbolic::CASContext& ctx, ExprPtr lhs, ExprPtr rhs) {
    return ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Sub, lhs, rhs));
}

[[nodiscard, maybe_unused]] bool is_zero_expr(ExprPtr expr) {
    if (const auto* il = expr_cast<IntegerLit>(expr)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(expr)) return rl->numerator.is_zero();
    return false;
}

[[nodiscard, maybe_unused]] Result<MatrixExpr> mat_sub(const MatrixExpr& A, const MatrixExpr& B, symbolic::CASContext& ctx) {
    return subtract(A, B, ctx);
}

[[nodiscard]] Result<MatrixExpr> mat_mul(const MatrixExpr& A, const MatrixExpr& B, symbolic::CASContext& ctx) {
    return multiply(A, B, ctx);
}

[[nodiscard]] Result<MatrixExpr> matrix_power(const MatrixExpr& A, unsigned int p, symbolic::CASContext& ctx) {
    auto res = identity(A.rows(), ctx);
    if (res.is_error()) return res;
    if (p == 0) return res;
    MatrixExpr result = res.value();
    MatrixExpr base = A;
    while (p > 0) {
        if (p % 2 == 1) {
            auto next = mat_mul(result, base, ctx);
            if (next.is_error()) return next;
            result = next.value();
        }
        if (p > 1) {
            auto next = mat_mul(base, base, ctx);
            if (next.is_error()) return next;
            base = next.value();
        }
        p /= 2;
    }
    return ok(result);
}

// Check if a vector is in the span of a set of vectors.
// We do this by putting basis as rows of a matrix, adding v as a row, and checking rank.
[[nodiscard]] bool is_in_span(const std::vector<std::vector<ExprPtr>>& basis, const std::vector<ExprPtr>& v, symbolic::CASContext& ctx) {
    if (basis.empty()) return false;
    const std::size_t n = v.size();
    MatrixExpr M(basis.size() + 1, n);
    for (std::size_t i = 0; i < basis.size(); ++i) {
        for (std::size_t j = 0; j < n; ++j) M(i, j) = basis[i][j];
    }
    for (std::size_t j = 0; j < n; ++j) M(basis.size(), j) = v[j];
    
    auto r = rank(M, ctx);
    if (r.is_error()) return false;
    return r.value() == basis.size();
}

// Extends a basis for space V_small to a basis for V_large.
// V_small is assumed to be a subspace of V_large.
[[nodiscard, maybe_unused]] Result<std::vector<std::vector<ExprPtr>>> extend_basis(
    const std::vector<std::vector<ExprPtr>>& basis_small,
    const std::vector<std::vector<ExprPtr>>& basis_large,
    symbolic::CASContext& ctx) {
    
    std::vector<std::vector<ExprPtr>> result = basis_small;
    for (const auto& v : basis_large) {
        if (!is_in_span(result, v, ctx)) {
            result.push_back(v);
        }
    }
    return ok(result);
}

[[nodiscard]] Result<std::vector<ExprPtr>> apply_matrix(const MatrixExpr& A, const std::vector<ExprPtr>& v, symbolic::CASContext& ctx) {
    const std::size_t rows = A.rows();
    const std::size_t cols = A.cols();
    if (v.size() != cols) return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument, "Vector dimension mismatch"));
    
    std::vector<ExprPtr> res(rows);
    for (std::size_t i = 0; i < rows; ++i) {
        ExprPtr sum = integer(ctx, 0);
        for (std::size_t j = 0; j < cols; ++j) {
            auto prod = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Mul, A(i, j), v[j]));
            if (prod.is_error()) return fail<std::vector<ExprPtr>>(prod.error());
            auto next = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Add, sum, prod.value()));
            if (next.is_error()) return fail<std::vector<ExprPtr>>(next.error());
            sum = next.value();
        }
        res[i] = sum;
    }
    return ok(res);
}

} // namespace

Result<JordanDecomposition> jordan_normal_form(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    if (matrix.rows() != matrix.cols()) {
        return fail<JordanDecomposition>(make_error(CASErrorKind::InvalidArgument, "Jordan Normal Form requires a square matrix"));
    }
    const std::size_t n = matrix.rows();
    if (n == 0) return ok(JordanDecomposition{MatrixExpr(0, 0), MatrixExpr(0, 0)});

    const Symbol lambda_var("_lambda_");
    auto char_poly_res = characteristic_polynomial(matrix, lambda_var, ctx);
    if (char_poly_res.is_error()) return fail<JordanDecomposition>(char_poly_res.error());
    
    auto factors_res = algebra::factor_over_integers(char_poly_res.value(), lambda_var, ctx);
    if (factors_res.is_error()) return fail<JordanDecomposition>(factors_res.error());
    
    const auto& factorization = factors_res.value();
    
    // We need to find roots of each factor.
    // For simplicity, we assume we can solve for eigenvalues.
    std::vector<std::pair<ExprPtr, unsigned int>> eigenvalues_with_multiplicity;
    for (const auto& f : factorization.factors) {
        auto roots = algebra::solve_polynomial(f.factor, lambda_var, ctx);
        if (roots.is_error()) return fail<JordanDecomposition>(roots.error());
        for (auto root : roots.value()) {
            eigenvalues_with_multiplicity.push_back({root, f.multiplicity});
        }
    }

    MatrixExpr J(n, n);
    for(size_t i=0; i<n; ++i) for(size_t j=0; j<n; ++j) J(i,j) = integer(ctx, 0);
    MatrixExpr P(n, n);
    std::size_t current_col = 0;

    for (const auto& [eigenval, mult] : eigenvalues_with_multiplicity) {
        // (A - lambda*I)
        MatrixExpr A_minus_lambdaI(n, n);
        for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t c = 0; c < n; ++c) {
                if (r == c) {
                    auto val = sub_expr(ctx, matrix(r, c), eigenval);
                    if (val.is_error()) return fail<JordanDecomposition>(val.error());
                    A_minus_lambdaI(r, c) = val.value();
                } else {
                    A_minus_lambdaI(r, c) = matrix(r, c);
                }
            }
        }

        // Compute kernels of (A - lambda*I)^k
        std::vector<std::vector<std::vector<ExprPtr>>> kernels; // kernels[k] is kernel of S^k
        std::vector<std::size_t> dims;
        dims.push_back(0);
        
        unsigned int k = 1;
        while (dims.back() < mult && k <= mult) {
            auto Sk = matrix_power(A_minus_lambdaI, k, ctx);
            if (Sk.is_error()) return fail<JordanDecomposition>(Sk.error());
            auto ker = null_space(Sk.value(), ctx);
            if (ker.is_error()) return fail<JordanDecomposition>(ker.error());
            kernels.push_back(ker.value());
            dims.push_back(ker.value().size());
            k++;
        }
        
        unsigned int max_k = k - 1;
        
        // Construct Jordan chains
        // We work from max_k down to 1.
        // For each k, we need vectors in Ker(S^k) that are not in Ker(S^{k-1}) + S(Ker(S^{k+1}))
        // Actually, the vectors we chose at k+1, their images S(v) are in Ker(S^k).
        
        std::vector<std::vector<std::vector<ExprPtr>>> chosen_vectors(max_k + 1); // vectors chosen AT level k (not their images)
        
        for (int i = (int)max_k; i >= 1; --i) {
            // Space already covered by images of vectors from i+1
            std::vector<std::vector<ExprPtr>> covered_basis;
            if (i < (int)max_k) {
                for (const auto& v : chosen_vectors[i + 1]) {
                    auto img = apply_matrix(A_minus_lambdaI, v, ctx);
                    if (img.is_error()) return fail<JordanDecomposition>(img.error());
                    covered_basis.push_back(img.value());
                }
            }
            // And also by Ker(S^{i-1})
            std::vector<std::vector<ExprPtr>> full_small_basis;
            if (i > 1) {
                full_small_basis = kernels[i - 2]; // kernels is 0-indexed, so ker(S^1) is kernels[0]
            }
            
            // Combine them
            std::vector<std::vector<ExprPtr>> basis_to_extend = full_small_basis;
            for (const auto& v : covered_basis) {
                if (!is_in_span(basis_to_extend, v, ctx)) {
                    basis_to_extend.push_back(v);
                }
            }
            
            // Extend to Ker(S^i)
            const auto& ker_i = kernels[i - 1];
            for (const auto& v : ker_i) {
                if (!is_in_span(basis_to_extend, v, ctx)) {
                    basis_to_extend.push_back(v);
                    chosen_vectors[i].push_back(v);
                }
            }
        }
        
        // Now build the blocks and P columns
        for (unsigned int i = 1; i <= max_k; ++i) {
            for (const auto& v : chosen_vectors[i]) {
                // This vector generates a block of size i.
                // The chain is S^{i-1}v, S^{i-2}v, ..., v
                std::vector<std::vector<ExprPtr>> chain;
                std::vector<ExprPtr> current_v = v;
                chain.push_back(current_v);
                for (unsigned int step = 1; step < i; ++step) {
                    auto next_v = apply_matrix(A_minus_lambdaI, current_v, ctx);
                    if (next_v.is_error()) return fail<JordanDecomposition>(next_v.error());
                    current_v = next_v.value();
                    chain.push_back(current_v);
                }
                
                // Add to P in reverse order to get 1s on the superdiagonal
                std::reverse(chain.begin(), chain.end());
                for (const auto& vec : chain) {
                    for (std::size_t r = 0; r < n; ++r) {
                        P(r, current_col) = vec[r];
                    }
                    current_col++;
                }
                
                // Add block to J
                std::size_t block_start = current_col - i;
                for (std::size_t b = 0; b < i; ++b) {
                    J(block_start + b, block_start + b) = eigenval;
                    if (b < i - 1) {
                        J(block_start + b, block_start + b + 1) = integer(ctx, 1);
                    }
                }
            }
        }
    }

    return ok(JordanDecomposition{std::move(J), std::move(P)});
}

} // namespace cas::linalg
