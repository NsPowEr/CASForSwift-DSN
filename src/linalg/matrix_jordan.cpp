#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include <vector>
#include <string>

namespace cas::linalg {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

[[nodiscard]] Result<MatrixExpr> matrix_power(const MatrixExpr& A, unsigned int k, symbolic::CASContext& ctx) {
    const std::size_t n = A.rows();
    if (k == 0) return identity(n, ctx);
    if (k == 1) return ok(A);

    auto res = identity(n, ctx);
    if (res.is_error()) return res;
    MatrixExpr result = std::move(res.value());

    for (unsigned int i = 0; i < k; ++i) {
        auto next = multiply(result, A, ctx);
        if (next.is_error()) return next;
        result = std::move(next.value());
    }
    return ok(std::move(result));
}

[[nodiscard]] Result<std::vector<ExprPtr>> apply_matrix(const MatrixExpr& A, const std::vector<ExprPtr>& v, symbolic::CASContext& ctx) {
    const std::size_t rows = A.rows();
    const std::size_t cols = A.cols();
    if (v.size() != cols) return fail<std::vector<ExprPtr>>(make_error(CASErrorKind::InvalidArgument, "Vector dimension mismatch"));

    std::vector<ExprPtr> res(rows);
    for (std::size_t i = 0; i < rows; ++i) {
        ExprPtr sum = integer(ctx, 0);
        for (std::size_t j = 0; j < cols; ++j) {
            auto prod = mul_expr(ctx, A(i, j), v[j]);
            if (prod.is_error()) return fail<std::vector<ExprPtr>>(prod.error());
            auto next_sum = add_expr(ctx, sum, prod.value());
            if (next_sum.is_error()) return fail<std::vector<ExprPtr>>(next_sum.error());
            sum = next_sum.value();
        }
        res[i] = sum;
    }
    return ok(res);
}

// Build MatrixExpr (n × cols) whose columns are the given vectors.
[[nodiscard]] MatrixExpr stack_columns(const std::vector<std::vector<ExprPtr>>& cols, std::size_t n, symbolic::CASContext& ctx) {
    MatrixExpr M(n, cols.size());
    for (std::size_t j = 0; j < cols.size(); ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            M(i, j) = (i < cols[j].size()) ? cols[j][i] : integer(ctx, 0);
        }
    }
    return M;
}

// Returns true if vector v belongs to span(basis) over the symbolic field, by
// comparing rank([basis]) with rank([basis | v]). Uses rank() helper which
// dispatches Bareiss for symbolic matrices.
[[nodiscard]] Result<bool> is_in_span(const std::vector<ExprPtr>& v,
                                       const std::vector<std::vector<ExprPtr>>& basis,
                                       symbolic::CASContext& ctx) {
    if (basis.empty()) {
        for (const auto& e : v) if (!is_zero_expr(e)) return ok(false);
        return ok(true);
    }
    const std::size_t n = v.size();
    MatrixExpr M_basis = stack_columns(basis, n, ctx);
    auto rank_basis = rank(M_basis, ctx);
    if (rank_basis.is_error()) return fail<bool>(rank_basis.error());

    std::vector<std::vector<ExprPtr>> with_v = basis;
    with_v.push_back(v);
    MatrixExpr M_full = stack_columns(with_v, n, ctx);
    auto rank_full = rank(M_full, ctx);
    if (rank_full.is_error()) return fail<bool>(rank_full.error());

    return ok(rank_full.value() == rank_basis.value());
}

} // namespace

Result<JordanDecomposition> jordan_normal_form(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n != matrix.cols()) {
        return fail<JordanDecomposition>(make_error(CASErrorKind::InvalidArgument, "Jordan form requires a square matrix"));
    }

    Symbol lambda{"lambda"};
    auto char_poly_res = characteristic_polynomial(matrix, lambda, ctx);
    if (char_poly_res.is_error()) return fail<JordanDecomposition>(char_poly_res.error());
    ExprPtr char_poly = char_poly_res.value();

    auto factors_res = algebra::factor_polynomial(char_poly, lambda, ctx);
    if (factors_res.is_error()) return fail<JordanDecomposition>(factors_res.error());
    const auto& factorization = factors_res.value();

    MatrixExpr J(n, n);
    J.fill(integer(ctx, 0));
    MatrixExpr P(n, n);
    P.fill(integer(ctx, 0));
    std::size_t p_col = 0;
    std::size_t j_block_start = 0;

    for (const auto& fact : factorization.factors) {
        ExprPtr root = fact.factor;
        auto roots = algebra::solve_polynomial(root, lambda, ctx);
        if (roots.is_error() || roots.value().empty()) continue;
        ExprPtr val = roots.value()[0];
        const unsigned int m = fact.multiplicity;

        // A - val*I
        MatrixExpr root_i(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (i == j) {
                    auto diff = sub_expr(ctx, matrix(i, j), val);
                    if (diff.is_error()) return fail<JordanDecomposition>(diff.error());
                    root_i(i, j) = diff.value();
                } else {
                    root_i(i, j) = matrix(i, j);
                }
            }
        }

        // kernels[k] = basis of null_space((A-val*I)^k) for k = 0..m+1.
        // kernels[0] is conventionally {} (dim 0). kernels[m+1] stabilizes to kernels[m].
        std::vector<std::vector<std::vector<ExprPtr>>> kernels(m + 2);
        kernels[0] = {};
        for (unsigned int k = 1; k <= m; ++k) {
            auto pow_res = matrix_power(root_i, k, ctx);
            if (pow_res.is_error()) return fail<JordanDecomposition>(pow_res.error());
            auto ns_res = null_space(pow_res.value(), ctx);
            if (ns_res.is_error()) return fail<JordanDecomposition>(ns_res.error());
            kernels[k] = ns_res.value();
        }
        kernels[m + 1] = kernels[m];

        // d[k] = dim ker((A-val*I)^k). Used in Filippov's formula
        //   n_k = 2·d[k] − d[k−1] − d[k+1]
        // which gives the number of Jordan blocks of size exactly k for this eigenvalue.
        std::vector<std::size_t> d(m + 2, 0);
        for (unsigned int k = 0; k <= m + 1; ++k) d[k] = kernels[k].size();

        std::vector<unsigned int> n_blocks(m + 1, 0);
        for (unsigned int k = 1; k <= m; ++k) {
            const long long v = 2LL * static_cast<long long>(d[k])
                              - static_cast<long long>(d[k - 1])
                              - static_cast<long long>(d[k + 1]);
            n_blocks[k] = (v > 0) ? static_cast<unsigned int>(v) : 0U;
        }

        // Build chains top-down. Each chain c of length s contributes:
        //   c[s-1] = top vector v_s ∈ ker((A-val*I)^s) \ ker((A-val*I)^{s-1})
        //   c[i]   = (A-val*I) · c[i+1] for i = s-2..0  (chain[0] is eigenvector)
        std::vector<std::vector<std::vector<ExprPtr>>> chains;

        for (int kk = static_cast<int>(m); kk >= 1; --kk) {
            const unsigned int k = static_cast<unsigned int>(kk);
            if (n_blocks[k] == 0) continue;

            // Vectors already used at level k from larger blocks: each chain of size s ≥ k
            // owns one vector at level k, namely chain[k-1].
            std::vector<std::vector<ExprPtr>> tops_at_k;
            for (const auto& chain : chains) {
                if (chain.size() >= k) tops_at_k.push_back(chain[k - 1]);
            }

            for (unsigned int b = 0; b < n_blocks[k]; ++b) {
                // span_basis = ker((A-val*I)^{k-1}) ∪ tops already placed at level k
                std::vector<std::vector<ExprPtr>> span_basis;
                for (const auto& bv : kernels[k - 1]) span_basis.push_back(bv);
                for (const auto& tv : tops_at_k) span_basis.push_back(tv);

                std::vector<ExprPtr> chosen;
                for (const auto& cand : kernels[k]) {
                    auto in_span = is_in_span(cand, span_basis, ctx);
                    if (in_span.is_error()) return fail<JordanDecomposition>(in_span.error());
                    if (!in_span.value()) { chosen = cand; break; }
                }
                if (chosen.empty()) {
                    return make_unimplemented<JordanDecomposition>(
                        "linalg", "jordan_normal_form",
                        "could not find generalized eigenvector outside span",
                        error::reason_codes::LINALG_LINEAR_DEPENDENT,
                        "Symbolic rank determination may be incomplete for this matrix",
                        "F4-Jordan");
                }

                // Build chain v_1..v_k (chain[i] = v_{i+1}); chain[k-1] = chosen top.
                std::vector<std::vector<ExprPtr>> chain(k);
                chain[k - 1] = chosen;
                for (int i = static_cast<int>(k) - 2; i >= 0; --i) {
                    auto next = apply_matrix(root_i, chain[i + 1], ctx);
                    if (next.is_error()) return fail<JordanDecomposition>(next.error());
                    chain[i] = next.value();
                }

                chains.push_back(chain);
                tops_at_k.push_back(chosen);

                // Place this block into J and P. Convention: J(j+i, j+i) = val,
                // J(j+i, j+i+1) = 1 for i < k-1. Columns of P are chain[0]..chain[k-1].
                for (unsigned int i = 0; i < k; ++i) {
                    J(j_block_start + i, j_block_start + i) = val;
                    if (i + 1 < k) J(j_block_start + i, j_block_start + i + 1) = integer(ctx, 1);
                    for (std::size_t r = 0; r < n; ++r) P(r, p_col) = chain[i][r];
                    p_col++;
                }
                j_block_start += k;
            }
        }
    }

    return ok(JordanDecomposition{std::move(J), std::move(P)});
}

} // namespace cas::linalg
