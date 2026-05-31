#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"
#include "cas/numtheory.hpp"
#include "cas/error_helpers.hpp"
#include "matrix_determinant.hpp"

#include <cmath>
#include <vector>
#include <algorithm>

namespace cas::linalg {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{.kind = kind, .message = std::move(message), .hint = std::nullopt};
}

// Determinant over Z_p using Gaussian elimination.
[[nodiscard]] Result<BigInt> determinant_z_p(const MatrixExpr& matrix, const BigInt& p) {
    const std::size_t n = matrix.rows();
    std::vector<std::vector<BigInt>> mat(n, std::vector<BigInt>(n));
    
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            auto val = try_get_bigint(matrix(i, j));
            if (!val) return fail<BigInt>(make_error(CASErrorKind::InvalidArgument, "determinant_z_p: non-integer entry"));
            mat[i][j] = (*val % p + p) % p;
        }
    }

    BigInt det(1);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t pivot = i;
        while (pivot < n && (mat[pivot][i] % p).is_zero()) ++pivot;
        if (pivot == n) return ok(BigInt(0));

        if (pivot != i) {
            std::swap(mat[i], mat[pivot]);
            det = (p - det) % p;
        }

        det = (det * mat[i][i]) % p;
        auto inv_res = numtheory::modular_inverse(mat[i][i], p);
        if (inv_res.is_error()) return fail<BigInt>(inv_res.error());
        BigInt inv = inv_res.value();

        for (std::size_t row = i + 1; row < n; ++row) {
            BigInt factor = (mat[row][i] * inv) % p;
            for (std::size_t col = i; col < n; ++col) {
                mat[row][col] = (mat[row][col] - factor * mat[i][col] % p + p) % p;
            }
        }
    }
    return ok(det);
}

} // namespace

Result<ExprPtr> determinant_modular(const MatrixExpr& matrix, symbolic::CASContext& ctx) {
    const std::size_t n = matrix.rows();
    if (n == 0) return ok(integer(ctx, 1));
    if (n != matrix.cols()) return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Matrix must be square"));

    // Check if all entries are integers
    bool all_integers = true;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (!try_get_bigint(matrix(i, j)).has_value()) {
                all_integers = false;
                break;
            }
        }
        if (!all_integers) break;
    }

    if (!all_integers) {
        return detail::bareiss_determinant(matrix, ctx);
    }

    // Hadamard bound for determinant: |det| ≤ ∏_i ‖row_i‖₂.
    // log₂(‖row_i‖₂²) = log₂(Σ_j val[i][j]²) ≤ log₂(n) + 2·max_j log₂|val[i][j]|
    // Use BigInt::bit_length() to avoid stod overflow on giant integers
    // (val.to_double() blows up for entries >~10³⁰⁸).
    double log_det_bound = 0;
    const double log2_n = std::log2(static_cast<double>(std::max<std::size_t>(n, 1U)));
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t max_bits = 0;
        for (std::size_t j = 0; j < n; ++j) {
            auto val = try_get_bigint(matrix(i, j));
            if (val) max_bits = std::max(max_bits, val->bit_length());
        }
        // log₂(‖row_i‖₂) ≤ (1/2)·(log₂(n) + 2·max_bits) = (log₂(n))/2 + max_bits
        log_det_bound += 0.5 * log2_n + static_cast<double>(max_bits);
    }
    
    BigInt bound(1);
    for (std::size_t i = 0; i < static_cast<std::size_t>(std::ceil(log_det_bound) + 1); ++i) {
        bound = bound * BigInt(2);
    }
    
    std::vector<BigInt> remainders;
    std::vector<BigInt> moduli;
    BigInt current_m(1);
    
    BigInt p(1000000007); 
    while (current_m < bound) {
        auto next_p_res = numtheory::next_prime(p);
        if (next_p_res.is_error()) return fail<ExprPtr>(next_p_res.error());
        p = next_p_res.value();

        auto det_p = determinant_z_p(matrix, p);
        if (det_p.is_error()) {
            // If modular inverse failed, it might be that p divides a denominator.
            // But we checked they are integers. 
            // Most likely p is not prime (if next_prime failed) or matrix is singular mod p.
            return fail<ExprPtr>(det_p.error());
        }
        
        remainders.push_back(det_p.value());
        moduli.push_back(p);
        current_m = current_m * p;
        
        if (remainders.size() == 1) {
             BigInt test_p(997);
             auto test_det = determinant_z_p(matrix, test_p);
             if (test_det.is_ok()) {
                 if ((det_p.value() % test_p) != test_det.value()) {
                     // Potential collision or bug
                 }
             }
        }
    }
    
    auto result = numtheory::chinese_remainder_theorem(remainders, moduli);
    if (result.is_error()) return fail<ExprPtr>(result.error());
    
    BigInt det = result.value();
    if (det > (current_m / BigInt(2))) det = det - current_m;
    
    return ok(ctx.arena().make<IntegerLit>(det));
}

} // namespace cas::linalg
