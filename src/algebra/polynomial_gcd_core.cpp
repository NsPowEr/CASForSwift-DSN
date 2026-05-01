#include "polynomial_internal.hpp"

#include <algorithm>
#include <utility>

namespace cas {
namespace algebra {

struct IntegerSubresultantExecution {
    IntPoly gcd;
    IntegerGcdPath path{IntegerGcdPath::Subresultant};
};

[[nodiscard]] IntegerSubresultantExecution run_integer_subresultant(IntPoly A, IntPoly B) {
    if (A.is_zero()) {
        return IntegerSubresultantExecution{
            .gcd = std::move(B),
            .path = IntegerGcdPath::Subresultant,
        };
    }
    if (B.is_zero()) {
        return IntegerSubresultantExecution{
            .gcd = std::move(A),
            .path = IntegerGcdPath::Subresultant,
        };
    }

    const BigInt contA = integer_content(A);
    const BigInt contB = integer_content(B);
    const BigInt common_content = gcd(contA, contB);

    IntPoly H1 = primitive_integer_poly(std::move(A));
    IntPoly H2 = primitive_integer_poly(std::move(B));

    if (H1.degree() < H2.degree()) {
        std::swap(H1, H2);
    }

    std::size_t d1 = H1.degree();
    std::size_t d2 = H2.degree();
    std::size_t delta = d1 - d2;

    // Brown/Collins Subresultant PRS
    // beta_2 = (-1)^(delta_1 + 1)
    BigInt beta = (delta % 2 == 0) ? BigInt(-1) : BigInt(1);
    BigInt psi = BigInt(-1);

    while (true) {
        IntPoly R = pseudo_remainder_integer_poly(H1, H2);
        if (R.is_zero()) {
            break;
        }

        if (!try_divide_integer_coefficients_by_scalar(R, beta)) {
            return IntegerSubresultantExecution{
                .gcd = IntPoly{},
                .path = IntegerGcdPath::PrimitiveFallback,
            };
        }

        H1 = std::move(H2);
        H2 = std::move(R);

        const std::size_t prev_delta = delta;
        d1 = H1.degree();
        d2 = H2.degree();
        delta = d1 - d2;

        const BigInt lcH1 = H1.leading_coeff();

        // psi_{i+1} = (-lc(H_i))^delta_{i-1} * psi_i^(1-delta_{i-1})
        if (prev_delta > 0) {
            psi = bigint_pow_nonnegative(-lcH1, prev_delta) / bigint_pow_nonnegative(psi, prev_delta - 1);
        } else if (prev_delta == 0) {
            // psi remains unchanged: psi = (-lcH1)^0 * psi^1
        }

        // beta_{i+1} = -lc(H_i) * psi_{i+1}^delta_i
        beta = -lcH1 * bigint_pow_nonnegative(psi, delta);
    }

    IntPoly result = primitive_integer_poly(std::move(H2));
    multiply_integer_coefficients_by_scalar(result, common_content);

    return IntegerSubresultantExecution{
        .gcd = std::move(result),
        .path = IntegerGcdPath::Subresultant,
    };
}

bool is_zero_integer_poly(const IntPoly& coefficients) {
    return coefficients.is_zero();
}

IntPoly gcd_integer_poly_primitive(IntPoly lhs, IntPoly rhs) {
    lhs = primitive_integer_poly(std::move(lhs));
    rhs = primitive_integer_poly(std::move(rhs));

    if (lhs.empty()) {
        return rhs;
    }
    if (rhs.empty()) {
        return lhs;
    }

    while (!rhs.empty()) {
        IntPoly remainder = pseudo_remainder_integer_poly(lhs, rhs);
        remainder = primitive_integer_poly(std::move(remainder));
        lhs = std::move(rhs);
        rhs = std::move(remainder);
    }

    return primitive_integer_poly(std::move(lhs));
}

std::optional<IntPoly> gcd_integer_poly_subresultant(IntPoly lhs, IntPoly rhs) {
    IntegerSubresultantExecution execution = run_integer_subresultant(std::move(lhs), std::move(rhs));
    if (execution.path != IntegerGcdPath::Subresultant) {
        return std::nullopt;
    }
    return std::move(execution.gcd);
}

IntegerGcdResult gcd_integer_poly_with_subresultant(IntPoly lhs, IntPoly rhs) {
    BigInt c_lhs = integer_content(lhs);
    BigInt c_rhs = integer_content(rhs);
    BigInt c_gcd = gcd(c_lhs, c_rhs);

    lhs = primitive_integer_poly(std::move(lhs));
    rhs = primitive_integer_poly(std::move(rhs));

    IntegerSubresultantExecution execution = run_integer_subresultant(lhs, rhs);
    
    IntPoly result_gcd;
    if (execution.path == IntegerGcdPath::Subresultant) {
        result_gcd = std::move(execution.gcd);
    } else {
        result_gcd = gcd_integer_poly_primitive(std::move(lhs), std::move(rhs));
    }
    
    if (!c_gcd.is_zero()) {
        multiply_integer_coefficients_by_scalar(result_gcd, c_gcd);
    }
    
    if (!result_gcd.is_zero() && result_gcd.leading_coeff().is_negative()) {
        multiply_integer_coefficients_by_scalar(result_gcd, BigInt(-1));
    }

    return IntegerGcdResult{
        .gcd = std::move(result_gcd),
        .path = execution.path,
    };
}

}  // namespace algebra
}  // namespace cas
