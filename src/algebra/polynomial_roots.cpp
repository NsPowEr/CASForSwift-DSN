#include "polynomial_internal.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {

std::optional<RationalRootCandidate> find_rational_root_candidate(const IntPoly& coefficients) {
    if (coefficients.size() <= 1U) {
        return std::nullopt;
    }
    if (coefficients.constant_term().is_zero()) {
        return RationalRootCandidate{
            .numerator = BigInt(0),
            .denominator = BigInt(1),
        };
    }

    const std::vector<BigInt> numerator_divisors = positive_divisors_or_one(coefficients.constant_term());
    const std::vector<BigInt> denominator_divisors = positive_divisors_or_one(coefficients.leading_coeff());

    for (const BigInt& denominator_divisor : denominator_divisors) {
        for (const BigInt& numerator_divisor : numerator_divisors) {
            for (int sign : {1, -1}) {
                BigInt numerator = sign > 0 ? numerator_divisor : -numerator_divisor;
                BigInt denominator = denominator_divisor;
                const BigInt common = gcd(numerator.abs(), denominator);
                if (common != BigInt(1)) {
                    numerator /= common;
                    denominator /= common;
                }

                if (evaluate_integer_polynomial_at(coefficients, Rational(numerator, denominator)).numerator().is_zero()) {
                    return RationalRootCandidate{
                        .numerator = std::move(numerator),
                        .denominator = std::move(denominator),
                    };
                }
            }
        }
    }

    return std::nullopt;
}

}  // namespace algebra
}  // namespace cas
