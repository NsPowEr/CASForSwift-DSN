#include "cas/numtheory.hpp"

#include "cas/error.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cas::numtheory {
namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    return CASError{
        .kind = kind,
        .message = std::move(message),
        .hint = std::nullopt,
    };
}

[[nodiscard]] Result<Integer> positive_mod(const Integer& value, const Integer& modulus) {
    auto reduced = checked_mod(value, modulus);
    if (reduced.is_error()) {
        return fail<Integer>(reduced.error());
    }
    if (reduced.value().is_negative()) {
        return ok(reduced.value() + modulus);
    }
    return ok(reduced.value());
}

[[nodiscard]] Integer absolute_difference(const Integer& lhs, const Integer& rhs) {
    return lhs >= rhs ? (lhs - rhs) : (rhs - lhs);
}

[[nodiscard]] Result<bool> is_probably_prime(const Integer& n) {
    static const Integer max_u64 = Integer::from_u64(std::numeric_limits<std::uint64_t>::max());
    if (!n.is_negative() && n <= max_u64) {
        return is_prime(n);
    }
    return is_prime_miller_rabin(n);
}

[[nodiscard]] Result<Integer> rho_next(const Integer& x, const Integer& c, const Integer& modulus) {
    return positive_mod(x * x + c, modulus);
}

[[nodiscard]] Result<void> collect_factor(const Integer& n, std::vector<Integer>& factors, std::size_t pollard_max_iter);

[[nodiscard]] Result<void> factor_with_trial_division(Integer& n, std::vector<Integer>& factors) {
    static const std::array<std::int64_t, 25> small_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97,
    };

    for (std::int64_t prime_value : small_primes) {
        const Integer prime(prime_value);
        while ((n % prime).is_zero()) {
            factors.push_back(prime);
            n /= prime;
        }
    }
    return ok();
}

[[nodiscard]] Result<void> collect_factor(const Integer& n, std::vector<Integer>& factors, std::size_t pollard_max_iter) {
    if (n == Integer(1)) {
        return ok();
    }

    auto primality = is_probably_prime(n);
    if (primality.is_error()) {
        return fail<void>(primality.error());
    }
    if (primality.value()) {
        factors.push_back(n);
        return ok();
    }

    auto factor = pollards_rho_factor(n, pollard_max_iter);
    if (factor.is_error()) {
        return fail<void>(factor.error());
    }

    auto left = collect_factor(factor.value(), factors, pollard_max_iter);
    if (left.is_error()) {
        return left;
    }
    return collect_factor(n / factor.value(), factors, pollard_max_iter);
}

[[nodiscard]] std::vector<std::pair<Integer, unsigned int>> compress_factors(std::vector<Integer> factors) {
    std::sort(factors.begin(), factors.end());

    std::vector<std::pair<Integer, unsigned int>> compressed;
    for (const Integer& factor : factors) {
        if (!compressed.empty() && compressed.back().first == factor) {
            ++compressed.back().second;
        } else {
            compressed.emplace_back(factor, 1U);
        }
    }
    return compressed;
}

}  // namespace

Result<Integer> pollards_rho_factor(const Integer& n, std::size_t max_iter) {
    static const Integer zero(0);
    static const Integer one(1);
    static const Integer two(2);

    if (n <= one) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "Pollard Rho richiede un intero maggiore di 1"));
    }
    if ((n % two).is_zero()) {
        return ok(two);
    }

    auto primality = is_probably_prime(n);
    if (primality.is_error()) {
        return fail<Integer>(primality.error());
    }
    if (primality.value()) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "Pollard Rho richiede un intero composto"));
    }

    static const std::array<std::int64_t, 8> seeds = {2, 3, 5, 7, 11, 13, 17, 19};
    static const std::array<std::int64_t, 8> constants = {1, 3, 5, 7, 11, 13, 17, 19};
    // HPP-021 CLOSED: max_iter exposed via overload param (default 4096).
    // Caller con CASContext: pollards_rho_factor(n, ctx.pollard_rho_max_iter()).
    const std::size_t max_iterations = max_iter;

    for (std::int64_t seed_value : seeds) {
        for (std::int64_t constant_value : constants) {
            Integer x(seed_value);
            Integer y(seed_value);
            const Integer c(constant_value);
            Integer d(one);

            for (std::size_t iteration = 0; iteration < max_iterations && d == one; ++iteration) {
                auto next_x = rho_next(x, c, n);
                if (next_x.is_error()) {
                    return next_x;
                }
                x = next_x.value();

                auto next_y = rho_next(y, c, n);
                if (next_y.is_error()) {
                    return next_y;
                }
                auto next_next_y = rho_next(next_y.value(), c, n);
                if (next_next_y.is_error()) {
                    return next_next_y;
                }
                y = next_next_y.value();

                d = gcd(absolute_difference(x, y), n);
            }

            if (d != one && d != n && d != zero) {
                return ok(d);
            }
        }
    }

    return fail<Integer>(make_error(
        CASErrorKind::Unimplemented,
        "Pollard Rho non ha trovato un fattore con i parametri deterministici correnti"));
}

Result<IntegerFactorization> factor_integer(const Integer& n, std::size_t pollard_max_iter) {
    if (n.is_zero()) {
        return fail<IntegerFactorization>(make_error(
            CASErrorKind::InvalidArgument,
            "factor_integer non e' definita per zero"));
    }

    IntegerFactorization result;
    Integer remaining = n;
    if (remaining.is_negative()) {
        result.sign = Integer(-1);
        remaining = -remaining;
    }

    if (remaining == Integer(1)) {
        return ok(result);
    }

    std::vector<Integer> factors;
    auto trial_division = factor_with_trial_division(remaining, factors);
    if (trial_division.is_error()) {
        return fail<IntegerFactorization>(trial_division.error());
    }

    if (remaining != Integer(1)) {
        auto recursive = collect_factor(remaining, factors, pollard_max_iter);
        if (recursive.is_error()) {
            return fail<IntegerFactorization>(recursive.error());
        }
    }

    result.prime_factors = compress_factors(std::move(factors));
    return ok(result);
}

Result<Integer> euler_phi(const Integer& n) {
    if (n <= Integer(0)) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "Euler phi richiede un intero strettamente positivo"));
    }

    auto factorization = factor_integer(n);
    if (factorization.is_error()) {
        return fail<Integer>(factorization.error());
    }

    Integer result = n;
    for (const auto& [prime, multiplicity] : factorization.value().prime_factors) {
        static_cast<void>(multiplicity);
        result /= prime;
        result *= (prime - Integer(1));
    }
    return ok(result);
}

Result<int> moebius_mu(const Integer& n) {
    if (n <= Integer(0)) {
        return fail<int>(make_error(
            CASErrorKind::InvalidArgument,
            "La funzione di Moebius richiede un intero strettamente positivo"));
    }

    auto factorization = factor_integer(n);
    if (factorization.is_error()) {
        return fail<int>(factorization.error());
    }

    for (const auto& [prime, multiplicity] : factorization.value().prime_factors) {
        static_cast<void>(prime);
        if (multiplicity > 1U) {
            return ok(0);
        }
    }

    const std::size_t distinct_prime_count = factorization.value().prime_factors.size();
    return ok((distinct_prime_count % 2U) == 0U ? 1 : -1);
}

Result<std::vector<Integer>> divisors(const Integer& n) {
    if (n.is_zero()) {
        return fail<std::vector<Integer>>(make_error(
            CASErrorKind::InvalidArgument,
            "L'enumerazione dei divisori non e' definita per zero"));
    }

    auto factorization = factor_integer(n.abs());
    if (factorization.is_error()) {
        return fail<std::vector<Integer>>(factorization.error());
    }

    std::vector<Integer> result{Integer(1)};
    for (const auto& [prime, multiplicity] : factorization.value().prime_factors) {
        const std::size_t base_size = result.size();
        Integer prime_power(1);
        for (unsigned int exponent = 1U; exponent <= multiplicity; ++exponent) {
            prime_power *= prime;
            for (std::size_t index = 0; index < base_size; ++index) {
                result.push_back(result[index] * prime_power);
            }
        }
    }

    std::sort(result.begin(), result.end());
    return ok(result);
}

Result<Integer> binomial(const Integer& n, const Integer& k) {
    if (n.is_negative() || k.is_negative()) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "binomial richiede argomenti interi non negativi"));
    }
    if (k > n) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "binomial richiede k <= n"));
    }

    Integer effective_k = k;
    const Integer mirrored = n - k;
    if (mirrored < effective_k) {
        effective_k = mirrored;
    }

    Integer result(1);
    Integer i(1);
    while (i <= effective_k) {
        Integer numerator = n - effective_k + i;
        Integer denominator = i;

        const Integer g1 = gcd(numerator, denominator);
        numerator /= g1;
        denominator /= g1;

        const Integer g2 = gcd(result, denominator);
        result /= g2;
        denominator /= g2;

        result *= numerator;
        if (denominator != Integer(1)) {
            result /= denominator;
        }

        i += Integer(1);
    }

    return ok(result);
}

}  // namespace cas::numtheory
