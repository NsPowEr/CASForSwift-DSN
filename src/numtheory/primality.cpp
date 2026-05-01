#include "cas/numtheory.hpp"

#include "cas/error.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
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

[[nodiscard]] std::vector<Integer> miller_rabin_bases(int k) {
    static constexpr std::array<std::int64_t, 40> bases = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
        127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    };

    const std::size_t count = k <= 0 ? 0U : std::min<std::size_t>(static_cast<std::size_t>(k), bases.size());
    std::vector<Integer> result;
    result.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        result.emplace_back(bases[index]);
    }
    return result;
}

[[nodiscard]] Result<bool> run_miller_rabin(const Integer& n, const std::vector<Integer>& bases) {
    static const Integer zero(0);
    static const Integer one(1);
    static const Integer two(2);

    if (n < two) {
        return ok(false);
    }
    if (n == two) {
        return ok(true);
    }
    if (n.to_u64() % 2 == 0) {
        return ok(false);
    }

    static const std::array<std::int64_t, 12> small_primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (std::int64_t prime_value : small_primes) {
        const Integer prime(prime_value);
        if (n == prime) {
            return ok(true);
        }
        if ((n % prime).is_zero()) {
            return ok(false);
        }
    }

    Integer d = n - one;
    std::size_t s = 0U;
    while (d.to_u64() % 2 == 0 && !d.is_zero()) {
        d = d.shift_right_bits(1);
        ++s;
    }

    const Integer n_minus_one = n - one;

    for (const Integer& raw_base : bases) {
        const Integer base = raw_base % n;
        if (base.is_zero()) {
            continue;
        }

        auto witness = power_mod(base, d, n);
        if (witness.is_error()) {
            return fail<bool>(witness.error());
        }

        if (witness.value() == one || witness.value() == n_minus_one) {
            continue;
        }

        bool composite = true;
        Integer value = witness.value();
        for (std::size_t round = 1U; round < s; ++round) {
            auto squared = power_mod(value, two, n);
            if (squared.is_error()) {
                return fail<bool>(squared.error());
            }
            value = squared.value();
            if (value == n_minus_one) {
                composite = false;
                break;
            }
        }

        if (composite) {
            return ok(false);
        }
    }

    return ok(true);
}

[[nodiscard]] bool fits_in_u64(const Integer& n) {
    static const Integer max_u64 = Integer::from_u64(std::numeric_limits<std::uint64_t>::max());
    return !n.is_negative() && n <= max_u64;
}

// Miller-Rabin test implementato in is_prime.

[[nodiscard]] Result<std::size_t> integer_to_size_t(const Integer& value) {
    if (value.is_negative()) {
        return fail<std::size_t>(make_error(
            CASErrorKind::InvalidArgument,
            "Il valore deve essere non negativo"));
    }

    if (!fits_in_u64(value)) {
        return fail<std::size_t>(make_error(
            CASErrorKind::Unimplemented,
            "Il valore eccede il range gestibile dal crivello locale corrente"));
    }

    const std::uint64_t u64_val = value.to_u64();
    if (u64_val > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        return fail<std::size_t>(make_error(
            CASErrorKind::Unimplemented,
            "Il valore eccede la capacita' del size_t nativo"));
    }

    return ok(static_cast<std::size_t>(u64_val));
}

[[nodiscard]] Result<Integer> nth_prime_with_sieve(std::size_t target_index) {
    if (target_index == 0U) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "nth_prime richiede un indice strettamente positivo"));
    }

    std::size_t upper_bound = 16U;
    while (true) {
        std::vector<bool> is_composite(upper_bound + 1U, false);
        std::vector<std::size_t> primes;
        primes.reserve(target_index);

        for (std::size_t candidate = 2U; candidate <= upper_bound; ++candidate) {
            if (is_composite[candidate]) {
                continue;
            }

            primes.push_back(candidate);
            if (primes.size() == target_index) {
                return ok(Integer::from_u64(candidate));
            }

            if (candidate > upper_bound / candidate) {
                continue;
            }

            for (std::size_t multiple = candidate * candidate; multiple <= upper_bound; multiple += candidate) {
                is_composite[multiple] = true;
            }
        }

        if (upper_bound > std::numeric_limits<std::size_t>::max() / 2U) {
            return fail<Integer>(make_error(
                CASErrorKind::Unimplemented,
                "nth_prime supera il range gestibile dal crivello locale corrente"));
        }
        upper_bound *= 2U;
    }
}

[[nodiscard]] Integer normalize_prime_candidate(Integer candidate) {
    static const Integer six(6);

    if (candidate.to_u64() % 2 == 0) {
        candidate += Integer(1);
    }

    const Integer remainder = candidate % six;
    if (remainder == Integer(1) || remainder == Integer(5)) {
        return candidate;
    }
    if (remainder == Integer(3)) {
        return candidate + Integer(2);
    }
    return candidate + (Integer(5) - remainder);
}

}  // namespace

Result<bool> is_prime_miller_rabin(const Integer& n, int k) {
    if (k <= 0) {
        return fail<bool>(make_error(
            CASErrorKind::InvalidArgument,
            "Miller-Rabin richiede almeno una base di test"));
    }
    return run_miller_rabin(n, miller_rabin_bases(k));
}

Result<bool> is_prime(const Integer& n) {
    if (n < Integer(2)) return ok(false);
    if (n == Integer(2) || n == Integer(3)) return ok(true);
    if (n.to_u64() % 2 == 0) return ok(false);

    std::uint64_t n_u64 = 0;
    bool small = false;
    if (fits_in_u64(n)) {
        n_u64 = n.to_u64();
        small = true;
    }

    std::vector<Integer> bases;
    if (small) {
        if (n_u64 < 2047ULL) {
            bases = {Integer(2)};
        } else if (n_u64 < 1373653ULL) {
            bases = {Integer(2), Integer(3)};
        } else if (n_u64 < 9080191ULL) {
            bases = {Integer(31), Integer(73)};
        } else if (n_u64 < 25326001ULL) {
            bases = {Integer(2), Integer(3), Integer(5)};
        } else if (n_u64 < 3215031751ULL) {
            bases = {Integer(2), Integer(3), Integer(5), Integer(7)};
        } else if (n_u64 < 4759123141ULL) {
            bases = {Integer(2), Integer(7), Integer(61)};
        } else if (n_u64 < 1122004669633ULL) {
            bases = {Integer(2), Integer(13), Integer(23), Integer(1662803)};
        } else if (n_u64 < 3474749660383ULL) {
            bases = {Integer(2), Integer(3), Integer(5), Integer(7), Integer(11), Integer(13)};
        } else if (n_u64 < 341550071728321ULL) {
            bases = {Integer(2), Integer(3), Integer(5), Integer(7), Integer(11), Integer(13), Integer(17)};
        } else {
            // Per n < 2^64
            bases = {Integer(2), Integer(3), Integer(5), Integer(7), Integer(11), Integer(13), Integer(17), Integer(19), Integer(23)};
        }
    } else {
        // Probabilistico per numeri grandi
        return run_miller_rabin(n, miller_rabin_bases(40));
    }

    return run_miller_rabin(n, bases);
}

Result<Integer> next_prime(const Integer& n) {
    static const Integer two(2);
    static const Integer three(3);
    static const Integer six(6);

    if (n < two) {
        return ok(two);
    }

    Integer candidate = n + Integer(1);
    if (candidate <= three) {
        return ok(three);
    }
    candidate = normalize_prime_candidate(std::move(candidate));
    Integer step = (candidate % six) == Integer(5) ? Integer(2) : Integer(4);

    while (true) {
        auto prime = is_prime(candidate);
        if (prime.is_error()) {
            return fail<Integer>(prime.error());
        }
        if (prime.value()) {
            return ok(candidate);
        }
        candidate += step;
        step = six - step;
    }
}

Result<Integer> nth_prime(const Integer& n) {
    if (n <= Integer(0)) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "nth_prime richiede un indice strettamente positivo"));
    }

    auto target_index = integer_to_size_t(n);
    if (target_index.is_ok()) {
        return nth_prime_with_sieve(target_index.value());
    }

    Integer count(0);
    Integer candidate(1);
    while (count < n) {
        auto next = next_prime(candidate);
        if (next.is_error()) {
            return next;
        }
        candidate = next.value();
        count += Integer(1);
    }
    return ok(candidate);
}

}  // namespace cas::numtheory
