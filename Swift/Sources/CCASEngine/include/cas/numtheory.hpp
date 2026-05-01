#pragma once

#include "cas/bigint.hpp"
#include "cas/result.hpp"

#include <tuple>
#include <utility>
#include <vector>

namespace cas::numtheory {

using Integer = BigInt;

struct IntegerFactorization {
    Integer sign{Integer(1)};
    std::vector<std::pair<Integer, unsigned int>> prime_factors;
};

[[nodiscard]] Result<Integer> power_mod(const Integer& base, const Integer& exp, const Integer& modulus);
[[nodiscard]] Result<std::tuple<Integer, Integer, Integer>> extended_gcd(const Integer& a, const Integer& b);
[[nodiscard]] Result<Integer> modular_inverse(const Integer& a, const Integer& modulus);
[[nodiscard]] Result<Integer> chinese_remainder_theorem(
    const std::vector<Integer>& remainders,
    const std::vector<Integer>& moduli);

[[nodiscard]] Result<bool> is_prime_miller_rabin(const Integer& n, int k = 40);
[[nodiscard]] Result<bool> is_prime(const Integer& n);
[[nodiscard]] Result<Integer> next_prime(const Integer& n);
[[nodiscard]] Result<Integer> nth_prime(const Integer& n);

[[nodiscard]] Result<Integer> pollards_rho_factor(const Integer& n);
[[nodiscard]] Result<IntegerFactorization> factor_integer(const Integer& n);

[[nodiscard]] Result<Integer> euler_phi(const Integer& n);
[[nodiscard]] Result<int> moebius_mu(const Integer& n);
[[nodiscard]] Result<std::vector<Integer>> divisors(const Integer& n);
[[nodiscard]] Result<Integer> binomial(const Integer& n, const Integer& k);

[[nodiscard]] Result<std::pair<Integer, Integer>> solve_linear_diophantine(
    const Integer& a,
    const Integer& b,
    const Integer& c);

}  // namespace cas::numtheory
