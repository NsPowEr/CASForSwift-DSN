#include "cas/numtheory.hpp"

#include <cstdint>

namespace cas::numtheory {

std::vector<Rational> bernoulli_numbers(unsigned int max_index) {
    std::vector<Rational> a(max_index + 1U);
    std::vector<Rational> bernoulli(max_index + 1U);

    for (unsigned int m = 0U; m <= max_index; ++m) {
        a[m] = Rational(BigInt(1), BigInt(static_cast<std::int64_t>(m + 1U)));
        for (unsigned int j = m; j > 0U; --j) {
            a[j - 1U] = Rational(BigInt(static_cast<std::int64_t>(j))) * (a[j - 1U] - a[j]);
        }
        bernoulli[m] = a[0];
    }

    return bernoulli;
}

Rational bernoulli_number(unsigned int n) {
    auto v = bernoulli_numbers(n);
    return v[n];
}

}  // namespace cas::numtheory
