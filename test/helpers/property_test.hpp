#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace cas::test {

class DeterministicRng {
public:
    explicit DeterministicRng(std::uint64_t seed) noexcept : state_(seed) {}

    [[nodiscard]] std::uint64_t next() noexcept {
        state_ = state_ * 6364136223846793005ULL + 1442695040888963407ULL;
        return state_;
    }

    [[nodiscard]] std::size_t next_index(std::size_t upper_exclusive) noexcept {
        return upper_exclusive == 0U ? 0U : static_cast<std::size_t>(next() % upper_exclusive);
    }

    [[nodiscard]] int next_int(int min_inclusive, int max_inclusive) noexcept {
        const auto span = static_cast<std::uint64_t>(max_inclusive - min_inclusive + 1);
        return min_inclusive + static_cast<int>(next() % span);
    }

    [[nodiscard]] bool coin_flip() noexcept {
        return (next() & 1ULL) != 0ULL;
    }

private:
    std::uint64_t state_;
};

[[nodiscard]] inline std::string sample_symbol(DeterministicRng& rng) {
    static constexpr std::string_view kSymbols[] = {"x", "y", "z", "a", "b", "c"};
    return std::string(kSymbols[rng.next_index(sizeof(kSymbols) / sizeof(kSymbols[0]))]);
}

[[nodiscard]] inline std::string sample_exact_atom(DeterministicRng& rng) {
    if (rng.coin_flip()) {
        return sample_symbol(rng);
    }
    return std::to_string(rng.next_int(1, 5));
}

[[nodiscard]] inline std::string maybe_parenthesize(std::string text) {
    return "(" + std::move(text) + ")";
}

[[nodiscard]] inline std::string generate_polynomial_expr(DeterministicRng& rng, int depth) {
    if (depth <= 0) {
        return sample_exact_atom(rng);
    }

    switch (rng.next_int(0, 4)) {
        case 0:
            return sample_exact_atom(rng);
        case 1:
            return maybe_parenthesize(generate_polynomial_expr(rng, depth - 1)) + " + " +
                maybe_parenthesize(generate_polynomial_expr(rng, depth - 1));
        case 2:
            return maybe_parenthesize(generate_polynomial_expr(rng, depth - 1)) + " * " +
                maybe_parenthesize(generate_polynomial_expr(rng, depth - 1));
        case 3:
            return sample_symbol(rng) + "^" + std::to_string(rng.next_int(1, 3));
        case 4:
        default:
            return std::to_string(rng.next_int(2, 4)) + " * " + sample_symbol(rng);
    }
}

template <typename Fn>
void run_seeded_cases(std::uint64_t seed, std::size_t iterations, Fn&& fn) {
    DeterministicRng rng(seed);
    for (std::size_t index = 0; index < iterations; ++index) {
        fn(rng, index);
    }
}

}  // namespace cas::test
