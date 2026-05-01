#include "cas/numtheory.hpp"

#include "cas/error.hpp"

#include <string>
#include <tuple>
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

[[nodiscard]] Result<void> require_positive_modulus(const Integer& modulus) {
    if (modulus <= Integer(0)) {
        return fail<void>(make_error(
            CASErrorKind::InvalidArgument,
            "Il modulo deve essere strettamente positivo"));
    }
    return ok();
}

[[nodiscard]] std::tuple<Integer, Integer, Integer> extended_gcd_nonnegative(
    const Integer& a,
    const Integer& b) {
    Integer old_r = a;
    Integer r = b;
    Integer old_s = Integer(1);
    Integer s = Integer(0);
    Integer old_t = Integer(0);
    Integer t = Integer(1);

    while (!r.is_zero()) {
        const Integer q = old_r / r;

        Integer next_r = old_r - q * r;
        old_r = r;
        r = std::move(next_r);

        Integer next_s = old_s - q * s;
        old_s = s;
        s = std::move(next_s);

        Integer next_t = old_t - q * t;
        old_t = t;
        t = std::move(next_t);
    }

    return {old_r, old_s, old_t};
}

}  // namespace

Result<Integer> power_mod(const Integer& base, const Integer& exp, const Integer& modulus) {
    auto modulus_check = require_positive_modulus(modulus);
    if (modulus_check.is_error()) {
        return fail<Integer>(modulus_check.error());
    }
    if (exp.is_negative()) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "L'esponente modulare non puo' essere negativo"));
    }
    if (modulus == Integer(1)) {
        return ok(Integer(0));
    }

    auto reduced_base = positive_mod(base, modulus);
    if (reduced_base.is_error()) {
        return reduced_base;
    }

    Integer result(1);
    Integer factor = reduced_base.value();
    Integer exponent = exp;
    const Integer zero(0);
    const Integer one(1);
    const Integer two(2);

    while (exponent > zero) {
        if ((exponent % two) == one) {
            auto reduced = positive_mod(result * factor, modulus);
            if (reduced.is_error()) {
                return reduced;
            }
            result = reduced.value();
        }

        exponent /= two;
        if (exponent > zero) {
            auto squared = positive_mod(factor * factor, modulus);
            if (squared.is_error()) {
                return squared;
            }
            factor = squared.value();
        }
    }

    return ok(result);
}

Result<std::tuple<Integer, Integer, Integer>> extended_gcd(const Integer& a, const Integer& b) {
    auto [g, x, y] = extended_gcd_nonnegative(a.abs(), b.abs());
    if (a.is_negative()) {
        x = -x;
    }
    if (b.is_negative()) {
        y = -y;
    }
    return ok(std::tuple<Integer, Integer, Integer>{std::move(g), std::move(x), std::move(y)});
}

Result<Integer> modular_inverse(const Integer& a, const Integer& modulus) {
    auto modulus_check = require_positive_modulus(modulus);
    if (modulus_check.is_error()) {
        return fail<Integer>(modulus_check.error());
    }

    auto gcd_result = extended_gcd(a, modulus);
    if (gcd_result.is_error()) {
        return fail<Integer>(gcd_result.error());
    }

    const auto& [g, x, y] = gcd_result.value();
    static_cast<void>(y);
    if (g != Integer(1)) {
        return fail<Integer>(make_error(
            CASErrorKind::Undefined,
            "L'inverso modulare non esiste quando gcd(a, m) != 1"));
    }

    return positive_mod(x, modulus);
}

Result<Integer> chinese_remainder_theorem(
    const std::vector<Integer>& remainders,
    const std::vector<Integer>& moduli) {
    if (remainders.empty() || moduli.empty()) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "CRT richiede almeno una congruenza"));
    }
    if (remainders.size() != moduli.size()) {
        return fail<Integer>(make_error(
            CASErrorKind::InvalidArgument,
            "CRT richiede array di resti e moduli della stessa dimensione"));
    }

    Integer solution(0);
    Integer combined_modulus(1);

    for (std::size_t index = 0; index < moduli.size(); ++index) {
        const Integer& modulus = moduli[index];
        auto modulus_check = require_positive_modulus(modulus);
        if (modulus_check.is_error()) {
            return fail<Integer>(modulus_check.error());
        }

        if (gcd(combined_modulus, modulus) != Integer(1)) {
            return fail<Integer>(make_error(
                CASErrorKind::InvalidArgument,
                "CRT richiede moduli coprimi a due a due"));
        }

        auto target_remainder = positive_mod(remainders[index], modulus);
        if (target_remainder.is_error()) {
            return target_remainder;
        }

        auto current_remainder = positive_mod(solution, modulus);
        if (current_remainder.is_error()) {
            return current_remainder;
        }

        auto inverse = modular_inverse(combined_modulus, modulus);
        if (inverse.is_error()) {
            return inverse;
        }

        auto delta = positive_mod(target_remainder.value() - current_remainder.value(), modulus);
        if (delta.is_error()) {
            return delta;
        }

        auto step = positive_mod(delta.value() * inverse.value(), modulus);
        if (step.is_error()) {
            return step;
        }

        solution += combined_modulus * step.value();
        combined_modulus *= modulus;

        auto normalized = positive_mod(solution, combined_modulus);
        if (normalized.is_error()) {
            return normalized;
        }
        solution = normalized.value();
    }

    return ok(solution);
}

Result<std::pair<Integer, Integer>> solve_linear_diophantine(
    const Integer& a,
    const Integer& b,
    const Integer& c) {
    auto gcd_result = extended_gcd(a, b);
    if (gcd_result.is_error()) {
        return fail<std::pair<Integer, Integer>>(gcd_result.error());
    }

    const auto& [g, x0, y0] = gcd_result.value();
    if ((c % g) != Integer(0)) {
        return fail<std::pair<Integer, Integer>>(make_error(
            CASErrorKind::Undefined,
            "L'equazione diofantea lineare non ha soluzioni intere"));
    }

    const Integer scale = c / g;
    return ok(std::pair<Integer, Integer>{x0 * scale, y0 * scale});
}

}  // namespace cas::numtheory
