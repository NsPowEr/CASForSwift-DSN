#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <algorithm>
#include <array>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cas::algebra {

namespace {

using FactorKey = std::vector<std::pair<std::string, unsigned int>>;

struct FactorKeyLess {
    [[nodiscard]] bool operator()(const FactorKey& lhs, const FactorKey& rhs) const noexcept {
        return lhs < rhs;
    }
};

using CoeffMap = std::map<FactorKey, BigInt, FactorKeyLess>;
using Monomial = std::vector<unsigned int>;
using SparsePoly = std::map<Monomial, BigInt>;

[[nodiscard]] FactorKey make_factor_key(const std::vector<std::pair<Symbol, unsigned int>>& factors,
                                        const std::optional<std::string>& omit_var = std::nullopt) {
    FactorKey key;
    key.reserve(factors.size());
    for (const auto& [symbol, exponent] : factors) {
        if (exponent == 0U) {
            continue;
        }
        if (omit_var.has_value() && symbol.name == omit_var.value()) {
            continue;
        }
        key.emplace_back(symbol.name, exponent);
    }

    std::sort(key.begin(), key.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });

    FactorKey merged;
    merged.reserve(key.size());
    for (const auto& [name, exponent] : key) {
        if (!merged.empty() && merged.back().first == name) {
            merged.back().second += exponent;
        } else {
            merged.emplace_back(name, exponent);
        }
    }
    return merged;
}

[[nodiscard]] std::vector<std::pair<Symbol, unsigned int>> key_to_factors(const FactorKey& key) {
    std::vector<std::pair<Symbol, unsigned int>> factors;
    factors.reserve(key.size());
    for (const auto& [name, exponent] : key) {
        if (exponent > 0U) {
            factors.emplace_back(Symbol(name), exponent);
        }
    }
    return factors;
}

[[nodiscard]] std::vector<Symbol> collect_all_variables(const MultivariatePolynomial& p,
                                                        const MultivariatePolynomial& q) {
    auto vars_p = p.variables();
    auto vars_q = q.variables();

    std::vector<Symbol> all_vars = vars_p;
    for (const auto& symbol : vars_q) {
        if (std::find_if(all_vars.begin(), all_vars.end(),
                         [&](const Symbol& candidate) { return candidate.name == symbol.name; }) == all_vars.end()) {
            all_vars.push_back(symbol);
        }
    }

    std::sort(all_vars.begin(), all_vars.end(), [](const Symbol& lhs, const Symbol& rhs) {
        return lhs.name < rhs.name;
    });
    return all_vars;
}

[[nodiscard]] std::size_t degree_in_var(const MultivariatePolynomial& poly, const Symbol& var) {
    std::size_t degree = 0U;
    for (const auto& term : poly.terms()) {
        for (const auto& [symbol, exponent] : term.factors) {
            if (symbol.name == var.name) {
                degree = std::max(degree, static_cast<std::size_t>(exponent));
            }
        }
    }
    return degree;
}

[[nodiscard]] CoeffMap to_coeff_map(const MultivariatePolynomial& poly,
                                    const std::optional<std::string>& omit_var = std::nullopt) {
    CoeffMap coefficients;
    for (const auto& term : poly.terms()) {
        coefficients[make_factor_key(term.factors, omit_var)] += term.coefficient;
    }

    for (auto it = coefficients.begin(); it != coefficients.end();) {
        if (it->second.is_zero()) {
            it = coefficients.erase(it);
        } else {
            ++it;
        }
    }
    return coefficients;
}

[[nodiscard]] bool same_polynomial(const MultivariatePolynomial& lhs, const MultivariatePolynomial& rhs) {
    return to_coeff_map(lhs) == to_coeff_map(rhs);
}

[[nodiscard]] MultivariatePolynomial multiply_by_scalar(const MultivariatePolynomial& poly, const BigInt& scalar) {
    if (poly.is_zero() || scalar.is_zero()) {
        return MultivariatePolynomial{};
    }

    std::vector<MultivariateTerm> terms;
    terms.reserve(poly.terms().size());
    for (const auto& term : poly.terms()) {
        terms.push_back(MultivariateTerm{
            .coefficient = term.coefficient * scalar,
            .factors = term.factors,
        });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] MultivariatePolynomial normalize_multivariate_gcd(const MultivariatePolynomial& poly) {
    if (poly.is_zero()) {
        return poly;
    }

    BigInt content(0);
    for (const auto& term : poly.terms()) {
        content = gcd(content, term.coefficient.abs());
    }

    std::vector<MultivariateTerm> normalized;
    normalized.reserve(poly.terms().size());
    for (const auto& term : poly.terms()) {
        MultivariateTerm reduced = term;
        if (content > BigInt(1)) {
            reduced.coefficient /= content;
        }
        normalized.push_back(std::move(reduced));
    }

    MultivariatePolynomial result(std::move(normalized));
    CoeffMap coeffs = to_coeff_map(result);
    if (!coeffs.empty() && std::prev(coeffs.end())->second.is_negative()) {
        result = multiply_by_scalar(result, BigInt(-1));
    }

    return result;
}

[[nodiscard]] bool is_unit_polynomial(const MultivariatePolynomial& poly) {
    if (poly.terms().size() != 1U) {
        return false;
    }
    const auto& term = poly.terms().front();
    return term.factors.empty() && term.coefficient.abs() == BigInt(1);
}

[[nodiscard]] Result<IntPoly> multivariate_single_var_to_intpoly(const MultivariatePolynomial& poly, const Symbol& var) {
    if (poly.is_zero()) {
        return ok(IntPoly{});
    }

    std::map<std::size_t, BigInt> coefficient_map;
    std::size_t max_degree = 0U;

    for (const auto& term : poly.terms()) {
        std::size_t degree = 0U;
        for (const auto& [symbol, exponent] : term.factors) {
            if (symbol.name == var.name) {
                degree = exponent;
            } else {
                return fail<IntPoly>(make_error(
                    CASErrorKind::Unimplemented,
                    "polynomial_gcd_multivariate: expected univariate polynomial after specialization"));
            }
        }
        coefficient_map[degree] += term.coefficient;
        max_degree = std::max(max_degree, degree);
    }

    IntPoly result;
    result.resize(max_degree + 1U, BigInt(0));
    for (const auto& [degree, coefficient] : coefficient_map) {
        result[degree] = coefficient;
    }
    result.normalize([](const BigInt& value) { return value.is_zero(); });
    return ok(std::move(result));
}

[[nodiscard]] MultivariatePolynomial intpoly_to_multivariate(const IntPoly& poly, const Symbol& var) {
    std::vector<MultivariateTerm> terms;
    terms.reserve(poly.size());
    for (std::size_t degree = 0; degree < poly.size(); ++degree) {
        if (poly[degree].is_zero()) {
            continue;
        }
        std::vector<std::pair<Symbol, unsigned int>> factors;
        if (degree > 0U) {
            factors.emplace_back(var, static_cast<unsigned int>(degree));
        }
        terms.push_back(MultivariateTerm{
            .coefficient = poly[degree],
            .factors = std::move(factors),
        });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] Result<std::vector<Rational>> lagrange_interpolate(
    const std::vector<BigInt>& values,
    const std::vector<BigInt>& points) {
    if (values.size() != points.size()) {
        return fail<std::vector<Rational>>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd_multivariate: interpolation points/value size mismatch"));
    }
    if (values.empty()) {
        return ok(std::vector<Rational>{});
    }

    const std::size_t n = values.size();
    std::vector<Rational> result(n, Rational(BigInt(0)));

    for (std::size_t i = 0U; i < n; ++i) {
        if (values[i].is_zero()) {
            continue;
        }

        std::vector<Rational> basis(1U, Rational(BigInt(1)));
        for (std::size_t j = 0U; j < n; ++j) {
            if (j == i) {
                continue;
            }

            const BigInt denominator = points[i] - points[j];
            if (denominator.is_zero()) {
                return fail<std::vector<Rational>>(make_error(
                    CASErrorKind::InvalidArgument,
                    "polynomial_gcd_multivariate: duplicate interpolation points"));
            }

            const Rational scale(BigInt(1), denominator);
            std::vector<Rational> next(basis.size() + 1U, Rational(BigInt(0)));
            for (std::size_t k = 0U; k < basis.size(); ++k) {
                next[k + 1U] = next[k + 1U] + basis[k] * scale;
                next[k] = next[k] - basis[k] * scale * Rational(points[j]);
            }
            basis = std::move(next);
        }

        const Rational value(values[i]);
        if (basis.size() > result.size()) {
            result.resize(basis.size(), Rational(BigInt(0)));
        }
        for (std::size_t k = 0U; k < basis.size(); ++k) {
            result[k] = result[k] + value * basis[k];
        }
    }

    while (!result.empty() && result.back().numerator().is_zero()) {
        result.pop_back();
    }
    return ok(std::move(result));
}

[[nodiscard]] MultivariatePolynomial coefficient_poly_in_var(
    const MultivariatePolynomial& poly,
    const Symbol& var,
    std::size_t degree) {
    std::vector<MultivariateTerm> terms;
    for (const auto& term : poly.terms()) {
        unsigned int exponent = 0U;
        for (const auto& [symbol, factor_exponent] : term.factors) {
            if (symbol.name == var.name) {
                exponent = factor_exponent;
                break;
            }
        }

        if (static_cast<std::size_t>(exponent) != degree) {
            continue;
        }

        std::vector<std::pair<Symbol, unsigned int>> reduced_factors;
        reduced_factors.reserve(term.factors.size());
        for (const auto& [symbol, factor_exponent] : term.factors) {
            if (symbol.name != var.name) {
                reduced_factors.emplace_back(symbol, factor_exponent);
            }
        }

        terms.push_back(MultivariateTerm{
            .coefficient = term.coefficient,
            .factors = std::move(reduced_factors),
        });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] MultivariatePolynomial multiply_by_variable_power(
    const MultivariatePolynomial& poly,
    const Symbol& var,
    std::size_t power) {
    if (power == 0U || poly.is_zero()) {
        return poly;
    }

    std::vector<MultivariateTerm> terms;
    terms.reserve(poly.terms().size());
    for (const auto& term : poly.terms()) {
        auto factors = term.factors;
        factors.emplace_back(var, static_cast<unsigned int>(power));
        terms.push_back(MultivariateTerm{
            .coefficient = term.coefficient,
            .factors = std::move(factors),
        });
    }
    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] Result<MultivariatePolynomial> interpolate_polynomial_values(
    const std::vector<MultivariatePolynomial>& values,
    const std::vector<BigInt>& points,
    const Symbol& interpolation_var) {
    if (values.size() != points.size()) {
        return fail<MultivariatePolynomial>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd_multivariate: interpolation polynomial points/value size mismatch"));
    }

    std::set<FactorKey, FactorKeyLess> all_keys;
    std::vector<CoeffMap> value_maps;
    value_maps.reserve(values.size());

    for (const auto& value_poly : values) {
        auto map = to_coeff_map(value_poly);
        for (const auto& [key, coefficient] : map) {
            static_cast<void>(coefficient);
            all_keys.insert(key);
        }
        value_maps.push_back(std::move(map));
    }

    std::vector<MultivariateTerm> terms;

    for (const auto& key : all_keys) {
        std::vector<BigInt> sequence(values.size(), BigInt(0));
        for (std::size_t i = 0U; i < value_maps.size(); ++i) {
            const auto it = value_maps[i].find(key);
            if (it != value_maps[i].end()) {
                sequence[i] = it->second;
            }
        }

        auto interpolation = lagrange_interpolate(sequence, points);
        if (interpolation.is_error()) {
            return fail<MultivariatePolynomial>(interpolation.error());
        }

        const auto& coeffs = interpolation.value();
        for (std::size_t degree = 0U; degree < coeffs.size(); ++degree) {
            const Rational& coefficient = coeffs[degree];
            if (coefficient.numerator().is_zero()) {
                continue;
            }
            if (coefficient.denominator() != BigInt(1)) {
                return fail<MultivariatePolynomial>(make_error(
                    CASErrorKind::InternalError,
                    "polynomial_gcd_multivariate: non-integer coefficient after interpolation"));
            }

            auto factors = key_to_factors(key);
            if (degree > 0U) {
                factors.emplace_back(interpolation_var, static_cast<unsigned int>(degree));
            }

            terms.push_back(MultivariateTerm{
                .coefficient = coefficient.numerator(),
                .factors = std::move(factors),
            });
        }
    }

    return ok(MultivariatePolynomial(std::move(terms)));
}

[[nodiscard]] SparsePoly to_sparse(const MultivariatePolynomial& poly, const std::vector<Symbol>& vars) {
    SparsePoly sparse;
    std::map<std::string, std::size_t> indices;
    for (std::size_t i = 0U; i < vars.size(); ++i) {
        indices[vars[i].name] = i;
    }

    for (const auto& term : poly.terms()) {
        Monomial monomial(vars.size(), 0U);
        for (const auto& [symbol, exponent] : term.factors) {
            const auto it = indices.find(symbol.name);
            if (it != indices.end()) {
                monomial[it->second] += exponent;
            }
        }
        sparse[monomial] += term.coefficient;
    }

    for (auto it = sparse.begin(); it != sparse.end();) {
        if (it->second.is_zero()) {
            it = sparse.erase(it);
        } else {
            ++it;
        }
    }

    return sparse;
}

[[nodiscard]] bool monomial_divides(const Monomial& divisor, const Monomial& dividend) {
    for (std::size_t i = 0U; i < divisor.size(); ++i) {
        if (divisor[i] > dividend[i]) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] Monomial monomial_quotient(const Monomial& dividend, const Monomial& divisor) {
    Monomial quotient(dividend.size(), 0U);
    for (std::size_t i = 0U; i < dividend.size(); ++i) {
        quotient[i] = dividend[i] - divisor[i];
    }
    return quotient;
}

void sparse_add_term(SparsePoly& poly, const Monomial& monomial, const BigInt& coefficient) {
    if (coefficient.is_zero()) {
        return;
    }

    poly[monomial] += coefficient;
    if (poly[monomial].is_zero()) {
        poly.erase(monomial);
    }
}

[[nodiscard]] SparsePoly multiply_sparse_by_term(
    const SparsePoly& poly,
    const Monomial& monomial,
    const BigInt& coefficient) {
    SparsePoly result;
    if (coefficient.is_zero()) {
        return result;
    }

    for (const auto& [base_monomial, base_coefficient] : poly) {
        Monomial product(base_monomial.size(), 0U);
        for (std::size_t i = 0U; i < base_monomial.size(); ++i) {
            product[i] = base_monomial[i] + monomial[i];
        }
        sparse_add_term(result, product, base_coefficient * coefficient);
    }
    return result;
}

void sparse_subtract(SparsePoly& lhs, const SparsePoly& rhs) {
    for (const auto& [monomial, coefficient] : rhs) {
        sparse_add_term(lhs, monomial, -coefficient);
    }
}

[[nodiscard]] MultivariatePolynomial sparse_to_multivariate(const SparsePoly& sparse, const std::vector<Symbol>& vars) {
    std::vector<MultivariateTerm> terms;
    terms.reserve(sparse.size());

    for (const auto& [monomial, coefficient] : sparse) {
        if (coefficient.is_zero()) {
            continue;
        }

        std::vector<std::pair<Symbol, unsigned int>> factors;
        factors.reserve(vars.size());
        for (std::size_t i = 0U; i < vars.size(); ++i) {
            if (monomial[i] > 0U) {
                factors.emplace_back(vars[i], monomial[i]);
            }
        }

        terms.push_back(MultivariateTerm{
            .coefficient = coefficient,
            .factors = std::move(factors),
        });
    }

    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] Result<std::optional<MultivariatePolynomial>> exact_quotient(
    const MultivariatePolynomial& dividend,
    const MultivariatePolynomial& divisor,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx) {
    if (divisor.is_zero()) {
        return fail<std::optional<MultivariatePolynomial>>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd_multivariate: exact division by zero divisor"));
    }

    SparsePoly remainder = to_sparse(dividend, vars);
    SparsePoly divisor_sparse = to_sparse(divisor, vars);
    SparsePoly quotient;

    if (divisor_sparse.empty()) {
        return fail<std::optional<MultivariatePolynomial>>(make_error(
            CASErrorKind::InvalidArgument,
            "polynomial_gcd_multivariate: zero sparse divisor"));
    }

    const auto [divisor_lm, divisor_lc] = *std::prev(divisor_sparse.end());
    const std::size_t max_steps =
        std::max(ctx.min_gcd_division_steps(), (remainder.size() + 1U) * (vars.size() + 1U) * 16U);
    std::size_t steps = 0U;

    while (!remainder.empty()) {
        if (++steps > max_steps) {
            return fail<std::optional<MultivariatePolynomial>>(make_error(
                CASErrorKind::Unimplemented,
                "polynomial_gcd_multivariate: exact division budget exceeded during certification"));
        }

        const auto [remainder_lm, remainder_lc] = *std::prev(remainder.end());

        if (!monomial_divides(divisor_lm, remainder_lm)) {
            return ok(std::optional<MultivariatePolynomial>{std::nullopt});
        }
        if ((remainder_lc % divisor_lc) != BigInt(0)) {
            return ok(std::optional<MultivariatePolynomial>{std::nullopt});
        }

        const BigInt term_coefficient = remainder_lc / divisor_lc;
        const Monomial term_monomial = monomial_quotient(remainder_lm, divisor_lm);

        sparse_add_term(quotient, term_monomial, term_coefficient);

        SparsePoly subtraction = multiply_sparse_by_term(divisor_sparse, term_monomial, term_coefficient);
        sparse_subtract(remainder, subtraction);

        if (!remainder.empty() && !(std::prev(remainder.end())->first < remainder_lm)) {
            return fail<std::optional<MultivariatePolynomial>>(make_error(
                CASErrorKind::Unimplemented,
                "polynomial_gcd_multivariate: exact division did not reduce leading monomial"));
        }
    }

    return ok(std::optional<MultivariatePolynomial>{sparse_to_multivariate(quotient, vars)});
}

[[nodiscard]] BigInt gcd_abs_values(const std::vector<int>& values) {
    BigInt result(0);
    for (int value : values) {
        result = gcd(result, BigInt(value).abs());
    }
    return result;
}

[[nodiscard]] std::optional<MultivariatePolynomial> make_primitive_linear_candidate(
    const std::vector<Symbol>& vars,
    const std::vector<int>& coefficients,
    int constant) {
    bool has_variable_term = false;
    std::vector<int> all_values = coefficients;
    all_values.push_back(constant);
    const BigInt content = gcd_abs_values(all_values);
    if (content != BigInt(1)) {
        return std::nullopt;
    }

    int sign = 1;
    for (int coefficient : coefficients) {
        if (coefficient != 0) {
            sign = coefficient < 0 ? -1 : 1;
            break;
        }
    }

    std::vector<MultivariateTerm> terms;
    for (std::size_t i = 0U; i < vars.size(); ++i) {
        int coefficient = coefficients[i] * sign;
        if (coefficient == 0) {
            continue;
        }
        has_variable_term = true;
        terms.push_back(MultivariateTerm{
            .coefficient = BigInt(coefficient),
            .factors = {{vars[i], 1U}},
        });
    }

    const int normalized_constant = constant * sign;
    if (normalized_constant != 0) {
        terms.push_back(MultivariateTerm{
            .coefficient = BigInt(normalized_constant),
            .factors = {},
        });
    }

    if (!has_variable_term) {
        return std::nullopt;
    }

    return MultivariatePolynomial(std::move(terms));
}

[[nodiscard]] std::vector<MultivariatePolynomial> primitive_linear_candidates(const std::vector<Symbol>& vars) {
    if (vars.empty() || vars.size() > 4U) {
        return {};
    }

    std::vector<MultivariatePolynomial> candidates;
    std::set<CoeffMap> seen;
    std::vector<int> coefficients(vars.size(), 0);

    auto emit = [&](int constant) {
        auto candidate = make_primitive_linear_candidate(vars, coefficients, constant);
        if (!candidate.has_value()) {
            return;
        }
        CoeffMap key = to_coeff_map(candidate.value());
        if (seen.insert(key).second) {
            candidates.push_back(std::move(candidate.value()));
        }
    };

    auto walk = [&](auto&& self, std::size_t index) -> void {
        if (index == vars.size()) {
            for (int constant = -3; constant <= 3; ++constant) {
                emit(constant);
            }
            return;
        }
        for (int coefficient = -2; coefficient <= 2; ++coefficient) {
            coefficients[index] = coefficient;
            self(self, index + 1U);
        }
    };

    walk(walk, 0U);
    std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.terms().size() != rhs.terms().size()) {
            return lhs.terms().size() < rhs.terms().size();
        }
        return to_coeff_map(lhs) < to_coeff_map(rhs);
    });
    return candidates;
}

[[nodiscard]] Result<std::optional<MultivariatePolynomial>> try_certified_linear_gcd(
    const MultivariatePolynomial& p,
    const MultivariatePolynomial& q,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx,
    std::size_t depth);

[[nodiscard]] Result<MultivariatePolynomial> gcd_multivariate_recursive(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx,
    std::size_t depth) {
    if (depth > ctx.max_gcd_recursion_depth()) {
        return fail<MultivariatePolynomial>(make_error(
            CASErrorKind::Unimplemented,
            "polynomial_gcd_multivariate: recursion depth limit reached while certifying gcd"));
    }

    MultivariatePolynomial p = normalize_multivariate_gcd(P);
    MultivariatePolynomial q = normalize_multivariate_gcd(Q);

    if (p.is_zero()) {
        return ok(q);
    }
    if (q.is_zero()) {
        return ok(p);
    }

    const std::vector<Symbol> vars = collect_all_variables(p, q);
    if (vars.empty()) {
        BigInt gcd_value = gcd(p.terms().front().coefficient.abs(), q.terms().front().coefficient.abs());
        return ok(MultivariatePolynomial({MultivariateTerm{.coefficient = gcd_value, .factors = {}}}));
    }

    if (vars.size() == 1U) {
        auto p_int = multivariate_single_var_to_intpoly(p, vars[0]);
        if (p_int.is_error()) {
            return fail<MultivariatePolynomial>(p_int.error());
        }
        auto q_int = multivariate_single_var_to_intpoly(q, vars[0]);
        if (q_int.is_error()) {
            return fail<MultivariatePolynomial>(q_int.error());
        }

        IntPoly gcd_int = gcd_integer_poly_primitive(std::move(p_int.value()), std::move(q_int.value()));
        if (!gcd_int.empty() && gcd_int.leading_coeff().is_negative()) {
            for (auto& coefficient : gcd_int.coefficients()) {
                coefficient = -coefficient;
            }
        }
        return ok(intpoly_to_multivariate(gcd_int, vars[0]));
    }

    // Identical-polynomial shortcut: covers gcd(p,p)=p for all variable counts,
    // including linear trivariate cases where the degree≤1 guard below is absent.
    if (same_polynomial(p, q)) {
        return ok(normalize_multivariate_gcd(p));
    }

    // Linear-candidate fast-path: only worthwhile for ≤2 variables.
    // For 3+ variables the candidate count grows as 5^n which makes the scan
    // more expensive than the eval-interpolation path it is meant to accelerate.
    if (vars.size() <= 2U) {
        auto linear_gcd = try_certified_linear_gcd(p, q, vars, ctx, depth);
        if (linear_gcd.is_error()) {
            return fail<MultivariatePolynomial>(linear_gcd.error());
        }
        if (linear_gcd.value().has_value()) {
            return ok(linear_gcd.value().value());
        }
        // Bivariate termination: try_certified_linear_gcd calls gcd_multivariate_recursive
        // for cofactor GCDs, which would loop forever if we fell through to interpolation
        // on degree-1 coprime inputs (e.g., gcd(x-y, x+y) = 1).  The identical-polynomial
        // check above already handled p==q; here p≠q so returning 1 is correct.
        if (std::min(p.total_degree(), q.total_degree()) <= 1U) {
            return ok(MultivariatePolynomial({MultivariateTerm{.coefficient = BigInt(1), .factors = {}}}));
        }
    }

    const Symbol& main_var = vars.front();
    const Symbol& interpolation_var = vars.back();

    const std::size_t interpolation_degree_bound =
        std::min(degree_in_var(p, interpolation_var), degree_in_var(q, interpolation_var));
    const std::size_t required_samples = std::max<std::size_t>(interpolation_degree_bound + 1U, 2U);
    // Schwartz-Zippel bound: lc(GCD, interp_var) vanishes at at most D = interpolation_degree_bound
    // values in any evaluation domain.  Therefore 2D+3 samples always contain at least D+2 lucky
    // ones — enough for unique interpolation + one extra for the bucket selection heuristic.
    // This replaces the previous O(log(1/δ)) safety margin which caused O(N^k) fan-out blowup for
    // polynomials in k≥3 variables.
    const std::size_t max_samples = 2U * interpolation_degree_bound + 3U;

    struct SamplePoint {
        BigInt value;
        MultivariatePolynomial gcd;
        std::size_t main_degree{0U};
    };

    const std::array<long long, 3> offsets = {0LL, 17LL, 43LL};

    for (long long offset : offsets) {
        std::vector<SamplePoint> samples;
        samples.reserve(max_samples);

        bool sampling_ok = true;
        for (std::size_t i = 0U; i < max_samples; ++i) {
            const BigInt point(static_cast<long long>(i) + 1LL + offset);
            ExprPtr value_expr = ctx.arena().make<IntegerLit>(point);

            auto p_eval = p.evaluate_at(interpolation_var, value_expr);
            if (p_eval.is_error()) {
                sampling_ok = false;
                break;
            }
            auto q_eval = q.evaluate_at(interpolation_var, value_expr);
            if (q_eval.is_error()) {
                sampling_ok = false;
                break;
            }

            auto gcd_eval = gcd_multivariate_recursive(p_eval.value(), q_eval.value(), ctx, depth + 1U);
            if (gcd_eval.is_error()) {
                sampling_ok = false;
                break;
            }

            MultivariatePolynomial normalized_eval = normalize_multivariate_gcd(gcd_eval.value());
            samples.push_back(SamplePoint{
                .value = point,
                .gcd = std::move(normalized_eval),
                .main_degree = degree_in_var(samples.empty() ? MultivariatePolynomial{} : samples.back().gcd, main_var),
            });
            samples.back().main_degree = degree_in_var(samples.back().gcd, main_var);
        }

        if (!sampling_ok) {
            continue;
        }

        std::map<std::size_t, std::vector<SamplePoint>> buckets;
        for (const auto& sample : samples) {
            buckets[sample.main_degree].push_back(sample);
        }

        std::optional<std::size_t> target_degree;
        for (const auto& [degree, bucket] : buckets) {
            if (bucket.size() >= required_samples) {
                if (!target_degree.has_value() || degree > target_degree.value()) {
                    target_degree = degree;
                }
            }
        }

        if (!target_degree.has_value()) {
            continue;
        }

        const auto& chosen_bucket = buckets[target_degree.value()];
        std::vector<BigInt> chosen_points;
        std::vector<MultivariatePolynomial> chosen_polys;
        chosen_points.reserve(required_samples);
        chosen_polys.reserve(required_samples);

        for (std::size_t i = 0U; i < required_samples; ++i) {
            chosen_points.push_back(chosen_bucket[i].value);
            chosen_polys.push_back(chosen_bucket[i].gcd);
        }

        std::vector<MultivariateTerm> candidate_terms;

        bool interpolation_ok = true;
        for (std::size_t degree = 0U; degree <= target_degree.value(); ++degree) {
            std::vector<MultivariatePolynomial> coefficient_values;
            coefficient_values.reserve(chosen_polys.size());
            for (const auto& sample_poly : chosen_polys) {
                coefficient_values.push_back(coefficient_poly_in_var(sample_poly, main_var, degree));
            }

            auto interpolation = interpolate_polynomial_values(coefficient_values, chosen_points, interpolation_var);
            if (interpolation.is_error()) {
                interpolation_ok = false;
                break;
            }

            MultivariatePolynomial lifted = multiply_by_variable_power(interpolation.value(), main_var, degree);
            for (const auto& term : lifted.terms()) {
                candidate_terms.push_back(term);
            }
        }

        if (!interpolation_ok) {
            continue;
        }

        MultivariatePolynomial candidate = normalize_multivariate_gcd(MultivariatePolynomial(std::move(candidate_terms)));
        if (candidate.is_zero()) {
            continue;
        }

        auto quotient_p = exact_quotient(p, candidate, vars, ctx);
        if (quotient_p.is_error() || !quotient_p.value().has_value()) {
            continue;
        }
        auto quotient_q = exact_quotient(q, candidate, vars, ctx);
        if (quotient_q.is_error() || !quotient_q.value().has_value()) {
            continue;
        }

        MultivariatePolynomial certified_gcd = candidate;
        MultivariatePolynomial cofactor_p = quotient_p.value().value();
        MultivariatePolynomial cofactor_q = quotient_q.value().value();

        bool certified = false;
        for (std::size_t refine = 0U; refine < 3U; ++refine) {
            auto cofactor_gcd = gcd_multivariate_recursive(cofactor_p, cofactor_q, ctx, depth + 1U);
            if (cofactor_gcd.is_error()) {
                certified = false;
                break;
            }

            MultivariatePolynomial extra_factor = normalize_multivariate_gcd(cofactor_gcd.value());
            if (is_unit_polynomial(extra_factor)) {
                certified = true;
                break;
            }

            MultivariatePolynomial next_gcd = normalize_multivariate_gcd(certified_gcd * extra_factor);
            auto next_quotient_p = exact_quotient(p, next_gcd, vars, ctx);
            auto next_quotient_q = exact_quotient(q, next_gcd, vars, ctx);
            if (next_quotient_p.is_error() || next_quotient_q.is_error() ||
                !next_quotient_p.value().has_value() || !next_quotient_q.value().has_value()) {
                certified = false;
                break;
            }

            if (same_polynomial(cofactor_p, next_quotient_p.value().value()) &&
                same_polynomial(cofactor_q, next_quotient_q.value().value())) {
                certified = false;
                break;
            }

            certified_gcd = std::move(next_gcd);
            cofactor_p = std::move(next_quotient_p.value().value());
            cofactor_q = std::move(next_quotient_q.value().value());
        }

        if (certified) {
            return ok(normalize_multivariate_gcd(certified_gcd));
        }
    }

    return fail<MultivariatePolynomial>(make_error(
        CASErrorKind::Unimplemented,
        "polynomial_gcd_multivariate: unable to certify multivariate gcd for this input"));
}

[[nodiscard]] Result<std::optional<MultivariatePolynomial>> try_certified_linear_gcd(
    const MultivariatePolynomial& p,
    const MultivariatePolynomial& q,
    const std::vector<Symbol>& vars,
    symbolic::CASContext& ctx,
    std::size_t depth) {
    for (const auto& candidate : primitive_linear_candidates(vars)) {
        auto quotient_p = exact_quotient(p, candidate, vars, ctx);
        if (quotient_p.is_error()) {
            continue;
        }
        if (!quotient_p.value().has_value()) {
            continue;
        }

        auto quotient_q = exact_quotient(q, candidate, vars, ctx);
        if (quotient_q.is_error()) {
            continue;
        }
        if (!quotient_q.value().has_value()) {
            continue;
        }

        auto cofactor_gcd = gcd_multivariate_recursive(
            quotient_p.value().value(),
            quotient_q.value().value(),
            ctx,
            depth + 1U);
        if (cofactor_gcd.is_error()) {
            return fail<std::optional<MultivariatePolynomial>>(cofactor_gcd.error());
        }

        MultivariatePolynomial certified = normalize_multivariate_gcd(candidate * cofactor_gcd.value());
        return ok(std::optional<MultivariatePolynomial>{std::move(certified)});
    }

    return ok(std::optional<MultivariatePolynomial>{std::nullopt});
}

} // namespace

Result<MultivariatePolynomial> gcd_multivariate_eval_interp(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    symbolic::CASContext& ctx) {
    return gcd_multivariate_recursive(P, Q, ctx, 0U);
}

Result<ExprPtr> polynomial_gcd_multivariate(ExprPtr p, ExprPtr q, symbolic::CASContext& ctx) {
    auto expanded_p = expand_expr_impl(p, ctx);
    if (expanded_p.is_error()) {
        return fail<ExprPtr>(expanded_p.error());
    }

    auto expanded_q = expand_expr_impl(q, ctx);
    if (expanded_q.is_error()) {
        return fail<ExprPtr>(expanded_q.error());
    }

    auto P = parse_multivariate_polynomial(expanded_p.value(), ctx);
    if (P.is_error()) {
        return fail<ExprPtr>(P.error());
    }

    auto Q = parse_multivariate_polynomial(expanded_q.value(), ctx);
    if (Q.is_error()) {
        return fail<ExprPtr>(Q.error());
    }

    auto G = gcd_multivariate_eval_interp(P.value(), Q.value(), ctx);
    if (G.is_error()) {
        return fail<ExprPtr>(G.error());
    }

    auto gcd_expr = multivariate_to_expr(G.value(), ctx);
    if (gcd_expr.is_error()) {
        return fail<ExprPtr>(gcd_expr.error());
    }

    return ok(gcd_expr.value());
}

} // namespace cas::algebra
