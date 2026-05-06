#pragma once

#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include <map>
#include <vector>

namespace cas {
namespace symbolic {

struct Monomial {
    std::vector<std::pair<ExprPtr, unsigned int>> factors;

    bool operator<(const Monomial& other) const {
        if (factors.size() != other.factors.size()) return factors.size() < other.factors.size();
        for (std::size_t i = 0; i < factors.size(); ++i) {
            int cmp = canonical_compare(factors[i].first, other.factors[i].first);
            if (cmp != 0) return cmp < 0;
            if (factors[i].second != other.factors[i].second) {
                return factors[i].second < other.factors[i].second;
            }
        }
        return false;
    }

    bool operator==(const Monomial& other) const {
        if (factors.size() != other.factors.size()) return false;
        for (std::size_t i = 0; i < factors.size(); ++i) {
            if (canonical_compare(factors[i].first, other.factors[i].first) != 0) return false;
            if (factors[i].second != other.factors[i].second) return false;
        }
        return true;
    }
};

[[nodiscard]] Result<std::map<Monomial, Rational>> collect_polynomial_terms(ExprPtr expr, CASContext& ctx);

[[nodiscard]] Result<ExprPtr> polynomial_normal_form(ExprPtr expr, CASContext& ctx);

[[nodiscard]] Result<ExprPtr> transcendental_normal_form(ExprPtr expr, CASContext& ctx);

} // namespace symbolic
} // namespace cas
