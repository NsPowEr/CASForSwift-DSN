#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"

namespace cas::symbolic {

[[nodiscard]] bool range_is_exact_zero(ExprPtr lower, ExprPtr upper) {
    const auto exact_lower = exact_scalar_from_expr(lower);
    const auto exact_upper = exact_scalar_from_expr(upper);
    if (!exact_lower.has_value() || !exact_upper.has_value()) {
        return false;
    }

    return exact_lower->numerator().is_zero() &&
        exact_upper->numerator().is_zero();
}

[[nodiscard]] bool exact_range_excludes_zero(ExprPtr lower, ExprPtr upper) {
    const auto exact_lower = exact_scalar_from_expr(lower);
    const auto exact_upper = exact_scalar_from_expr(upper);
    if (!exact_lower.has_value() || !exact_upper.has_value()) {
        return false;
    }

    return compare_exact_scalars(*exact_lower, Rational(BigInt(0))) > 0 ||
        compare_exact_scalars(*exact_upper, Rational(BigInt(0))) < 0;
}

void Assumptions::assume_positive(const Symbol& symbol) {
    real_symbols_.insert(symbol.name);
    positive_symbols_.insert(symbol.name);
    nonzero_symbols_.insert(symbol.name);
}

Result<void> Assumptions::check_consistency() const {
    for (const auto& [name, range] : range_symbols_) {
        if (positive_symbols_.contains(name)) {
            const auto exact_upper = exact_scalar_from_expr(range.upper);
            if (exact_upper.has_value() && 
                compare_exact_scalars(*exact_upper, Rational(BigInt(0))) <= 0) {
                return fail<void>(make_error(
                    CASErrorKind::InvalidArgument,
                    "Positive assumption for '" + name + "' conflicts with its upper bound"));
            }
        }
        
        if (nonzero_symbols_.contains(name)) {
            if (range_is_exact_zero(range.lower, range.upper)) {
                return fail<void>(make_error(
                    CASErrorKind::InvalidArgument,
                    "Nonzero assumption for '" + name + "' conflicts with a range fixed at zero"));
            }
        }
    }
    return ok();
}

void Assumptions::assume_real(const Symbol& symbol) {
    real_symbols_.insert(symbol.name);
}

void Assumptions::assume_integer(const Symbol& symbol) {
    real_symbols_.insert(symbol.name);
    integer_symbols_.insert(symbol.name);
}

void Assumptions::assume_nonzero(const Symbol& symbol) {
    nonzero_symbols_.insert(symbol.name);
}

void Assumptions::assume_in_range(const Symbol& symbol, ExprPtr lower, ExprPtr upper) {
    real_symbols_.insert(symbol.name);
    range_symbols_[symbol.name] = RangeAssumption{
        .lower = lower,
        .upper = upper,
    };
}

bool Assumptions::is_real(const Symbol& symbol) const {
    return real_symbols_.contains(symbol.name);
}

bool Assumptions::is_positive(const Symbol& symbol) const {
    return positive_symbols_.contains(symbol.name);
}

bool Assumptions::is_nonzero(const Symbol& symbol) const {
    return nonzero_symbols_.contains(symbol.name);
}

bool Assumptions::could_be_zero(const Symbol& symbol) const {
    if (is_nonzero(symbol)) {
        return false;
    }

    const auto found = range_symbols_.find(symbol.name);
    if (found == range_symbols_.end()) {
        return true;
    }

    return !exact_range_excludes_zero(found->second.lower, found->second.upper);
}

bool Assumptions::is_integer(const Symbol& symbol) const {
    return integer_symbols_.contains(symbol.name);
}

std::optional<RangeAssumption> Assumptions::get_range(const Symbol& symbol) const {
    const auto found = range_symbols_.find(symbol.name);
    if (found == range_symbols_.end()) {
        return std::nullopt;
    }
    return found->second;
}

} // namespace cas::symbolic
