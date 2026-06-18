// assumptions_core.cpp — Assumptions mutators, consistency, range/domain getters.
// Split from assumptions.cpp (W9.3 anti-monolith). Predicate/prover logic lives
// in assumptions_predicates.cpp. Shared scalar helpers are declared in
// symbolic_internal.hpp and defined here (single definition for the whole engine).
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"
#include <unordered_set>

namespace cas::symbolic {

[[nodiscard]] bool is_zero_expr(ExprPtr expr) {
    if (!expr) return false;
    auto scalar = exact_scalar_from_expr(expr);
    return scalar.has_value() && scalar->numerator().is_zero();
}

[[nodiscard]] bool is_positive_scalar(ExprPtr expr) {
    auto scalar = exact_scalar_from_expr(expr);
    return scalar.has_value() && scalar->numerator() > BigInt(0);
}

[[nodiscard]] bool is_negative_scalar(ExprPtr expr) {
    auto scalar = exact_scalar_from_expr(expr);
    return scalar.has_value() && scalar->numerator() < BigInt(0);
}

void Assumptions::assume_domain(const Symbol& symbol, Domain domain) {
    ++revision_;  // F7.0-A4.1: signal change for cache invalidation
    symbol_domains_[symbol.name] = domain;
    switch (domain) {
        case Domain::Positive:
            assume_positive(symbol);
            break;
        case Domain::Negative:
            real_symbols_.insert(symbol.name);
            negative_symbols_.insert(symbol.name);
            nonzero_symbols_.insert(symbol.name);
            break;
        case Domain::Natural:
            real_symbols_.insert(symbol.name);
            integer_symbols_.insert(symbol.name);
            assume_greater_equal(ExprPtr(&symbol), nullptr);
            break;
        case Domain::NonZero:
            assume_nonzero(symbol);
            break;
        case Domain::Integer:
            assume_integer(symbol);
            break;
        case Domain::Real:
            assume_real(symbol);
            break;
        case Domain::Rational:
            real_symbols_.insert(symbol.name);
            break;
        default:
            break;
    }
}

Domain Assumptions::get_domain(const Symbol& symbol) const {
    auto it = symbol_domains_.find(symbol.name);
    if (it != symbol_domains_.end()) return it->second;

    if (positive_symbols_.contains(symbol.name)) return Domain::Positive;
    if (negative_symbols_.contains(symbol.name)) return Domain::Negative;
    if (integer_symbols_.contains(symbol.name)) return Domain::Integer;
    if (real_symbols_.contains(symbol.name)) return Domain::Real;
    if (nonzero_symbols_.contains(symbol.name)) return Domain::NonZero;

    return Domain::Complex;
}

void Assumptions::assume_positive(const Symbol& symbol) {
    ++revision_;  // F7.0-A4.1
    real_symbols_.insert(symbol.name);
    positive_symbols_.insert(symbol.name);
    nonzero_symbols_.insert(symbol.name);
}

Result<void> Assumptions::check_consistency() const {
    // 1. Check Range vs Domain
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

    // 2. Check Positive vs Negative
    for (const auto& name : positive_symbols_) {
        if (negative_symbols_.contains(name)) {
            return fail<void>(make_error(
                CASErrorKind::InvalidArgument,
                "Symbol '" + name + "' cannot be both positive and negative"));
        }
    }

    // 3. Check for x > x or x < x in relation graph
    for (auto const& [node, rels] : relations_) {
        std::unordered_set<const ExprNode*> visited;
        if (prove_relation(node, node, true, visited)) {
            return fail<void>(make_error(
                CASErrorKind::InvalidArgument,
                "Contradiction detected in relation graph (strict self-loop)"));
        }
    }

    // 4. Check x > y and y > x
    for (auto const& [lhs, rels] : relations_) {
        for (const auto& rel : rels) {
            if (rel.type == RelType::Less) { // lhs < rel.target
                std::unordered_set<const ExprNode*> visited;
                if (prove_relation(rel.target, lhs, false, visited)) {
                    return fail<void>(make_error(
                        CASErrorKind::InvalidArgument,
                        "Contradiction: cyclic strict inequality"));
                }
            }
        }
    }

    return ok();
}

void Assumptions::assume_real(const Symbol& symbol) {
    ++revision_;  // F7.0-A4.1
    real_symbols_.insert(symbol.name);
}

void Assumptions::assume_integer(const Symbol& symbol) {
    ++revision_;  // F7.0-A4.1
    real_symbols_.insert(symbol.name);
    integer_symbols_.insert(symbol.name);
}

void Assumptions::assume_nonzero(const Symbol& symbol) {
    ++revision_;  // F7.0-A4.1
    nonzero_symbols_.insert(symbol.name);
}

void Assumptions::assume_in_range(const Symbol& symbol, ExprPtr lower, ExprPtr upper) {
    ++revision_;  // F7.0-A4.1
    real_symbols_.insert(symbol.name);
    range_symbols_[symbol.name] = RangeAssumption{
        .lower = lower,
        .upper = upper,
    };
}

void Assumptions::assume_greater(ExprPtr lhs, ExprPtr rhs) {
    ++revision_;  // F7.0-A4.1
    // lhs > rhs  =>  rhs < lhs
    ExprPtr r = is_zero_expr(rhs) ? nullptr : rhs;
    ExprPtr l = is_zero_expr(lhs) ? nullptr : lhs;
    relations_[r].push_back({l, RelType::Less});
}

void Assumptions::assume_greater_equal(ExprPtr lhs, ExprPtr rhs) {
    ++revision_;  // F7.0-A4.1
    // lhs >= rhs  =>  rhs <= lhs
    ExprPtr r = is_zero_expr(rhs) ? nullptr : rhs;
    ExprPtr l = is_zero_expr(lhs) ? nullptr : lhs;
    relations_[r].push_back({l, RelType::LessEqual});
}

void Assumptions::assume(ExprPtr condition) {
    if (!condition) return;
    // Note: revision_ is bumped by the called assume_*/assume_greater_* mutators,
    // so we don't double-bump here.

    if (const auto* fc = expr_cast<FuncCall>(condition)) {
        const std::string& name = fc->name;
        if (name == "greater" && fc->args.size() == 2U) {
            assume_greater(fc->args[0], fc->args[1]);
        } else if (name == "greater_equal" && fc->args.size() == 2U) {
            assume_greater_equal(fc->args[0], fc->args[1]);
        } else if (name == "less" && fc->args.size() == 2U) {
            assume_greater(fc->args[1], fc->args[0]);
        } else if (name == "less_equal" && fc->args.size() == 2U) {
            assume_greater_equal(fc->args[1], fc->args[0]);
        } else if (name == "positive" && fc->args.size() == 1U) {
            if (const auto* sym = expr_cast<Symbol>(fc->args[0])) {
                assume_positive(*sym);
            } else {
                assume_greater(fc->args[0], nullptr);
            }
        } else if (name == "real" && fc->args.size() == 1U) {
            if (const auto* sym = expr_cast<Symbol>(fc->args[0]))
                assume_real(*sym);
        } else if (name == "nonzero" && fc->args.size() == 1U) {
            if (const auto* sym = expr_cast<Symbol>(fc->args[0]))
                assume_nonzero(*sym);
        } else if (name == "integer" && fc->args.size() == 1U) {
            if (const auto* sym = expr_cast<Symbol>(fc->args[0]))
                assume_integer(*sym);
        } else if ((name == "nonneg" || name == "nonnegative") && fc->args.size() == 1U) {
            assume_greater_equal(fc->args[0], nullptr);
        } else if (name == "negative" && fc->args.size() == 1U) {
            if (const auto* sym = expr_cast<Symbol>(fc->args[0]))
                assume_domain(*sym, Domain::Negative);
            else
                assume_greater(nullptr, fc->args[0]);
        }
        return;
    }

    if (const auto* bin = expr_cast<Binary>(condition)) {
        if (bin->op == BinaryOp::Equal) {
            assume_greater_equal(bin->left, bin->right);
            assume_greater_equal(bin->right, bin->left);
        } else if (bin->op == BinaryOp::Less) {
            assume_greater(bin->right, bin->left);
        } else if (bin->op == BinaryOp::Greater) {
            assume_greater(bin->left, bin->right);
        } else if (bin->op == BinaryOp::LessEqual) {
            assume_greater_equal(bin->right, bin->left);
        } else if (bin->op == BinaryOp::GreaterEqual) {
            assume_greater_equal(bin->left, bin->right);
        }
        return;
    }

    if (const auto* sym = expr_cast<Symbol>(condition)) {
        assume_nonzero(*sym);
    }
}

std::optional<RangeAssumption> Assumptions::get_range(const Symbol& symbol) const {
    const auto found = range_symbols_.find(symbol.name);
    if (found == range_symbols_.end()) return std::nullopt;
    return found->second;
}

void Assumptions::update_roots(AstArena& target, std::unordered_map<ExprPtr, ExprPtr>& cache) {
    for (auto& [name, range] : range_symbols_) {
        range.lower = clone_into_arena(range.lower, target, cache);
        range.upper = clone_into_arena(range.upper, target, cache);
    }

    std::unordered_map<ExprPtr, std::vector<Relation>, ExprHash, ExprEqual> new_relations;
    for (auto& [key, rels] : relations_) {
        ExprPtr new_key = key ? clone_into_arena(key, target, cache) : key;
        std::vector<Relation> new_rels;
        for (auto& rel : rels) {
            new_rels.push_back({clone_into_arena(rel.target, target, cache), rel.type});
        }
        new_relations[new_key] = std::move(new_rels);
    }
    relations_ = std::move(new_relations);
}

} // namespace cas::symbolic
