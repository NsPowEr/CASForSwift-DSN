#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"
#include <unordered_set>

namespace cas::symbolic {

bool Assumptions::is_greater_equal(ExprPtr lhs, ExprPtr rhs) const {
    if (lhs == rhs) return true;

    if (is_zero_expr(rhs)) {
        if (const auto* lhs_symbol = expr_cast<Symbol>(lhs)) {
            if (positive_symbols_.contains(lhs_symbol->name)) return true;
        }
    }

    if (is_zero_expr(lhs)) {
        if (const auto* rhs_symbol = expr_cast<Symbol>(rhs)) {
            if (negative_symbols_.contains(rhs_symbol->name)) return true;
        }
    }

    // 1. Scalar comparison
    auto l_scalar = exact_scalar_from_expr(lhs);
    auto r_scalar = exact_scalar_from_expr(rhs);
    if (l_scalar && r_scalar) return compare_exact_scalars(*l_scalar, *r_scalar) >= 0;

    // 2. Direct graph proof
    ExprPtr l = is_zero_expr(lhs) ? nullptr : lhs;
    ExprPtr r = is_zero_expr(rhs) ? nullptr : rhs;
    if (l == r) return true;
    std::unordered_set<const ExprNode*> visited;
    if (prove_relation(r, l, false, visited)) return true;

    auto get_terms = [](ExprPtr e) -> std::vector<ExprPtr> {
        if (const auto* s = expr_cast<Sum>(e)) return s->terms;
        if (const auto* b = expr_cast<Binary>(e)) {
            if (b->op == BinaryOp::Add) return {b->left, b->right};
        }
        return {e};
    };

    // 3. Decompose sums: (a + b) >= (c + d) if a >= c and b >= d
    auto l_terms = get_terms(lhs);
    auto r_terms = get_terms(rhs);

    if (l_terms.size() > 1 || r_terms.size() > 1) {
        if (l_terms.size() == r_terms.size()) {
            bool all_ge = true;
            for (size_t i = 0; i < l_terms.size(); ++i) {
                if (!is_greater_equal(l_terms[i], r_terms[i])) {
                    all_ge = false;
                    break;
                }
            }
            if (all_ge) return true;
        }

        if (l_terms.size() > 1 && r_terms.size() == 1) {
            for (size_t i = 0; i < l_terms.size(); ++i) {
                bool others_nonnegative = true;
                for (size_t j = 0; j < l_terms.size(); ++j) {
                    if (i == j) continue;
                    if (!is_nonnegative(l_terms[j])) { others_nonnegative = false; break; }
                }
                if (others_nonnegative && is_greater_equal(l_terms[i], rhs)) return true;
            }
        }
    }

    if (is_greater(lhs, rhs)) return true;

    return false;
}

bool Assumptions::is_nonzero(const Symbol& symbol) const {
    if (nonzero_symbols_.contains(symbol.name)) return true;
    if (positive_symbols_.contains(symbol.name)) return true;
    if (negative_symbols_.contains(symbol.name)) return true;
    return false;
}

bool Assumptions::is_nonzero(ExprPtr expr) const {
    if (!expr) return false;
    if (is_positive_scalar(expr) || is_negative_scalar(expr)) return true;
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        if (is_nonzero(*sym)) return true;
        // derive nonzero from relational graph: x>0 or 0>x stored as relations
        return is_greater(expr, ExprPtr()) || is_greater(ExprPtr(), expr);
    }
    if (is_positive(expr) || is_negative(expr)) return true;
    return is_greater(expr, ExprPtr()) || is_greater(ExprPtr(), expr);
}

bool Assumptions::could_be_zero(const Symbol& symbol) const {
    if (is_nonzero(symbol)) return false;
    const auto found = range_symbols_.find(symbol.name);
    if (found == range_symbols_.end()) return true;
    return !exact_range_excludes_zero(found->second.lower, found->second.upper);
}

bool Assumptions::could_be_zero(ExprPtr expr) const {
    if (!expr) return false;
    if (is_zero_expr(expr)) return true;
    if (is_nonzero(expr)) return false;
    if (const auto* sym = expr_cast<Symbol>(expr)) return could_be_zero(*sym);
    return true;
}

bool Assumptions::is_integer(const Symbol& symbol) const {
    return integer_symbols_.contains(symbol.name);
}

bool Assumptions::is_integer(ExprPtr expr) const {
    if (!expr) return false;
    if (expr_is<IntegerLit>(expr)) return true;
    if (const auto* sym = expr_cast<Symbol>(expr)) return is_integer(*sym);
    return false;
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

bool Assumptions::prove_relation(ExprPtr current, ExprPtr target, bool strict_needed, std::unordered_set<const ExprNode*>& visited) const {
    struct State {
        ExprPtr current;
        bool strict_met;
    };

    std::vector<State> stack;
    stack.push_back({current, false});

    auto target_scalar = exact_scalar_from_expr(target);

    while (!stack.empty()) {
        auto [curr, strict_so_far] = stack.back();
        stack.pop_back();

        if (curr == target) {
            if (!strict_needed || strict_so_far) return true;
        }

        // If we reached a scalar and target is a scalar, we can bridge the gap
        if (target_scalar) {
            if (auto curr_scalar = exact_scalar_from_expr(curr)) {
                int cmp = compare_exact_scalars(*curr_scalar, *target_scalar);
                if (cmp < 0) return true; // strict gap found
                if (cmp == 0 && (!strict_needed || strict_so_far)) return true;
            }
        }

        if (visited.contains(curr.get())) continue;
        visited.insert(curr.get());

        // If target is 0, any node proven to be strictly negative satisfies the strict relation < 0.
        if (target == nullptr) {
            if (const auto* sym = expr_cast<Symbol>(curr)) {
                if (negative_symbols_.contains(sym->name)) return true;
            }
        }

        auto it = relations_.find(curr);
        if (it != relations_.end()) {
            for (const auto& rel : it->second) {
                bool new_strict = strict_so_far || (rel.type == RelType::Less);
                stack.push_back({rel.target, new_strict});
            }
        }

        // F2.x-C1: Assumptions transitive closure.
        // Implicitly bridge 0 to all positive symbols that act as sources in the relation graph.
        if (curr == nullptr) {
            for (const auto& [node, rels] : relations_) {
                if (!node) continue;
                if (const auto* sym = expr_cast<Symbol>(node)) {
                    if (positive_symbols_.contains(sym->name)) {
                        stack.push_back({node, true});
                    }
                }
            }
        }
    }

    return false;
}

bool Assumptions::prove_positive_product(const Product& prod) const {
    if (prod.factors.empty()) return false;
    int negative_count = 0;
    for (ExprPtr factor : prod.factors) {
        if (is_positive(factor)) continue;
        if (is_negative(factor)) {
            negative_count++;
            continue;
        }
        return false;
    }
    return (negative_count % 2 == 0);
}

bool Assumptions::prove_positive_linear(ExprPtr expr) const {
    if (const auto* sum = expr_cast<Sum>(expr)) {
        if (sum->terms.size() == 2U) {
            auto check_sub = [&](ExprPtr pos, ExprPtr neg) -> bool {
                const auto* neg_node = expr_cast<Unary>(neg);
                if (neg_node != nullptr && neg_node->op == UnaryOp::Neg)
                    return is_greater(pos, neg_node->operand);
                return false;
            };
            if (check_sub(sum->terms[0], sum->terms[1]) || check_sub(sum->terms[1], sum->terms[0]))
                return true;
        }
        bool has_strict_positive = false;
        for (ExprPtr term : sum->terms) {
            if (is_positive(term)) {
                has_strict_positive = true;
            } else if (!is_nonnegative(term)) {
                return false;
            }
        }
        return has_strict_positive;
    }

    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Sub) return is_greater(bin->left, bin->right);
        if (bin->op == BinaryOp::Add) {
            return (is_positive(bin->left) && is_nonnegative(bin->right)) ||
                   (is_nonnegative(bin->left) && is_positive(bin->right));
        }
        if (bin->op == BinaryOp::Mul) {
            return (is_positive(bin->left) && is_positive(bin->right)) ||
                   (is_negative(bin->left) && is_negative(bin->right));
        }
    }

    return false;
}

} // namespace cas::symbolic
