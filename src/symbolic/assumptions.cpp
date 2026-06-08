#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "symbolic_internal.hpp"
#include <unordered_set>
#include <algorithm>

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

void Assumptions::assume_greater(ExprPtr lhs, ExprPtr rhs) {
    // lhs > rhs  =>  rhs < lhs
    ExprPtr r = is_zero_expr(rhs) ? nullptr : rhs;
    ExprPtr l = is_zero_expr(lhs) ? nullptr : lhs;
    relations_[r].push_back({l, RelType::Less});
}

void Assumptions::assume_greater_equal(ExprPtr lhs, ExprPtr rhs) {
    // lhs >= rhs  =>  rhs <= lhs
    ExprPtr r = is_zero_expr(rhs) ? nullptr : rhs;
    ExprPtr l = is_zero_expr(lhs) ? nullptr : lhs;
    relations_[r].push_back({l, RelType::LessEqual});
}

void Assumptions::assume(ExprPtr condition) {
    if (!condition) return;

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

bool Assumptions::is_real(const Symbol& symbol) const {
    return real_symbols_.contains(symbol.name);
}

bool Assumptions::is_real(ExprPtr expr) const {
    if (!expr) return false;
    if (expr_is<IntegerLit>(expr) || expr_is<RationalLit>(expr) || expr_is<DecimalLit>(expr)) return true;
    if (const auto* sym = expr_cast<Symbol>(expr)) return is_real(*sym);
    if (const auto* con = expr_cast<Constant>(expr)) {
        return con->value == MathConstant::Pi || con->value == MathConstant::E;
    }
    
    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr term : sum->terms) if (!is_real(term)) return false;
        return true;
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr factor : prod->factors) if (!is_real(factor)) return false;
        return true;
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) {
            if (is_positive(bin->left) && is_real(bin->right)) return true;
            if (is_negative(bin->left) && is_integer(bin->right)) return true;
            return false;
        }
        return is_real(bin->left) && is_real(bin->right);
    }
    if (const auto* unary = expr_cast<Unary>(expr)) {
        return is_real(unary->operand);
    }
    if (const auto* func = expr_cast<FuncCall>(expr)) {
        if (func->func_id == BuiltinOp::Exp || func->func_id == BuiltinOp::Sin || 
            func->func_id == BuiltinOp::Cos || func->func_id == BuiltinOp::Tan ||
            func->func_id == BuiltinOp::Abs) {
            return is_real(func->args[0]);
        }
        if (func->func_id == BuiltinOp::Sqrt) {
            return is_nonnegative(func->args[0]);
        }
        if (func->func_id == BuiltinOp::Ln) {
            return is_positive(func->args[0]);
        }
    }
    return false;
}

bool Assumptions::is_positive(const Symbol& symbol) const {
    if (positive_symbols_.contains(symbol.name)) return true;
    return false;
}

bool Assumptions::is_positive(ExprPtr expr) const {
    if (!expr) return false;
    if (is_positive_scalar(expr)) return true;
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        if (positive_symbols_.contains(sym->name)) return true;
    }
    if (is_greater(expr, ExprPtr())) return true;

    if (const auto* func = expr_cast<FuncCall>(expr)) {
        if (func->func_id == BuiltinOp::Exp) return is_real(func->args[0]);
        if (func->func_id == BuiltinOp::Abs) return is_nonzero(func->args[0]);
        if (func->func_id == BuiltinOp::Sqrt) return is_positive(func->args[0]);
    }

    if (const auto* prod = expr_cast<Product>(expr)) return prove_positive_product(*prod);
    if (prove_positive_linear(expr)) return true;

    return false;
}

bool Assumptions::is_nonnegative(ExprPtr expr) const {
    if (!expr) return false;
    if (is_zero_expr(expr)) return true;
    if (is_positive(expr)) return true;
    if (is_greater_equal(expr, ExprPtr())) return true;
    
    if (const auto* func = expr_cast<FuncCall>(expr)) {
        if (func->func_id == BuiltinOp::Abs) return true;
        if (func->func_id == BuiltinOp::Sqrt) return true;
    }

    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) {
            if (auto scalar = exact_scalar_from_expr(bin->right)) {
                if (scalar->is_integer() && scalar->numerator() > BigInt(0)) {
                    if (scalar->numerator() % BigInt(2) == BigInt(0)) {
                        return is_real(bin->left);
                    }
                }
            }
        }
        if (bin->op == BinaryOp::Add) {
            return is_nonnegative(bin->left) && is_nonnegative(bin->right);
        }
        if (bin->op == BinaryOp::Mul) {
            if (is_positive(bin->left) && is_positive(bin->right)) return true;
            if (is_negative(bin->left) && is_negative(bin->right)) return true;
            if (is_zero_expr(bin->left) || is_zero_expr(bin->right)) return true;
            return false;
        }
    }

    if (const auto* sum = expr_cast<Sum>(expr)) {
        for (ExprPtr term : sum->terms) {
            if (!is_nonnegative(term)) return false;
        }
        return true;
    }
    
    if (const auto* prod = expr_cast<Product>(expr)) {
        int negative_count = 0;
        for (ExprPtr factor : prod->factors) {
            if (is_zero_expr(factor)) return true;

            if (is_positive(factor)) continue;
            if (is_negative(factor)) {
                negative_count++;
                continue;
            }
            return false; 
        }
        return (negative_count % 2 == 0);
    }

    return false;
}

bool Assumptions::is_negative(ExprPtr expr) const {
    if (!expr) return false;
    if (is_negative_scalar(expr)) return true;
    if (is_greater(ExprPtr(), expr)) return true;
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        if (negative_symbols_.contains(sym->name)) return true;
    }
    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) {
        if (is_positive(unary->operand)) return true;
    }

    return false;
}

bool Assumptions::is_greater(ExprPtr lhs, ExprPtr rhs) const {
    if (lhs == rhs) return false;

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
    if (l_scalar && r_scalar) return compare_exact_scalars(*l_scalar, *r_scalar) > 0;

    // 2. Direct graph proof
    ExprPtr l = is_zero_expr(lhs) ? nullptr : lhs;
    ExprPtr r = is_zero_expr(rhs) ? nullptr : rhs;
    if (l != r) {
        std::unordered_set<const ExprNode*> visited;
        if (prove_relation(r, l, true, visited)) return true;
    }

    // Helper to extract terms from Sum or Binary Add/Sub
    auto get_terms = [](ExprPtr e) -> std::vector<ExprPtr> {
        if (const auto* s = expr_cast<Sum>(e)) return s->terms;
        if (const auto* b = expr_cast<Binary>(e)) {
            if (b->op == BinaryOp::Add) return {b->left, b->right};
            if (b->op == BinaryOp::Sub) return {b->left, ExprPtr()}; // simplified
        }
        return {e};
    };

    // 3. Decompose sums: (a + b) > (c + d) if a > c and b >= d
    auto l_terms = get_terms(lhs);
    auto r_terms = get_terms(rhs);
    
    if (l_terms.size() > 1 || r_terms.size() > 1) {
        if (l_terms.size() == r_terms.size()) {
            bool all_ge = true;
            bool one_gt = false;
            for (size_t i = 0; i < l_terms.size(); ++i) {
                if (is_greater(l_terms[i], r_terms[i])) {
                    one_gt = true;
                } else if (!is_greater_equal(l_terms[i], r_terms[i])) {
                    all_ge = false;
                    break;
                }
            }
            if (all_ge && one_gt) return true;
        }
        
        // (a + b) > c if a > c and b >= 0
        if (l_terms.size() > 1 && r_terms.size() == 1) {
            for (size_t i = 0; i < l_terms.size(); ++i) {
                bool others_nonnegative = true;
                for (size_t j = 0; j < l_terms.size(); ++j) {
                    if (i == j) continue;
                    if (!is_nonnegative(l_terms[j])) { others_nonnegative = false; break; }
                }
                if (others_nonnegative && is_greater(l_terms[i], rhs)) return true;
            }
        }
    }

    // 4. Linear combination: x - y > 0
    if ((expr_is<Sum>(lhs) || expr_is<Binary>(lhs)) && is_zero_expr(rhs)) {
        if (prove_positive_linear(lhs)) return true;
    }

    return false;
}

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
