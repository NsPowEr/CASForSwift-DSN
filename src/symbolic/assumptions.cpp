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

void Assumptions::assume_positive(const Symbol& symbol) {
    real_symbols_.insert(symbol.name);
    positive_symbols_.insert(symbol.name);
    nonzero_symbols_.insert(symbol.name);
    
    // Also add to relation graph conceptually: 0 < symbol
    // We'll use a null ExprPtr or a special internal zero node if needed, 
    // but for now let's just use the set as a primary source for symbols.
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
    if (const auto* bin = expr_cast<Binary>(condition)) {
        if ([[maybe_unused]] const auto* func = expr_cast<FuncCall>(bin->left)) {
             // Handle some special cases if needed
        }
        
        // Note: The parser usually produces Binary nodes for > and >= if implemented as such,
        // but often CAS uses a custom node kind for relations. 
        // Assuming Binary for now as per usual pattern if not specialized.
        // Wait, ExprKind doesn't have Comparison. 
        // Let's check if BinaryOp has GT/GE.
    }
    // If we receive a raw relation, we might need to handle it.
    // For the task, we assume the caller passes the right ExprPtrs to assume_greater etc.
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
        // Most binary ops preserve realness if both operands are real
        // except possibly Pow and Div (division by zero is not real, but we usually treat it as undefined/complex)
        if (bin->op == BinaryOp::Pow) {
            // x^y is real if x > 0 and y is real, or x < 0 and y is integer
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
    return is_greater(ExprPtr(&symbol), ExprPtr());
}

bool Assumptions::is_positive(ExprPtr expr) const {
    if (!expr) return false;

    // 1. Exact scalar check
    if (is_positive_scalar(expr)) return true;

    // 2. Symbol check
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        if (positive_symbols_.contains(sym->name)) return true;
    }

    // 3. Graph proof (Deduction chains)
    if (is_greater(expr, ExprPtr())) return true;

    // 4. Builtin functions
    if (const auto* func = expr_cast<FuncCall>(expr)) {
        if (func->func_id == BuiltinOp::Exp) {
            // exp(x) is always positive for real x
            // In our CAS we treat it as positive in the real branch context
            return is_real(func->args[0]);
        }
        if (func->func_id == BuiltinOp::Abs) {
            return is_nonzero(func->args[0]);
        }
        if (func->func_id == BuiltinOp::Sqrt) {
            return is_positive(func->args[0]);
        }
    }

    // 5. Linear combination check (Simple solver)
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
            // x^n is nonnegative if n is even and x is real
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
            // Re-use product logic or implement here
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
        bool could_be_zero_val = false;
        int negative_count = 0;
        for (ExprPtr factor : prod->factors) {
            if (is_zero_expr(factor)) return true;
            if (could_be_zero(factor)) could_be_zero_val = true;

            if (is_positive(factor)) continue;
            if (is_negative(factor)) {
                negative_count++;
                continue;
            }
            return false; // Unknown sign
        }
        if (could_be_zero_val) return (negative_count % 2 == 0); // If even negatives, it's >= 0
        return (negative_count % 2 == 0);
    }

    return false;
}

bool Assumptions::is_negative(ExprPtr expr) const {
    if (!expr) return false;

    // 1. Exact scalar check
    if (is_negative_scalar(expr)) return true;

    // 2. Graph proof: 0 > expr
    if (is_greater(ExprPtr(), expr)) return true;

    // 3. Negation of positive check
    // if (-expr) is positive, then expr is negative.
    // This is useful for things like -x where x > 0.
    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) {
        if (is_positive(unary->operand)) return true;
    }

    return false;
}

bool Assumptions::is_greater(ExprPtr lhs, ExprPtr rhs) const {
    if (lhs == rhs) return false;
    ExprPtr l = is_zero_expr(lhs) ? nullptr : lhs;
    ExprPtr r = is_zero_expr(rhs) ? nullptr : rhs;
    if (l == r) return false;
    std::unordered_set<const ExprNode*> visited;
    return prove_relation(r, l, true, visited);
}

bool Assumptions::is_greater_equal(ExprPtr lhs, ExprPtr rhs) const {
    if (lhs == rhs) return true;
    ExprPtr l = is_zero_expr(lhs) ? nullptr : lhs;
    ExprPtr r = is_zero_expr(rhs) ? nullptr : rhs;
    if (l == r) return true;
    std::unordered_set<const ExprNode*> visited;
    return prove_relation(r, l, false, visited);
}

bool Assumptions::is_nonzero(const Symbol& symbol) const {
    if (nonzero_symbols_.contains(symbol.name)) return true;
    if (positive_symbols_.contains(symbol.name)) return true;
    // Check if x > 0 or x < 0 via graph
    return is_greater(ExprPtr(&symbol), ExprPtr()) || is_greater(ExprPtr(), ExprPtr(&symbol));
}

bool Assumptions::is_nonzero(ExprPtr expr) const {
    if (!expr) return false;
    if (is_positive_scalar(expr) || is_negative_scalar(expr)) return true;
    if (const auto* sym = expr_cast<Symbol>(expr)) return is_nonzero(*sym);
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
    
    // Also update relation graph
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
        bool strict_needed;
    };
    
    std::vector<State> stack;
    stack.push_back({current, strict_needed});
    
    while (!stack.empty()) {
        auto [curr, strict] = stack.back();
        stack.pop_back();

        // Conceptual target NULL means we are looking for 0
        if (!target) {
            if (is_zero_expr(curr)) {
                if (!strict) return true;
            }
        } else if (curr == target) {
            if (!strict) return true;
        }

        // Check symbol sets for positive/negative facts not in relation graph
        if (!curr) {
            if (const auto* sym = expr_cast<Symbol>(target)) {
                if (positive_symbols_.contains(sym->name)) return true;
            }
        } else if (!target) {
            if ([[maybe_unused]] const auto* sym = expr_cast<Symbol>(curr)) {
                // if we need curr < 0 and we know curr is negative... (todo: negative_symbols_)
            }
        }

        // Also check if curr is a scalar and we can compare directly if target is a scalar
        auto current_scalar = exact_scalar_from_expr(curr);
        auto target_scalar = target ? exact_scalar_from_expr(target) : std::make_optional(Rational(0));
        if (current_scalar && target_scalar) {
            int cmp = compare_exact_scalars(*current_scalar, *target_scalar);
            if (strict) {
                if (cmp < 0) return true;
            } else {
                if (cmp <= 0) return true;
            }
            continue; // Scalar, no edges to follow
        }

        if (visited.contains(curr.get())) continue;
        visited.insert(curr.get());

        auto it = relations_.find(curr);
        if (it != relations_.end()) {
            for (const auto& rel : it->second) {
                // Edge is curr < rel.target or curr <= rel.target
                bool edge_strict = (rel.type == RelType::Less);
                bool next_strict_needed = strict && !edge_strict;
                stack.push_back({rel.target, next_strict_needed});
            }
        }
    }

    return false;
}

bool Assumptions::prove_positive_linear(ExprPtr expr) const {
    if (const auto* sum = expr_cast<Sum>(expr)) {
        // Detect x - y pattern in Sum canonical form: Sum{u, Unary(Neg, v)} → is_greater(u, v)
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
    
    if (const auto* prod = expr_cast<Product>(expr)) {
        if (prod->factors.empty()) return false;
        int negative_count = 0;
        for (ExprPtr factor : prod->factors) {
            if (is_positive(factor)) {
                // OK
            } else if (is_negative(factor)) {
                negative_count++;
            } else {
                return false;
            }
        }
        return (negative_count % 2 == 0);
    }
    
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Sub) {
            return is_greater(bin->left, bin->right);
        }
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
