#pragma once
// runner_format.hpp — output size cap for FAIL diagnostics.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Runner_Robustness.md (F7.5.A3).
//
// integrate(sin(x)^N, x) for N ≥ 5 can produce expressions whose
// textual rendering is hundreds of KB. Without a cap the shell stdout
// buffer truncates the runner's summary table. The cap applies only to
// diagnostic dumps — mathematically_equal still consumes the full AST.

#include "cas/ast.hpp"
#include "cas/formatter.hpp"

#include <cstddef>
#include <string>

namespace cas::golden {

inline constexpr std::size_t kDefaultFormatCapBytes = 4096;
inline constexpr std::size_t kDefaultMaxAstNodes = 2000;

// Recursively count AST nodes up to `budget`; returns budget+1 if exceeded.
// Bounded recursion: every visit decrements budget; early exit on overflow.
inline std::size_t count_ast_nodes_capped(cas::ExprPtr expr, std::size_t budget) {
    if (!expr) return 0;
    if (budget == 0) return 1;
    std::size_t total = 1;
    auto add = [&](std::size_t k) {
        if (total > budget - k) total = budget + 1;
        else total += k;
    };
    if (total > budget) return total;
    switch (expr->kind) {
        case cas::ExprKind::Sum: {
            const auto* s = cas::expr_cast<cas::Sum>(expr);
            if (s) for (auto& t : s->terms) {
                if (total > budget) return total;
                add(count_ast_nodes_capped(t, budget - total));
            }
            break;
        }
        case cas::ExprKind::Product: {
            const auto* p = cas::expr_cast<cas::Product>(expr);
            if (p) for (auto& f : p->factors) {
                if (total > budget) return total;
                add(count_ast_nodes_capped(f, budget - total));
            }
            break;
        }
        case cas::ExprKind::Binary: {
            const auto* b = cas::expr_cast<cas::Binary>(expr);
            if (b) {
                add(count_ast_nodes_capped(b->left, budget - total));
                if (total > budget) return total;
                add(count_ast_nodes_capped(b->right, budget - total));
            }
            break;
        }
        case cas::ExprKind::Unary: {
            const auto* u = cas::expr_cast<cas::Unary>(expr);
            if (u) add(count_ast_nodes_capped(u->operand, budget - total));
            break;
        }
        case cas::ExprKind::FuncCall: {
            const auto* fc = cas::expr_cast<cas::FuncCall>(expr);
            if (fc) for (auto& a : fc->args) {
                if (total > budget) return total;
                add(count_ast_nodes_capped(a, budget - total));
            }
            break;
        }
        default:
            break;  // atomic
    }
    return total;
}

// Format expression as text, but bypass formatter entirely if the AST
// is huge (formatter recursion alone may exceed sane time). Then cap
// resulting string to `cap_bytes`.
inline std::string format_expr_capped(cas::ExprPtr expr,
                                       std::size_t cap_bytes = kDefaultFormatCapBytes,
                                       std::size_t max_nodes = kDefaultMaxAstNodes) {
    if (!expr) return "<null>";
    std::size_t n = count_ast_nodes_capped(expr, max_nodes);
    if (n > max_nodes) {
        return "<expr too large: > " + std::to_string(max_nodes) + " AST nodes>";
    }
    cas::formatter::TextFormatter fmt;
    std::string s = fmt.format(expr);
    if (s.size() <= cap_bytes) return s;
    std::size_t trimmed = s.size() - cap_bytes;
    s.resize(cap_bytes > 64 ? cap_bytes - 48 : cap_bytes);
    s += " ... <truncated " + std::to_string(trimmed) + " bytes>";
    return s;
}

}  // namespace cas::golden
