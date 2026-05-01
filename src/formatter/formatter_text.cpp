#include "cas/formatter.hpp"
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <vector>

namespace cas::formatter {

namespace {

// Pass 1: Conta le occorrenze di ogni nodo nel DAG
void count_occurrences(ExprPtr expr, std::unordered_map<ExprPtr, int>& counts) {
    if (!expr) return;
    auto& c = counts[expr];
    c++;
    if (c > 1) return; // Già visitato questo ramo condiviso

    visit_expr(expr, [&](const auto& node) {
        using NodeT = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<NodeT, Unary>) {
            count_occurrences(node.operand, counts);
        } else if constexpr (std::is_same_v<NodeT, Binary>) {
            count_occurrences(node.left, counts);
            count_occurrences(node.right, counts);
        } else if constexpr (std::is_same_v<NodeT, FuncCall>) {
            for (auto arg : node.args) count_occurrences(arg, counts);
        } else if constexpr (std::is_same_v<NodeT, Sum>) {
            for (auto term : node.terms) count_occurrences(term, counts);
        } else if constexpr (std::is_same_v<NodeT, Product>) {
            for (auto factor : node.factors) count_occurrences(factor, counts);
        } else if constexpr (std::is_same_v<NodeT, RootOf>) {
            count_occurrences(node.polynomial, counts);
        } else if constexpr (std::is_same_v<NodeT, Matrix>) {
            for (auto elem : node.elements) count_occurrences(elem, counts);
        }
    });
}

// Verifica se un nodo è "complesso" abbastanza da giustificare il CSE
bool is_complex(ExprPtr expr) {
    if (!expr) return false;
    return visit_expr(expr, [](const auto& node) -> bool {
        using NodeT = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<NodeT, IntegerLit> || std::is_same_v<NodeT, Symbol> || std::is_same_v<NodeT, Constant>) {
            return false;
        }
        return true;
    });
}

} // namespace

int TextFormatter::precedence(ExprPtr expr) {
    if (!expr) return 0;
    return visit_expr(expr, [](const auto& node) -> int {
        using NodeT = std::decay_t<decltype(node)>;
        if constexpr (std::is_same_v<NodeT, IntegerLit> || std::is_same_v<NodeT, DecimalLit> || 
                      std::is_same_v<NodeT, Symbol> || std::is_same_v<NodeT, Constant> || 
                      std::is_same_v<NodeT, FuncCall> || std::is_same_v<NodeT, Matrix>) {
            return 100; // Atomici
        } else if constexpr (std::is_same_v<NodeT, Unary>) {
            return 80;
        } else if constexpr (std::is_same_v<NodeT, Binary>) {
            switch (node.op) {
                case BinaryOp::Pow: return 90;
                case BinaryOp::Mul:
                case BinaryOp::Div: return 70;
                case BinaryOp::Add:
                case BinaryOp::Sub: return 60;
                default: return 0;
            }
        } else if constexpr (std::is_same_v<NodeT, Sum>) {
            return 60;
        } else if constexpr (std::is_same_v<NodeT, Product>) {
            return 70;
        }
        return 0;
    });
}

std::string TextFormatter::format_child(ExprPtr parent, ExprPtr child, bool is_right) {
    int p_prec = precedence(parent);
    int c_prec = precedence(child);
    
    std::string s = format(child);
    
    if (c_prec < p_prec) return "(" + s + ")";
    
    if (c_prec == p_prec && is_right) {
        if (const auto* b = expr_cast<Binary>(parent)) {
            if (b->op == BinaryOp::Sub || b->op == BinaryOp::Div) return "(" + s + ")";
        }
    }
    
    return s;
}

std::string TextFormatter::format_internal(ExprPtr expr, std::unordered_map<ExprPtr, std::string>& cse_map) {
    if (!expr) return "";

    // Se il nodo è nel cse_map, usiamo la variabile
    if (auto it = cse_map.find(expr); it != cse_map.end()) {
        return it->second;
    }

    return visit_expr(expr, [this, expr, &cse_map](const auto& node) -> std::string {
        using NodeT = std::decay_t<decltype(node)>;

        if constexpr (std::is_same_v<NodeT, IntegerLit>) {
            return node.value.decimal();
        } else if constexpr (std::is_same_v<NodeT, RationalLit>) {
            return node.numerator.decimal() + "/" + node.denominator.decimal();
        } else if constexpr (std::is_same_v<NodeT, DecimalLit>) {
            return node.text;
        } else if constexpr (std::is_same_v<NodeT, Symbol>) {
            return node.name;
        } else if constexpr (std::is_same_v<NodeT, Constant>) {
            switch (node.value) {
                case MathConstant::Pi: return "pi";
                case MathConstant::E: return "e";
                case MathConstant::I: return "I";
                case MathConstant::Infinity: return "inf";
                default: return "const";
            }
        } else if constexpr (std::is_same_v<NodeT, Unary>) {
            return "-" + format_child(expr, node.operand);
        } else if constexpr (std::is_same_v<NodeT, Binary>) {
            std::string op_str;
            switch (node.op) {
                case BinaryOp::Add: op_str = " + "; break;
                case BinaryOp::Sub: op_str = " - "; break;
                case BinaryOp::Mul: op_str = " * "; break;
                case BinaryOp::Div: op_str = " / "; break;
                case BinaryOp::Mod: op_str = " % "; break;
                case BinaryOp::Pow: return format_child(expr, node.left) + "^" + format_child(expr, node.right, true);
                case BinaryOp::Equal: op_str = " = "; break;
            }
            return format_child(expr, node.left) + op_str + format_child(expr, node.right, true);
        } else if constexpr (std::is_same_v<NodeT, Sum>) {
            std::string s;
            for (size_t i = 0; i < node.terms.size(); ++i) {
                ExprPtr term = node.terms[i];
                const auto* unary = expr_cast<Unary>(term);
                bool is_neg = (unary != nullptr && unary->op == UnaryOp::Neg);
                
                if (i == 0) {
                    s += format_child(expr, term);
                } else {
                    if (is_neg) {
                        s += " - " + format_child(expr, unary->operand);
                    } else {
                        s += " + " + format_child(expr, term);
                    }
                }
            }
            return s;
        } else if constexpr (std::is_same_v<NodeT, Product>) {
            std::string s;
            for (size_t i = 0; i < node.factors.size(); ++i) {
                if (i > 0) s += " * ";
                s += format_child(expr, node.factors[i]);
            }
            return s;
        } else if constexpr (std::is_same_v<NodeT, RootOf>) {
            std::string s = "RootOf(" + format_internal(node.polynomial, cse_map) + ", " + node.variable.name;
            if (node.root_index.has_value()) {
                s += ", " + std::to_string(*node.root_index);
            }
            return s + ")";
        } else if constexpr (std::is_same_v<NodeT, FuncCall>) {
            std::string s = node.name + "(";
            for (size_t i = 0; i < node.args.size(); ++i) {
                if (i > 0) s += ", ";
                s += format_internal(node.args[i], cse_map);
            }
            return s + ")";
        } else if constexpr (std::is_same_v<NodeT, Matrix>) {
            std::string s = "[";
            for (size_t r = 0; r < node.rows; ++r) {
                s += "[";
                for (size_t c = 0; c < node.cols; ++c) {
                    if (c > 0) s += ", ";
                    s += format_internal(node.elements[r * node.cols + c], cse_map);
                }
                s += "]";
                if (r < node.rows - 1) s += ", ";
            }
            return s + "]";
        }
        return "?";
    });
}

std::string TextFormatter::format(ExprPtr expr) {
    if (!expr) return "";

    std::unordered_map<ExprPtr, int> counts;
    count_occurrences(expr, counts);

    std::unordered_map<ExprPtr, std::string> cse_map;
    std::vector<std::pair<std::string, std::string>> definitions;
    int var_idx = 1;

    // Identifichiamo i nodi da estrarre (complessi e ripetuti)
    // Usiamo una visita BFS o DFS ordinata per definire le variabili in ordine corretto
    // Per semplicità qui usiamo un approccio greedy ma corretto per la visualizzazione
    std::vector<ExprPtr> to_extract;
    for (auto const& [node, count] : counts) {
        if (count > 1 && is_complex(node) && node != expr) {
            to_extract.push_back(node);
        }
    }

    // Sort per profondità (opzionale, per ordine naturale v1, v2...)
    // Qui definiamo le variabili
    for (auto node : to_extract) {
        std::string var_name = "v" + std::to_string(var_idx++);
        cse_map[node] = var_name;
    }

    // Formattiamo le definizioni
    for (auto node : to_extract) {
        std::string var_name = cse_map[node];
        // Rimuoviamo temporaneamente se stesso dal map per formattare la definizione
        std::string val = var_name;
        auto it = cse_map.find(node);
        cse_map.erase(it);
        val = format_internal(node, cse_map);
        cse_map[node] = var_name;
        definitions.push_back({var_name, val});
    }

    std::string main_expr = format_internal(expr, cse_map);

    if (definitions.empty()) return main_expr;

    std::string result = main_expr + " where { ";
    for (size_t i = 0; i < definitions.size(); ++i) {
        if (i > 0) result += ", ";
        result += definitions[i].first + " = " + definitions[i].second;
    }
    result += " }";
    return result;
}

} // namespace cas::formatter
