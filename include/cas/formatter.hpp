#pragma once

#include "cas/ast.hpp"
#include <string>
#include <vector>
#include <unordered_map>

namespace cas::formatter {

/**
 * Formatter 1D: Genera testo parser-safe.
 */
class TextFormatter {
public:
    [[nodiscard]] std::string format(ExprPtr expr);

private:
    [[nodiscard]] int precedence(ExprPtr expr);
    [[nodiscard]] std::string format_child(ExprPtr parent, ExprPtr child, bool is_right = false);
    [[nodiscard]] std::string format_internal(ExprPtr expr, std::unordered_map<ExprPtr, std::string>& cse_map);
};

/**
 * Formatter LaTeX.
 */
class LaTeXFormatter {
public:
    [[nodiscard]] std::string format(ExprPtr expr);
};

/**
 * Box Model per il rendering 2D (ASCII/Grafico).
 */
struct LayoutBox {
    std::vector<std::string> lines;
    int width{0};
    int height{0};
    int baseline{0}; // Riga dell'asse principale
};

class Ascii2DFormatter {
public:
    [[nodiscard]] std::string format(ExprPtr expr);
    [[nodiscard]] LayoutBox layout(ExprPtr expr);
};

} // namespace cas::formatter
