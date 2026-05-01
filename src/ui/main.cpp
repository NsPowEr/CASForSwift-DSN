#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include "cas/formatter.hpp"
#include <iostream>
#include <string>

using namespace cas;

void print_help() {
    std::cout << "REAL CAS ENGINE REPL\n";
    std::cout << "Comandi speciali:\n";
    std::cout << "  :q, :quit   Esci\n";
    std::cout << "  :h, :help   Mostra questo aiuto\n";
}

int main() {
    symbolic::CASContext ctx;
    formatter::TextFormatter text_fmt;
    formatter::LaTeXFormatter latex_fmt;
    formatter::Ascii2DFormatter ascii_fmt;

    std::string line;
    print_help();

    while (true) {
        std::cout << "\ncas> ";
        if (!std::getline(std::cin, line)) break;
        if (line == ":q" || line == ":quit") break;
        if (line == ":h" || line == ":help") { print_help(); continue; }
        if (line.empty()) continue;

        auto tokens = Lexer(line).tokenize();
        if (tokens.is_error()) {
            std::cerr << "Errore Lexer: " << tokens.error().message << "\n";
            continue;
        }

        Parser parser(tokens.value(), ctx.arena());
        auto expr = parser.parse();
        if (expr.is_error()) {
            std::cerr << "Errore Parser: " << expr.error().message << "\n";
            continue;
        }

        auto simplified = ctx.simplify(expr.value());
        if (simplified.is_error()) {
            std::cerr << "Errore Semplificatore: " << simplified.error().message << "\n";
            continue;
        }

        ExprPtr res = simplified.value();

        std::cout << "\n[RESULT]\n";
        std::cout << "Text:  " << text_fmt.format(res) << "\n";
        std::cout << "LaTeX: " << latex_fmt.format(res) << "\n";
        std::cout << "\n[ASCII 2D]\n";
        std::cout << ascii_fmt.format(res) << "\n";
    }

    return 0;
}
