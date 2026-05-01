#include "cas_api.h"
#include "cas/symbolic.hpp"
#include "cas/parser.hpp"
#include "cas/lexer.hpp"
#include "cas/formatter.hpp"

#include <cstring>
#include <string>
#include <new>
#include <exception>
#include <cstdlib>

// Internal implementation of handles
struct CASContext {
    cas::symbolic::CASContext impl;
};

struct CASExpr {
    cas::ExprPtr ptr;
};

struct CASError {
    std::string message;
};

extern "C" {

CASContextRef cas_context_create() {
    try {
        return new CASContext();
    } catch (...) {
        return nullptr;
    }
}

void cas_context_destroy(CASContextRef ctx) {
    if (ctx) {
        delete ctx;
    }
}

const char* cas_error_message(CASErrorRef err) {
    if (!err) {
        return "No error";
    }
    return err->message.c_str();
}

void cas_error_destroy(CASErrorRef err) {
    if (err) {
        delete err;
    }
}

void cas_expr_destroy(CASExprRef expr) {
    if (expr) {
        delete expr;
    }
}

static char* cas_internal_strdup(const char* s) {
    if (!s) return nullptr;
    size_t len = std::strlen(s) + 1;
    char* d = static_cast<char*>(std::malloc(len));
    if (d) {
        std::memcpy(d, s, len);
    }
    return d;
}

void cas_string_destroy(char* str) {
    if (str) {
        std::free(str);
    }
}

CASErrorRef cas_parse(CASContextRef ctx, const char* input, CASExprRef* out_expr) {
    if (!ctx || !input || !out_expr) {
        return new CASError{"Invalid arguments to cas_parse"};
    }

    try {
        cas::Lexer lexer(input);
        auto tokens_res = lexer.tokenize();
        if (tokens_res.is_error()) {
            return new CASError{tokens_res.error().message};
        }

        cas::Parser parser(tokens_res.value(), ctx->impl.arena());
        auto expr_res = parser.parse();
        if (expr_res.is_error()) {
            return new CASError{expr_res.error().message};
        }

        *out_expr = new CASExpr{expr_res.value()};
        return nullptr;
    } catch (const std::exception& e) {
        return new CASError{std::string("C++ Exception in cas_parse: ") + e.what()};
    } catch (...) {
        return new CASError{"Unknown C++ Exception in cas_parse"};
    }
}

CASErrorRef cas_simplify(CASContextRef ctx, CASExprRef expr, CASExprRef* out_expr) {
    if (!ctx || !expr || !out_expr) {
        return new CASError{"Invalid arguments to cas_simplify"};
    }

    try {
        auto res = ctx->impl.simplify(expr->ptr);
        if (res.is_error()) {
            return new CASError{res.error().message};
        }

        *out_expr = new CASExpr{res.value()};
        return nullptr;
    } catch (const std::exception& e) {
        return new CASError{std::string("C++ Exception in cas_simplify: ") + e.what()};
    } catch (...) {
        return new CASError{"Unknown C++ Exception in cas_simplify"};
    }
}

CASErrorRef cas_format_text(CASContextRef ctx, CASExprRef expr, char** out_text) {
    (void)ctx;
    if (!expr || !out_text) {
        return new CASError{"Invalid arguments to cas_format_text"};
    }

    try {
        cas::formatter::TextFormatter formatter;
        std::string s = formatter.format(expr->ptr);
        *out_text = cas_internal_strdup(s.c_str());
        if (!*out_text) {
            return new CASError{"Memory allocation failed for formatted string"};
        }
        return nullptr;
    } catch (const std::exception& e) {
        return new CASError{std::string("C++ Exception in cas_format_text: ") + e.what()};
    } catch (...) {
        return new CASError{"Unknown C++ Exception in cas_format_text"};
    }
}

CASErrorRef cas_format_latex(CASContextRef ctx, CASExprRef expr, char** out_text) {
    (void)ctx;
    if (!expr || !out_text) {
        return new CASError{"Invalid arguments to cas_format_latex"};
    }

    try {
        cas::formatter::LaTeXFormatter formatter;
        std::string s = formatter.format(expr->ptr);
        *out_text = cas_internal_strdup(s.c_str());
        if (!*out_text) {
            return new CASError{"Memory allocation failed for LaTeX string"};
        }
        return nullptr;
    } catch (const std::exception& e) {
        return new CASError{std::string("C++ Exception in cas_format_latex: ") + e.what()};
    } catch (...) {
        return new CASError{"Unknown C++ Exception in cas_format_latex"};
    }
}

} // extern "C"
