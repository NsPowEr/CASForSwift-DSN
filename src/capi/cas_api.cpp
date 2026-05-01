#include "cas_api.h"
#include "cas/symbolic.hpp"
#include "cas/parser.hpp"
#include "cas/lexer.hpp"
#include "cas/formatter.hpp"
#include "cas/numeric.hpp"

#include <string>

struct CASContext_s {
    cas::symbolic::CASContext impl;
};

struct CASError_s {
    std::string message;
};

extern "C" {

CASContextRef cas_context_create(void) {
    try {
        return new CASContext_s();
    } catch (...) {
        return nullptr;
    }
}

void cas_context_destroy(CASContextRef ctx) {
    delete ctx;
}

void cas_error_destroy(CASErrorRef error) {
    delete error;
}

const char* cas_error_get_message(CASErrorRef error) {
    return error ? error->message.c_str() : nullptr;
}

CASExprRef cas_parse(CASContextRef ctx, const char* input, CASErrorRef* out_error) {
    if (!ctx || !input) return nullptr;
    
    try {
        cas::Lexer lexer(input);
        auto tokens_res = lexer.tokenize();
        if (tokens_res.is_error()) {
            if (out_error) *out_error = new CASError_s{tokens_res.error().message};
            return nullptr;
        }
        
        cas::Parser parser(tokens_res.value(), ctx->impl.arena());
        auto expr_res = parser.parse();
        if (expr_res.is_error()) {
            if (out_error) *out_error = new CASError_s{expr_res.error().message};
            return nullptr;
        }
        
        return reinterpret_cast<CASExprRef>(const_cast<cas::ExprNode*>(expr_res.value().get()));
    } catch (const std::exception& e) {
        if (out_error) *out_error = new CASError_s{e.what()};
        return nullptr;
    } catch (...) {
        if (out_error) *out_error = new CASError_s{"Unknown error during parsing"};
        return nullptr;
    }
}

CASExprRef cas_simplify(CASContextRef ctx, CASExprRef expr, CASErrorRef* out_error) {
    if (!ctx || !expr) return nullptr;
    
    try {
        cas::ExprPtr e(reinterpret_cast<const cas::ExprNode*>(expr));
        auto res = ctx->impl.simplify(e);
        if (res.is_error()) {
            if (out_error) *out_error = new CASError_s{res.error().message};
            return nullptr;
        }
        
        return reinterpret_cast<CASExprRef>(const_cast<cas::ExprNode*>(res.value().get()));
    } catch (const std::exception& e) {
        if (out_error) *out_error = new CASError_s{e.what()};
        return nullptr;
    } catch (...) {
        if (out_error) *out_error = new CASError_s{"Unknown error during simplification"};
        return nullptr;
    }
}

CASExprRef cas_evaluate_numeric(CASContextRef ctx, CASExprRef expr, CASErrorRef* out_error) {
    if (!ctx || !expr) return nullptr;
    
    try {
        cas::ExprPtr e(reinterpret_cast<const cas::ExprNode*>(expr));
        cas::numeric::NumericEvaluator evaluator;
        auto res = evaluator.evaluate(e);
        if (res.is_error()) {
            if (out_error) *out_error = new CASError_s{res.error().message};
            return nullptr;
        }
        
        // Wrap the double in a DecimalLit node in the arena
        auto result_expr = ctx->impl.arena().make<cas::DecimalLit>(res.value());
        return reinterpret_cast<CASExprRef>(const_cast<cas::ExprNode*>(result_expr.get()));
    } catch (const std::exception& e) {
        if (out_error) *out_error = new CASError_s{e.what()};
        return nullptr;
    } catch (...) {
        if (out_error) *out_error = new CASError_s{"Unknown error during numeric evaluation"};
        return nullptr;
    }
}

const char* cas_format_latex(CASExprRef expr) {
    if (!expr) return "";
    try {
        cas::ExprPtr e(reinterpret_cast<const cas::ExprNode*>(expr));
        cas::formatter::LaTeXFormatter formatter;
        static thread_local std::string buf;
        buf = formatter.format(e);
        return buf.c_str();
    } catch (...) {
        return "Error";
    }
}

const char* cas_format_text(CASExprRef expr) {
    if (!expr) return "";
    try {
        cas::ExprPtr e(reinterpret_cast<const cas::ExprNode*>(expr));
        cas::formatter::TextFormatter formatter;
        static thread_local std::string buf;
        buf = formatter.format(e);
        return buf.c_str();
    } catch (...) {
        return "Error";
    }
}

}
