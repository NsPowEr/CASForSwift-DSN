#include "cas_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    printf("Testing CAS C API...\n");

    CASContextRef ctx = cas_context_create();
    if (!ctx) {
        fprintf(stderr, "Failed to create context\n");
        return 1;
    }

    const char* input = "1 + 2 * 3";
    CASExprRef expr = NULL;
    CASErrorRef err = cas_parse(ctx, input, &expr);
    if (err) {
        fprintf(stderr, "Parse error: %s\n", cas_error_get_message(err));
        cas_error_destroy(err);
        cas_context_destroy(ctx);
        return 1;
    }

    char* text = NULL;
    err = cas_format_text(ctx, expr, &text);
    if (err) {
        fprintf(stderr, "Format error: %s\n", cas_error_get_message(err));
        cas_error_destroy(err);
        cas_expr_destroy(expr);
        cas_context_destroy(ctx);
        return 1;
    }

    printf("Parsed expression: %s\n", text);

    CASExprRef simplified = NULL;
    err = cas_simplify(ctx, expr, &simplified);
    if (err) {
        fprintf(stderr, "Simplify error: %s\n", cas_error_get_message(err));
        cas_error_destroy(err);
    } else {
        char* simplified_text = NULL;
        cas_format_text(ctx, simplified, &simplified_text);
        printf("Simplified expression: %s\n", simplified_text);
        cas_string_destroy(simplified_text);
        cas_expr_destroy(simplified);
    }

    // Test numeric evaluation
    CASExprRef numeric = NULL;
    err = cas_evaluate_numeric(ctx, expr, &numeric);
    if (err) {
        fprintf(stderr, "Numeric eval error: %s\n", cas_error_get_message(err));
        cas_error_destroy(err);
    } else {
        char* numeric_text = NULL;
        cas_format_text(ctx, numeric, &numeric_text);
        printf("Numeric evaluation: %s\n", numeric_text);
        cas_string_destroy(numeric_text);
        cas_expr_destroy(numeric);
    }

    cas_string_destroy(text);
    cas_expr_destroy(expr);
    cas_context_destroy(ctx);

    printf("C API test passed!\n");
    return 0;
}
