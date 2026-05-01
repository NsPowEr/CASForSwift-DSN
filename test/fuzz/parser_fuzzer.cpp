#include "cas_api.h"
#include <stddef.h>
#include <stdint.h>
#include <string>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    if (Size == 0) return 0;

    // Convert input data to a null-terminated string
    std::string input(reinterpret_cast<const char*>(Data), Size);
    
    // Create a context
    CASContextRef ctx = cas_context_create();
    if (!ctx) return 0;

    // Try to parse the input
    CASExprRef expr = nullptr;
    CASErrorRef err = cas_parse(ctx, input.c_str(), &expr);

    // If parsing was successful, we might want to do something with the expression
    // to fuzz further, but for now, just parsing is enough as requested.
    if (err) {
        cas_error_destroy(err);
    }

    if (expr) {
        // Optionally evaluate or format to fuzz other parts
        char* text = nullptr;
        CASErrorRef format_err = cas_format_text(ctx, expr, &text);
        if (format_err) {
            cas_error_destroy(format_err);
        }
        if (text) {
            cas_string_destroy(text);
        }
        
        cas_expr_destroy(expr);
    }

    // Cleanup
    cas_context_destroy(ctx);
    
    return 0;
}
