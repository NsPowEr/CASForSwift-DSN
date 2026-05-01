#ifndef CAS_API_H
#define CAS_API_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

// Opaque handles
typedef struct CASContext_s* CASContextRef;
typedef struct CASExpr_s* CASExprRef;
typedef struct CASError_s* CASErrorRef;

// Lifecycle
CASContextRef cas_context_create(void);
void cas_context_destroy(CASContextRef ctx);
void cas_error_destroy(CASErrorRef error);

// Error inspection
const char* cas_error_get_message(CASErrorRef error);

// Operations
CASExprRef cas_parse(CASContextRef ctx, const char* input, CASErrorRef* out_error);
CASExprRef cas_simplify(CASContextRef ctx, CASExprRef expr, CASErrorRef* out_error);
CASExprRef cas_evaluate_numeric(CASContextRef ctx, CASExprRef expr, CASErrorRef* out_error);

// Formatting
const char* cas_format_latex(CASExprRef expr);
const char* cas_format_text(CASExprRef expr);

#ifdef __cplusplus
}
#endif

#endif // CAS_API_H
