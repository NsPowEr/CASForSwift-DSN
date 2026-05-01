#ifndef CAS_API_H
#define CAS_API_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Opaque handle to a CAS context.
 */
typedef struct CASContext* CASContextRef;

/**
 * @brief Opaque handle to a CAS expression.
 */
typedef struct CASExpr* CASExprRef;

/**
 * @brief Opaque handle to a CAS error.
 */
typedef struct CASError* CASErrorRef;

// Context
CASContextRef cas_context_create(void);
void cas_context_destroy(CASContextRef ctx);

// Error
const char* cas_error_get_message(CASErrorRef err);
void cas_error_destroy(CASErrorRef err);

// Operations
CASErrorRef cas_parse(CASContextRef ctx, const char* input, CASExprRef* out_expr);
CASErrorRef cas_simplify(CASContextRef ctx, CASExprRef expr, CASExprRef* out_expr);
CASErrorRef cas_evaluate_numeric(CASContextRef ctx, CASExprRef expr, CASExprRef* out_expr);

// Formatting
CASErrorRef cas_format_text(CASContextRef ctx, CASExprRef expr, char** out_text);
CASErrorRef cas_format_latex(CASContextRef ctx, CASExprRef expr, char** out_text);

// Utility
void cas_string_destroy(char* str);
void cas_expr_destroy(CASExprRef expr);

#ifdef __cplusplus
}
#endif

#endif // CAS_API_H
