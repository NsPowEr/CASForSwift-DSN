#ifndef CAS_C_API_H
#define CAS_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdbool.h>

typedef struct cas_context_t cas_context_t;

// Context lifecycle
cas_context_t* cas_create_context(void);
void cas_destroy_context(cas_context_t* ctx);

// Computation
typedef struct {
    bool ok;
    const char* error;
    const char* result_latex;
    const char* result_ascii;
    double numeric_value;
    bool has_numeric;
} cas_compute_result_t;

cas_compute_result_t cas_simplify(cas_context_t* ctx, const char* input);

// Memory management for strings (owned by API)
void cas_free_result(cas_compute_result_t res);

#ifdef __cplusplus
}
#endif

#endif // CAS_C_API_H
