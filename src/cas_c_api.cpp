#include "cas_c_api.h"
#include "CasGuiSession.hpp"
#include <cstring>

struct cas_context_t {
    cas::gui::CasGuiSession session;
};

static char* cas_strdup(const std::string& s) {
    if (s.empty()) return nullptr;
    char* res = (char*)malloc(s.size() + 1);
    std::strcpy(res, s.c_str());
    return res;
}

extern "C" {

cas_context_t* cas_create_context() {
    return new cas_context_t();
}

void cas_destroy_context(cas_context_t* ctx) {
    delete ctx;
}

cas_compute_result_t cas_simplify(cas_context_t* ctx, const char* input) {
    auto res = ctx->session.simplify(input);
    cas_compute_result_t out;
    out.ok = res.ok;
    out.error = cas_strdup(res.error);
    out.result_latex = cas_strdup(res.latex);
    out.result_ascii = cas_strdup(res.ascii);
    out.numeric_value = res.numeric_value.value_or(0.0);
    out.has_numeric = res.numeric_value.has_value();
    return out;
}

void cas_free_result(cas_compute_result_t res) {
    free((void*)res.error);
    free((void*)res.result_latex);
    free((void*)res.result_ascii);
}

}
