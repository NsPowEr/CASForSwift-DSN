#pragma once
#include <string>
#include <string_view>
#include "cas/result.hpp"

namespace cas::gui {

/**
 * @brief Bridge nativo verso il kernel simbolico locale (Giac/embedded).
 * In un'implementazione reale, questo include <giac/giac.h>.
 */
class EmbeddedGiac {
public:
    EmbeddedGiac();
    ~EmbeddedGiac();

    [[nodiscard]] std::string evaluate(std::string_view input);
    
    // Altre funzioni Giac-specifiche (es. plot, integrali con step)
    [[nodiscard]] std::string version() const;

private:
    // struct GiacInternal; 
    // std::unique_ptr<GiacInternal> m_impl;
};

} // namespace cas::gui
