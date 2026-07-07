#pragma once

#include <string>

namespace cas {

/// Structured diagnostic payload for Unimplemented errors.
struct UnimplementedInfo {
    std::string module;       ///< e.g. "calculus", "algebra"
    std::string function;     ///< enclosing C++ function name
    std::string input_shape;  ///< human-readable descriptor of the input form
    std::string reason;       ///< reason string or reason code constant
    std::string suggestion;   ///< actionable next step for the user/developer
    std::string ticket;       ///< e.g. "F0.8", "L3-04" — empty string if no ticket

    // Backward compatibility aliases for existing code/tests
    [[nodiscard]] const std::string& reason_code() const noexcept { return reason; }
    [[nodiscard]] const std::string& ticket_id() const noexcept { return ticket; }
};

using UnimplementedPayload = UnimplementedInfo;

}  // namespace cas
