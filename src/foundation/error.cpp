// error.cpp — implementation of CASError::format_user_message()
// F0.8 Error diagnostic framework.

#include "cas/error.hpp"

#include <sstream>

namespace cas {

std::string CASError::format_user_message() const {
    if (!payload.has_value()) {
        // Backward-compatible fallback: just return the plain message.
        return message;
    }

    const UnimplementedInfo& p = *payload;
    std::ostringstream oss;
    oss << "[Unimplemented] module=" << p.module
        << " function=" << p.function << "\n"
        << "  Input shape: " << p.input_shape << "\n"
        << "  Reason: "      << p.reason << "\n"
        << "  Suggestion: "  << p.suggestion;
    if (!p.ticket.empty()) {
        oss << "\n  Ticket: " << p.ticket;
    }
    return oss.str();
}

}  // namespace cas
