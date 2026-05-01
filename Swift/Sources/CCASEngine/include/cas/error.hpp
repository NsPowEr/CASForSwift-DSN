#pragma once

#include <cstdint>
#include <optional>
#include <string>

namespace cas {

enum class CASErrorKind : std::uint8_t {
    Unsupported,
    Unimplemented,
    Undefined,
    Overflow,
    Timeout,
    ParseError,
    AssumptionConflict,
    InvalidArgument,
    InternalError,
};

struct CASError {
    static constexpr std::uint32_t API_VERSION = 1;

    CASErrorKind kind{CASErrorKind::InternalError};
    std::string message;
    std::optional<std::string> hint;
};

}  // namespace cas
