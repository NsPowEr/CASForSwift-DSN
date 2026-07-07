#pragma once
// error_helpers.hpp — typed Result<T> helper for make_unimplemented.
// Included at the bottom of error.hpp so it can use both CASError and Result<T>.

#include "cas/result.hpp"

namespace cas {

/// Create a Result<T> carrying a structured Unimplemented error.
/// This is the preferred call-site helper for F0.8-migrated code.
template <typename T>
[[nodiscard]] Result<T> make_unimplemented(
    UnimplementedInfo info,
    std::string brief_message = {})
{
    return fail<T>(make_unimplemented_error(std::move(info), std::move(brief_message)));
}

template <typename T>
[[nodiscard]] Result<T> make_unimplemented(
    std::string module,
    std::string function,
    std::string input_shape,
    std::string reason,
    std::string suggestion,
    std::string ticket = {},
    std::string brief_message = {})
{
    return fail<T>(make_unimplemented_error(
        std::move(module),
        std::move(function),
        std::move(input_shape),
        std::move(reason),
        std::move(suggestion),
        std::move(ticket),
        std::move(brief_message)));
}

}  // namespace cas
