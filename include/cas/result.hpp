#pragma once

#include "cas/error.hpp"

#include <optional>
#include <type_traits>
#include <utility>

namespace cas {

template <typename T>
class Result {
public:
    Result(const T& value) : value_(value) {}
    Result(T&& value) : value_(std::move(value)) {}
    Result(const CASError& error) : error_(error) {}
    Result(CASError&& error) : error_(std::move(error)) {}

    [[nodiscard]] bool is_ok() const noexcept { return value_.has_value(); }
    [[nodiscard]] bool is_error() const noexcept { return error_.has_value(); }

    [[nodiscard]] const T& value() const& { return *value_; }
    [[nodiscard]] T& value() & { return *value_; }
    [[nodiscard]] T&& value() && { return std::move(*value_); }

    [[nodiscard]] const CASError& error() const& { return *error_; }
    [[nodiscard]] CASError& error() & { return *error_; }

    [[nodiscard]] static Result<T> unimplemented(UnimplementedInfo info, std::string brief_message = {}) {
        return Result<T>(make_unimplemented_error(std::move(info), std::move(brief_message)));
    }

private:
    std::optional<T> value_;
    std::optional<CASError> error_;
};

template <>
class Result<void> {
public:
    Result() : ok_(true) {}
    Result(const CASError& error) : ok_(false), error_(error) {}
    Result(CASError&& error) : ok_(false), error_(std::move(error)) {}

    [[nodiscard]] bool is_ok() const noexcept { return ok_; }
    [[nodiscard]] bool is_error() const noexcept { return !ok_; }
    [[nodiscard]] const CASError& error() const& { return *error_; }

    [[nodiscard]] static Result<void> unimplemented(UnimplementedInfo info, std::string brief_message = {}) {
        return Result<void>(make_unimplemented_error(std::move(info), std::move(brief_message)));
    }

private:
    bool ok_{false};
    std::optional<CASError> error_;
};

template <typename T>
[[nodiscard]] Result<std::decay_t<T>> ok(T&& value) {
    return Result<std::decay_t<T>>(std::forward<T>(value));
}

inline Result<void> ok() {
    return Result<void>();
}

template <typename T>
[[nodiscard]] Result<T> fail(CASError error) {
    return Result<T>(std::move(error));
}

}  // namespace cas
