#pragma once

#include "cas/ast.hpp"

#include <string>

namespace cas {

[[nodiscard]] std::string debug_print(ExprPtr expr);

}  // namespace cas
