/// @file ast.hpp
/// @brief Backward-compatible umbrella header for the AST subsystem.
///
/// **F8.0 / Task 1.1 — Monolith Split**
/// The AST public interface has been split into three focused headers:
///
///   - `cas/ast_kinds.hpp`  — ExprKind, MathConstant, UnaryOp, BinaryOp, LimitDirection
///   - `cas/ast_nodes.hpp`  — ExprNode, ExprPtr, all concrete node types, cast helpers
///   - `cas/ast_arena.hpp`  — AstArena (interning allocator), visit_expr<Visitor>
///
/// Existing code that includes `cas/ast.hpp` continues to work unchanged because
/// this file simply includes all three sub-headers.  New code that only needs a
/// subset (e.g. `AlgebraicNumber` needs only `ast_kinds.hpp` + `ast_nodes.hpp`)
/// can include the narrower header to avoid pulling in the full Arena machinery.
#pragma once

#include "cas/ast_kinds.hpp"
#include "cas/ast_nodes.hpp"
#include "cas/ast_arena.hpp"
