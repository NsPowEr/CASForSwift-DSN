#pragma once
// matrix_adapter.hpp — corpus `[[…]]` notation + Maxima `matrix(…)`
// output mapped to cas::linalg::MatrixExpr and the standard linalg
// dispatch (det/trace/transpose/rank/inverse/eigenvalues).
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Matrix_Adapter.md (F7.5.A2).

#include "corpus_runner.hpp"
#include "maxima_parser.hpp"
#include "solve_set_equal.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cas::golden {

// ---------------------------------------------------------------------------
// Split balanced bracketed text "[a, b, [c, d], …]" into top-level
// elements, respecting nested (), [], {}. Strips outer brackets first.
// Returns empty vector on malformed input.
// ---------------------------------------------------------------------------
inline std::vector<std::string> split_brackets(const std::string& raw) {
    return cas::golden::split_maxima_list(raw);  // reuse F7.5.A1 helper
}

// ---------------------------------------------------------------------------
// Parse "[[a,b,…],[c,d,…],…]" → MatrixExpr. Each element parsed via
// existing CAS parser (parse_expr).
// ---------------------------------------------------------------------------
inline Result<cas::linalg::MatrixExpr> parse_matrix_lit(
    const std::string& raw, cas::symbolic::CASContext& ctx) {
    auto rows = split_brackets(raw);
    if (rows.empty()) {
        return CASError{CASErrorKind::Unimplemented,
                        "parse_matrix_lit: empty or malformed: " + raw,
                        std::nullopt};
    }
    std::vector<std::vector<ExprPtr>> grid;
    grid.reserve(rows.size());
    std::size_t cols = 0;
    for (const auto& row : rows) {
        auto elems = split_brackets(row);
        if (elems.empty()) {
            return CASError{CASErrorKind::Unimplemented,
                            "parse_matrix_lit: empty row: " + row,
                            std::nullopt};
        }
        if (cols == 0) cols = elems.size();
        else if (elems.size() != cols) {
            return CASError{CASErrorKind::Unimplemented,
                            "parse_matrix_lit: ragged rows",
                            std::nullopt};
        }
        std::vector<ExprPtr> row_vec;
        row_vec.reserve(elems.size());
        for (const auto& e : elems) {
            auto p = cas::golden::parse_expr(e, ctx);
            if (!p.is_ok()) return p.error();
            row_vec.push_back(p.value());
        }
        grid.push_back(std::move(row_vec));
    }
    cas::linalg::MatrixExpr m(grid.size(), cols);
    for (std::size_t i = 0; i < grid.size(); ++i)
        for (std::size_t j = 0; j < cols; ++j)
            m(i, j) = grid[i][j];
    return ok(std::move(m));
}

// ---------------------------------------------------------------------------
// Parse Maxima output `matrix([a,b],[c,d])` → MatrixExpr.
// ---------------------------------------------------------------------------
inline Result<cas::linalg::MatrixExpr> parse_maxima_matrix(
    const std::string& raw_in, cas::symbolic::CASContext& ctx) {
    std::string s = raw_in;
    // Strip leading/trailing whitespace.
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t' ||
                          s.front() == '\n' || s.front() == '\r'))
        s.erase(s.begin());
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' ||
                          s.back() == '\n' || s.back() == '\r' ||
                          s.back() == ';' || s.back() == '$'))
        s.pop_back();
    if (s.rfind("matrix(", 0) != 0 || s.empty() || s.back() != ')') {
        return CASError{CASErrorKind::Unimplemented,
                        "parse_maxima_matrix: not a matrix(...) form: " + s,
                        std::nullopt};
    }
    std::string inner = s.substr(7, s.size() - 8);
    // Top-level split on commas at depth 0 yields the [row] elements.
    int depth = 0;
    std::vector<std::string> rows;
    std::string cur;
    for (char c : inner) {
        if (c == '(' || c == '[' || c == '{') ++depth;
        else if (c == ')' || c == ']' || c == '}') --depth;
        if (c == ',' && depth == 0) {
            rows.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!cur.empty()) rows.push_back(cur);

    std::vector<std::vector<ExprPtr>> grid;
    grid.reserve(rows.size());
    std::size_t cols = 0;
    for (const auto& row : rows) {
        auto elems = split_brackets(row);
        if (elems.empty()) {
            return CASError{CASErrorKind::Unimplemented,
                            "parse_maxima_matrix: empty row",
                            std::nullopt};
        }
        if (cols == 0) cols = elems.size();
        else if (elems.size() != cols) {
            return CASError{CASErrorKind::Unimplemented,
                            "parse_maxima_matrix: ragged",
                            std::nullopt};
        }
        std::vector<ExprPtr> rv;
        rv.reserve(elems.size());
        for (const auto& e : elems) {
            auto normalized = normalize_maxima_output(e);
            if (normalized.empty()) {
                return CASError{CASErrorKind::Unimplemented,
                                "parse_maxima_matrix: element normalize empty",
                                std::nullopt};
            }
            auto p = parse_maxima_expr(normalized, ctx.arena());
            if (!p.is_ok()) return p.error();
            rv.push_back(p.value());
        }
        grid.push_back(std::move(rv));
    }
    cas::linalg::MatrixExpr m(grid.size(), cols);
    for (std::size_t i = 0; i < grid.size(); ++i)
        for (std::size_t j = 0; j < cols; ++j)
            m(i, j) = grid[i][j];
    return ok(std::move(m));
}

// ---------------------------------------------------------------------------
// Compare two matrices element-wise via mathematically_equal.
// ---------------------------------------------------------------------------
inline Result<bool> compare_matrices(
    const cas::linalg::MatrixExpr& a, const cas::linalg::MatrixExpr& b,
    cas::symbolic::CASContext& ctx) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return ok(false);
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t j = 0; j < a.cols(); ++j) {
            auto eq = cas::symbolic::mathematically_equal(a(i, j), b(i, j), ctx);
            if (!eq.is_ok()) return eq.error();
            if (!eq.value()) return ok(false);
        }
    }
    return ok(true);
}

// ---------------------------------------------------------------------------
// Output of a matrix CAS dispatch: tagged union (scalar | matrix |
// eigenvalues list). Default constructed = NONE.
// ---------------------------------------------------------------------------
struct MatrixDispatchResult {
    enum class Kind { None, Scalar, Matrix, Eigenvalues };
    Kind kind{Kind::None};
    ExprPtr scalar;
    cas::linalg::MatrixExpr matrix{0, 0};
    std::vector<ExprPtr> eigenvalues_list;
};

// ---------------------------------------------------------------------------
// HC-F75-A2-MATRIX-SCALAR-OP (F7.5.A2 closure):
//
// A recursive, precedence-aware top-level evaluator for expressions that
// mix matrix literals `[[…]]` with scalar subexpressions, using the
// existing CAS parser for scalar operands and `cas::linalg` for matrix
// operations. No closed pattern table: any well-formed +/-/* tree is
// handled uniformly, and unsupported operand combinations (e.g.
// scalar + matrix) yield an explicit `Unimplemented` diagnostic — never
// silent skips.
// ---------------------------------------------------------------------------
struct MatrixOrScalar {
    bool is_matrix{false};
    cas::linalg::MatrixExpr matrix{0, 0};
    ExprPtr scalar{nullptr};
};

// Trim ASCII whitespace at both ends and strip a single matching outer
// pair of parentheses if it brackets the whole string (depth-aware).
inline std::string strip_outer_parens(std::string s) {
    auto is_ws = [](char c) {
        return c == ' ' || c == '\t' || c == '\n' || c == '\r';
    };
    while (!s.empty() && is_ws(s.front())) s.erase(s.begin());
    while (!s.empty() && is_ws(s.back())) s.pop_back();
    while (s.size() >= 2 && s.front() == '(' && s.back() == ')') {
        int depth = 0;
        bool wraps_whole = true;
        for (std::size_t i = 0; i < s.size(); ++i) {
            char c = s[i];
            if (c == '(' || c == '[' || c == '{') ++depth;
            else if (c == ')' || c == ']' || c == '}') --depth;
            if (depth == 0 && i + 1 < s.size()) { wraps_whole = false; break; }
        }
        if (!wraps_whole) break;
        s = s.substr(1, s.size() - 2);
        while (!s.empty() && is_ws(s.front())) s.erase(s.begin());
        while (!s.empty() && is_ws(s.back())) s.pop_back();
    }
    return s;
}

// Scan `s` for the rightmost top-level operator of the lowest available
// precedence. Precedence levels: +,- (lowest); *,/. Unary +/- (operator
// at position 0 or immediately following another operator or an opening
// bracket) is skipped — it belongs to the operand.
//
// Returns {op, idx} or nullopt if no top-level binary operator exists.
struct TopBinop {
    char op;
    std::size_t pos;
};

inline std::optional<TopBinop> find_top_level_binop(const std::string& s) {
    auto is_ws = [](char c) {
        return c == ' ' || c == '\t' || c == '\n' || c == '\r';
    };
    auto is_op = [](char c) {
        return c == '+' || c == '-' || c == '*' || c == '/';
    };
    // Pass A: rightmost + or -.
    // Pass B: rightmost * or /.
    for (int precedence = 0; precedence < 2; ++precedence) {
        std::optional<TopBinop> best;
        int depth = 0;
        char prev_nonws = '\0';
        for (std::size_t i = 0; i < s.size(); ++i) {
            char c = s[i];
            if (c == '(' || c == '[' || c == '{') { ++depth; prev_nonws = c; continue; }
            if (c == ')' || c == ']' || c == '}') { --depth; prev_nonws = c; continue; }
            if (depth == 0 && is_op(c)) {
                const bool is_low = (c == '+' || c == '-');
                const bool match  = (precedence == 0 ? is_low : !is_low);
                if (match) {
                    // Unary check: if prev_nonws is empty, an operator, or
                    // an opening bracket, treat as unary.
                    const bool unary =
                        prev_nonws == '\0' || is_op(prev_nonws) ||
                        prev_nonws == '(' || prev_nonws == '[' ||
                        prev_nonws == '{' || prev_nonws == ',';
                    if (!unary) best = TopBinop{c, i};
                }
            }
            if (!is_ws(c)) prev_nonws = c;
        }
        if (best.has_value()) return best;
    }
    return std::nullopt;
}

[[nodiscard]] Result<MatrixOrScalar> evaluate_matrix_expression(
    const std::string& raw, cas::symbolic::CASContext& ctx);

// Build a scalar*matrix product element-wise via CAS simplifier.
inline Result<cas::linalg::MatrixExpr> matrix_scalar_multiply(
    const cas::linalg::MatrixExpr& m, ExprPtr s,
    cas::symbolic::CASContext& ctx) {
    cas::linalg::MatrixExpr out(m.rows(), m.cols());
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            ExprPtr prod = ctx.arena().make<Product>(
                std::vector<ExprPtr>{s, m(i, j)});
            auto r = ctx.simplify(prod);
            if (!r.is_ok()) return r.error();
            out(i, j) = r.value();
        }
    }
    return ok(std::move(out));
}

inline Result<MatrixOrScalar> evaluate_matrix_expression_impl(
    const std::string& raw, cas::symbolic::CASContext& ctx) {
    std::string s = strip_outer_parens(raw);
    if (s.empty()) {
        return CASError{CASErrorKind::Unimplemented,
                        "evaluate_matrix_expression: empty input",
                        std::nullopt};
    }
    // Binary split if a top-level binop exists.
    if (auto op = find_top_level_binop(s); op.has_value()) {
        std::string l_raw = s.substr(0, op->pos);
        std::string r_raw = s.substr(op->pos + 1);
        auto lhs = evaluate_matrix_expression(l_raw, ctx);
        if (!lhs.is_ok()) return lhs.error();
        auto rhs = evaluate_matrix_expression(r_raw, ctx);
        if (!rhs.is_ok()) return rhs.error();
        const auto& L = lhs.value();
        const auto& R = rhs.value();
        switch (op->op) {
            case '+':
            case '-': {
                if (!L.is_matrix || !R.is_matrix) {
                    if (!L.is_matrix && !R.is_matrix) {
                        // Pure scalar; combine via Sum/Difference and simplify.
                        ExprPtr neg_r = (op->op == '-')
                            ? ctx.arena().make<Unary>(UnaryOp::Neg, R.scalar)
                            : R.scalar;
                        ExprPtr sum = ctx.arena().make<Sum>(
                            std::vector<ExprPtr>{L.scalar, neg_r});
                        auto si = ctx.simplify(sum);
                        if (!si.is_ok()) return si.error();
                        MatrixOrScalar out;
                        out.is_matrix = false;
                        out.scalar = si.value();
                        return ok(std::move(out));
                    }
                    return CASError{CASErrorKind::Unimplemented,
                        "evaluate_matrix_expression: scalar±matrix is undefined",
                        std::nullopt};
                }
                auto r = (op->op == '+')
                    ? cas::linalg::add(L.matrix, R.matrix, ctx)
                    : cas::linalg::subtract(L.matrix, R.matrix, ctx);
                if (!r.is_ok()) return r.error();
                MatrixOrScalar out;
                out.is_matrix = true;
                out.matrix = r.value();
                return ok(std::move(out));
            }
            case '*': {
                if (L.is_matrix && R.is_matrix) {
                    auto r = cas::linalg::multiply(L.matrix, R.matrix, ctx);
                    if (!r.is_ok()) return r.error();
                    MatrixOrScalar out;
                    out.is_matrix = true;
                    out.matrix = r.value();
                    return ok(std::move(out));
                }
                if (L.is_matrix && !R.is_matrix) {
                    auto r = matrix_scalar_multiply(L.matrix, R.scalar, ctx);
                    if (!r.is_ok()) return r.error();
                    MatrixOrScalar out;
                    out.is_matrix = true;
                    out.matrix = r.value();
                    return ok(std::move(out));
                }
                if (!L.is_matrix && R.is_matrix) {
                    auto r = matrix_scalar_multiply(R.matrix, L.scalar, ctx);
                    if (!r.is_ok()) return r.error();
                    MatrixOrScalar out;
                    out.is_matrix = true;
                    out.matrix = r.value();
                    return ok(std::move(out));
                }
                // scalar * scalar — defer to CAS simplifier.
                ExprPtr prod = ctx.arena().make<Product>(
                    std::vector<ExprPtr>{L.scalar, R.scalar});
                auto si = ctx.simplify(prod);
                if (!si.is_ok()) return si.error();
                MatrixOrScalar out;
                out.is_matrix = false;
                out.scalar = si.value();
                return ok(std::move(out));
            }
            case '/': {
                if (R.is_matrix) {
                    return CASError{CASErrorKind::Unimplemented,
                        "evaluate_matrix_expression: division by matrix is not "
                        "supported in this adapter (use inverse() explicitly)",
                        std::nullopt};
                }
                if (L.is_matrix) {
                    ExprPtr inv = ctx.arena().make<Binary>(
                        BinaryOp::Pow, R.scalar,
                        ctx.arena().make<IntegerLit>(BigInt(-1)));
                    auto inv_s = ctx.simplify(inv);
                    if (!inv_s.is_ok()) return inv_s.error();
                    auto r = matrix_scalar_multiply(L.matrix, inv_s.value(), ctx);
                    if (!r.is_ok()) return r.error();
                    MatrixOrScalar out;
                    out.is_matrix = true;
                    out.matrix = r.value();
                    return ok(std::move(out));
                }
                // scalar / scalar.
                ExprPtr div = ctx.arena().make<Binary>(
                    BinaryOp::Div, L.scalar, R.scalar);
                auto si = ctx.simplify(div);
                if (!si.is_ok()) return si.error();
                MatrixOrScalar out;
                out.is_matrix = false;
                out.scalar = si.value();
                return ok(std::move(out));
            }
        }
    }
    // Unary +/- prefix (after strip_outer_parens leaves leading sign):
    // recursively evaluate the rest and apply the sign at the result level
    // so a `-[[…]]` form still produces a Matrix-tagged value.
    if (s.front() == '+') {
        return evaluate_matrix_expression(s.substr(1), ctx);
    }
    if (s.front() == '-') {
        auto inner = evaluate_matrix_expression(s.substr(1), ctx);
        if (!inner.is_ok()) return inner.error();
        ExprPtr neg_one = ctx.arena().make<IntegerLit>(BigInt(-1));
        if (inner.value().is_matrix) {
            auto neg = matrix_scalar_multiply(inner.value().matrix, neg_one, ctx);
            if (!neg.is_ok()) return neg.error();
            MatrixOrScalar out;
            out.is_matrix = true;
            out.matrix = neg.value();
            return ok(std::move(out));
        }
        ExprPtr negated = ctx.arena().make<Unary>(UnaryOp::Neg, inner.value().scalar);
        auto si = ctx.simplify(negated);
        if (!si.is_ok()) return si.error();
        MatrixOrScalar out;
        out.is_matrix = false;
        out.scalar = si.value();
        return ok(std::move(out));
    }
    // Leaf: either a `[[…]]` matrix literal or a scalar subexpression.
    if (s.size() >= 2 && s.front() == '[' && s[1] == '[') {
        auto m = parse_matrix_lit(s, ctx);
        if (!m.is_ok()) return m.error();
        MatrixOrScalar out;
        out.is_matrix = true;
        out.matrix = m.value();
        return ok(std::move(out));
    }
    auto p = cas::golden::parse_expr(s, ctx);
    if (!p.is_ok()) return p.error();
    MatrixOrScalar out;
    out.is_matrix = false;
    out.scalar = p.value();
    return ok(std::move(out));
}

inline Result<MatrixOrScalar> evaluate_matrix_expression(
    const std::string& raw, cas::symbolic::CASContext& ctx) {
    return evaluate_matrix_expression_impl(raw, ctx);
}

// ---------------------------------------------------------------------------
// Dispatch a matrix command `fn([[…]])`. Returns the appropriate tag.
// ---------------------------------------------------------------------------
inline Result<MatrixDispatchResult> evaluate_cas_matrix(
    const std::string& input_str, cas::symbolic::CASContext& ctx) {
    auto cmd = cas::golden::parse_command(input_str);
    // If parse_command did not recognize a function call wrapper, treat the
    // whole input as a top-level expression mixing matrix literals and
    // scalars (HC-F75-A2-MATRIX-SCALAR-OP).
    if (cmd.fn.empty() || cmd.arg_strs.empty()) {
        auto val = evaluate_matrix_expression(input_str, ctx);
        if (!val.is_ok()) return val.error();
        MatrixDispatchResult out;
        if (val.value().is_matrix) {
            out.kind = MatrixDispatchResult::Kind::Matrix;
            out.matrix = val.value().matrix;
        } else {
            out.kind = MatrixDispatchResult::Kind::Scalar;
            out.scalar = val.value().scalar;
        }
        return ok(std::move(out));
    }
    if (cmd.arg_strs.size() != 1U) {
        return CASError{CASErrorKind::Unimplemented,
                        "evaluate_cas_matrix: expected 1 arg",
                        std::nullopt};
    }
    auto m_res = parse_matrix_lit(cmd.arg_strs[0], ctx);
    if (!m_res.is_ok()) return m_res.error();
    const auto& m = m_res.value();
    MatrixDispatchResult out;
    if (cmd.fn == "det") {
        auto r = cas::linalg::determinant(m, ctx);
        if (!r.is_ok()) return r.error();
        out.kind = MatrixDispatchResult::Kind::Scalar;
        out.scalar = r.value();
        return ok(std::move(out));
    }
    if (cmd.fn == "trace") {
        auto r = cas::linalg::trace(m, ctx);
        if (!r.is_ok()) return r.error();
        out.kind = MatrixDispatchResult::Kind::Scalar;
        out.scalar = r.value();
        return ok(std::move(out));
    }
    if (cmd.fn == "rank") {
        auto r = cas::linalg::rank(m, ctx);
        if (!r.is_ok()) return r.error();
        out.kind = MatrixDispatchResult::Kind::Scalar;
        out.scalar = ctx.arena().make<IntegerLit>(BigInt(static_cast<long long>(r.value())));
        return ok(std::move(out));
    }
    if (cmd.fn == "transpose") {
        auto r = cas::linalg::transpose(m);
        if (!r.is_ok()) return r.error();
        out.kind = MatrixDispatchResult::Kind::Matrix;
        out.matrix = r.value();
        return ok(std::move(out));
    }
    if (cmd.fn == "inverse") {
        auto r = cas::linalg::inverse(m, ctx);
        if (!r.is_ok()) return r.error();
        out.kind = MatrixDispatchResult::Kind::Matrix;
        out.matrix = r.value();
        return ok(std::move(out));
    }
    if (cmd.fn == "eigenvalues") {
        auto r = cas::linalg::eigenvalues(m, ctx);
        if (!r.is_ok()) return r.error();
        out.kind = MatrixDispatchResult::Kind::Eigenvalues;
        out.eigenvalues_list = r.value();
        return ok(std::move(out));
    }
    return CASError{CASErrorKind::Unimplemented,
                    "evaluate_cas_matrix: unknown fn " + cmd.fn,
                    std::nullopt};
}

}  // namespace cas::golden
