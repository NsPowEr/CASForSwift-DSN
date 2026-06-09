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
// Dispatch a matrix command `fn([[…]])`. Returns the appropriate tag.
// ---------------------------------------------------------------------------
inline Result<MatrixDispatchResult> evaluate_cas_matrix(
    const std::string& input_str, cas::symbolic::CASContext& ctx) {
    auto cmd = cas::golden::parse_command(input_str);
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
