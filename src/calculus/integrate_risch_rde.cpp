// integrate_risch_rde.cpp — Risch Differential Equation solvers over Q(x).
// Bronstein "Symbolic Integration I", Chapters 5–6.
//
// Public API:
//   solve_risch_de_q()   — dispatcher (polynomial fast-path then rational)
//
// Private (anonymous namespace):
//   solve_risch_de_poly_q()     — polynomial-coefficient case (§5.6)
//   solve_risch_de_rational_q() — rational-coefficient case  (§6.1)
//   risch_as_rational / risch_extract_rational_coeffs / risch_build_poly_expr
//                                — coefficient helpers

#include "integrate_risch_internal.hpp"
#include "integrate_engine.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

// Extract a Rational from an expression that simplifies to a rational literal.
[[nodiscard]] std::optional<Rational> risch_as_rational(ExprPtr e) {
    if (const auto* lit = expr_cast<IntegerLit>(e)) return Rational(lit->value);
    if (const auto* lit = expr_cast<RationalLit>(e)) {
        return Rational(lit->numerator, lit->denominator);
    }
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (auto inner = risch_as_rational(un->operand)) return -*inner;
    }
    return std::nullopt;
}

// Convert a PolyExpr (ExprPtr coefficients) to a vector<Rational> when every
// coefficient is a literal rational.  Returns std::nullopt if any coefficient
// depends on a symbol, on transcendentals, or on extensions.
[[nodiscard]] std::optional<std::vector<Rational>>
risch_extract_rational_coeffs(const algebra::PolyExpr& poly) {
    std::vector<Rational> out;
    out.reserve(poly.size());
    for (ExprPtr c : poly.coefficients()) {
        if (auto r = risch_as_rational(c)) {
            out.push_back(*r);
        } else {
            return std::nullopt;
        }
    }
    return out;
}

// Build an ExprPtr representing Σ_k c_k · var^k from rational coefficients.
[[nodiscard]] ExprPtr risch_build_poly_expr(
    AstArena& arena, const std::vector<Rational>& coeffs, const Symbol& var) {
    std::vector<ExprPtr> terms;
    terms.reserve(coeffs.size());
    ExprPtr var_sym = arena.make<Symbol>(var);
    for (std::size_t k = 0; k < coeffs.size(); ++k) {
        const Rational& c = coeffs[k];
        if (c.numerator().is_zero()) continue;
        ExprPtr c_expr = (c.denominator() == BigInt(1))
            ? static_cast<ExprPtr>(arena.make<IntegerLit>(c.numerator()))
            : static_cast<ExprPtr>(arena.make<RationalLit>(c.numerator(), c.denominator()));
        if (k == 0U) { terms.push_back(c_expr); continue; }
        ExprPtr xk = (k == 1U)
            ? var_sym
            : static_cast<ExprPtr>(arena.make<Binary>(BinaryOp::Pow,
                  var_sym,
                  arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k)))));
        if (c.numerator() == BigInt(1) && c.denominator() == BigInt(1)) {
            terms.push_back(xk);
        } else {
            terms.push_back(arena.make<Product>(std::vector<ExprPtr>{c_expr, xk}));
        }
    }
    if (terms.empty()) return arena.make<IntegerLit>(BigInt(0));
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

// Solve  y' + f · y = g  for  y ∈ Q[var]  (polynomial-coefficient case).
// Bronstein §5.6 degree-bound + linear-system inversion.
// Degree analysis: if deg(f)≥1 then deg(y)=deg(g)-deg(f); if deg(f)=0
// then deg(y)=deg(g); if f≡0 then deg(y)=deg(g)+1 (pure quadrature).
// Gaussian elimination over Q; inconsistency or singular ⇒ Unimplemented.
[[nodiscard]] Result<ExprPtr> solve_risch_de_poly_q(
    ExprPtr f_expr, ExprPtr g_expr, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    // F0.8-MIGRATED
    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "solve_risch_de_poly_q",
            msg,
            cas::error::reason_codes::RISCH_NO_POLYNOMIAL_SOLUTION,
            "Risch DE poly: extend coefficient solver or fall back to exponential branch",
            "F0.8");
    };

    auto f_poly_res = algebra::parse_polynomial(f_expr, var, ctx);
    auto g_poly_res = algebra::parse_polynomial(g_expr, var, ctx);
    if (f_poly_res.is_error()) return fail_unimpl("Risch DE: f is not a polynomial in the base variable");
    if (g_poly_res.is_error()) return fail_unimpl("Risch DE: g is not a polynomial in the base variable");

    auto f_opt = risch_extract_rational_coeffs(f_poly_res.value());
    auto g_opt = risch_extract_rational_coeffs(g_poly_res.value());
    if (!f_opt.has_value()) return fail_unimpl("Risch DE: f has non-rational coefficients");
    if (!g_opt.has_value()) return fail_unimpl("Risch DE: g has non-rational coefficients");

    std::vector<Rational> f = std::move(*f_opt);
    std::vector<Rational> g = std::move(*g_opt);

    auto strip_trailing_zeros = [](std::vector<Rational>& v) {
        while (!v.empty() && v.back().numerator().is_zero()) v.pop_back();
    };
    strip_trailing_zeros(f);
    strip_trailing_zeros(g);

    if (g.empty()) {
        // g ≡ 0: y = 0 is a solution.
        return ok(arena.make<IntegerLit>(BigInt(0)));
    }

    const int deg_g = static_cast<int>(g.size()) - 1;
    int deg_f = static_cast<int>(f.size()) - 1; // -1 when f ≡ 0.
    int e;
    if (f.empty()) {
        // f ≡ 0:  y' = g  ⇒  y_{k+1} = g_k / (k+1) ; y_0 free, take 0.
        std::vector<Rational> y(g.size() + 1, Rational(BigInt(0)));
        for (std::size_t k = 0; k < g.size(); ++k) {
            y[k + 1] = g[k] / Rational(BigInt(static_cast<std::int64_t>(k + 1)));
        }
        return ok(risch_build_poly_expr(arena, y, var));
    }
    if (deg_f == 0) {
        e = deg_g;
    } else {
        if (deg_g < deg_f) return fail_unimpl("Risch DE: no polynomial solution (deg g < deg f)");
        e = deg_g - deg_f;
    }
    if (e < 0) return fail_unimpl("Risch DE: negative degree bound");

    // Build dense linear system M · y_vec = g_vec over Q.
    // Equation index j ∈ [0, deg_g] picks out coefficient of var^j of
    // y' + f·y - g.  Unknown vector y_vec = (y_0, …, y_e).
    const std::size_t n_unk = static_cast<std::size_t>(e) + 1U;
    const std::size_t n_eq = static_cast<std::size_t>(deg_g) + 1U;
    std::vector<std::vector<Rational>> M(n_eq,
        std::vector<Rational>(n_unk + 1U, Rational(BigInt(0))));
    for (std::size_t j = 0; j < n_eq; ++j) {
        // RHS: g[j].
        M[j][n_unk] = (j < g.size()) ? g[j] : Rational(BigInt(0));
        // y' contribution: coefficient of var^j of y' is (j+1)·y_{j+1}.
        if (j + 1U <= static_cast<std::size_t>(e)) {
            M[j][j + 1U] = M[j][j + 1U] + Rational(BigInt(static_cast<std::int64_t>(j + 1U)));
        }
        // f·y contribution: coefficient of var^j is Σ_k f_k · y_{j-k}.
        for (std::size_t k = 0; k <= static_cast<std::size_t>(deg_f); ++k) {
            if (k > j) break;
            std::size_t idx = j - k;
            if (idx > static_cast<std::size_t>(e)) continue;
            M[j][idx] = M[j][idx] + f[k];
        }
    }

    // Gaussian elimination with partial pivoting over Q.  Pure Rational/BigInt
    // arithmetic, no simplify()/substitute() in the body → not transitively
    // interruptible: poll per pivot column (A2 / HC-F70-A33, matrix_bareiss
    // precedent).  Wires the previously-dead timeout stub.
    const std::size_t cols = n_unk + 1U;
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < n_unk && pivot_row < n_eq; ++col) {
        if (auto chk = ctx.check_interrupt(); chk.is_error()) return fail<ExprPtr>(chk.error());
        std::size_t best = pivot_row;
        while (best < n_eq && M[best][col].numerator().is_zero()) ++best;
        if (best == n_eq) {
            // No pivot in this column — unknown stays free.  Force to zero
            // by recording no constraint; we still need a unique answer, so
            // declare Unimplemented if any equation in lower rows depends
            // on this unknown.  Detection happens by full RREF in a future
            // generalisation; here we treat singular as Unimplemented.
            return fail_unimpl("Risch DE: singular linear system (free variable in y)");
        }
        if (best != pivot_row) std::swap(M[best], M[pivot_row]);
        Rational pivot = M[pivot_row][col];
        for (std::size_t c = col; c < cols; ++c) M[pivot_row][c] = M[pivot_row][c] / pivot;
        for (std::size_t r = 0; r < n_eq; ++r) {
            if (r == pivot_row) continue;
            Rational factor = M[r][col];
            if (factor.numerator().is_zero()) continue;
            for (std::size_t c = col; c < cols; ++c) {
                M[r][c] = M[r][c] - factor * M[pivot_row][c];
            }
        }
        ++pivot_row;
    }

    // Check residual equations beyond pivot_row: must all be zero for
    // consistency.
    for (std::size_t r = pivot_row; r < n_eq; ++r) {
        if (!M[r][n_unk].numerator().is_zero()) {
            return fail_unimpl("Risch DE: inconsistent linear system (no polynomial solution)");
        }
    }

    // Read solution: row r has pivot in column r (after RREF on pivot rows).
    std::vector<Rational> y(n_unk, Rational(BigInt(0)));
    for (std::size_t r = 0; r < std::min(pivot_row, n_unk); ++r) {
        y[r] = M[r][n_unk];
    }

    return ok(risch_build_poly_expr(arena, y, var));
}

// Solve  y' + f · y = g  for  y ∈ Q[var]  with f, g RATIONAL in `var`.
//
// Strategy (Bronstein §6.1 simplified for the polynomial-y case):
//   1. Decompose f and g into num/den via `algebra::apart_num_den`.
//   2. Compute the common denominator D = lcm(den(f), den(g)) (here taken
//      as den(f)·den(g) divided by gcd; falls back to product when gcd is
//      a unit).  Multiply through:
//          D · y' + f_n · y = g_n
//      where  f_n = f · D    and   g_n = g · D    are polynomials in var.
//   3. Hypothesise y ∈ Q[var] of degree ≤ e_max.  The dominant-degree
//      analysis is now over the extended identity:
//          deg(D · y')    = deg(D) + e − 1   (if y' ≠ 0)
//          deg(f_n · y)   = deg(f_n) + e
//          deg(g_n)
//      e_max is chosen as deg(g_n) − min(deg(D)−1, deg(f_n)) when positive;
//      we then build a square Q-linear system and read off the solution.
//      The special case e = 0 (y constant) is checked separately to avoid
//      false negatives when the D·y' term vanishes identically.
//
// Returns the polynomial y as ExprPtr, or Unimplemented if no polynomial
// solution exists / the rational structure cannot be canonicalised.
[[nodiscard]] Result<ExprPtr> solve_risch_de_rational_q(
    ExprPtr f_expr, ExprPtr g_expr, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    // F0.8-MIGRATED
    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "solve_risch_de_rational_q",
            msg,
            cas::error::reason_codes::RISCH_SINGULAR_SYSTEM,
            "Risch DE rational: check for singular Hermite system or missing partial fraction coefficients",
            "F0.8");
    };

    // ----- Stage 1: decompose f and g into num/den and build LCM-scaled
    //                polynomial coefficients (D, F_n, G_n) for the
    //                multiplied form  D · y' + F_n · y = G_n.
    auto f_parts = algebra::apart_num_den(f_expr, ctx);
    auto g_parts = algebra::apart_num_den(g_expr, ctx);
    if (f_parts.is_error()) return fail_unimpl("Risch DE rational: cannot split f into num/den");
    if (g_parts.is_error()) return fail_unimpl("Risch DE rational: cannot split g into num/den");

    ExprPtr fn = f_parts.value().numerator;
    ExprPtr fd = f_parts.value().denominator;
    ExprPtr gn = g_parts.value().numerator;
    ExprPtr gd = g_parts.value().denominator;

    // D = lcm(fd, gd) = fd · gd / gcd(fd, gd).
    auto gcd_fd_gd = algebra::polynomial_gcd(fd, gd, var, ctx);
    ExprPtr D_expr;
    if (gcd_fd_gd.is_ok()) {
        auto prod_simp = ctx.simplify(arena.make<Binary>(BinaryOp::Mul, fd, gd));
        if (prod_simp.is_error()) return prod_simp;
        auto D_div = ctx.simplify(arena.make<Binary>(BinaryOp::Div, prod_simp.value(), gcd_fd_gd.value()));
        if (D_div.is_error()) return D_div;
        D_expr = D_div.value();
    } else {
        auto prod_simp = ctx.simplify(arena.make<Binary>(BinaryOp::Mul, fd, gd));
        if (prod_simp.is_error()) return prod_simp;
        D_expr = prod_simp.value();
    }

    auto D_over_fd = ctx.simplify(arena.make<Binary>(BinaryOp::Div, D_expr, fd));
    if (D_over_fd.is_error()) return D_over_fd;
    auto D_over_gd = ctx.simplify(arena.make<Binary>(BinaryOp::Div, D_expr, gd));
    if (D_over_gd.is_error()) return D_over_gd;
    auto fn_full = ctx.simplify(arena.make<Binary>(BinaryOp::Mul, fn, D_over_fd.value()));
    if (fn_full.is_error()) return fn_full;
    auto gn_full = ctx.simplify(arena.make<Binary>(BinaryOp::Mul, gn, D_over_gd.value()));
    if (gn_full.is_error()) return gn_full;

    auto D_poly_res = algebra::parse_polynomial(D_expr, var, ctx);
    auto fn_poly_res = algebra::parse_polynomial(fn_full.value(), var, ctx);
    auto gn_poly_res = algebra::parse_polynomial(gn_full.value(), var, ctx);
    if (D_poly_res.is_error()) return fail_unimpl("Risch DE rational: D not polynomial");
    if (fn_poly_res.is_error()) return fail_unimpl("Risch DE rational: f_n not polynomial");
    if (gn_poly_res.is_error()) return fail_unimpl("Risch DE rational: g_n not polynomial");

    auto D_opt = risch_extract_rational_coeffs(D_poly_res.value());
    auto fn_opt = risch_extract_rational_coeffs(fn_poly_res.value());
    auto gn_opt = risch_extract_rational_coeffs(gn_poly_res.value());
    if (!D_opt || !fn_opt || !gn_opt) return fail_unimpl("Risch DE rational: non-rational coefficients");

    std::vector<Rational> D = std::move(*D_opt);
    std::vector<Rational> Fn = std::move(*fn_opt);
    std::vector<Rational> Gn = std::move(*gn_opt);

    auto strip_trailing = [](std::vector<Rational>& v) {
        while (!v.empty() && v.back().numerator().is_zero()) v.pop_back();
    };
    strip_trailing(D); strip_trailing(Fn); strip_trailing(Gn);

    if (D.empty()) return fail_unimpl("Risch DE rational: zero denominator");

    if (Gn.empty()) return ok(arena.make<IntegerLit>(BigInt(0)));

    // ----- Stage 2: ansatz  y = P(x) / D(x)  with P polynomial.
    // Substituting y = P/D into  D·y' + F_n·y = G_n  gives the
    // polynomial identity
    //
    //     D · P' + (F_n − D') · P  =  G_n · D                       (*)
    //
    // (See Bronstein, Symbolic Integration I, §6.1; this is the
    // "denominator over-bound" reduction: the true denominator of y
    // always divides D, hence P = y·D is polynomial whenever y is
    // rational-with-denominator-dividing-D, which covers every
    // elementary integral whose Risch DE arises from rational f, g.)
    //
    // Build the new coefficient vectors H = F_n − D'  and  Gd = G_n · D
    // and search for polynomial P of degree ≤ M_bound via a Q-linear
    // system on its coefficients.

    // D'.
    std::vector<Rational> Dprime;
    if (D.size() >= 2) {
        Dprime.resize(D.size() - 1U, Rational(BigInt(0)));
        for (std::size_t k = 1; k < D.size(); ++k) {
            Dprime[k - 1U] = Rational(BigInt(static_cast<std::int64_t>(k))) * D[k];
        }
    }

    // H = F_n − D'.
    std::vector<Rational> H(std::max(Fn.size(), Dprime.size()), Rational(BigInt(0)));
    for (std::size_t i = 0; i < Fn.size(); ++i) H[i] = H[i] + Fn[i];
    for (std::size_t i = 0; i < Dprime.size(); ++i) H[i] = H[i] - Dprime[i];
    strip_trailing(H);

    // Gd = G_n · D.
    std::vector<Rational> Gd(Gn.size() + D.size() - 1U, Rational(BigInt(0)));
    for (std::size_t i = 0; i < Gn.size(); ++i) {
        for (std::size_t j = 0; j < D.size(); ++j) {
            Gd[i + j] = Gd[i + j] + Gn[i] * D[j];
        }
    }
    strip_trailing(Gd);

    if (Gd.empty()) {
        // Trivially y = 0.
        return ok(arena.make<IntegerLit>(BigInt(0)));
    }

    const int deg_D = static_cast<int>(D.size()) - 1;
    const int deg_H = H.empty() ? -1 : static_cast<int>(H.size()) - 1;
    const int deg_Gd = static_cast<int>(Gd.size()) - 1;

    // Degree bound for P from the dominant-degree analysis of (*):
    //   deg(D · P')     = deg_D + (deg_P − 1)   for deg_P ≥ 1
    //   deg(H · P)      = deg_H + deg_P
    //   deg(G_n · D)    = deg_Gd
    // M_dom = max attainable contribution → equate to deg_Gd.
    int M_bound;
    if (deg_H >= deg_D - 1) {
        // H · P term dominates (or ties).
        M_bound = deg_Gd - deg_H;
    } else {
        // D · P' term dominates.
        M_bound = deg_Gd - deg_D + 1;
    }
    // Cancellation case (M_bound = deg_D) loses leading coefficient → allow extra slack.
    if (M_bound < deg_D) M_bound = deg_D;
    // Always allow degenerate constant solutions.
    if (M_bound < 0) M_bound = 0;
    // Safety cap — configurable via CASContext (CLAUDE.md Cat. 1: no magic constants).
    const int ansatz_cap = static_cast<int>(ctx.max_risch_rational_ansatz_degree());
    if (M_bound > ansatz_cap) return fail_unimpl("Risch DE rational: degree bound for P exceeds ctx.max_risch_rational_ansatz_degree (BUG-HANG-001)");

    const std::size_t n_unk = static_cast<std::size_t>(M_bound) + 1U;
    int max_eq_deg = std::max({deg_Gd,
        deg_D + M_bound - 1,
        (deg_H >= 0 ? deg_H + M_bound : -1)});
    if (max_eq_deg < 0) max_eq_deg = 0;
    const std::size_t n_eq = static_cast<std::size_t>(max_eq_deg) + 1U;

    std::vector<std::vector<Rational>> M(n_eq,
        std::vector<Rational>(n_unk + 1U, Rational(BigInt(0))));
    for (std::size_t j = 0; j < n_eq; ++j) {
        // RHS column.
        if (j < Gd.size()) M[j][n_unk] = Gd[j];
        // D · P' contribution: Σ_{p+k-1=j, k≥1} k · D_p · P_k.
        for (std::size_t k = 1; k <= static_cast<std::size_t>(M_bound); ++k) {
            std::size_t p_needed = j + 1U;
            if (p_needed < k) continue;
            std::size_t p = p_needed - k;
            if (p >= D.size()) continue;
            Rational add = Rational(BigInt(static_cast<std::int64_t>(k))) * D[p];
            M[j][k] = M[j][k] + add;
        }
        // H · P contribution: Σ_{p+k=j} H_p · P_k.
        for (std::size_t k = 0; k <= static_cast<std::size_t>(M_bound); ++k) {
            if (k > j) break;
            std::size_t p = j - k;
            if (p >= H.size()) continue;
            M[j][k] = M[j][k] + H[p];
        }
    }

    // Gaussian elimination over Q with partial pivoting (full RREF on pivot rows).
    // Pure-Rational arithmetic, no simplify() in the loop body → poll the
    // interrupt flag per pivot column to keep the (degree-capped) rational
    // Risch-DE solver interruptible (A2 / HC-F70-A33).
    const std::size_t cols = n_unk + 1U;
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < n_unk && pivot_row < n_eq; ++col) {
        if (auto chk = ctx.check_interrupt(); chk.is_error()) return fail<ExprPtr>(chk.error());
        std::size_t best = pivot_row;
        while (best < n_eq && M[best][col].numerator().is_zero()) ++best;
        if (best == n_eq) continue;
        if (best != pivot_row) std::swap(M[best], M[pivot_row]);
        Rational pivot = M[pivot_row][col];
        for (std::size_t c = col; c < cols; ++c) M[pivot_row][c] = M[pivot_row][c] / pivot;
        for (std::size_t r = 0; r < n_eq; ++r) {
            if (r == pivot_row) continue;
            Rational factor = M[r][col];
            if (factor.numerator().is_zero()) continue;
            for (std::size_t c = col; c < cols; ++c) {
                M[r][c] = M[r][c] - factor * M[pivot_row][c];
            }
        }
        ++pivot_row;
    }

    for (std::size_t r = pivot_row; r < n_eq; ++r) {
        if (!M[r][n_unk].numerator().is_zero()) {
            return fail_unimpl("Risch DE rational: inconsistent linear system");
        }
    }

    std::vector<Rational> P(n_unk, Rational(BigInt(0)));
    for (std::size_t r = 0; r < pivot_row; ++r) {
        std::size_t leading = n_unk;
        for (std::size_t c = 0; c < n_unk; ++c) {
            if (!M[r][c].numerator().is_zero()) { leading = c; break; }
        }
        if (leading < n_unk) P[leading] = M[r][n_unk];
    }

    // Build y = P(x) / D(x) and simplify.
    ExprPtr P_expr = risch_build_poly_expr(arena, P, var);
    // Fast path: P ≡ 0 ⇒ y = 0.
    bool all_zero = true;
    for (const auto& c : P) if (!c.numerator().is_zero()) { all_zero = false; break; }
    if (all_zero) return ok(arena.make<IntegerLit>(BigInt(0)));

    auto y_expr = ctx.simplify(arena.make<Binary>(BinaryOp::Div, P_expr, D_expr));
    if (y_expr.is_error()) return y_expr;
    auto y_together = algebra::together(y_expr.value(), ctx);
    if (y_together.is_ok()) y_expr = ctx.simplify(y_together.value());
    return y_expr;
}

}  // namespace

// Dispatch: try polynomial-fast path first; if either f or g has a
// non-trivial denominator, fall back to the rational-coefficient solver.
// Public to cas::calculus per riuso da risch_logarithmic.cpp (cap.8) e
// risch_parametric.cpp (cap.7).
Result<ExprPtr> solve_risch_de_q(
    ExprPtr f_expr, ExprPtr g_expr, const Symbol& var, symbolic::CASContext& ctx) {
    // A2: cancellation must short-circuit the dispatcher and surface as Timeout.
    // Poll at entry — this building block is also called directly by the
    // cap.7/cap.8 tower extensions, bypassing the integrate() head poll.
    if (auto chk = ctx.check_interrupt(); chk.is_error()) return fail<ExprPtr>(chk.error());
    // Simplify failures are tolerated (use original expr), EXCEPT a cancellation.
    auto f_simp = ctx.simplify(f_expr);
    if (f_simp.is_error() && f_simp.error().kind == CASErrorKind::Timeout)
        return fail<ExprPtr>(f_simp.error());
    auto g_simp = ctx.simplify(g_expr);
    if (g_simp.is_error() && g_simp.error().kind == CASErrorKind::Timeout)
        return fail<ExprPtr>(g_simp.error());
    if (f_simp.is_ok()) f_expr = f_simp.value();
    if (g_simp.is_ok()) g_expr = g_simp.value();
    auto poly_attempt = solve_risch_de_poly_q(f_expr, g_expr, var, ctx);
    if (poly_attempt.is_ok()) return poly_attempt;
    // Only fall back to the rational branch for a genuine "no polynomial
    // solution" dead-end — a cancellation must not be retried as if algebraic.
    if (poly_attempt.error().kind == CASErrorKind::Timeout) return poly_attempt;
    return solve_risch_de_rational_q(f_expr, g_expr, var, ctx);
}

} // namespace cas::calculus
