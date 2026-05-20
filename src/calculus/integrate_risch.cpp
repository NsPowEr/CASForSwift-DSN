#include "calculus_internal.hpp"
#include "integrate_engine.hpp"
#include "cas/differential_algebra.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"
#include "../algebra/polynomial_internal.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
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

// Solve  y' + f · y = g  for  y ∈ Q[var]
// (Risch differential equation, polynomial-coefficient case).
//
// Given f, g polynomials in `var` with rational coefficients, returns the
// unique polynomial solution y (if one exists), or CASErrorKind::Unimplemented
// otherwise.  This generalises the previous "f must be a constant" inline
// solver to arbitrary polynomial f via degree-bound + linear-system inversion,
// following the standard Bronstein degree analysis (Symbolic Integration I,
// §5.6).
//
// Degree analysis:
//   deg(y')      = deg(y) − 1
//   deg(f · y)   = deg(f) + deg(y)
// The right-hand side g of degree d must therefore match the dominant
// term on the left:
//   • If deg(f) ≥ 1,   deg(f · y) > deg(y'),
//                      so  d = deg(f) + deg(y)  ⇒  deg(y) = d − deg(f).
//     If d < deg(f), no polynomial solution exists for nonzero g.
//   • If deg(f) = 0 with f ≠ 0,  deg(y) = d (constant f case).
//   • If f ≡ 0,        deg(y) = d + 1 (pure quadrature).
//
// With the degree e fixed, the e+1 unknown coefficients y_0..y_e are
// determined by a square linear system over Q (one equation per power of
// `var` in the identity y' + f·y = g).  Gaussian elimination over the
// rationals returns the solution exactly.  Inconsistency or singular
// system ⇒ Unimplemented (no polynomial y).
[[nodiscard]] Result<ExprPtr> solve_risch_de_poly_q(
    ExprPtr f_expr, ExprPtr g_expr, const Symbol& var, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto fail_unimpl = [&](const char* msg) {
        return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, msg, std::nullopt});
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

    // Gaussian elimination with partial pivoting over Q.
    auto sysop_timeout = ctx.timeout_check_interval(); (void)sysop_timeout;
    const std::size_t cols = n_unk + 1U;
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < n_unk && pivot_row < n_eq; ++col) {
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
    auto fail_unimpl = [&](const char* msg) {
        return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, msg, std::nullopt});
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
    // Safety cap (config knob; conservative for now).
    if (M_bound > 256) return fail_unimpl("Risch DE rational: degree bound for P too large");

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
    const std::size_t cols = n_unk + 1U;
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < n_unk && pivot_row < n_eq; ++col) {
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

// Dispatch: try polynomial-fast path first; if either f or g has a
// non-trivial denominator, fall back to the rational-coefficient solver.
[[nodiscard]] Result<ExprPtr> solve_risch_de_q(
    ExprPtr f_expr, ExprPtr g_expr, const Symbol& var, symbolic::CASContext& ctx) {
    auto f_simp = ctx.simplify(f_expr);
    auto g_simp = ctx.simplify(g_expr);
    if (f_simp.is_ok()) f_expr = f_simp.value();
    if (g_simp.is_ok()) g_expr = g_simp.value();
    auto poly_attempt = solve_risch_de_poly_q(f_expr, g_expr, var, ctx);
    if (poly_attempt.is_ok()) return poly_attempt;
    return solve_risch_de_rational_q(f_expr, g_expr, var, ctx);
}

// Integrate the polynomial-in-t part of a single logarithmic extension
// tower:  ∫ Σ_{k=0..n} a_k(x) * t^k dx,   where t = ln(u(x)),  Dt = u'/u.
//
// Standard ansatz:  the antiderivative is again a polynomial in t,
//   B(t) = Σ_{k=0..n} b_k(x) * t^k,   with
//
//   d/dx B(t)
//     = Σ b_k'(x) * t^k + Σ k * b_k(x) * (u'/u) * t^{k-1}
//     = b_n' * t^n + Σ_{k=0..n-1} [ b_k' + (k+1) * b_{k+1} * (u'/u) ] * t^k.
//
// Matching coefficients with a_k * t^k gives the descending recursion
//   b_n = ∫ a_n dx
//   b_k = ∫ [ a_k - (k+1) * b_{k+1} * (u'/u) ] dx,   for k = n-1, n-2, ..., 0.
//
// Each integration is carried out in the lower field (here Q(x) for a single
// log extension) via the existing integrate() routine.  Failures propagate
// as Unimplemented so that the caller can fall back to other strategies.
[[nodiscard]] Result<ExprPtr> integrate_log_polynomial_part(
    const algebra::PolyExpr& quot,
    ExprPtr u_arg,
    const Symbol& t_top,
    const Symbol& var,
    symbolic::CASContext& context) {
    AstArena& arena = context.arena();

    if (quot.empty()) {
        return ok(arena.make<IntegerLit>(BigInt(0)));
    }
    const std::size_t deg = quot.size() - 1U;

    // u'/u (simplified once up front)
    auto du_res = diff(u_arg, var, 1U, context);
    if (du_res.is_error()) return fail<ExprPtr>(du_res.error());
    ExprPtr du_over_u = arena.make<Binary>(BinaryOp::Div, du_res.value(), u_arg);
    if (auto s = context.simplify(du_over_u); s.is_ok()) du_over_u = s.value();

    std::vector<ExprPtr> b(deg + 1U, ExprPtr{});

    for (std::ptrdiff_t k = static_cast<std::ptrdiff_t>(deg); k >= 0; --k) {
        const std::size_t kz = static_cast<std::size_t>(k);
        ExprPtr a_k = (kz < quot.size()) ? quot[kz] : ExprPtr{};
        if (!a_k) a_k = arena.make<IntegerLit>(BigInt(0));

        // rhs = a_k  -  (k+1) * b_{k+1} * (u'/u)
        ExprPtr rhs = a_k;
        if (kz + 1U <= deg && b[kz + 1U]) {
            ExprPtr kp1 = arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(kz + 1U)));
            ExprPtr correction = arena.make<Product>(std::vector<ExprPtr>{kp1, b[kz + 1U], du_over_u});
            rhs = arena.make<Binary>(BinaryOp::Sub, rhs, correction);
        }
        if (auto s = context.simplify(rhs); s.is_ok()) rhs = s.value();

        auto b_k_res = integrate(rhs, var, context);
        if (b_k_res.is_error()) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Risch log: lower-field integration failed at degree k=" + std::to_string(kz),
                .hint = std::nullopt});
        }
        b[kz] = b_k_res.value();
    }

    // Build B(t) = Σ b_k * t^k.
    std::vector<ExprPtr> terms;
    terms.reserve(b.size());
    ExprPtr t_sym = arena.make<Symbol>(t_top.name);
    for (std::size_t k = 0; k < b.size(); ++k) {
        if (!b[k]) continue;
        if (const auto* il = expr_cast<IntegerLit>(b[k]); il && il->value.is_zero()) continue;
        ExprPtr term;
        if (k == 0U) {
            term = b[k];
        } else if (k == 1U) {
            term = arena.make<Binary>(BinaryOp::Mul, b[k], t_sym);
        } else {
            ExprPtr t_pow = arena.make<Binary>(
                BinaryOp::Pow,
                t_sym,
                arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k))));
            term = arena.make<Binary>(BinaryOp::Mul, b[k], t_pow);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return ok(arena.make<IntegerLit>(BigInt(0)));
    if (terms.size() == 1U) return ok(terms.front());
    ExprPtr raw = arena.make<Sum>(std::move(terms));
    if (auto s = context.simplify(raw); s.is_ok()) return ok(s.value());
    return ok(raw);
}

}  // namespace

Result<ExprPtr> integrate_risch(ExprPtr expr, const Symbol& var, symbolic::CASContext& context) {
    AstArena& arena = context.arena();

    // 0. Simple pattern matching for fundamental transcendental functions
    if (const auto* call = expr_cast<FuncCall>(expr)) {
        if (call->args.size() == 1 && integrate_detail::is_same_symbol(call->args[0], var)) {
            if (call->func_id == BuiltinOp::Exp) {
                return ok(expr);
            }
            if (call->func_id == BuiltinOp::Ln) {
                ExprPtr x = arena.make<Symbol>(var);
                return ok(arena.make<Sum>(std::vector<ExprPtr>{
                    arena.make<Binary>(BinaryOp::Mul, x, expr),
                    arena.make<Unary>(UnaryOp::Neg, x)
                }));
            }
        }
    }

    // 1. Build Differential Extension Tower
    auto field_res = DifferentialField::build(expr, var, context);
    if (field_res.is_error()) return fail<ExprPtr>(field_res.error());
    const auto& field = field_res.value();

    // 2. Map expression to the differential field (generators t_1, ..., t_n)
    auto gen_expr_res = field.to_field_generators(expr, context);
    if (gen_expr_res.is_error()) return fail<ExprPtr>(gen_expr_res.error());
    ExprPtr gen_expr = gen_expr_res.value();

    // 2b. Logarithmic-derivative recognition (Risch structure theorem).
    // If integrand == c · D(g)/g for some generator g in the field tower
    // and c constant, then ∫ = c · ln(g). Handles ∫ 1/(x·ln(x)) dx = ln(ln(x))
    // and similar nested-log cases that Hermite/Rothstein-Trager cannot
    // express (resultant root is non-constant rational function of x).
    //
    // Verification by differentiation roundtrip: candidate F = c·ln(g);
    // accept iff D(F) - gen_expr simplifies to 0. This bypasses fragile
    // shape inspection on the candidate constant c.
    // Try each extension as candidate g such that integrand = c · D(g)/g
    // for some constant c. For pragmatic detection, try small rational
    // constants (c ∈ {1, -1, 1/2, 2, -1/2, -2}) and verify by subtraction.
    auto try_constant = [&](ExprPtr cF, ExprPtr DF_val) -> Result<ExprPtr> {
        ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, DF_val, expr);
        // together() reduces rational sub-expressions over a common
        // denominator before simplify — required because raw simplify
        // does not commute the (1/x)·(1/y) → 1/(xy) transform that the
        // subtraction needs to reach zero.
        auto delta_tog = algebra::together(delta, context);
        if (delta_tog.is_error()) return fail<ExprPtr>(delta_tog.error());
        auto delta_simp = context.simplify(delta_tog.value());
        if (delta_simp.is_error()) return fail<ExprPtr>(delta_simp.error());
        if (expr_is<IntegerLit>(delta_simp.value())
            && expr_ref<IntegerLit>(delta_simp.value()).value.is_zero()) {
            return context.simplify(cF);
        }
        return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, "no match", std::nullopt});
    };
    const std::array<std::pair<long long, long long>, 6> trial_consts = {{
        {1, 1}, {-1, 1}, {1, 2}, {2, 1}, {-1, 2}, {-2, 1},
    }};
    for (const auto& ext : field.extensions()) {
        ExprPtr g;
        if (ext.type == ExtensionType::Logarithmic) {
            g = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{ext.argument});
        } else if (ext.type == ExtensionType::Exponential) {
            g = arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{ext.argument});
        } else {
            continue;
        }
        ExprPtr F_unit = arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{g});
        auto DF_res = diff(F_unit, var, 1U, context);
        if (DF_res.is_error()) continue;
        ExprPtr DF = DF_res.value();
        for (const auto& [num, den] : trial_consts) {
            ExprPtr c_const = (den == 1)
                ? static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(num)))
                : static_cast<ExprPtr>(arena.make<RationalLit>(BigInt(num), BigInt(den)));
            ExprPtr cF = arena.make<Binary>(BinaryOp::Mul, c_const, F_unit);
            ExprPtr cDF = arena.make<Binary>(BinaryOp::Mul, c_const, DF);
            auto res = try_constant(cF, cDF);
            if (res.is_ok()) return res;
        }
    }

    // 3. Decompose into P/Q with respect to the topmost generator t_n
    Symbol t_top = field.extensions().empty() ? var : field.extensions().back().t_var;

    auto rational_res = algebra::apart_num_den(gen_expr, context);
    if (rational_res.is_error()) return fail<ExprPtr>(rational_res.error());
    ExprPtr P = rational_res.value().numerator;
    ExprPtr Q = rational_res.value().denominator;

    // 3b. Polynomial part
    ExprPtr poly_integral_part = arena.make<IntegerLit>(BigInt(0));
    auto P_poly = algebra::parse_polynomial(P, t_top, context);
    auto Q_poly = algebra::parse_polynomial(Q, t_top, context);

    if (P_poly.is_ok() && Q_poly.is_ok() && !algebra::is_zero_poly(Q_poly.value())) {
        if (algebra::poly_degree(P_poly.value()) >= algebra::poly_degree(Q_poly.value())) {
            auto div_res = algebra::divide_poly_with_remainder(P_poly.value(), Q_poly.value(), context);
            if (div_res.is_ok()) {
                const auto& quot = div_res.value().quotient;
                const auto& rem  = div_res.value().remainder;
                
                std::vector<ExprPtr> int_terms;
                if (field.extensions().empty()) {
                    // Base case: t_top = x, D(t_top) = 1
                    for (std::size_t k = 0; k < quot.size(); ++k) {
                        ExprPtr coeff = quot[k];
                        if (!coeff || algebra::poly_is_zero_expr(coeff)) continue;
                        if (k == 0) {
                            int_terms.push_back(arena.make<Binary>(BinaryOp::Mul, coeff, arena.make<Symbol>(t_top)));
                        } else {
                            ExprPtr kp1 = arena.make<IntegerLit>(BigInt(static_cast<long long>(k + 1)));
                            ExprPtr t_pow = arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(t_top), kp1);
                            int_terms.push_back(arena.make<Binary>(BinaryOp::Div, arena.make<Binary>(BinaryOp::Mul, coeff, t_pow), kp1));
                        }
                    }
                } else {
                    // Transcendental case
                    const auto& ext = field.extensions().back();
                    const bool handle_log_polynomial = (ext.type == ExtensionType::Logarithmic);

                    // Logarithmic extension: descending recursion on B(t),
                    // see integrate_log_polynomial_part.  Performs the full
                    // polynomial-in-t integration in one structured pass.
                    if (handle_log_polynomial) {
                        auto B_res = integrate_log_polynomial_part(quot, ext.argument, t_top, var, context);
                        if (B_res.is_error()) return fail<ExprPtr>(B_res.error());
                        int_terms.push_back(B_res.value());
                    }

                    for (std::size_t k = 0; !handle_log_polynomial && k < quot.size(); ++k) {
                        ExprPtr coeff = quot[k];
                        if (!coeff) continue;
                        // Pre-simplify coefficient: divide_poly_with_remainder
                        // can leave Product([0, …]) un-collapsed when Q has
                        // x-dependent coefficients.  Without this collapse
                        // the zero-check below misses true zeroes and the
                        // recursive integrate path injects spurious terms
                        // (e.g. arbitrary x-multiples of 1/x).
                        {
                            auto coeff_tog = algebra::together(coeff, context);
                            if (coeff_tog.is_ok()) {
                                auto coeff_simp = context.simplify(coeff_tog.value());
                                if (coeff_simp.is_ok()) coeff = coeff_simp.value();
                            }
                        }
                        if (algebra::poly_is_zero_expr(coeff)) continue;

                        if (ext.type == ExtensionType::Exponential) {
                            // ∫ a_k * t^k where t = exp(u), Dt = u't
                            if (k == 0) {
                                auto base_int = integrate(coeff, var, context);
                                if (base_int.is_ok()) int_terms.push_back(base_int.value());
                                else return base_int;
                            } else {
                                // solve Dy + k*u'*y = a_k
                                auto du_res = diff(ext.argument, var, 1U, context);
                                if (du_res.is_error()) return fail<ExprPtr>(du_res.error());
                                ExprPtr du = du_res.value();
                                
                                // Risch DE for the exponential extension at level k:
                                //   y' + k · u' · y = a_k        (a_k = `coeff`)
                                // Solved exactly over Q[var] via the polynomial-coefficient
                                // engine `solve_risch_de_poly_q` (degree-bound + linear
                                // system).  Handles every polynomial u, not just linear.
                                ExprPtr f_expr = arena.make<Binary>(BinaryOp::Mul,
                                    arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(k))), du);
                                auto y_res = solve_risch_de_q(f_expr, coeff, var, context);
                                if (y_res.is_ok()) {
                                    ExprPtr t_pow = (k == 1) ? arena.make<Symbol>(t_top) : arena.make<Binary>(BinaryOp::Pow, arena.make<Symbol>(t_top), arena.make<IntegerLit>(BigInt(k)));
                                    int_terms.push_back(arena.make<Binary>(BinaryOp::Mul, y_res.value(), t_pow));
                                } else {
                                    return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, "Risch: could not solve DE for exponential term", std::nullopt});
                                }
                            }
                        } else if (ext.type == ExtensionType::Logarithmic) {
                            if (k == 1 && integrate_detail::is_same_symbol(ext.argument, var) && integrate_detail::is_one(coeff)) {
                                int_terms.push_back(arena.make<Sum>(std::vector<ExprPtr>{
                                    arena.make<Binary>(BinaryOp::Mul, arena.make<Symbol>(var), arena.make<Symbol>(t_top)),
                                    arena.make<Unary>(UnaryOp::Neg, arena.make<Symbol>(var))
                                }));
                            } else {
                                return fail<ExprPtr>(CASError{CASErrorKind::Unimplemented, "Risch: log extension integration not fully implemented", std::nullopt});
                            }
                        }
                    }
                }
                
                if (!int_terms.empty()) {
                    ExprPtr raw = int_terms.size() == 1 ? int_terms[0] : arena.make<Sum>(std::move(int_terms));
                    auto simp = context.simplify(raw);
                    if (simp.is_ok()) poly_integral_part = simp.value();
                }
                auto rem_expr = algebra::polynomial_to_expr(rem, t_top, context);
                if (rem_expr.is_ok()) P = rem_expr.value();
            }
        }
    }

    // Short-circuit: if the polynomial division consumed the whole gen_expr
    // (rem ≡ 0 ⇒ P was zeroed at line 664), the rational/log layers have
    // nothing left to do — Hermite/Trager would still try and may bail on
    // Q whose coefficients are x-dependent rational expressions.  Skip
    // straight to the back-mapping.
    ExprPtr rational_part;
    ExprPtr rem_P;
    ExprPtr rem_Q;
    // Aggressive re-simplification of P after polynomial division.  The
    // divide_poly_with_remainder routine produces correct but unreduced
    // coefficient algebra (e.g. x²·x⁻² is left in product form).  Without
    // this normalisation the downstream Hermite/Trager pass parses non-
    // polynomial coefficients and bails.  algebra::together folds rational
    // sub-expressions into a single fraction, then simplify reduces.
    {
        auto P_tog = algebra::together(P, context);
        if (P_tog.is_ok()) {
            auto P_simp = context.simplify(P_tog.value());
            if (P_simp.is_ok()) P = P_simp.value();
        }
    }
    bool P_is_zero = expr_is<IntegerLit>(P) && expr_ref<IntegerLit>(P).value.is_zero();
    Result<ExprPtr> log_part_res = ok(static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(0))));
    if (P_is_zero) {
        rational_part = arena.make<IntegerLit>(BigInt(0));
        rem_P = arena.make<IntegerLit>(BigInt(0));
        rem_Q = Q;
    } else {
        // 4. Hermite Reduction: ∫ P/Q dt = A/B + ∫ C/D dt where D is square-free
        auto hermite_res = hermite_reduce(P, Q, t_top, field, context);
        if (hermite_res.is_error()) return fail<ExprPtr>(hermite_res.error());
        rational_part = hermite_res.value().rational_part;
        rem_P = hermite_res.value().remaining_P;
        rem_Q = hermite_res.value().remaining_Q;
        // 5. Rothstein-Trager for the remaining square-free part ∫ C/D dt
        log_part_res = integrate_rothstein_trager(rem_P, rem_Q, t_top, field, context);
    }
    
    if (log_part_res.is_ok()) {
        // combine rational + log parts and map back from field generators
        ExprPtr total_gen;
        bool rat_zero = expr_is<IntegerLit>(rational_part) && expr_ref<IntegerLit>(rational_part).value.is_zero();
        bool log_zero = expr_is<IntegerLit>(log_part_res.value()) && expr_ref<IntegerLit>(log_part_res.value()).value.is_zero();
        if (rat_zero && log_zero) {
            total_gen = arena.make<IntegerLit>(BigInt(0));
        } else if (rat_zero) {
            total_gen = log_part_res.value();
        } else if (log_zero) {
            total_gen = rational_part;
        } else {
            total_gen = arena.make<Sum>(std::vector<ExprPtr>{rational_part, log_part_res.value()});
        }

        auto back_res = field.from_field_generators(total_gen, context);
        if (back_res.is_error()) return back_res;

        // The polynomial integral part is also expressed in field generators
        // (t_k symbols).  Mapping it back through from_field_generators is
        // mandatory; otherwise high-degree polynomial Risch DE solutions like
        //   ∫ (3x^3 + x) exp(x^2) dx  =  ((3/2) x^2 − 1) exp(x^2)
        // would leak the t_0 symbol into the final result.
        auto poly_back_res = field.from_field_generators(poly_integral_part, context);
        if (poly_back_res.is_error()) return poly_back_res;
        bool poly_zero = expr_is<IntegerLit>(poly_back_res.value())
            && expr_ref<IntegerLit>(poly_back_res.value()).value.is_zero();
        if (poly_zero) return back_res;
        bool back_zero = expr_is<IntegerLit>(back_res.value())
            && expr_ref<IntegerLit>(back_res.value()).value.is_zero();
        if (back_zero) return context.simplify(poly_back_res.value());
        return context.simplify(arena.make<Sum>(std::vector<ExprPtr>{
            poly_back_res.value(), back_res.value()}));
    }

    // Fallback to simpler cases if full Risch fails
    return fail<ExprPtr>(CASError{
        .kind = CASErrorKind::Unimplemented,
        .message = "Risch algorithm: integrability could not be decided for this transcendental extension",
        .hint = std::nullopt
    });
}

} // namespace cas::calculus
