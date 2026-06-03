// F5.5 — Puiseux series leading-term extractor via Newton polygon.
//
// Given an algebraic curve f(x, y) = 0 with f ∈ K[x, y] and K a field of
// characteristic 0, the branches at x = 0 are described by Puiseux series
//   y(x) = c · x^(p/q) + (higher-order terms in x^(1/q)).
//
// The Newton polygon of f is the lower convex hull of the support
//   S(f) = { (i, j) : a_{ij} ≠ 0 }    in the (x-exponent, y-exponent) plane,
// where f(x, y) = Σ a_{ij} x^i y^j.
//
// Each edge of the lower hull with negative slope -μ = (j2 - j1) / (i2 - i1)
// (with i1 < i2, j1 > j2) supplies the Puiseux exponent μ = (j1 - j2) /
// (i2 - i1).  Substituting the ansatz y = c · x^μ into f, the dominant terms
// — those for which i + j·μ attains its minimum over the support — are
// exactly the lattice points on the edge.  Cancellation of the dominant
// power x^L (L = i1 + j1·μ = i2 + j2·μ) forces the characteristic
// polynomial
//     Φ_edge(c) = Σ_{(i, j) on edge} a_{ij} · c^j
// to vanish.  The roots of Φ_edge in an algebraic closure of K are the
// leading coefficients of the branches on this edge; their multiplicities
// in Φ_edge match the multiplicity tag returned to the caller.
//
// References
//   * R. Walker, "Algebraic Curves" (1950) §IV.2-3 — Newton polygon
//     algorithm, lower-hull characterisation, recursive branch refinement.
//   * D. Duval, "Rational Puiseux expansions" (1989) — characteristic
//     polynomial / leading-coefficient construction with algebraic
//     coefficients.

#include "cas/algebra.hpp"
#include "cas/calculus.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

struct Monomial {
    long long i;          // x-exponent (≥ 0; we work at x = 0)
    long long j;          // y-exponent (≥ 0)
    ExprPtr   coefficient;
};

[[nodiscard]] bool is_literal_zero(ExprPtr e) {
    if (!e) return true;
    if (const auto* il = expr_cast<IntegerLit>(e)) return il->value.is_zero();
    if (const auto* rl = expr_cast<RationalLit>(e)) return rl->numerator.is_zero();
    return false;
}

// Reduce f to monomial list (i, j, a_{ij}) with a_{ij} ≠ 0 by two univariate
// expansions: first f as a polynomial in y, then each y-coefficient as a
// polynomial in x.  Coefficients that survive both passes are constants in
// {x, y} but may contain other symbols (the only requirement Puiseux makes
// is that they live in K).
[[nodiscard]] Result<std::vector<Monomial>> extract_monomials(
    ExprPtr f,
    const Symbol& x,
    const Symbol& y,
    symbolic::CASContext& ctx) {
    auto expanded_res = algebra::expand(f, ctx);
    if (expanded_res.is_error()) return fail<std::vector<Monomial>>(expanded_res.error());

    auto y_coeffs_res = algebra::univariate_coefficients(expanded_res.value(), y, ctx);
    if (y_coeffs_res.is_error()) return fail<std::vector<Monomial>>(y_coeffs_res.error());
    const std::vector<ExprPtr>& y_coeffs = y_coeffs_res.value();

    std::vector<Monomial> result;
    result.reserve(y_coeffs.size() * 2U);
    for (std::size_t j = 0; j < y_coeffs.size(); ++j) {
        ExprPtr cj = y_coeffs[j];
        if (is_literal_zero(cj)) continue;
        auto x_coeffs_res = algebra::univariate_coefficients(cj, x, ctx);
        if (x_coeffs_res.is_error()) return fail<std::vector<Monomial>>(x_coeffs_res.error());
        const std::vector<ExprPtr>& x_coeffs = x_coeffs_res.value();
        for (std::size_t i = 0; i < x_coeffs.size(); ++i) {
            auto a_simp = ctx.simplify(x_coeffs[i]);
            if (a_simp.is_error()) return fail<std::vector<Monomial>>(a_simp.error());
            if (is_literal_zero(a_simp.value())) continue;
            result.push_back(Monomial{
                static_cast<long long>(i),
                static_cast<long long>(j),
                a_simp.value(),
            });
        }
    }
    return ok(std::move(result));
}

// Lower convex hull of a finite point set in 2-D.  Andrew's monotone chain
// returns the vertices of the lower hull in left-to-right order (increasing
// i).  Points are deduplicated by (i, j); ties on i keep the minimum j.
//
// We are interested in the "lower-right" part of the hull (the part touching
// the support from the side of the origin), so the standard lower-hull
// algorithm is exactly what we want.
[[nodiscard]] std::vector<std::pair<long long, long long>> lower_hull(
    const std::vector<Monomial>& monomials) {
    std::vector<std::pair<long long, long long>> pts;
    pts.reserve(monomials.size());
    for (const auto& m : monomials) pts.emplace_back(m.i, m.j);
    std::sort(pts.begin(), pts.end(), [](const auto& a, const auto& b) {
        return a.first < b.first || (a.first == b.first && a.second < b.second);
    });
    // Per x-column keep only the minimum j (the bottom of that column is
    // what defines the hull edge from below).
    std::vector<std::pair<long long, long long>> compact;
    for (const auto& p : pts) {
        if (!compact.empty() && compact.back().first == p.first) continue;
        compact.push_back(p);
    }
    // Andrew's monotone chain — lower hull only.
    std::vector<std::pair<long long, long long>> hull;
    auto cross = [](const std::pair<long long, long long>& O,
                    const std::pair<long long, long long>& A,
                    const std::pair<long long, long long>& B) -> long long {
        return (A.first - O.first) * (B.second - O.second) -
               (A.second - O.second) * (B.first - O.first);
    };
    for (const auto& p : compact) {
        while (hull.size() >= 2 && cross(hull[hull.size() - 2], hull.back(), p) <= 0) {
            hull.pop_back();
        }
        hull.push_back(p);
    }
    return hull;
}

[[nodiscard]] long long gcd_ll(long long a, long long b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        long long t = a % b;
        a = b;
        b = t;
    }
    return a;
}

// Build the characteristic polynomial Φ(c) = Σ_{(i,j) on edge} a_{ij} · c^j
// in the auxiliary variable `c`, then return it as ExprPtr ready for
// solve_polynomial.  The exponents are shifted so the polynomial has c-zero
// term at j = j_min (degree = j_max - j_min) — this is harmless because we
// only care about non-zero roots in c, but it keeps the degree minimal and
// avoids spurious c = 0 roots that would correspond to non-leading terms.
[[nodiscard]] ExprPtr characteristic_polynomial(
    const std::vector<Monomial>& edge_points,
    const Symbol& c,
    AstArena& arena) {
    long long j_min = edge_points.front().j;
    for (const auto& m : edge_points) j_min = std::min(j_min, m.j);

    std::vector<ExprPtr> terms;
    terms.reserve(edge_points.size());
    ExprPtr c_sym = arena.make<Symbol>(c);
    for (const auto& m : edge_points) {
        long long e = m.j - j_min;
        ExprPtr term;
        if (e == 0) {
            term = m.coefficient;
        } else if (e == 1) {
            term = arena.make<Binary>(BinaryOp::Mul, m.coefficient, c_sym);
        } else {
            ExprPtr power = arena.make<Binary>(BinaryOp::Pow, c_sym,
                arena.make<IntegerLit>(BigInt(e)));
            term = arena.make<Binary>(BinaryOp::Mul, m.coefficient, power);
        }
        terms.push_back(term);
    }
    if (terms.size() == 1U) return terms.front();
    return arena.make<Sum>(std::move(terms));
}

// Decide multiplicity of `root` in the characteristic polynomial by
// repeated synthetic division — equivalent to evaluating Φ, Φ', Φ'', ... at
// the root.  Returns 1 unless we can prove a higher multiplicity via the
// simplifier (so the answer is always a sound lower bound).
[[nodiscard]] unsigned int decide_multiplicity(
    ExprPtr char_poly,
    const Symbol& c_var,
    ExprPtr root,
    symbolic::CASContext& ctx) {
    unsigned int mult = 0U;
    ExprPtr current = char_poly;
    for (unsigned int k = 0; k < 32U; ++k) {
        auto sub = ctx.substitute(current, c_var, root);
        if (sub.is_error()) return std::max(mult, 1U);
        auto val = ctx.simplify(sub.value());
        if (val.is_error()) return std::max(mult, 1U);
        if (!is_literal_zero(val.value())) break;
        ++mult;
        auto next = diff(current, c_var, 1U, ctx);
        if (next.is_error()) break;
        current = next.value();
    }
    return std::max(mult, 1U);
}

}  // namespace

Result<std::vector<PuiseuxBranch>> puiseux_leading_terms(
    ExprPtr f,
    const Symbol& x,
    const Symbol& y,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();

    auto monomials_res = extract_monomials(f, x, y, ctx);
    if (monomials_res.is_error()) return fail<std::vector<PuiseuxBranch>>(monomials_res.error());
    const std::vector<Monomial>& monomials = monomials_res.value();

    if (monomials.empty()) {
        return fail<std::vector<PuiseuxBranch>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "Puiseux: f vanishes identically; no branches to compute",
        });
    }

    auto hull = lower_hull(monomials);
    if (hull.size() < 2U) {
        // A single hull vertex means f has the form a·x^i·y^j (monomial
        // times unit) or coincident column; no Puiseux branches with
        // positive exponent at x = 0.
        return fail<std::vector<PuiseuxBranch>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Puiseux: Newton polygon degenerates to a single "
                       "vertex; no branches with positive exponent at x = 0",
        });
    }

    // For coefficient lookup along the edges.
    auto coeff_at = [&](long long i, long long j) -> ExprPtr {
        for (const auto& m : monomials) {
            if (m.i == i && m.j == j) return m.coefficient;
        }
        return arena.make<IntegerLit>(BigInt(0));
    };

    Symbol c_var = ctx.make_fresh_symbol("c");

    std::vector<PuiseuxBranch> branches;

    for (std::size_t e = 0; e + 1 < hull.size(); ++e) {
        long long i1 = hull[e].first;
        long long j1 = hull[e].second;
        long long i2 = hull[e + 1].first;
        long long j2 = hull[e + 1].second;
        // Only edges going right-and-down (negative slope) carry positive
        // Puiseux exponents.  Skip flat or rising edges.
        if (j1 <= j2) continue;

        // μ = (i2 − i1) / (j1 − j2)  reduced p/q with gcd(p, q) = 1, q > 0.
        long long dj = j1 - j2;     // > 0 (drop in y-exponent across edge)
        long long di = i2 - i1;     // > 0 (rise in x-exponent across edge)
        long long g = gcd_ll(dj, di);
        long long mu_num = di / g;        // numerator of μ (reduced)
        long long mu_den = dj / g;        // denominator of μ (reduced)
        long long step_i = di / g;        // lattice step on the edge in i
        long long step_j = dj / g;        // lattice step on the edge in −j

        // Collect every monomial whose (i, j) sits on the affine line
        //   (i − i1) · dj + (j − j1) · di = 0,  i ∈ [i1, i2].
        std::vector<Monomial> edge_points;
        for (long long s = 0; s <= g; ++s) {
            long long i_s = i1 + s * step_i;
            long long j_s = j1 - s * step_j;
            ExprPtr a = coeff_at(i_s, j_s);
            if (is_literal_zero(a)) continue;
            edge_points.push_back(Monomial{i_s, j_s, a});
        }
        if (edge_points.size() < 2U) continue;

        Rational exponent(BigInt(static_cast<std::int64_t>(mu_num)),
                          BigInt(static_cast<std::int64_t>(mu_den)));
        ExprPtr char_poly = characteristic_polynomial(edge_points, c_var, arena);

        auto roots_res = algebra::solve_polynomial(char_poly, c_var, ctx);
        if (roots_res.is_error()) {
            // If the characteristic polynomial cannot be solved over Q (or
            // RootOf-extended), surface a diagnostic that names the edge —
            // never a silent omission.
            return fail<std::vector<PuiseuxBranch>>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Puiseux: characteristic polynomial on Newton "
                           "polygon edge cannot be solved by the polynomial "
                           "solver; closure requires extending the "
                           "coefficient field",
            });
        }
        for (const ExprPtr& root : roots_res.value()) {
            if (is_literal_zero(root)) continue;  // spurious c = 0
            unsigned int mult = decide_multiplicity(char_poly, c_var, root, ctx);
            branches.push_back(PuiseuxBranch{
                exponent,
                root,
                mult,
            });
        }
    }

    if (branches.empty()) {
        return fail<std::vector<PuiseuxBranch>>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Puiseux: no negative-slope edge with solvable "
                       "characteristic polynomial; the curve may be smooth at "
                       "the origin or require an extended coefficient field",
        });
    }

    return ok(std::move(branches));
}

}  // namespace cas::calculus
