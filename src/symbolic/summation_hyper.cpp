// F5.7 — Petkovšek's Hyper algorithm (rational characteristic root z ∈ ℚ).
//
// Solves  Σ_{i=0}^J p[i](n)·y(n+i) = 0  for a hypergeometric term ratio
// ρ(n) = y(n+1)/y(n) ∈ ℚ(n), following Petkovšek 1992 ("A=B" ch.8):
//
//   1. a ranges over monic divisors of p[0](n);
//      b ranges over monic divisors of p[J](n−J+1).
//   2. For each (a,b) build  P_i(n) = p[i](n)·∏_{j<i}a(n+j)·∏_{j≥i}b(n+j);
//      the characteristic polynomial χ(z) = Σ_{deg P_i = δ} lc(P_i)·z^i gives
//      the candidate constants z (we keep rational roots).
//   3. For each z, the transformed recurrence Σ_i z^i P_i(n)·c(n+i)=0 is linear
//      in the unknown polynomial c(n); solve it (degree ≤ ctx bound) and form
//      ρ(n) = z·(a(n)/b(n))·(c(n+1)/c(n)).
//
// Every returned ρ is verified to annihilate the recurrence (soundness gate),
// so spurious (a,b,z) candidates can never produce a wrong closed form.

#include "summation_hyper.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas::symbolic {
namespace {

ExprPtr int_lit(AstArena& a, long long v) { return a.make<IntegerLit>(BigInt(v)); }

ExprPtr rat_lit(AstArena& a, const Rational& r) {
    if (r.denominator() == BigInt(1)) return a.make<IntegerLit>(r.numerator());
    return a.make<RationalLit>(r.numerator(), r.denominator());
}

ExprPtr expand_simplify(ExprPtr e, symbolic::CASContext& ctx) {
    auto ex = algebra::expand(e, ctx);
    ExprPtr r = ex.is_ok() ? ex.value() : e;
    auto s = ctx.simplify(r);
    return s.is_ok() ? s.value() : r;
}

// n → n + h
ExprPtr shift_poly(ExprPtr poly, const Symbol& n, long long h, symbolic::CASContext& ctx) {
    if (h == 0) return poly;
    AstArena& a = ctx.arena();
    ExprPtr arg = a.make<Binary>(BinaryOp::Add, a.make<Symbol>(n), int_lit(a, h));
    auto s = ctx.substitute(poly, n, arg);
    return s.is_ok() ? expand_simplify(s.value(), ctx) : poly;
}

// Coefficients of poly as a univariate polynomial in n, as exact rationals.
// Returns nullopt if poly is not a rational-coefficient polynomial in n.
std::optional<std::vector<Rational>> rational_coeffs(
    ExprPtr poly, const Symbol& n, symbolic::CASContext& ctx) {
    auto cs = algebra::univariate_coefficients(poly, n, ctx);
    if (cs.is_error()) return std::nullopt;
    std::vector<Rational> out;
    for (ExprPtr c : cs.value()) {
        auto sc = ctx.simplify(c);
        ExprPtr cc = sc.is_ok() ? sc.value() : c;
        if (const auto* il = expr_cast<IntegerLit>(cc)) { out.emplace_back(il->value); continue; }
        if (const auto* rl = expr_cast<RationalLit>(cc)) { out.emplace_back(rl->numerator, rl->denominator); continue; }
        if (const auto* un = expr_cast<Unary>(cc); un && un->op == UnaryOp::Neg) {
            if (const auto* il2 = expr_cast<IntegerLit>(un->operand)) { out.emplace_back(-il2->value); continue; }
            if (const auto* rl2 = expr_cast<RationalLit>(un->operand)) { out.emplace_back(-rl2->numerator, rl2->denominator); continue; }
        }
        return std::nullopt;
    }
    while (out.size() > 1U && out.back().numerator().is_zero()) out.pop_back();
    return out;
}

// Make poly monic (leading coefficient 1) over ℚ[n].
ExprPtr make_monic(ExprPtr poly, const Symbol& n, symbolic::CASContext& ctx) {
    auto cv = rational_coeffs(poly, n, ctx);
    if (!cv || cv->empty()) return poly;
    Rational lc = cv->back();
    if (lc.numerator().is_zero() || lc == Rational(1)) return poly;
    AstArena& a = ctx.arena();
    Rational inv = Rational(1) / lc;
    return expand_simplify(a.make<Binary>(BinaryOp::Mul, rat_lit(a, inv), poly), ctx);
}

// All monic divisors of poly over ℚ[n] (1 included).  nullopt if the divisor
// count would exceed the configured cap.
std::optional<std::vector<ExprPtr>> monic_divisors(
    ExprPtr poly, const Symbol& n, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    auto fac = algebra::factor_polynomial(poly, n, ctx);
    if (fac.is_error()) return std::vector<ExprPtr>{int_lit(a, 1)};

    // Keep only the non-constant irreducible factors, made monic.
    std::vector<ExprPtr> irr;
    std::vector<unsigned int> mult;
    for (const auto& f : fac.value().factors) {
        auto dg = algebra::polynomial_degree(f.factor, n, ctx);
        if (dg.is_error() || dg.value() == 0U) continue;
        irr.push_back(make_monic(f.factor, n, ctx));
        mult.push_back(f.multiplicity);
    }

    std::size_t count = 1U;
    for (unsigned int m : mult) count *= static_cast<std::size_t>(m) + 1U;
    if (count > ctx.max_hyper_divisors()) return std::nullopt;

    std::vector<ExprPtr> divisors{int_lit(a, 1)};
    for (std::size_t i = 0U; i < irr.size(); ++i) {
        std::vector<ExprPtr> next;
        for (ExprPtr d : divisors) {
            ExprPtr acc = d;
            for (unsigned int e = 0U; e <= mult[i]; ++e) {
                next.push_back(acc);
                if (e < mult[i]) acc = expand_simplify(a.make<Binary>(BinaryOp::Mul, acc, irr[i]), ctx);
            }
        }
        divisors = std::move(next);
    }
    return divisors;
}

// ∏_{j=lo}^{hi-1} base(n+j)   (empty product = 1)
ExprPtr shifted_product(ExprPtr base, const Symbol& n, long long lo, long long hi,
                        symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    ExprPtr acc = int_lit(a, 1);
    for (long long j = lo; j < hi; ++j)
        acc = a.make<Binary>(BinaryOp::Mul, acc, shift_poly(base, n, j, ctx));
    return expand_simplify(acc, ctx);
}

bool is_rational_literal(ExprPtr e, Rational& out) {
    if (const auto* il = expr_cast<IntegerLit>(e)) { out = Rational(il->value); return true; }
    if (const auto* rl = expr_cast<RationalLit>(e)) { out = Rational(rl->numerator, rl->denominator); return true; }
    if (const auto* un = expr_cast<Unary>(e); un && un->op == UnaryOp::Neg) {
        if (const auto* il2 = expr_cast<IntegerLit>(un->operand)) { out = Rational(-il2->value); return true; }
        if (const auto* rl2 = expr_cast<RationalLit>(un->operand)) { out = Rational(-rl2->numerator, rl2->denominator); return true; }
    }
    return false;
}

bool is_zero_poly(ExprPtr e, const Symbol& n, symbolic::CASContext& ctx) {
    auto cv = rational_coeffs(expand_simplify(e, ctx), n, ctx);
    if (!cv) return false;
    for (const auto& c : *cv) if (!c.numerator().is_zero()) return false;
    return true;
}

// ρ(n) verified: Σ_{i=0}^J p_i(n)·∏_{j<i} ρ(n+j) ≡ 0 ?
bool verify_ratio(const std::vector<ExprPtr>& p, ExprPtr rho,
                  const Symbol& n, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    ExprPtr sum = int_lit(a, 0);
    ExprPtr prod = int_lit(a, 1);  // ∏_{j<i} ρ(n+j); i=0 → 1
    for (std::size_t i = 0U; i < p.size(); ++i) {
        if (i > 0U) {
            ExprPtr rho_shift = shift_poly(rho, n, static_cast<long long>(i) - 1, ctx);
            prod = a.make<Binary>(BinaryOp::Mul, prod, rho_shift);
        }
        sum = a.make<Binary>(BinaryOp::Add, sum, a.make<Binary>(BinaryOp::Mul, p[i], prod));
    }
    auto tog = algebra::together(sum, ctx);
    ExprPtr combined = tog.is_ok() ? tog.value() : sum;
    auto parts = algebra::apart_num_den(combined, ctx);
    ExprPtr num = parts.is_ok() ? parts.value().numerator : combined;
    return is_zero_poly(num, n, ctx);
}

// Solve  Σ_i Q[i](n)·c(n+i) = 0  for a nonzero polynomial c of degree ≤ d_max.
// Returns the lowest-degree solution found, else nullopt.
std::optional<ExprPtr> solve_poly_recurrence(
    const std::vector<ExprPtr>& Q, const Symbol& n,
    unsigned int d_max, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();

    for (unsigned int d = 0U; d <= d_max; ++d) {
        // c(n) = Σ_{t=0}^d g_t n^t  with fresh unknowns g_t.
        std::vector<ExprPtr> g;
        ExprPtr c = int_lit(a, 0);
        for (unsigned int t = 0U; t <= d; ++t) {
            ExprPtr gt = a.make<Symbol>(ctx.make_fresh_symbol("hc"));
            g.push_back(gt);
            ExprPtr mon = (t == 0U) ? gt
                : a.make<Binary>(BinaryOp::Mul, gt,
                    a.make<Binary>(BinaryOp::Pow, a.make<Symbol>(n), int_lit(a, t)));
            c = a.make<Binary>(BinaryOp::Add, c, mon);
        }

        // E(n) = Σ_i Q[i]·c(n+i)
        ExprPtr E = int_lit(a, 0);
        for (std::size_t i = 0U; i < Q.size(); ++i) {
            ExprPtr c_shift = (i == 0U) ? c
                : [&] {
                    ExprPtr arg = a.make<Binary>(BinaryOp::Add, a.make<Symbol>(n),
                        int_lit(a, static_cast<long long>(i)));
                    auto s = ctx.substitute(c, n, arg);
                    return s.is_ok() ? s.value() : c;
                  }();
            E = a.make<Binary>(BinaryOp::Add, E, a.make<Binary>(BinaryOp::Mul, Q[i], c_shift));
        }

        auto Eexp = algebra::expand(E, ctx);
        ExprPtr Es = Eexp.is_ok() ? Eexp.value() : E;
        auto eqs_res = algebra::univariate_coefficients(Es, n, ctx);
        if (eqs_res.is_error()) continue;
        std::vector<ExprPtr> eqs = eqs_res.value();
        if (eqs.empty()) continue;

        // Homogeneous linear system in g_t.  Probe g_t = 1 (highest t first) and
        // solve the rest; the verification gate rejects bad guesses.
        for (std::size_t probe = g.size(); probe-- > 0U;) {
            ExprPtr one = int_lit(a, 1);
            std::vector<ExprPtr> sub_eqs;
            for (ExprPtr eq : eqs) {
                auto se = ctx.substitute(eq, *expr_cast<Symbol>(g[probe]), one);
                sub_eqs.push_back(se.is_ok() ? se.value() : eq);
            }
            std::vector<ExprPtr> sub_vars;
            for (std::size_t i = 0U; i < g.size(); ++i) if (i != probe) sub_vars.push_back(g[i]);

            ExprPtr c_built = c;
            bool solved = true;
            if (!sub_vars.empty()) {
                ExprPtr eqs_mat = a.make<Matrix>(sub_eqs.size(), 1U, sub_eqs);
                ExprPtr vars_mat = a.make<Matrix>(sub_vars.size(), 1U, sub_vars);
                auto sol = algebra::csolve(eqs_mat, vars_mat, ctx);
                if (sol.is_error()) { solved = false; }
                else {
                    const auto* m = expr_cast<Matrix>(sol.value());
                    if (!m || m->elements.size() != sub_vars.size()) { solved = false; }
                    else {
                        std::size_t k = 0U;
                        for (std::size_t i = 0U; i < g.size(); ++i) {
                            ExprPtr val = (i == probe) ? one : m->elements[k++];
                            auto su = ctx.substitute(c_built, *expr_cast<Symbol>(g[i]), val);
                            if (su.is_ok()) c_built = su.value();
                        }
                    }
                }
            } else {
                // single unknown: c = g_0, probe sets it to 1 ⇒ all eqs must vanish.
                for (ExprPtr eq : sub_eqs) if (!is_zero_poly(eq, n, ctx)) { solved = false; break; }
                auto su = ctx.substitute(c_built, *expr_cast<Symbol>(g[probe]), one);
                if (su.is_ok()) c_built = su.value();
            }
            if (!solved) continue;

            c_built = expand_simplify(c_built, ctx);
            if (is_zero_poly(c_built, n, ctx)) continue;
            return c_built;
        }
    }
    return std::nullopt;
}

}  // namespace

Result<HyperResult> hyper_term_ratio(
    const std::vector<ExprPtr>& p, const Symbol& n, symbolic::CASContext& ctx) {
    AstArena& a = ctx.arena();
    HyperResult result;

    if (p.size() < 2U) return ok(result);
    const std::size_t J = p.size() - 1U;
    if (is_zero_poly(p.front(), n, ctx) || is_zero_poly(p.back(), n, ctx)) {
        return fail<HyperResult>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "hyper_term_ratio: p[0] and p[J] must be nonzero."});
    }

    auto A = monic_divisors(p.front(), n, ctx);
    ExprPtr pJ_shift = shift_poly(p.back(), n, -static_cast<long long>(J) + 1, ctx);
    auto B = monic_divisors(pJ_shift, n, ctx);
    if (!A || !B) {
        return fail<HyperResult>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "hyper_term_ratio: monic-divisor search exceeds "
                       "ctx.max_hyper_divisors(); raise the limit."});
    }

    const unsigned int d_max = ctx.max_zeilberger_poly_degree();
    Symbol z_sym(ctx.make_fresh_symbol("hz"));

    for (ExprPtr aa : *A) {
        for (ExprPtr bb : *B) {
            // P_i(n) = p[i]·∏_{j<i} a(n+j)·∏_{j≥i} b(n+j)
            std::vector<ExprPtr> P(J + 1U);
            std::vector<std::optional<std::vector<Rational>>> Pcoeffs(J + 1U);
            std::size_t delta = 0U;
            for (std::size_t i = 0U; i <= J; ++i) {
                ExprPtr prod_a = shifted_product(aa, n, 0, static_cast<long long>(i), ctx);
                ExprPtr prod_b = shifted_product(bb, n, static_cast<long long>(i), static_cast<long long>(J), ctx);
                ExprPtr Pi = expand_simplify(
                    a.make<Binary>(BinaryOp::Mul, p[i],
                        a.make<Binary>(BinaryOp::Mul, prod_a, prod_b)), ctx);
                P[i] = Pi;
                Pcoeffs[i] = rational_coeffs(Pi, n, ctx);
                if (Pcoeffs[i] && !Pcoeffs[i]->empty())
                    delta = std::max(delta, Pcoeffs[i]->size() - 1U);
            }

            // Characteristic polynomial χ(z) = Σ_{deg P_i = δ} lc(P_i)·z^i.
            ExprPtr chi = int_lit(a, 0);
            bool chi_ok = true;
            for (std::size_t i = 0U; i <= J; ++i) {
                if (!Pcoeffs[i]) { chi_ok = false; break; }
                if (Pcoeffs[i]->size() != delta + 1U) continue;  // degree < δ
                Rational lc = Pcoeffs[i]->back();
                ExprPtr term = (i == 0U) ? rat_lit(a, lc)
                    : a.make<Binary>(BinaryOp::Mul, rat_lit(a, lc),
                        a.make<Binary>(BinaryOp::Pow, a.make<Symbol>(z_sym), int_lit(a, static_cast<long long>(i))));
                chi = a.make<Binary>(BinaryOp::Add, chi, term);
            }
            if (!chi_ok) continue;
            chi = expand_simplify(chi, ctx);

            auto roots = algebra::solve_polynomial(chi, z_sym, ctx);
            if (roots.is_error()) continue;

            for (ExprPtr root : roots.value()) {
                Rational z;
                if (!is_rational_literal(expand_simplify(root, ctx), z)) { result.needs_algebraic = true; continue; }
                if (z.numerator().is_zero()) continue;

                // Q_i(n) = z^i · P_i(n)
                std::vector<ExprPtr> Q(J + 1U);
                Rational zi(1);
                for (std::size_t i = 0U; i <= J; ++i) {
                    Q[i] = expand_simplify(a.make<Binary>(BinaryOp::Mul, rat_lit(a, zi), P[i]), ctx);
                    zi = zi * z;
                }

                auto c_opt = solve_poly_recurrence(Q, n, d_max, ctx);
                if (!c_opt) continue;
                ExprPtr c = *c_opt;

                // ρ(n) = z · a(n)/b(n) · c(n+1)/c(n)
                ExprPtr c1 = shift_poly(c, n, 1, ctx);
                ExprPtr rho = a.make<Binary>(BinaryOp::Mul, rat_lit(a, z),
                    a.make<Binary>(BinaryOp::Mul,
                        a.make<Binary>(BinaryOp::Div, aa, bb),
                        a.make<Binary>(BinaryOp::Div, c1, c)));
                auto tog = algebra::together(rho, ctx);
                ExprPtr rho_s = tog.is_ok() ? tog.value() : rho;
                auto rho_simp = ctx.simplify(rho_s);
                if (rho_simp.is_ok()) rho_s = rho_simp.value();

                if (verify_ratio(p, rho_s, n, ctx)) {
                    result.ratio = rho_s;
                    return ok(result);
                }
            }
        }
    }

    return ok(result);
}

}  // namespace cas::symbolic
