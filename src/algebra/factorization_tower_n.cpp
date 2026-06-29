// F3.5: factor_polynomial_tower_n — single-extension Trager factorisation
// over Q(α₁,…,α_n) via F3.4 primitive element collapse + ring-GCD in Q(θ)[x].
// Refs: Trager 1976; Cohen §3.6.4; GCL §8.7.

#include "cas/algebra.hpp"
#include "cas/algebraic_number.hpp"
#include "cas/algebraic_number_bridge.hpp"
#include "cas/algebraic_tower.hpp"
#include "cas/algebraic_tower_bridge.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include "polynomial_internal.hpp"
#include "factorization_tower_internal.hpp"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {

namespace {

namespace fti = factorization_tower_internal;

CASError nfact_error(CASErrorKind kind, std::string msg) {
    return CASError{.kind = kind, .message = std::move(msg), .hint = std::nullopt};
}
bool is_fatal_inner(const CASError& err) {
    return err.kind == CASErrorKind::Timeout
        || err.kind == CASErrorKind::InternalError;
}
ExprPtr rational_lit(AstArena& arena, const Rational& r) {
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}
AlgebraicNumber an_zero(const AlgebraicNumber::CoeffVec& q) { return AlgebraicNumber({}, q); }
AlgebraicNumber an_one(const AlgebraicNumber::CoeffVec& q) { return AlgebraicNumber({Rational(BigInt(1))}, q); }
AlgebraicNumber an_from_rational(const Rational& r, const AlgebraicNumber::CoeffVec& q) {
    if (r.numerator().is_zero()) return an_zero(q);
    return AlgebraicNumber({r}, q);
}

// Compute N(x_j) = Res_y(q_θ(y), f(x_j − s·y)) — single evaluation point.
[[nodiscard]] Result<Rational> norm_at_point(
    const RatPoly& q_theta,
    const RatPoly& f,
    const Rational& s_rat,
    const Rational& x_j) {
    // Build Q_j(y) = f(x_j − s_rat · y) as RatPoly in y.
    //   f(t) = Σ_i f_i · t^i
    //   t = x_j − s·y  →  t^i = Σ_{k} C(i,k) · x_j^{i-k} · (−s)^k · y^k
    const std::size_t df = f.size();
    if (df == 0U) {
        return ok(Rational(BigInt(0)));
    }
    std::vector<Rational> Q_j(df, Rational(BigInt(0)));
    // Precompute powers of x_j and -s.
    std::vector<Rational> x_pows(df, Rational(BigInt(0)));
    std::vector<Rational> neg_s_pows(df, Rational(BigInt(0)));
    x_pows[0] = Rational(BigInt(1));
    neg_s_pows[0] = Rational(BigInt(1));
    const Rational neg_s = Rational(BigInt(0)) - s_rat;
    for (std::size_t p = 1U; p < df; ++p) {
        x_pows[p] = x_pows[p - 1U] * x_j;
        neg_s_pows[p] = neg_s_pows[p - 1U] * neg_s;
    }
    // Binomial table.
    std::vector<std::vector<BigInt>> binom(df, std::vector<BigInt>(df, BigInt(0)));
    for (std::size_t r = 0U; r < df; ++r) {
        binom[r][0] = BigInt(1);
        for (std::size_t c = 1U; c <= r; ++c) {
            binom[r][c] = binom[r - 1U][c - 1U] + binom[r - 1U][c];
        }
    }
    for (std::size_t i = 0U; i < df; ++i) {
        const Rational& fi = f[i];
        if (fi.numerator().is_zero()) continue;
        for (std::size_t k = 0U; k <= i; ++k) {
            const Rational contrib = fi * Rational(binom[i][k])
                                    * x_pows[i - k] * neg_s_pows[k];
            Q_j[k] = Q_j[k] + contrib;
        }
    }
    RatPoly Qj(Q_j);
    Qj.normalize([](const Rational& r) { return r.numerator().is_zero(); });

    auto res = resultant_generic<Rational>(
        q_theta.coefficients(), Qj.coefficients());
    if (res.is_error()) return fail<Rational>(res.error());
    return ok(res.value());
}

// Newton interpolate values into monomial coefficients (ascending).
[[nodiscard]] Result<std::vector<Rational>> newton_interpolate(
    const std::vector<Rational>& pts,
    const std::vector<Rational>& vals) {
    const std::size_t n = pts.size();
    if (n != vals.size() || n == 0U) {
        return fail<std::vector<Rational>>(nfact_error(
            CASErrorKind::InternalError,
            "newton_interpolate: mismatched point/value vectors"));
    }
    std::vector<std::vector<Rational>> dd(n, std::vector<Rational>(n, Rational(BigInt(0))));
    for (std::size_t i = 0U; i < n; ++i) dd[i][0] = vals[i];
    for (std::size_t j = 1U; j < n; ++j) {
        for (std::size_t i = j; i < n; ++i) {
            const Rational d = pts[i] - pts[i - j];
            if (d.numerator().is_zero()) {
                return fail<std::vector<Rational>>(nfact_error(
                    CASErrorKind::InternalError,
                    "newton_interpolate: duplicate evaluation points"));
            }
            const Rational d_inv{d.denominator(), d.numerator()};
            dd[i][j] = (dd[i][j - 1U] - dd[i - 1U][j - 1U]) * d_inv;
        }
    }
    std::vector<Rational> coeffs(n, Rational(BigInt(0)));
    std::vector<Rational> basis = {Rational(BigInt(1))};
    for (std::size_t j = 0U; j < n; ++j) {
        const Rational cj = dd[j][j];
        if (!cj.numerator().is_zero()) {
            for (std::size_t d = 0U; d < basis.size(); ++d) {
                coeffs[d] = coeffs[d] + cj * basis[d];
            }
        }
        if (j + 1U < n) {
            std::vector<Rational> nb(basis.size() + 1U, Rational(BigInt(0)));
            for (std::size_t d = 0U; d < basis.size(); ++d) {
                nb[d + 1U] = nb[d + 1U] + basis[d];
                nb[d] = nb[d] - pts[j] * basis[d];
            }
            while (!nb.empty() && nb.back().numerator().is_zero()) nb.pop_back();
            basis = std::move(nb);
        }
    }
    return ok(std::move(coeffs));
}

// Compute the full Trager norm  N(x) = Res_y(q_θ(y), f(x − s·y)) ∈ Q[x].
[[nodiscard]] Result<RatPoly> compute_trager_norm(
    const RatPoly& q_theta,
    const RatPoly& f,
    const Rational& s_rat) {
    const std::size_t deg_N = q_theta.degree() * f.degree();
    const std::size_t npts = deg_N + 1U;
    std::vector<Rational> pts;
    std::vector<Rational> vals;
    pts.reserve(npts);
    vals.reserve(npts);
    for (std::size_t j = 0U; j < npts; ++j) {
        pts.push_back(Rational(BigInt(static_cast<std::int64_t>(j))));
    }
    for (const Rational& xj : pts) {
        auto v = norm_at_point(q_theta, f, s_rat, xj);
        if (v.is_error()) return fail<RatPoly>(v.error());
        vals.push_back(v.value());
    }
    auto cf = newton_interpolate(pts, vals);
    if (cf.is_error()) return fail<RatPoly>(cf.error());
    RatPoly N(std::move(cf.value()));
    N.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return ok(std::move(N));
}

[[nodiscard]] bool ratpoly_squarefree(const RatPoly& f) {
    if (f.degree() == 0U) return true;
    // f' :
    RatPoly df;
    df.reserve(f.size());
    for (std::size_t i = 1U; i < f.size(); ++i) {
        df.push_back(f[i] * Rational(BigInt(static_cast<std::int64_t>(i))));
    }
    df.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    if (df.is_zero()) return false;
    auto [g, sg, tg] = extended_gcd_rational_poly(f, df);
    return g.degree() == 0U;
}

// Build the AlgebraicNumber-coefficient lift of f(x) (constant coefficients
// in Q viewed as elements of Q(θ)).
[[nodiscard]] std::vector<AlgebraicNumber> lift_f_to_qtheta(
    const RatPoly& f,
    const AlgebraicNumber::CoeffVec& q_theta) {
    std::vector<AlgebraicNumber> out;
    out.reserve(f.size());
    for (const Rational& c : f.coefficients()) {
        out.push_back(an_from_rational(c, q_theta));
    }
    return out;
}

// g_i(x + s·θ) as poly in x over Q(θ): coeff of x^j is Σ_{k≥j} a_k·C(k,j)·(sθ)^{k-j}.
[[nodiscard]] std::vector<AlgebraicNumber> shift_qfactor_in_qtheta(
    const RatPoly& g_i,
    const Rational& s_rat,
    const AlgebraicNumber::CoeffVec& q_theta) {
    if (g_i.is_zero()) return {};
    const std::size_t n = g_i.degree();
    // θ as AlgebraicNumber = [0,1].
    AlgebraicNumber theta({Rational(BigInt(0)), Rational(BigInt(1))}, q_theta);
    AlgebraicNumber s_theta = an_from_rational(s_rat, q_theta) * theta;
    // Powers (s·θ)^p for p = 0..n.
    std::vector<AlgebraicNumber> sp;
    sp.reserve(n + 1U);
    sp.push_back(an_one(q_theta));
    for (std::size_t p = 1U; p <= n; ++p) sp.push_back(sp.back() * s_theta);
    // Binomials.
    std::vector<std::vector<BigInt>> binom(n + 1U,
        std::vector<BigInt>(n + 1U, BigInt(0)));
    for (std::size_t r = 0U; r <= n; ++r) {
        binom[r][0] = BigInt(1);
        for (std::size_t c = 1U; c <= r; ++c) {
            binom[r][c] = binom[r - 1U][c - 1U] + binom[r - 1U][c];
        }
    }
    std::vector<AlgebraicNumber> out(n + 1U, an_zero(q_theta));
    for (std::size_t j = 0U; j <= n; ++j) {
        AlgebraicNumber acc = an_zero(q_theta);
        for (std::size_t k = j; k <= n; ++k) {
            const Rational& ak = g_i[k];
            if (ak.numerator().is_zero()) continue;
            const Rational scaled = ak * Rational(binom[k][j]);
            if (scaled.numerator().is_zero()) continue;
            AlgebraicNumber term = an_from_rational(scaled, q_theta) * sp[k - j];
            acc = acc + term;
        }
        out[j] = std::move(acc);
    }
    return out;
}

// Render an AlgebraicNumber (= polynomial of degree < D in y modulo q_θ)
// as an ExprPtr by substituting y → theta_expr.
[[nodiscard]] Result<ExprPtr> render_qtheta_elem(
    const AlgebraicNumber& a,
    ExprPtr theta_expr,
    symbolic::CASContext& ctx) {
    ExprPtr raw = algebraic_number_to_expr_raw(a, theta_expr, ctx.arena());
    auto simp = ctx.simplify(raw);
    if (simp.is_error()) return fail<ExprPtr>(simp.error());
    return ok(simp.value());
}

// Convert a vector of AlgebraicNumber polynomial coefficients (ascending) into
// an ExprPtr polynomial in `var`, with each coefficient rendered via theta_expr.
[[nodiscard]] Result<ExprPtr> render_qtheta_poly(
    const std::vector<AlgebraicNumber>& coeffs,
    const Symbol& var,
    ExprPtr theta_expr,
    symbolic::CASContext& ctx) {
    PolyExpr pe;
    pe.reserve(coeffs.size());
    for (const AlgebraicNumber& c : coeffs) {
        if (c.is_zero()) { pe.push_back(ExprPtr{}); continue; }
        auto e = render_qtheta_elem(c, theta_expr, ctx);
        if (e.is_error()) return fail<ExprPtr>(e.error());
        pe.push_back(e.value());
    }
    normalize_poly(pe);
    auto pexpr = polynomial_to_expr(pe, var, ctx);
    if (pexpr.is_error()) return fail<ExprPtr>(pexpr.error());
    return ctx.simplify(pexpr.value());
}

// Make poly coefficients monic by dividing all by leading.
[[nodiscard]] Result<std::vector<AlgebraicNumber>> monic_qtheta(
    std::vector<AlgebraicNumber> poly) {
    tower_detail::strip_trailing(poly);
    if (poly.empty()) return ok(std::move(poly));
    auto inv = poly.back().inverse();
    if (inv.is_error()) return fail<std::vector<AlgebraicNumber>>(inv.error());
    for (auto& c : poly) c = c * inv.value();
    tower_detail::strip_trailing(poly);
    return ok(std::move(poly));
}

[[nodiscard]] std::size_t qtheta_poly_degree(
    const std::vector<AlgebraicNumber>& p) {
    std::size_t i = p.size();
    while (i > 0U) {
        if (!p[i - 1U].is_zero()) return i - 1U;
        --i;
    }
    return 0U;
}

}  // namespace

// ── public entry point ───────────────────────────────────────────────────

Result<Factorization> factor_polynomial_tower_n(
    ExprPtr poly,
    const Symbol& var,
    const TowerGeneratorsN& gens,
    symbolic::CASContext& ctx) {
    if (!poly) {
        return fail<Factorization>(nfact_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower_n: null polynomial"));
    }
    if (gens.alphas.empty() || gens.min_polys.empty()
        || gens.alphas.size() != gens.min_polys.size()) {
        return fail<Factorization>(nfact_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower_n: tower generators malformed"));
    }
    for (const auto& mp : gens.min_polys) {
        if (mp.size() < 2U) {
            return fail<Factorization>(nfact_error(
                CASErrorKind::InvalidArgument,
                "factor_polynomial_tower_n: generator with degree < 1"));
        }
    }

    // Parse f into RatPoly.
    auto parsed = parse_polynomial(poly, var, ctx);
    if (parsed.is_error()) return fail<Factorization>(parsed.error());
    auto f_rat_res = poly_to_rational_poly(parsed.value());
    if (f_rat_res.is_error()) {
        return fail<Factorization>(nfact_error(
            CASErrorKind::Unimplemented,
            "factor_polynomial_tower_n: input polynomial must lie in Q[var]"));
    }
    RatPoly f_poly = std::move(f_rat_res.value());
    normalize_rational_coefficients(f_poly);
    if (f_poly.is_zero()) {
        return fail<Factorization>(nfact_error(
            CASErrorKind::InvalidArgument,
            "factor_polynomial_tower_n: zero polynomial"));
    }
    const std::size_t deg_f = f_poly.degree();
    if (deg_f == 0U) {
        Factorization triv;
        triv.content = rational_lit(ctx.arena(), f_poly.constant_term());
        return ok(std::move(triv));
    }

    // 1. Collapse to single extension Q(θ).
    auto pe_res = compute_primitive_element(gens.alphas, gens.min_polys, ctx);
    if (pe_res.is_error()) return fail<Factorization>(pe_res.error());
    PrimitiveElementResult pe = std::move(pe_res.value());
    const AlgebraicNumber::CoeffVec& q_theta_coeffs = pe.min_poly_theta;
    RatPoly q_theta(q_theta_coeffs);
    q_theta.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    const std::size_t D = q_theta.degree();
    if (D == 0U) {
        return fail<Factorization>(nfact_error(
            CASErrorKind::InternalError,
            "factor_polynomial_tower_n: primitive element produced trivial min-poly"));
    }

    // 2. Shift search.
    const std::size_t user_bound = ctx.max_trager_tower_shift_attempts();
    const std::size_t default_bound = 2U * deg_f * D + 1U;
    const std::size_t max_attempts = (user_bound > 0U) ? user_bound : default_bound;
    const auto deadline = std::chrono::steady_clock::now() + ctx.timeout();

    // Opt-in hard deadline so inner Kronecker honours the tower budget, restored on exit (HC-F8-FACTORIZATIONTOWER-PERF).
    struct DeadlineGuard {
        symbolic::CASContext& c;
        std::chrono::steady_clock::time_point prev;
        ~DeadlineGuard() { c.set_hard_deadline(prev); }
    } deadline_guard{ctx, ctx.hard_deadline()};
    ctx.set_hard_deadline(std::min(ctx.hard_deadline(), deadline));

    for (std::size_t s_val = 0U; s_val <= max_attempts; ++s_val) {
        if (std::chrono::steady_clock::now() >= deadline) {
            return fail<Factorization>(nfact_error(
                CASErrorKind::Unimplemented,
                "factor_polynomial_tower_n [reason=WALL_BUDGET, ticket=F3.5]: "
                "ctx.timeout() exceeded during Trager shift search; increase "
                "ctx.set_timeout() or check separability of the extension."));
        }
        const Rational s_rat(BigInt(static_cast<std::int64_t>(s_val)));

        auto N_res = compute_trager_norm(q_theta, f_poly, s_rat);
        if (N_res.is_error()) {
            if (is_fatal_inner(N_res.error())) return fail<Factorization>(N_res.error());
            continue;
        }
        RatPoly N = std::move(N_res.value());
        normalize_rational_coefficients(N);
        if (N.is_zero() || N.degree() != deg_f * D) continue;
        if (!ratpoly_squarefree(N)) continue;

        // 3. Factor N over Q.  Render N as ExprPtr in var.
        PolyExpr N_pe;
        N_pe.reserve(N.size());
        for (const Rational& c : N.coefficients()) {
            N_pe.push_back(c.numerator().is_zero()
                           ? ExprPtr{}
                           : rational_lit(ctx.arena(), c));
        }
        normalize_poly(N_pe);
        auto N_expr_res = polynomial_to_expr(N_pe, var, ctx);
        if (N_expr_res.is_error()) return fail<Factorization>(N_expr_res.error());
        auto N_int_res = fti::clear_denominators_to_integer_poly(N, var, ctx);
        if (N_int_res.is_error()) return fail<Factorization>(N_int_res.error());

        auto Nfact = factor_over_integers(N_int_res.value(), var, ctx);
        if (Nfact.is_error()) {
            if (is_fatal_inner(Nfact.error())) return fail<Factorization>(Nfact.error());
            continue;
        }

        // 4. For each Q-factor g_i, compute f_i = gcd_{Q(θ)[x]}(f, g_i(x+sθ)).
        Factorization out;
        out.content = rational_lit(ctx.arena(), f_poly.leading_coeff());

        std::vector<AlgebraicNumber> f_qt = lift_f_to_qtheta(f_poly, q_theta_coeffs);

        std::size_t cumulative_deg = 0U;

        for (const auto& pf : Nfact.value().factors) {
            auto gi_parsed = parse_polynomial(pf.factor, var, ctx);
            if (gi_parsed.is_error()) {
                if (is_fatal_inner(gi_parsed.error()))
                    return fail<Factorization>(gi_parsed.error());
                continue;
            }
            auto gi_rat_res = poly_to_rational_poly(gi_parsed.value());
            if (gi_rat_res.is_error()) return fail<Factorization>(gi_rat_res.error());
            RatPoly gi_rat = std::move(gi_rat_res.value());
            normalize_rational_coefficients(gi_rat);
            if (gi_rat.is_zero()) continue;

            std::vector<AlgebraicNumber> gi_shifted =
                shift_qfactor_in_qtheta(gi_rat, s_rat, q_theta_coeffs);
            auto monic_gi = monic_qtheta(std::move(gi_shifted));
            if (monic_gi.is_error()) {
                if (is_fatal_inner(monic_gi.error()))
                    return fail<Factorization>(monic_gi.error());
                continue;
            }
            std::vector<AlgebraicNumber> f_copy = f_qt;
            auto monic_f = monic_qtheta(std::move(f_copy));
            if (monic_f.is_error()) return fail<Factorization>(monic_f.error());

            auto gcd_res = tower_detail::poly_extended_gcd<AlgebraicNumber>(
                monic_f.value(), monic_gi.value());
            if (gcd_res.is_error()) {
                if (is_fatal_inner(gcd_res.error()))
                    return fail<Factorization>(gcd_res.error());
                continue;
            }
            auto [g_coeffs, s_co, t_co] = std::move(gcd_res.value());
            (void)s_co; (void)t_co;
            tower_detail::strip_trailing(g_coeffs);
            if (g_coeffs.size() <= 1U) continue;
            auto factor_monic = monic_qtheta(std::move(g_coeffs));
            if (factor_monic.is_error()) return fail<Factorization>(factor_monic.error());
            cumulative_deg += qtheta_poly_degree(factor_monic.value());

            auto fexpr = render_qtheta_poly(factor_monic.value(), var,
                                            pe.theta_expr, ctx);
            if (fexpr.is_error()) return fail<Factorization>(fexpr.error());
            out.factors.push_back(PolynomialFactor{fexpr.value(), 1U});
        }

        if (out.factors.empty()) {
            // f is irreducible over Q(θ): return monic representative.
            PolyExpr monic_pe;
            monic_pe.reserve(f_poly.size());
            const Rational inv_lc = Rational(BigInt(1)) / f_poly.leading_coeff();
            for (const Rational& c : f_poly.coefficients()) {
                const Rational sc = c * inv_lc;
                monic_pe.push_back(sc.numerator().is_zero()
                                   ? ExprPtr{} : rational_lit(ctx.arena(), sc));
            }
            normalize_poly(monic_pe);
            auto m_expr = polynomial_to_expr(monic_pe, var, ctx);
            if (m_expr.is_error()) return fail<Factorization>(m_expr.error());
            out.factors.push_back(PolynomialFactor{m_expr.value(), 1U});
            return ok(std::move(out));
        }

        if (cumulative_deg != deg_f) {
            return fail<Factorization>(nfact_error(
                CASErrorKind::InternalError,
                "factor_polynomial_tower_n: Σ deg(f_i)=" +
                std::to_string(cumulative_deg) + " != deg(f)=" +
                std::to_string(deg_f) + " (Trager invariant violated)"));
        }
        return ok(std::move(out));
    }

    return fail<Factorization>(nfact_error(
        CASErrorKind::Unimplemented,
        "factor_polynomial_tower_n: no square-free Trager norm found within "
        "ctx.max_trager_tower_shift_attempts=" + std::to_string(max_attempts)));
}

}  // namespace algebra
}  // namespace cas
