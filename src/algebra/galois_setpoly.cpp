// A6 / CAS-L3-18 — Exact set-resolvent machinery for Galois deg ≥ 6.
//
// For monic f ∈ Q[x] of degree n with roots α_1..α_n, the k-set resolvent
//   R_k(y) = ∏_{|T|=k} (y − Σ_{i∈T} α_i)          (degree C(n,k))
// is computed exactly through power sums, with NO resultants:
//
//   1. Newton's identities give the power sums p_m = Σ_i α_i^m of f.
//   2. In the variables x_i = e^{s·α_i}, the elementary symmetric
//      function e_k(x) = Σ_{|T|=k} e^{s·Σ_T α} is the exponential
//      generating function of the k-set-sum power sums
//      q_m = Σ_{|T|=k} (Σ_T α)^m, and Newton's identities for e_k(x)
//      give the truncated-series recurrence
//          E_k(s) = (1/k) · Σ_{j=1..k} (−1)^{j−1} · P(js) · E_{k−j}(s),
//      where P(s) = Σ_i e^{s·α_i} has EGF coefficients p_m·j^m at P(js).
//   3. Inverse Newton on q_1..q_M (M = C(n,k)) recovers the monic
//      coefficients of R_k.
//
// Everything is exact rational arithmetic; soundness is by algebraic
// identity (Newton), certified structurally by q_0 = C(n,k).
//
// When Gal(f) acts on the roots, the degrees of the irreducible factors
// of a *squarefree* R_k over Q are exactly the orbit lengths of Gal(f)
// on k-subsets (Soicher-McKay 1985). Root-sum collisions make R_k
// non-squarefree; the caller then applies tschirnhaus_general (power
// basis β = P(α)), which preserves the permutation action of the Galois
// group on the roots whenever it is injective on them (⇔ the transformed
// polynomial is squarefree). Completeness of the power basis: the map
// S ↦ (Σ_{α∈S} α^m)_{m=1..n−1} is injective on root subsets (Vandermonde
// on distinct roots), so every pair of distinct k-subsets is separated by
// some coefficient vector.
//
// Spec: MISSING_FEATURES_SPECS/Galois_Groups.md (§Tschirnhaus fallback).

#include "galois_setpoly_internal.hpp"

#include "algebraic_tower_resultant.hpp"
#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::algebra::galois_setpoly {

namespace {

using primitive_internal::Deadline;

[[nodiscard]] Rational rat_zero() { return Rational(BigInt(0)); }
[[nodiscard]] Rational rat_one() { return Rational(BigInt(1)); }
[[nodiscard]] bool rat_is_zero(const Rational& r) {
    return r.numerator().is_zero();
}

// Newton's identities: power sums p_0..p_M of the roots of monic f.
// For f = x^n + c_{n−1}x^{n−1} + … + c_0:
//   p_m = −m·c_{n−m} − Σ_{i=1..m−1} c_{n−i}·p_{m−i}       (m ≤ n)
//   p_m = −Σ_{i=1..n} c_{n−i}·p_{m−i}                     (m > n)
[[nodiscard]] std::vector<Rational> power_sums_of_monic(const RatPoly& f,
                                                        std::size_t M) {
    const std::size_t n = f.degree();
    std::vector<Rational> p(M + 1U, rat_zero());
    p[0] = Rational(BigInt(static_cast<std::int64_t>(n)));
    for (std::size_t m = 1U; m <= M; ++m) {
        Rational acc = rat_zero();
        const std::size_t upper = (m <= n) ? (m - 1U) : n;
        for (std::size_t i = 1U; i <= upper; ++i) {
            acc = acc + f[n - i] * p[m - i];
        }
        if (m <= n) {
            acc = acc + f[n - m] * Rational(BigInt(static_cast<std::int64_t>(m)));
        }
        p[m] = rat_zero() - acc;
    }
    return p;
}

// EGF-basis product: (A·B)_m = Σ_{r=0..m} C(m,r)·a_r·b_{m−r}, where the
// stored coefficient u_m denotes m!·[s^m] of the series.
[[nodiscard]] std::vector<Rational> egf_multiply(
    const std::vector<Rational>& a, const std::vector<Rational>& b,
    const std::vector<std::vector<BigInt>>& binom) {
    const std::size_t M = a.size() - 1U;
    std::vector<Rational> out(M + 1U, rat_zero());
    for (std::size_t m = 0U; m <= M; ++m) {
        Rational acc = rat_zero();
        for (std::size_t r = 0U; r <= m; ++r) {
            if (rat_is_zero(a[r]) || rat_is_zero(b[m - r])) continue;
            acc = acc + a[r] * b[m - r] * Rational(binom[m][r]);
        }
        out[m] = acc;
    }
    return out;
}

}  // namespace

Result<RatPoly> kset_resolvent(const RatPoly& f_monic, std::size_t k,
                               const Deadline& deadline) {
    if (f_monic.is_zero() || f_monic.degree() < k || k == 0U ||
        !(f_monic.leading_coeff() == rat_one())) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_setpoly: kset_resolvent requires a monic "
                       "polynomial of degree >= k >= 1"});
    }
    const std::size_t n = f_monic.degree();
    // M = C(n,k), computed exactly.
    BigInt m_big(1);
    for (std::size_t i = 0U; i < k; ++i) {
        m_big = m_big * BigInt(static_cast<std::int64_t>(n - i));
        m_big = m_big / BigInt(static_cast<std::int64_t>(i + 1U));
    }
    const std::size_t M = static_cast<std::size_t>(m_big.to_u64());

    // Pascal triangle up to M (exact binomials for the EGF products).
    std::vector<std::vector<BigInt>> binom(M + 1U);
    for (std::size_t m = 0U; m <= M; ++m) {
        binom[m].assign(m + 1U, BigInt(1));
        for (std::size_t r = 1U; r < m; ++r) {
            binom[m][r] = binom[m - 1U][r - 1U] + binom[m - 1U][r];
        }
    }

    // p_m of f, then P_j (EGF coefficients p_m·j^m) on demand.
    const std::vector<Rational> p = power_sums_of_monic(f_monic, M);
    auto p_scaled = [&](std::size_t j) {
        std::vector<Rational> out(M + 1U, rat_zero());
        BigInt jp(1);
        const BigInt jb(static_cast<std::int64_t>(j));
        for (std::size_t m = 0U; m <= M; ++m) {
            out[m] = p[m] * Rational(jp);
            jp = jp * jb;
        }
        return out;
    };

    // E_0 = 1; E_k = (1/k)·Σ_{j=1..k} (−1)^{j−1}·P_j·E_{k−j}.
    std::vector<std::vector<Rational>> E;
    E.reserve(k + 1U);
    {
        std::vector<Rational> e0(M + 1U, rat_zero());
        e0[0] = rat_one();
        E.push_back(std::move(e0));
    }
    for (std::size_t kk = 1U; kk <= k; ++kk) {
        if (primitive_internal::deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "galois_setpoly: deadline exceeded during "
                           "k-set power-sum recurrence"});
        }
        std::vector<Rational> acc(M + 1U, rat_zero());
        for (std::size_t j = 1U; j <= kk; ++j) {
            const auto term = egf_multiply(p_scaled(j), E[kk - j], binom);
            const bool negative = (j % 2U) == 0U;
            for (std::size_t m = 0U; m <= M; ++m) {
                acc[m] = negative ? (acc[m] - term[m]) : (acc[m] + term[m]);
            }
        }
        const Rational inv_k(BigInt(1), BigInt(static_cast<std::int64_t>(kk)));
        for (auto& v : acc) v = v * inv_k;
        E.push_back(std::move(acc));
    }
    const std::vector<Rational>& q = E[k];

    // Structural certificate: q_0 must equal C(n,k).
    if (!(q[0] == Rational(m_big))) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "galois_setpoly: k-set power-sum recurrence violated "
                       "q_0 = C(n,k)"});
    }

    // Inverse Newton: m·e_m = Σ_{i=1..m} (−1)^{i−1}·q_i·e_{m−i}.
    std::vector<Rational> e(M + 1U, rat_zero());
    e[0] = rat_one();
    for (std::size_t m = 1U; m <= M; ++m) {
        Rational acc = rat_zero();
        for (std::size_t i = 1U; i <= m; ++i) {
            const Rational term = q[i] * e[m - i];
            acc = ((i % 2U) == 1U) ? (acc + term) : (acc - term);
        }
        e[m] = acc / Rational(BigInt(static_cast<std::int64_t>(m)));
    }

    // R(y) = Σ_{m=0..M} (−1)^m·e_m·y^{M−m}.
    std::vector<Rational> coeffs(M + 1U, rat_zero());
    for (std::size_t m = 0U; m <= M; ++m) {
        coeffs[M - m] = ((m % 2U) == 0U) ? e[m] : (rat_zero() - e[m]);
    }
    RatPoly r(std::move(coeffs));
    r.normalize(rat_is_zero);
    if (r.is_zero() || r.degree() != M) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "galois_setpoly: k-set resolvent has unexpected "
                       "degree"});
    }
    return ok(std::move(r));
}

Result<RatPoly> two_set_resolvent(const RatPoly& f_monic,
                                  const Deadline& deadline) {
    if (f_monic.is_zero() || f_monic.degree() < 2U ||
        !(f_monic.leading_coeff() == rat_one())) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_setpoly: two_set_resolvent requires a monic "
                       "polynomial of degree >= 2"});
    }
    return kset_resolvent(f_monic, 2U, deadline);
}

Result<RatPoly> three_set_resolvent(const RatPoly& f_monic,
                                    const Deadline& deadline) {
    if (f_monic.is_zero() || f_monic.degree() < 3U ||
        !(f_monic.leading_coeff() == rat_one())) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_setpoly: three_set_resolvent requires a "
                       "monic polynomial of degree >= 3"});
    }
    return kset_resolvent(f_monic, 3U, deadline);
}

Result<bool> is_squarefree_q(const RatPoly& p) {
    if (p.is_zero()) {
        return fail<bool>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_setpoly: squarefree test on zero polynomial"});
    }
    if (p.degree() <= 1U) return ok(true);
    // p' (formal derivative).
    std::vector<Rational> dcoeffs;
    dcoeffs.reserve(p.degree());
    for (std::size_t i = 1U; i < p.size(); ++i) {
        dcoeffs.push_back(p[i] * Rational(BigInt(static_cast<std::int64_t>(i))));
    }
    RatPoly dp(std::move(dcoeffs));
    dp.normalize(rat_is_zero);
    // Euclidean gcd in Q[y]; squarefree ⇔ gcd(p, p') constant.
    RatPoly a = p;
    RatPoly b = dp;
    while (!b.is_zero() && b.degree() > 0U) {
        // a mod b.
        RatPoly rem = a;
        const Rational& blc = b.leading_coeff();
        while (!rem.is_zero() && rem.degree() >= b.degree()) {
            const std::size_t shift = rem.degree() - b.degree();
            const Rational c = rem.leading_coeff() / blc;
            for (std::size_t i = 0U; i < b.size(); ++i) {
                rem[shift + i] = rem[shift + i] - c * b[i];
            }
            rem.normalize(rat_is_zero);
        }
        a = std::move(b);
        b = std::move(rem);
    }
    // Loop exits with b either zero (last nonzero remainder had degree ≥ 1
    // ⇒ gcd(p,p′) non-constant ⇒ NOT squarefree) or a nonzero constant
    // (⇒ gcd is a unit ⇒ squarefree).
    return ok(!b.is_zero());
}

Result<RatPoly> tschirnhaus_general(const RatPoly& f_monic,
                                    const std::vector<BigInt>& c,
                                    const Deadline& deadline) {
    if (f_monic.is_zero() || f_monic.degree() < 2U ||
        !(f_monic.leading_coeff() == rat_one())) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_setpoly: tschirnhaus_general requires a "
                       "monic polynomial of degree >= 2"});
    }
    if (c.empty()) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_setpoly: tschirnhaus_general requires "
                       "deg P >= 1"});
    }
    bool all_zero = true;
    for (const auto& cm : c) {
        if (!cm.is_zero()) all_zero = false;
    }
    if (all_zero) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "galois_setpoly: tschirnhaus_general requires a "
                       "nonzero transform P"});
    }
    const std::size_t n = f_monic.degree();
    // g(y) = ∏_i (y − P(α_i)) = Res_x(f(x), y − P(x)) with
    // P(x) = Σ_{m=1..|c|} c_m·x^m, degree n in y (f monic) — eval-interp
    // with n+1 integer sample points.
    const std::size_t num_pts = n + 1U;
    std::vector<Rational> xs, ys;
    xs.reserve(num_pts);
    ys.reserve(num_pts);
    for (std::size_t j = 1U; j <= num_pts; ++j) {
        if (primitive_internal::deadline_exceeded(deadline)) {
            return fail<RatPoly>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "galois_setpoly: deadline exceeded during "
                           "Tschirnhaus eval-interp"});
        }
        const Rational y_j(BigInt(static_cast<std::int64_t>(j)));
        // q_j(x) = y_j − Σ_m c_m·x^m.
        std::vector<Rational> qc;
        qc.reserve(c.size() + 1U);
        qc.push_back(y_j);
        for (const auto& cm : c) {
            qc.push_back(Rational(BigInt(0)) - Rational(cm));
        }
        while (qc.size() > 1U && rat_is_zero(qc.back())) qc.pop_back();
        auto res = resultant_generic<Rational>(
            f_monic.coefficients(), qc, nullptr, deadline);
        if (res.is_error()) return fail<RatPoly>(res.error());
        xs.push_back(y_j);
        ys.push_back(res.value());
    }
    // Newton divided-difference interpolation (degree ≤ n).
    std::vector<Rational> coef = ys;
    for (std::size_t j = 1U; j < num_pts; ++j) {
        for (std::size_t i = num_pts - 1U; i >= j; --i) {
            coef[i] = (coef[i] - coef[i - 1U]) / (xs[i] - xs[i - j]);
            if (i == j) break;
        }
    }
    // Expand Newton form to monomial coefficients.
    std::vector<Rational> poly{rat_zero()};
    std::vector<Rational> basis{rat_one()};  // ∏_{k<j}(y − x_k)
    for (std::size_t j = 0U; j < num_pts; ++j) {
        if (poly.size() < basis.size()) poly.resize(basis.size(), rat_zero());
        for (std::size_t i = 0U; i < basis.size(); ++i) {
            poly[i] = poly[i] + coef[j] * basis[i];
        }
        // basis *= (y − x_j).
        std::vector<Rational> nb(basis.size() + 1U, rat_zero());
        for (std::size_t i = 0U; i < basis.size(); ++i) {
            nb[i + 1U] = nb[i + 1U] + basis[i];
            nb[i] = nb[i] - basis[i] * xs[j];
        }
        basis = std::move(nb);
    }
    RatPoly g(std::move(poly));
    g.normalize(rat_is_zero);
    if (g.is_zero() || g.degree() != n) {
        return fail<RatPoly>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "galois_setpoly: Tschirnhaus minimal polynomial has "
                       "unexpected degree"});
    }
    // Res_x(f, y−P(x)) for monic f can carry a unit sign — normalize to
    // monic (the roots are what matter).
    if (!(g.leading_coeff() == rat_one())) {
        const Rational lc = g.leading_coeff();
        for (auto& cc : g.coefficients()) cc = cc / lc;
    }
    return ok(std::move(g));
}

Result<RatPoly> tschirnhaus_quadratic(const RatPoly& f_monic,
                                      const BigInt& c,
                                      const Deadline& deadline) {
    // β = α² + c·α: the degree-2 special case of the general transform.
    return tschirnhaus_general(
        f_monic, std::vector<BigInt>{c, BigInt(1)}, deadline);
}

}  // namespace cas::algebra::galois_setpoly
