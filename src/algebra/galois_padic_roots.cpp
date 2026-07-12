// A6 Brick 3a — the mod-p layer of the p-adic splitting engine: modulus Φ
// (a degree-L factor of f mod p when one exists, else exhaustive search
// certified by Rabin) and the roots of f in the residue field GF(p^L)
// (linear read-off / Frobenius orbit of t / Cantor-Zassenhaus split).
// Lifting, Frobenius and certification live in galois_padic_split.cpp.

#include "galois_padic_detail.hpp"

#include "algebraic_tower_resultant.hpp"
#include "cas/bigint.hpp"
#include "cas/error.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"
#include "polynomial_internal.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <utility>
#include <vector>

namespace cas::algebra::galois_padic::detail {

namespace {

using primitive_internal::Deadline;
using primitive_internal::deadline_exceeded;

// The prime divisors of L by trial division (L is a tiny lcm of degrees).
[[nodiscard]] std::vector<std::size_t> prime_divisors(std::size_t L) {
    std::vector<std::size_t> out;
    std::size_t m = L;
    for (std::size_t d = 2U; d * d <= m; ++d) {
        if (m % d == 0U) {
            out.push_back(d);
            while (m % d == 0U) m /= d;
        }
    }
    if (m > 1U) out.push_back(m);
    return out;
}

// Rabin irreducibility over F_p for a monic candidate of degree L ≥ 2:
// X^{p^L} ≡ X mod (cand, p), and gcd(X^{p^{L/r}} − X, cand) = 1 for every
// prime r | L. Arithmetic runs over the trivial field presentation
// GR(p, 1) so the same RingPoly kernels serve every layer.
[[nodiscard]] Result<bool> rabin_irreducible_fp(const IntPoly& cand,
                                                const BigInt& p) {
    IntPoly t_poly(std::vector<BigInt>{BigInt(0), BigInt(1)});
    auto fp = PadicRing::make(p, 1U, t_poly);
    if (fp.is_error()) return fail<bool>(fp.error());
    const PadicRing& F = fp.value();
    auto cnd = to_field_poly(F, cand);
    if (cnd.is_error()) return fail<bool>(cnd.error());
    const RingPoly& C = cnd.value();
    const std::size_t L = cand.degree();
    const RingPoly X{F.zero(), F.one()};
    auto xq = rp_powmod(F, X, bigint_pow_nonnegative(p, L), C);
    if (xq.is_error()) return fail<bool>(xq.error());
    // X^{p^L} − X must vanish mod cand.
    RingPoly diff = xq.value();
    if (diff.size() < 2U) diff.resize(2U, F.zero());
    diff[1] = F.sub(diff[1], F.one());
    rp_normalize(F, diff);
    if (!diff.empty()) return ok(false);
    for (const std::size_t r : prime_divisors(L)) {
        auto xr = rp_powmod(F, X, bigint_pow_nonnegative(p, L / r), C);
        if (xr.is_error()) return fail<bool>(xr.error());
        RingPoly d2 = xr.value();
        if (d2.size() < 2U) d2.resize(2U, F.zero());
        d2[1] = F.sub(d2[1], F.one());
        rp_normalize(F, d2);
        auto g = rp_gcd_monic(F, C, std::move(d2));
        if (g.is_error()) return fail<bool>(g.error());
        if (g.value().size() != 1U) return ok(false);
    }
    return ok(true);
}

}  // namespace

Result<void> poll(symbolic::CASContext* ctx, const Deadline& dl) {
    if (ctx) {
        if (auto chk = ctx->check_interrupt(); chk.is_error()) {
            return chk;
        }
    }
    if (deadline_exceeded(dl)) {
        return fail<void>(
            CASError{.kind = CASErrorKind::Unimplemented,
                     .message = "galois_padic: deadline exceeded"});
    }
    return ok();
}

Result<RingPoly> to_field_poly(const PadicRing& K, const IntPoly& g) {
    RingPoly G;
    G.reserve(g.size());
    for (const auto& c : g.coefficients()) G.push_back(K.from_int(c));
    rp_normalize(K, G);
    if (G.empty()) {
        return fail<RingPoly>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "galois_padic: polynomial vanishes mod p"});
    }
    auto linv = K.inv(G.back());
    if (linv.is_error()) return fail<RingPoly>(linv.error());
    for (auto& c : G) c = K.mul(c, linv.value());
    return ok(std::move(G));
}

Result<IntPoly> find_phi(const BigInt& p, const std::vector<IntPoly>& factors,
                         std::size_t L, symbolic::CASContext* ctx,
                         const Deadline& dl) {
    for (const auto& g : factors) {
        if (g.degree() == L) return ok(g);
    }
    if (L == 1U) {
        return ok(IntPoly(std::vector<BigInt>{BigInt(0), BigInt(1)}));
    }
    const BigInt total = bigint_pow_nonnegative(p, L);
    for (BigInt m(0); m < total; m += BigInt(1)) {
        if (auto pr = poll(ctx, dl); pr.is_error()) {
            return fail<IntPoly>(pr.error());
        }
        std::vector<BigInt> coeffs(L + 1U, BigInt(0));
        BigInt rest = m;
        for (std::size_t i = 0U; i < L; ++i) {
            coeffs[i] = rest % p;
            rest /= p;
        }
        coeffs[L] = BigInt(1);
        IntPoly cand(std::move(coeffs));
        auto irr = rabin_irreducible_fp(cand, p);
        if (irr.is_error()) return fail<IntPoly>(irr.error());
        if (irr.value()) return ok(std::move(cand));
    }
    return fail<IntPoly>(CASError{
        .kind = CASErrorKind::InternalError,
        .message = "galois_padic: no irreducible of the required degree "
                   "found — impossible over a prime field"});
}

Result<std::vector<RingElem>> roots_in_field(const PadicRing& K,
                                             const IntPoly& g,
                                             std::mt19937& rng,
                                             symbolic::CASContext* ctx,
                                             const Deadline& dl) {
    const std::size_t d = g.degree();
    std::vector<RingElem> roots;
    if (d == 1U) {
        // Monic x + c (up to a unit): root −c/lead over F_p.
        auto gf = to_field_poly(K, g);
        if (gf.is_error()) return fail<std::vector<RingElem>>(gf.error());
        roots.push_back(K.neg(gf.value()[0]));
        return ok(std::move(roots));
    }
    // g equal to the modulus (mod p): roots are the Frobenius orbit of t.
    if (g.degree() == K.ext_degree()) {
        bool same = true;
        const auto& gc = g.coefficients();
        const auto& pc = K.phi().coefficients();
        for (std::size_t i = 0U; i <= d && same; ++i) {
            BigInt diff = (gc[i] - pc[i]) % K.prime();
            if (diff.is_negative()) diff += K.prime();
            same = diff.is_zero();
        }
        if (same) {
            RingElem r = K.basis_power(1U);
            for (std::size_t j = 0U; j < d; ++j) {
                roots.push_back(r);
                if (j + 1U < d) r = K.pow(r, K.prime());
            }
            return ok(std::move(roots));
        }
    }
    // Cantor-Zassenhaus split down to linear factors over K.
    auto gf = to_field_poly(K, g);
    if (gf.is_error()) return fail<std::vector<RingElem>>(gf.error());
    const BigInt q = bigint_pow_nonnegative(K.prime(), K.ext_degree());
    const std::uint64_t p64 = K.prime().to_u64();
    const auto random_elem = [&]() {
        RingElem a = K.zero();
        for (std::size_t i = 0U; i < K.ext_degree(); ++i) {
            a[i] = BigInt(static_cast<std::int64_t>(rng() % p64));
        }
        return a;
    };
    // Exact quotient by a monic divisor (remainder must vanish: both are
    // products of distinct linear factors of g).
    const auto div_exact = [&](const RingPoly& num,
                               const RingPoly& den) -> Result<RingPoly> {
        RingPoly rem = num;
        RingPoly quot(num.size() - den.size() + 1U, K.zero());
        while (rem.size() >= den.size() && !rem.empty()) {
            const std::size_t shift = rem.size() - den.size();
            const RingElem c = rem.back();  // den monic
            quot[shift] = c;
            for (std::size_t j = 0U; j < den.size(); ++j) {
                rem[shift + j] = K.sub(rem[shift + j], K.mul(c, den[j]));
            }
            rp_normalize(K, rem);
        }
        if (!rem.empty()) {
            return fail<RingPoly>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "galois_padic: non-exact division in CZ split"});
        }
        return ok(std::move(quot));
    };
    const std::function<Result<void>(RingPoly)> split =
        [&](RingPoly h) -> Result<void> {
        if (h.size() == 2U) {
            roots.push_back(K.neg(h[0]));  // h monic: X + c
            return ok();
        }
        while (true) {
            if (auto pr = poll(ctx, dl); pr.is_error()) return pr;
            const RingElem a = random_elem();
            RingPoly splitter;
            if (p64 == 2ULL) {
                // Trace splitter over GF(2^L): T(X) = Σ_{j<L} (aX)^{2^j}.
                const RingPoly ax{K.zero(), a};
                RingPoly cur = ax;
                RingPoly tr = ax;
                for (std::size_t j = 1U; j < K.ext_degree(); ++j) {
                    auto sq = rp_rem(K, rp_mul(K, cur, cur), h);
                    if (sq.is_error()) return fail<void>(sq.error());
                    cur = std::move(sq.value());
                    if (tr.size() < cur.size()) {
                        tr.resize(cur.size(), K.zero());
                    }
                    for (std::size_t i = 0U; i < cur.size(); ++i) {
                        tr[i] = K.add(tr[i], cur[i]);
                    }
                }
                rp_normalize(K, tr);
                splitter = std::move(tr);
            } else {
                // ((X + a)^{(q−1)/2} − 1) splits the roots by quadratic
                // character of (root + a).
                RingPoly xa{a, K.one()};
                auto pw = rp_powmod(K, std::move(xa),
                                    (q - BigInt(1)) / BigInt(2), h);
                if (pw.is_error()) return fail<void>(pw.error());
                splitter = std::move(pw.value());
                if (splitter.empty()) splitter.push_back(K.zero());
                splitter[0] = K.sub(splitter[0], K.one());
                rp_normalize(K, splitter);
            }
            if (splitter.empty()) continue;
            auto gc = rp_gcd_monic(K, h, std::move(splitter));
            if (gc.is_error()) return fail<void>(gc.error());
            RingPoly g1 = std::move(gc.value());
            if (g1.size() <= 1U || g1.size() >= h.size()) continue;
            auto g2 = div_exact(h, g1);
            if (g2.is_error()) return fail<void>(g2.error());
            auto r1 = split(std::move(g1));
            if (r1.is_error()) return r1;
            return split(std::move(g2.value()));
        }
    };
    auto sr = split(std::move(gf.value()));
    if (sr.is_error()) return fail<std::vector<RingElem>>(sr.error());
    return ok(std::move(roots));
}

IntPoly derivative(const IntPoly& f) {
    std::vector<BigInt> d;
    if (f.size() > 1U) {
        d.reserve(f.size() - 1U);
        for (std::size_t i = 1U; i < f.size(); ++i) {
            d.push_back(f.coefficients()[i] *
                        BigInt(static_cast<std::int64_t>(i)));
        }
    }
    IntPoly out(std::move(d));
    out.normalize([](const BigInt& v) { return v.is_zero(); });
    return out;
}

Result<RingElem> newton_lift(const IntPoly& g, const IntPoly& gprime,
                             const PadicRing& base_ring, RingElem x,
                             std::size_t m0, std::size_t k) {
    std::size_t m = m0;
    while (m < k) {
        const std::size_t m2 = std::min(2U * m, k);
        auto ring2 = base_ring.with_precision(m2);
        if (ring2.is_error()) return fail<RingElem>(ring2.error());
        const PadicRing& R2 = ring2.value();
        const RingElem gx = R2.eval_int_poly(g, x);
        auto dinv = R2.inv(R2.eval_int_poly(gprime, x));
        if (dinv.is_error()) {
            return fail<RingElem>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "galois_padic: derivative not a unit during "
                           "Newton lift (root not simple mod p)"});
        }
        x = R2.sub(x, R2.mul(gx, dinv.value()));
        m = m2;
    }
    return ok(std::move(x));
}

}  // namespace cas::algebra::galois_padic::detail
