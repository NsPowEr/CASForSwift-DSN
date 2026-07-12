// A6 Brick 3a — arithmetic of the unramified Galois ring GR(p^k, L) and of
// dense polynomials over it. See galois_padic_internal.hpp for the contract;
// root lifting and certification live in galois_padic_roots.cpp.
//
// Everything is exact BigInt arithmetic modulo p^k; the extension is
// presented by a fixed integer monic Φ irreducible mod p, so coordinates
// are stable across precision changes (Z[t]/(p^k, Φ) for growing k).

#include "galois_padic_internal.hpp"

#include "cas/bigint.hpp"
#include "cas/numtheory.hpp"
#include "cas/result.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace cas::algebra::galois_padic {

namespace {

// v mod m into [0, m).
[[nodiscard]] BigInt reduce_mod(const BigInt& v, const BigInt& m) {
    BigInt r = v % m;
    if (r.is_negative()) r += m;
    return r;
}

// Extended Euclid in F_p[t]: inverse of a modulo (p, phi). Dense vectors,
// index = degree. Returns empty vector iff a ≡ 0 mod (p, phi) is not a unit
// (phi irreducible mod p makes every nonzero residue invertible).
[[nodiscard]] std::vector<BigInt> fp_poly_inverse(std::vector<BigInt> a,
                                                  const IntPoly& phi,
                                                  const BigInt& p) {
    const auto norm = [&](std::vector<BigInt>& v) {
        for (auto& c : v) c = reduce_mod(c, p);
        while (!v.empty() && v.back().is_zero()) v.pop_back();
    };
    // r0 = phi, r1 = a; s0 = 0, s1 = 1; invariant s_i·a ≡ r_i (mod phi).
    std::vector<BigInt> r0(phi.coefficients());
    std::vector<BigInt> r1 = std::move(a);
    norm(r0);
    norm(r1);
    std::vector<BigInt> s0;
    std::vector<BigInt> s1{BigInt(1)};
    while (!r1.empty()) {
        // Divide r0 by r1 over F_p.
        std::vector<BigInt> q;
        auto linv_res =
            numtheory::modular_inverse(r1.back(), p);
        if (linv_res.is_error()) return {};  // p not prime — caller bug
        const BigInt& linv = linv_res.value();
        std::vector<BigInt> rem = r0;
        while (rem.size() >= r1.size() && !rem.empty()) {
            const std::size_t shift = rem.size() - r1.size();
            const BigInt c = reduce_mod(rem.back() * linv, p);
            if (q.size() < shift + 1U) q.resize(shift + 1U, BigInt(0));
            q[shift] = c;
            for (std::size_t j = 0U; j < r1.size(); ++j) {
                rem[shift + j] = reduce_mod(rem[shift + j] - c * r1[j], p);
            }
            while (!rem.empty() && rem.back().is_zero()) rem.pop_back();
        }
        // (r0, r1) ← (r1, rem); (s0, s1) ← (s1, s0 − q·s1) over F_p.
        std::vector<BigInt> s2 = s0;
        if (s2.size() < q.size() + s1.size()) {
            s2.resize(q.size() + s1.size(), BigInt(0));
        }
        for (std::size_t i = 0U; i < q.size(); ++i) {
            if (q[i].is_zero()) continue;
            for (std::size_t j = 0U; j < s1.size(); ++j) {
                s2[i + j] = reduce_mod(s2[i + j] - q[i] * s1[j], p);
            }
        }
        norm(s2);
        r0 = std::move(r1);
        r1 = std::move(rem);
        s0 = std::move(s1);
        s1 = std::move(s2);
    }
    if (r0.size() != 1U) return {};  // gcd not constant: not a unit
    auto ginv = numtheory::modular_inverse(r0[0], p);
    if (ginv.is_error()) return {};
    for (auto& c : s0) c = reduce_mod(c * ginv.value(), p);
    while (!s0.empty() && s0.back().is_zero()) s0.pop_back();
    return s0.empty() ? std::vector<BigInt>{BigInt(0)} : s0;
}

}  // namespace

Result<PadicRing> PadicRing::make(BigInt p, std::size_t k, IntPoly phi) {
    if (k == 0U || phi.is_zero() || phi.degree() < 1U ||
        !(phi.leading_coeff() == BigInt(1))) {
        return fail<PadicRing>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "PadicRing: need k >= 1 and a monic phi of degree "
                       ">= 1"});
    }
    if (p < BigInt(2)) {
        return fail<PadicRing>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "PadicRing: p must be a prime >= 2"});
    }
    PadicRing r;
    r.k_ = k;
    r.pk_ = bigint_pow_nonnegative(p, k);
    r.p_ = std::move(p);
    r.ext_ = phi.degree();
    // Reduce phi's coefficients into [0, p^k) once.
    std::vector<BigInt> pc = phi.coefficients();
    for (auto& c : pc) c = reduce_mod(c, r.pk_);
    pc.back() = BigInt(1);
    r.phi_ = IntPoly(std::move(pc));
    return ok(std::move(r));
}

RingElem PadicRing::zero() const { return RingElem(ext_, BigInt(0)); }

RingElem PadicRing::one() const {
    RingElem e = zero();
    e[0] = reduce_mod(BigInt(1), pk_);
    return e;
}

RingElem PadicRing::from_int(const BigInt& v) const {
    RingElem e = zero();
    e[0] = reduce_mod(v, pk_);
    return e;
}

RingElem PadicRing::basis_power(std::size_t j) const {
    RingElem e = zero();
    if (j < ext_) {
        e[j] = BigInt(1);
        return e;
    }
    // t^j for j ≥ L via repeated multiplication (rarely needed). For L = 1
    // the class of t is −Φ₀ (t ≡ −Φ₀ mod Φ = t + Φ₀).
    RingElem t = zero();
    if (ext_ > 1U) {
        t[1] = BigInt(1);
    } else {
        t = from_int(-phi_.coefficients()[0]);
    }
    RingElem acc = one();
    for (std::size_t i = 0U; i < j; ++i) acc = mul(acc, t);
    return acc;
}

RingElem PadicRing::add(const RingElem& a, const RingElem& b) const {
    RingElem c(ext_);
    for (std::size_t i = 0U; i < ext_; ++i) c[i] = reduce_mod(a[i] + b[i], pk_);
    return c;
}

RingElem PadicRing::sub(const RingElem& a, const RingElem& b) const {
    RingElem c(ext_);
    for (std::size_t i = 0U; i < ext_; ++i) c[i] = reduce_mod(a[i] - b[i], pk_);
    return c;
}

RingElem PadicRing::neg(const RingElem& a) const {
    RingElem c(ext_);
    for (std::size_t i = 0U; i < ext_; ++i) c[i] = reduce_mod(-a[i], pk_);
    return c;
}

RingElem PadicRing::mul(const RingElem& a, const RingElem& b) const {
    // Schoolbook convolution, then reduction by the monic Φ:
    // t^L ≡ −Σ_{j<L} Φ_j t^j.
    std::vector<BigInt> conv(2U * ext_ - 1U, BigInt(0));
    for (std::size_t i = 0U; i < ext_; ++i) {
        if (a[i].is_zero()) continue;
        for (std::size_t j = 0U; j < ext_; ++j) {
            if (b[j].is_zero()) continue;
            conv[i + j] += a[i] * b[j];
        }
    }
    const auto& pc = phi_.coefficients();
    for (std::size_t i = conv.size(); i-- > ext_;) {
        if (conv[i].is_zero()) continue;
        const BigInt c = reduce_mod(conv[i], pk_);
        conv[i] = BigInt(0);
        for (std::size_t j = 0U; j < ext_; ++j) {
            conv[i - ext_ + j] -= c * pc[j];
        }
    }
    RingElem out(ext_);
    for (std::size_t i = 0U; i < ext_; ++i) out[i] = reduce_mod(conv[i], pk_);
    return out;
}

RingElem PadicRing::pow(const RingElem& a, const BigInt& e) const {
    RingElem base = a;
    RingElem acc = one();
    BigInt exp = e;
    while (!exp.is_zero() && !exp.is_negative()) {
        if (!(exp % BigInt(2)).is_zero()) acc = mul(acc, base);
        exp /= BigInt(2);
        if (!exp.is_zero()) base = mul(base, base);
    }
    return acc;
}

bool PadicRing::is_zero(const RingElem& a) const {
    for (std::size_t i = 0U; i < ext_; ++i) {
        if (!a[i].is_zero()) return false;
    }
    return true;
}

bool PadicRing::equal(const RingElem& a, const RingElem& b) const {
    for (std::size_t i = 0U; i < ext_; ++i) {
        if (!(a[i] == b[i])) return false;
    }
    return true;
}

bool PadicRing::equal_mod_p(const RingElem& a, const RingElem& b) const {
    for (std::size_t i = 0U; i < ext_; ++i) {
        if (!reduce_mod(a[i] - b[i], p_).is_zero()) return false;
    }
    return true;
}

Result<RingElem> PadicRing::inv(const RingElem& a) const {
    // Inverse mod p by extended Euclid in F_p[t]/(phi), then Newton
    // z ← z(2 − a z): z correct mod p^{2^j} after j steps.
    std::vector<BigInt> a_p(a.begin(), a.end());
    std::vector<BigInt> z0 = fp_poly_inverse(std::move(a_p), phi_, p_);
    if (z0.empty()) {
        return fail<RingElem>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "PadicRing::inv: not a unit (zero mod p)"});
    }
    RingElem z = zero();
    for (std::size_t i = 0U; i < z0.size() && i < ext_; ++i) z[i] = z0[i];
    const RingElem two = from_int(BigInt(2));
    for (std::size_t m = 1U; m < k_; m *= 2U) {
        z = mul(z, sub(two, mul(a, z)));
    }
    return ok(std::move(z));
}

RingElem PadicRing::eval_int_poly(const IntPoly& f, const RingElem& x) const {
    RingElem acc = zero();
    const auto& fc = f.coefficients();
    for (std::size_t i = fc.size(); i-- > 0U;) {
        acc = mul(acc, x);
        acc[0] = reduce_mod(acc[0] + fc[i], pk_);
    }
    return acc;
}

std::optional<BigInt> PadicRing::integer_residue(const RingElem& a) const {
    for (std::size_t i = 1U; i < ext_; ++i) {
        if (!a[i].is_zero()) return std::nullopt;
    }
    BigInt half = pk_ / BigInt(2);
    if (a[0] > half) return a[0] - pk_;
    return a[0];
}

Result<PadicRing> PadicRing::with_precision(std::size_t k2) const {
    return make(p_, k2, phi_);
}

// ── dense polynomials over the ring ─────────────────────────────────────────

void rp_normalize(const PadicRing& R, RingPoly& a) {
    while (!a.empty() && R.is_zero(a.back())) a.pop_back();
}

RingPoly rp_mul(const PadicRing& R, const RingPoly& a, const RingPoly& b) {
    if (a.empty() || b.empty()) return {};
    RingPoly c(a.size() + b.size() - 1U, R.zero());
    for (std::size_t i = 0U; i < a.size(); ++i) {
        if (R.is_zero(a[i])) continue;
        for (std::size_t j = 0U; j < b.size(); ++j) {
            c[i + j] = R.add(c[i + j], R.mul(a[i], b[j]));
        }
    }
    rp_normalize(R, c);
    return c;
}

Result<RingPoly> rp_rem(const PadicRing& R, RingPoly a, const RingPoly& b) {
    if (b.empty()) {
        return fail<RingPoly>(
            CASError{.kind = CASErrorKind::InvalidArgument,
                     .message = "rp_rem: division by the zero polynomial"});
    }
    auto linv = R.inv(b.back());
    if (linv.is_error()) return fail<RingPoly>(linv.error());
    rp_normalize(R, a);
    while (a.size() >= b.size()) {
        const std::size_t shift = a.size() - b.size();
        const RingElem c = R.mul(a.back(), linv.value());
        for (std::size_t j = 0U; j < b.size(); ++j) {
            a[shift + j] = R.sub(a[shift + j], R.mul(c, b[j]));
        }
        rp_normalize(R, a);
    }
    return ok(std::move(a));
}

Result<RingPoly> rp_gcd_monic(const PadicRing& R, RingPoly a, RingPoly b) {
    if (R.precision() != 1U) {
        return fail<RingPoly>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "rp_gcd_monic: requires the residue field (k == 1)"});
    }
    rp_normalize(R, a);
    rp_normalize(R, b);
    while (!b.empty()) {
        auto r = rp_rem(R, std::move(a), b);
        if (r.is_error()) return r;
        a = std::move(b);
        b = std::move(r.value());
    }
    if (!a.empty()) {
        auto linv = R.inv(a.back());
        if (linv.is_error()) return fail<RingPoly>(linv.error());
        for (auto& c : a) c = R.mul(c, linv.value());
    }
    return ok(std::move(a));
}

Result<RingPoly> rp_powmod(const PadicRing& R, RingPoly base, const BigInt& e,
                           const RingPoly& m) {
    auto b0 = rp_rem(R, std::move(base), m);
    if (b0.is_error()) return b0;
    RingPoly b = std::move(b0.value());
    RingPoly acc{R.one()};
    BigInt exp = e;
    while (!exp.is_zero() && !exp.is_negative()) {
        if (!(exp % BigInt(2)).is_zero()) {
            auto t = rp_rem(R, rp_mul(R, acc, b), m);
            if (t.is_error()) return t;
            acc = std::move(t.value());
        }
        exp /= BigInt(2);
        if (!exp.is_zero()) {
            auto t = rp_rem(R, rp_mul(R, b, b), m);
            if (t.is_error()) return t;
            b = std::move(t.value());
        }
    }
    return ok(std::move(acc));
}

RingElem rp_eval(const PadicRing& R, const RingPoly& a, const RingElem& x) {
    RingElem acc = R.zero();
    for (std::size_t i = a.size(); i-- > 0U;) {
        acc = R.add(R.mul(acc, x), a[i]);
    }
    return acc;
}

}  // namespace cas::algebra::galois_padic
