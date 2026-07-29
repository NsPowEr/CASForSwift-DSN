// A6 — GF(q) construction and the Möbius/semilinear point permutations on
// P¹(F_q) (see perm_construct_fields_internal.hpp). F_q is built as
// F_p[x]/(f) with f found by exhaustive search and certified irreducible by
// Rabin's test (x^{p^e} ≡ x mod f, and gcd(x^{p^{e/r}} − x, f) = 1 for every
// prime r | e) — fully algorithmic, nothing transcribed.

#include "perm_construct_fields_internal.hpp"

#include "cas/error.hpp"
#include "cas/result.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::algebra::permgrp {

namespace {

[[nodiscard]] bool is_prime_u32(std::uint32_t p) {
    if (p < 2U) return false;
    for (std::uint32_t d = 2U; d * d <= p; ++d) {
        if (p % d == 0U) return false;
    }
    return true;
}

// Dense F_p[x] polynomial, coeffs[i] = coefficient of x^i (may carry
// trailing zeros; degree is tracked by inspection).
using Poly = std::vector<std::uint32_t>;

[[nodiscard]] std::size_t pdeg(const Poly& a) {
    std::size_t d = a.size();
    while (d > 0U && a[d - 1U] == 0U) --d;
    return d;  // number of significant coeffs; 0 means the zero polynomial
}

[[nodiscard]] Poly pmul(const Poly& a, const Poly& b, std::uint32_t p) {
    const std::size_t da = pdeg(a);
    const std::size_t db = pdeg(b);
    if (da == 0U || db == 0U) return {};
    Poly c(da + db - 1U, 0U);
    for (std::size_t i = 0U; i < da; ++i) {
        if (a[i] == 0U) continue;
        for (std::size_t j = 0U; j < db; ++j) {
            c[i + j] = (c[i + j] + a[i] * b[j]) % p;
        }
    }
    return c;
}

// a mod f, with f monic of significant length lf (degree lf−1 ≥ 1).
[[nodiscard]] Poly pmod(Poly a, const Poly& f, std::uint32_t p) {
    const std::size_t lf = pdeg(f);
    for (std::size_t la = pdeg(a); la >= lf; la = pdeg(a)) {
        const std::uint32_t lead = a[la - 1U];
        const std::size_t shift = la - lf;
        for (std::size_t i = 0U; i < lf; ++i) {
            const std::uint32_t sub = (lead * f[i]) % p;
            a[shift + i] = (a[shift + i] + p - sub) % p;
        }
    }
    return a;
}

[[nodiscard]] Poly pgcd(Poly a, Poly b, std::uint32_t p) {
    // Make the divisor monic each round so pmod applies.
    while (pdeg(b) != 0U) {
        const std::size_t lb = pdeg(b);
        std::uint32_t inv_lead = 1U;
        for (std::uint32_t y = 1U; y < p; ++y) {
            if ((b[lb - 1U] * y) % p == 1U) { inv_lead = y; break; }
        }
        Poly bm = b;
        for (auto& c : bm) c = (c * inv_lead) % p;
        Poly r = pmod(std::move(a), bm, p);
        a = std::move(b);
        b = std::move(r);
    }
    return a;
}

// x^(p^j) mod f by j successive p-th powers (each by repeated multiplication;
// p ≤ 15 here, exact and tiny).
[[nodiscard]] Poly x_pow_p_tower(const Poly& f, std::uint32_t p,
                                 std::size_t j) {
    Poly acc{0U, 1U};  // x
    for (std::size_t step = 0U; step < j; ++step) {
        Poly base = acc;
        Poly res{1U};
        for (std::uint32_t k = 0U; k < p; ++k) {
            res = pmod(pmul(res, base, p), f, p);
        }
        acc = std::move(res);
    }
    return acc;
}

// Rabin irreducibility for monic f of degree e ≥ 2 over F_p.
[[nodiscard]] bool rabin_irreducible(const Poly& f, std::uint32_t p,
                                     std::size_t e) {
    // x^{p^e} ≡ x (mod f)?
    Poly diff = x_pow_p_tower(f, p, e);
    diff.resize(std::max<std::size_t>(diff.size(), 2U), 0U);
    diff[1] = (diff[1] + p - 1U) % p;
    if (pdeg(pmod(std::move(diff), f, p)) != 0U) return false;
    // For each prime r | e: gcd(x^{p^{e/r}} − x, f) must be constant.
    std::size_t m = e;
    for (std::size_t r = 2U; r <= m; ++r) {
        if (m % r != 0U) continue;
        while (m % r == 0U) m /= r;
        Poly xr = x_pow_p_tower(f, p, e / r);
        xr.resize(std::max<std::size_t>(xr.size(), 2U), 0U);
        xr[1] = (xr[1] + p - 1U) % p;
        if (pdeg(pgcd(xr, f, p)) != 1U) return false;
    }
    return true;
}

[[nodiscard]] Poly decode_poly(std::uint32_t code, std::uint32_t p,
                               std::size_t e) {
    Poly c(e, 0U);
    for (std::size_t i = 0U; i < e; ++i) {
        c[i] = code % p;
        code /= p;
    }
    return c;
}

[[nodiscard]] std::uint32_t encode_poly(const Poly& c, std::uint32_t p) {
    std::uint32_t code = 0U;
    for (std::size_t i = c.size(); i > 0U; --i) code = code * p + c[i - 1U];
    return code;
}

}  // namespace

std::uint32_t Gf::add(std::uint32_t a, std::uint32_t b) const {
    std::uint32_t out = 0U;
    std::uint32_t mult = 1U;
    for (std::size_t i = 0U; i < e; ++i) {
        out += ((a % p + b % p) % p) * mult;
        a /= p;
        b /= p;
        mult *= p;
    }
    return out;
}

Result<Gf> build_gf(std::uint32_t p, std::size_t e) {
    if (!is_prime_u32(p) || e < 1U) {
        return fail<Gf>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "build_gf: p must be prime and e >= 1"});
    }
    std::uint64_t q64 = 1U;
    for (std::size_t i = 0U; i < e; ++i) q64 *= p;
    if (q64 + 1U > 255U) {
        return fail<Gf>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "build_gf: q + 1 exceeds the Perm image bound (255)"});
    }
    Gf gf;
    gf.p = p;
    gf.e = e;
    gf.q = static_cast<std::uint32_t>(q64);

    // Find a monic irreducible f = x^e + Σ c_i x^i by exhaustive search.
    // (Existence is classical; the search space p^e is tiny here.)
    Poly f;
    if (e == 1U) {
        f = {0U, 1U};  // x — GF(p) itself, reduction is plain mod p
    } else {
        for (std::uint32_t code = 0U; code < gf.q; ++code) {
            Poly cand = decode_poly(code, p, e);
            cand.push_back(1U);  // monic x^e
            if (rabin_irreducible(cand, p, e)) { f = std::move(cand); break; }
        }
        if (f.empty()) {
            // Unreachable: an irreducible of every degree exists over F_p.
            return fail<Gf>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = "build_gf: no irreducible polynomial found"});
        }
    }

    // Multiplication table (q ≤ 254 ⇒ at most ~64k entries).
    gf.mul.assign(static_cast<std::size_t>(gf.q) * gf.q, 0U);
    for (std::uint32_t a = 0U; a < gf.q; ++a) {
        const Poly pa = decode_poly(a, p, e);
        for (std::uint32_t b = 0U; b <= a; ++b) {
            const Poly pb = decode_poly(b, p, e);
            Poly prod = pmod(pmul(pa, pb, p), f, p);
            prod.resize(e, 0U);
            const std::uint32_t v = encode_poly(prod, p);
            gf.mul[static_cast<std::size_t>(a) * gf.q + b] = v;
            gf.mul[static_cast<std::size_t>(b) * gf.q + a] = v;
        }
    }
    // Negation, inverse, Frobenius — all from the tables just built.
    gf.neg.assign(gf.q, 0U);
    for (std::uint32_t a = 0U; a < gf.q; ++a) {
        Poly pa = decode_poly(a, p, e);
        for (auto& c : pa) c = (p - c) % p;
        gf.neg[a] = encode_poly(pa, p);
    }
    gf.inv.assign(gf.q, 0U);
    for (std::uint32_t a = 1U; a < gf.q; ++a) {
        for (std::uint32_t b = 1U; b < gf.q; ++b) {
            if (gf.mul[static_cast<std::size_t>(a) * gf.q + b] == 1U) {
                gf.inv[a] = b;
                break;
            }
        }
    }
    gf.frob.assign(gf.q, 0U);
    for (std::uint32_t a = 0U; a < gf.q; ++a) {
        std::uint32_t acc = 1U;
        for (std::uint32_t k = 0U; k < p; ++k) {
            acc = gf.mul[static_cast<std::size_t>(acc) * gf.q + a];
        }
        gf.frob[a] = acc;
    }
    // Generator of F_q^*: exhaustive exact order test.
    gf.gamma = 0U;
    for (std::uint32_t g = 1U; g < gf.q && gf.gamma == 0U; ++g) {
        std::uint32_t x = g;
        std::uint32_t ord = 1U;
        while (x != 1U) {
            x = gf.mul[static_cast<std::size_t>(x) * gf.q + g];
            ++ord;
        }
        if (ord == gf.q - 1U) gf.gamma = g;
    }
    return ok(std::move(gf));
}

namespace {

// Assemble a P¹(F_q) permutation from its action on field points and the
// image of ∞ (point q).
template <typename FieldMap>
[[nodiscard]] Perm p1_perm(const Gf& gf, FieldMap&& on_field,
                           std::uint32_t infinity_image) {
    Perm perm(gf.q + 1U);
    for (std::uint32_t x = 0U; x < gf.q; ++x) {
        perm[x] = static_cast<std::uint8_t>(on_field(x));
    }
    perm[gf.q] = static_cast<std::uint8_t>(infinity_image);
    return perm;
}

}  // namespace

std::vector<Perm> pgl2_point_gens(const Gf& gf) {
    std::vector<Perm> out;
    // x ↦ x + 1 (∞ fixed).
    out.push_back(p1_perm(
        gf, [&gf](std::uint32_t x) { return gf.add(x, 1U); }, gf.q));
    // x ↦ γx (0, ∞ fixed).
    out.push_back(gamma_mult_point_perm(gf));
    // x ↦ 1/x (0 ↔ ∞).
    out.push_back(p1_perm(
        gf,
        [&gf](std::uint32_t x) { return x == 0U ? gf.q : gf.inv[x]; }, 0U));
    return out;
}

std::vector<Perm> psl2_point_gens(const Gf& gf) {
    std::vector<Perm> out;
    out.push_back(p1_perm(
        gf, [&gf](std::uint32_t x) { return gf.add(x, 1U); }, gf.q));
    // x ↦ γ²x — the image of diag(γ, γ⁻¹) ∈ SL(2,q).
    const std::uint32_t g2 =
        gf.mul[static_cast<std::size_t>(gf.gamma) * gf.q + gf.gamma];
    out.push_back(p1_perm(
        gf,
        [&gf, g2](std::uint32_t x) {
            return gf.mul[static_cast<std::size_t>(g2) * gf.q + x];
        },
        gf.q));
    // x ↦ −1/x — the image of [[0,−1],[1,0]] ∈ SL(2,q).
    out.push_back(p1_perm(
        gf,
        [&gf](std::uint32_t x) {
            return x == 0U ? gf.q : gf.neg[gf.inv[x]];
        },
        0U));
    return out;
}

Perm frobenius_point_perm(const Gf& gf) {
    return p1_perm(
        gf, [&gf](std::uint32_t x) { return gf.frob[x]; }, gf.q);
}

Perm gamma_mult_point_perm(const Gf& gf) {
    return p1_perm(
        gf,
        [&gf](std::uint32_t x) {
            return gf.mul[static_cast<std::size_t>(gf.gamma) * gf.q + x];
        },
        gf.q);
}

}  // namespace cas::algebra::permgrp
