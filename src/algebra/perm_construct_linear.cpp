// A6 — Prime-field linear/affine/projective families (see
// perm_construct_fields_internal.hpp): AGL(d,p) on F_p^d and the projective
// action of GL(d,p) on P^{d−1}(F_p). Exact modular arithmetic throughout.

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

// Multiplicative order of g modulo prime p (g ∈ [1,p)).
[[nodiscard]] std::uint32_t order_mod_p(std::uint32_t g, std::uint32_t p) {
    std::uint32_t x = g % p;
    std::uint32_t ord = 1U;
    while (x != 1U) {
        x = (x * g) % p;
        ++ord;
    }
    return ord;
}

// Smallest primitive root modulo prime p ≥ 3 (exists by classical theory;
// found by exhaustive order test — exact, no heuristics).
[[nodiscard]] std::uint32_t primitive_root_mod_p(std::uint32_t p) {
    for (std::uint32_t g = 2U; g < p; ++g) {
        if (order_mod_p(g, p) == p - 1U) return g;
    }
    return 0U;  // unreachable for prime p ≥ 3
}

// Row-major d×d matrix over F_p.
using MatFp = std::vector<std::uint32_t>;

[[nodiscard]] MatFp mat_identity(std::size_t d) {
    MatFp m(d * d, 0U);
    for (std::size_t i = 0U; i < d; ++i) m[i * d + i] = 1U;
    return m;
}

// w = M·v over F_p.
[[nodiscard]] std::vector<std::uint32_t> mat_vec(std::size_t d,
                                                 std::uint32_t p,
                                                 const MatFp& m,
                                                 const std::vector<std::uint32_t>& v) {
    std::vector<std::uint32_t> w(d, 0U);
    for (std::size_t i = 0U; i < d; ++i) {
        std::uint32_t acc = 0U;
        for (std::size_t j = 0U; j < d; ++j) acc += m[i * d + j] * v[j];
        w[i] = acc % p;
    }
    return w;
}

[[nodiscard]] std::vector<std::uint32_t> decode_vec(std::size_t d,
                                                    std::uint32_t p,
                                                    std::size_t idx) {
    std::vector<std::uint32_t> v(d);
    for (std::size_t i = 0U; i < d; ++i) {
        v[i] = static_cast<std::uint32_t>(idx % p);
        idx /= p;
    }
    return v;
}

[[nodiscard]] std::size_t encode_vec(std::uint32_t p,
                                     const std::vector<std::uint32_t>& v) {
    std::size_t idx = 0U;
    for (std::size_t i = v.size(); i > 0U; --i) idx = idx * p + v[i - 1U];
    return idx;
}

// The GL(d,p) generating matrices ⟨T, C, D⟩ (see header): transvection
// I + E_{01}, coordinate d-cycle, torus diag(γ,1,…,1).
[[nodiscard]] std::vector<MatFp> gl_gen_mats(std::size_t d, std::uint32_t p) {
    std::vector<MatFp> mats;
    if (d >= 2U) {
        MatFp t = mat_identity(d);
        t[0U * d + 1U] = 1U;
        mats.push_back(std::move(t));
        MatFp c(d * d, 0U);
        for (std::size_t i = 0U; i < d; ++i) c[((i + 1U) % d) * d + i] = 1U;
        mats.push_back(std::move(c));
    }
    if (p >= 3U) {
        MatFp dm = mat_identity(d);
        dm[0] = primitive_root_mod_p(p);
        mats.push_back(std::move(dm));
    }
    return mats;
}

[[nodiscard]] Result<void> validate_field_degree(std::size_t d,
                                                 std::uint32_t p,
                                                 std::size_t min_d) {
    if (!is_prime_u32(p) || d < min_d) {
        return fail<void>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "linear family: p must be prime and the dimension "
                       "must meet the family's minimum"});
    }
    std::size_t points = 1U;
    for (std::size_t i = 0U; i < d; ++i) {
        points *= p;
        if (points > 255U) {
            return fail<void>(CASError{
                .kind = CASErrorKind::InvalidArgument,
                .message = "linear family: p^d exceeds the Perm image "
                           "bound (255)"});
        }
    }
    return ok();
}

}  // namespace

Result<std::vector<Perm>> agl_gens(std::size_t d, std::uint32_t p) {
    auto valid = validate_field_degree(d, p, 1U);
    if (valid.is_error()) return fail<std::vector<Perm>>(valid.error());
    std::size_t npoints = 1U;
    for (std::size_t i = 0U; i < d; ++i) npoints *= p;

    std::vector<Perm> out;
    for (const MatFp& m : gl_gen_mats(d, p)) {
        Perm perm(npoints);
        for (std::size_t idx = 0U; idx < npoints; ++idx) {
            perm[idx] = static_cast<std::uint8_t>(
                encode_vec(p, mat_vec(d, p, m, decode_vec(d, p, idx))));
        }
        out.push_back(std::move(perm));
    }
    // Translation v ↦ v + e_0 (digit-0 increment mod p): with GL transitive
    // on the nonzero vectors, its conjugates give every translation.
    Perm trans(npoints);
    for (std::size_t idx = 0U; idx < npoints; ++idx) {
        std::vector<std::uint32_t> v = decode_vec(d, p, idx);
        v[0] = (v[0] + 1U) % p;
        trans[idx] = static_cast<std::uint8_t>(encode_vec(p, v));
    }
    out.push_back(std::move(trans));
    return ok(std::move(out));
}

Result<std::vector<Perm>> projective_gl_gens(std::size_t d, std::uint32_t p) {
    auto valid = validate_field_degree(d, p, 2U);
    if (valid.is_error()) return fail<std::vector<Perm>>(valid.error());
    std::size_t nvecs = 1U;
    for (std::size_t i = 0U; i < d; ++i) nvecs *= p;

    // Modular inverse in F_p by exhaustive scan (p ≤ 255 — exact).
    auto inv_mod = [p](std::uint32_t x) -> std::uint32_t {
        for (std::uint32_t y = 1U; y < p; ++y) {
            if ((x * y) % p == 1U) return y;
        }
        return 0U;  // unreachable for x ≢ 0
    };
    // Canonical representative: scale so the lowest-index nonzero
    // coordinate is 1.
    auto normalize = [&](std::vector<std::uint32_t> v) {
        std::size_t j = 0U;
        while (j < v.size() && v[j] == 0U) ++j;
        const std::uint32_t s = inv_mod(v[j]);
        for (auto& c : v) c = (c * s) % p;
        return v;
    };

    // Enumerate the normalised representatives in increasing encoding order.
    std::vector<std::size_t> reps;
    for (std::size_t idx = 1U; idx < nvecs; ++idx) {
        const std::vector<std::uint32_t> v = decode_vec(d, p, idx);
        if (encode_vec(p, normalize(v)) == idx) reps.push_back(idx);
    }
    if (reps.size() > 255U) {
        return fail<std::vector<Perm>>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "projective family: point count exceeds the Perm "
                       "image bound (255)"});
    }

    std::vector<Perm> out;
    for (const MatFp& m : gl_gen_mats(d, p)) {
        Perm perm(reps.size());
        for (std::size_t i = 0U; i < reps.size(); ++i) {
            const std::size_t image = encode_vec(
                p, normalize(mat_vec(d, p, m, decode_vec(d, p, reps[i]))));
            const auto it =
                std::lower_bound(reps.begin(), reps.end(), image);
            perm[i] = static_cast<std::uint8_t>(
                static_cast<std::size_t>(it - reps.begin()));
        }
        out.push_back(std::move(perm));
    }
    return ok(std::move(out));
}

}  // namespace cas::algebra::permgrp
