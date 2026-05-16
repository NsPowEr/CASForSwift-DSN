#pragma once

// AlgebraicElement<Coeff> — generic simple algebraic extension over an
// arbitrary coefficient ring (a "tower" rung).  Specialisations of
// interest:
//
//     AlgebraicElement<Rational>                    = Q(α)
//     AlgebraicElement<AlgebraicElement<Rational>>  = Q(α₁)(α₂) = Q(α₁,α₂)
//     AlgebraicElement<AlgebraicElement<...>>       = arbitrary depth
//
// The element is represented by its canonical polynomial value modulo a
// monic minimal polynomial m(x) over the coefficient ring.  All
// arithmetic is performed by template polynomial helpers parameterised
// over the coefficient ring; the only ring-specific dispatch lives in
// the `CoeffOps<Coeff>` trait (zero/one/inverse/zero-predicate).
//
// Design rationale (CLAUDE.md REGOLA ZERO compliance):
//   * No hardcoded set of "supported field types" — any Coeff with a
//     CoeffOps specialisation is admissible.
//   * No fixed maximum tower depth — recursion is type-level.
//   * The minimal polynomial is stored monic-normalised at construction
//     to guarantee a canonical representation; equality therefore tests
//     value coefficients only (after assuming a shared min_poly).
//
// File counts toward the 500-line limit; helpers split into
// algebraic_tower.cpp keep this header focused on the public template.

#include "cas/algebraic_number.hpp"
#include "cas/error.hpp"
#include "cas/rational.hpp"
#include "cas/result.hpp"

#include <cstddef>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace cas {
namespace algebra {

// Forward declaration of the public template.
template <typename Coeff>
class AlgebraicElement;

// === Ring trait dispatch =====================================================
//
// CoeffOps must expose:
//   - static Coeff zero_like(const Coeff& sample)
//   - static Coeff one_like (const Coeff& sample)
//   - static bool  is_zero  (const Coeff& c)
//   - static Result<Coeff> inverse(const Coeff& c)
//
// `_like` variants allow contexts (e.g. minimal polynomial of an inner
// extension) to be inherited from an existing element when constructing
// neutral elements.
template <typename Coeff>
struct CoeffOps;

template <>
struct CoeffOps<Rational> {
    [[nodiscard]] static Rational zero_like(const Rational& /*sample*/) {
        return Rational(BigInt(0));
    }
    [[nodiscard]] static Rational one_like(const Rational& /*sample*/) {
        return Rational(BigInt(1));
    }
    [[nodiscard]] static bool is_zero(const Rational& c) {
        return c.numerator().is_zero();
    }
    [[nodiscard]] static Result<Rational> inverse(const Rational& c) {
        if (c.numerator().is_zero()) {
            return fail<Rational>(CASError{CASErrorKind::InvalidArgument,
                "Rational inverse: division by zero", std::nullopt});
        }
        return ok(Rational(BigInt(1)) / c);
    }
};

// Specialisation for AlgebraicNumber (= Q(α) bottom-most rung).  This
// pins the "Q(α)" case to the existing battle-tested AlgebraicNumber
// implementation; higher towers build on top via AlgebraicElement<AlgebraicNumber>.
template <>
struct CoeffOps<AlgebraicNumber> {
    [[nodiscard]] static AlgebraicNumber zero_like(const AlgebraicNumber& sample) {
        return AlgebraicNumber({}, sample.min_poly());
    }
    [[nodiscard]] static AlgebraicNumber one_like(const AlgebraicNumber& sample) {
        return AlgebraicNumber({Rational(BigInt(1))}, sample.min_poly());
    }
    [[nodiscard]] static bool is_zero(const AlgebraicNumber& c) {
        return c.is_zero();
    }
    [[nodiscard]] static Result<AlgebraicNumber> inverse(const AlgebraicNumber& c) {
        return c.inverse();
    }
};

// Recursive specialisation for further-nested towers
// (AlgebraicElement<AlgebraicElement<X>>, …).  Defined below the
// template class declaration.
template <typename C>
struct CoeffOps<AlgebraicElement<C>>;

// === Generic polynomial helpers (Coeff-templated) ===========================
namespace tower_detail {

template <typename Coeff>
inline void strip_trailing(std::vector<Coeff>& v) {
    while (!v.empty() && CoeffOps<Coeff>::is_zero(v.back())) v.pop_back();
}

// Find any non-zero element in (a, b) to use as a "context sample" for
// CoeffOps::zero_like.  Caller guarantees at least one of the arguments
// is non-empty.
template <typename Coeff>
[[nodiscard]] inline const Coeff& any_sample(const std::vector<Coeff>& a,
                                              const std::vector<Coeff>& b) {
    if (!a.empty()) return a.front();
    return b.front();
}

template <typename Coeff>
[[nodiscard]] std::vector<Coeff> poly_add(const std::vector<Coeff>& a,
                                          const std::vector<Coeff>& b) {
    if (a.empty()) return b;
    if (b.empty()) return a;
    const Coeff zero = CoeffOps<Coeff>::zero_like(any_sample(a, b));
    std::vector<Coeff> r(std::max(a.size(), b.size()), zero);
    for (std::size_t i = 0; i < a.size(); ++i) r[i] = r[i] + a[i];
    for (std::size_t i = 0; i < b.size(); ++i) r[i] = r[i] + b[i];
    strip_trailing(r);
    return r;
}

template <typename Coeff>
[[nodiscard]] std::vector<Coeff> poly_neg(const std::vector<Coeff>& a) {
    std::vector<Coeff> r;
    r.reserve(a.size());
    for (const auto& c : a) r.push_back(-c);
    return r;
}

template <typename Coeff>
[[nodiscard]] std::vector<Coeff> poly_sub(const std::vector<Coeff>& a,
                                          const std::vector<Coeff>& b) {
    return poly_add(a, poly_neg(b));
}

template <typename Coeff>
[[nodiscard]] std::vector<Coeff> poly_mul(const std::vector<Coeff>& a,
                                          const std::vector<Coeff>& b) {
    if (a.empty() || b.empty()) return {};
    const Coeff zero = CoeffOps<Coeff>::zero_like(a.front());
    std::vector<Coeff> r(a.size() + b.size() - 1U, zero);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] = r[i + j] + a[i] * b[j];
        }
    }
    strip_trailing(r);
    return r;
}

// Polynomial division with remainder over a field.  Returns (q, r) such
// that a = q·b + r with deg(r) < deg(b).  Fails if b is zero or its
// leading coefficient is not invertible.
template <typename Coeff>
[[nodiscard]] Result<std::pair<std::vector<Coeff>, std::vector<Coeff>>>
poly_divmod(std::vector<Coeff> a, const std::vector<Coeff>& b) {
    using PairT = std::pair<std::vector<Coeff>, std::vector<Coeff>>;
    if (b.empty()) {
        return fail<PairT>(CASError{CASErrorKind::InvalidArgument,
            "polynomial divmod by zero", std::nullopt});
    }
    strip_trailing(a);
    if (a.size() < b.size()) return ok(PairT{std::vector<Coeff>{}, std::move(a)});

    auto lead_inv_res = CoeffOps<Coeff>::inverse(b.back());
    if (lead_inv_res.is_error()) {
        return fail<PairT>(lead_inv_res.error());
    }
    const Coeff lead_inv = lead_inv_res.value();
    const Coeff zero = CoeffOps<Coeff>::zero_like(b.front());

    std::vector<Coeff> q(a.size() - b.size() + 1U, zero);
    while (a.size() >= b.size() && !a.empty()) {
        const Coeff factor = a.back() * lead_inv;
        const std::size_t shift = a.size() - b.size();
        q[shift] = factor;
        for (std::size_t i = 0; i < b.size(); ++i) {
            a[shift + i] = a[shift + i] - factor * b[i];
        }
        strip_trailing(a);
    }
    return ok(PairT{std::move(q), std::move(a)});
}

// Extended Euclidean algorithm: returns (g, s, t) with s·a + t·b = g.
template <typename Coeff>
[[nodiscard]] Result<std::tuple<std::vector<Coeff>, std::vector<Coeff>, std::vector<Coeff>>>
poly_extended_gcd(std::vector<Coeff> a, std::vector<Coeff> b) {
    using TripleT = std::tuple<std::vector<Coeff>, std::vector<Coeff>, std::vector<Coeff>>;
    strip_trailing(a);
    strip_trailing(b);
    if (a.empty() && b.empty()) return ok(TripleT{{}, {}, {}});

    const Coeff& sample = a.empty() ? b.front() : a.front();
    const Coeff one = CoeffOps<Coeff>::one_like(sample);

    std::vector<Coeff> s0{one}, s1{};
    std::vector<Coeff> t0{}, t1{one};
    while (!b.empty()) {
        auto dm = poly_divmod(a, b);
        if (dm.is_error()) return fail<TripleT>(dm.error());
        auto [q, r] = std::move(dm.value());
        std::vector<Coeff> s2 = poly_sub(s0, poly_mul(q, s1));
        std::vector<Coeff> t2 = poly_sub(t0, poly_mul(q, t1));
        a = std::move(b);
        b = std::move(r);
        s0 = std::move(s1); s1 = std::move(s2);
        t0 = std::move(t1); t1 = std::move(t2);
    }
    return ok(TripleT{std::move(a), std::move(s0), std::move(t0)});
}

// Monic-normalise: divide every coefficient by the leading coefficient.
// Required to canonicalise the minimal polynomial at construction.
template <typename Coeff>
[[nodiscard]] Result<std::vector<Coeff>> poly_monic(std::vector<Coeff> p) {
    strip_trailing(p);
    if (p.empty()) return ok(p);
    auto inv_res = CoeffOps<Coeff>::inverse(p.back());
    if (inv_res.is_error()) {
        return fail<std::vector<Coeff>>(inv_res.error());
    }
    const Coeff inv = inv_res.value();
    for (auto& c : p) c = c * inv;
    return ok(p);
}

} // namespace tower_detail

// === Public template ========================================================

template <typename Coeff>
class AlgebraicElement {
public:
    using CoeffVec = std::vector<Coeff>;

    AlgebraicElement(CoeffVec value, CoeffVec min_poly)
        : value_(std::move(value)), min_poly_(std::move(min_poly)) {
        auto monic = tower_detail::poly_monic(std::move(min_poly_));
        // Construction-time invariant: minimal polynomial must be
        // monic-normalisable.  If not (e.g. zero leading coefficient
        // after stripping trailing zeros), leave it as-is and assume the
        // caller knows what they're doing — arithmetic will surface
        // errors via Result returns from `inverse()`.
        if (monic.is_ok()) min_poly_ = std::move(monic.value());
        // Canonicalise value mod min_poly.
        if (!min_poly_.empty()) {
            auto dm = tower_detail::poly_divmod(std::move(value_), min_poly_);
            if (dm.is_ok()) value_ = std::move(dm.value().second);
        }
        tower_detail::strip_trailing(value_);
    }

    [[nodiscard]] const CoeffVec& value() const noexcept { return value_; }
    [[nodiscard]] const CoeffVec& min_poly() const noexcept { return min_poly_; }
    [[nodiscard]] bool is_zero() const noexcept { return value_.empty(); }

    [[nodiscard]] AlgebraicElement operator-() const {
        return AlgebraicElement(tower_detail::poly_neg(value_), min_poly_);
    }
    [[nodiscard]] AlgebraicElement operator+(const AlgebraicElement& other) const {
        return AlgebraicElement(tower_detail::poly_add(value_, other.value_), min_poly_);
    }
    [[nodiscard]] AlgebraicElement operator-(const AlgebraicElement& other) const {
        return AlgebraicElement(tower_detail::poly_sub(value_, other.value_), min_poly_);
    }
    [[nodiscard]] AlgebraicElement operator*(const AlgebraicElement& other) const {
        return AlgebraicElement(tower_detail::poly_mul(value_, other.value_), min_poly_);
    }

    [[nodiscard]] bool operator==(const AlgebraicElement& other) const {
        return value_ == other.value_ && min_poly_ == other.min_poly_;
    }
    [[nodiscard]] bool operator!=(const AlgebraicElement& other) const {
        return !(*this == other);
    }

    [[nodiscard]] Result<AlgebraicElement> inverse() const {
        if (is_zero()) {
            return fail<AlgebraicElement>(CASError{CASErrorKind::InvalidArgument,
                "AlgebraicElement inverse: zero is not invertible", std::nullopt});
        }
        auto egcd = tower_detail::poly_extended_gcd(value_, min_poly_);
        if (egcd.is_error()) return fail<AlgebraicElement>(egcd.error());
        auto [g, s, t] = std::move(egcd.value());
        (void)t;
        // For the inverse to exist, gcd(value, min_poly) must be a
        // (non-zero) constant — i.e. of degree zero.  Otherwise the
        // minimal polynomial was reducible at the operand's witness or
        // the operand has a non-trivial factor in common with it.
        if (g.empty() || g.size() != 1U) {
            return fail<AlgebraicElement>(CASError{CASErrorKind::InvalidArgument,
                "AlgebraicElement inverse: gcd with min_poly is not a unit", std::nullopt});
        }
        auto g_inv_res = CoeffOps<Coeff>::inverse(g.front());
        if (g_inv_res.is_error()) return fail<AlgebraicElement>(g_inv_res.error());
        const Coeff g_inv = g_inv_res.value();
        for (auto& c : s) c = c * g_inv;
        return ok(AlgebraicElement(std::move(s), min_poly_));
    }

    [[nodiscard]] Result<AlgebraicElement> div(const AlgebraicElement& other) const {
        auto inv = other.inverse();
        if (inv.is_error()) return fail<AlgebraicElement>(inv.error());
        return ok((*this) * inv.value());
    }

    [[nodiscard]] Result<AlgebraicElement> pow(std::size_t exponent) const {
        const Coeff one_coeff = min_poly_.empty()
            ? CoeffOps<Coeff>::one_like(value_.empty() ? Coeff{} : value_.front())
            : CoeffOps<Coeff>::one_like(min_poly_.front());
        AlgebraicElement result(CoeffVec{one_coeff}, min_poly_);
        if (exponent == 0U) return ok(result);
        AlgebraicElement base = *this;
        std::size_t power = exponent;
        while (power > 0U) {
            if ((power & 1U) != 0U) result = result * base;
            power >>= 1U;
            if (power > 0U) base = base * base;
        }
        return ok(result);
    }

private:
    CoeffVec value_;
    CoeffVec min_poly_;
};

// === Recursive CoeffOps specialisation (definition after class) ============

template <typename C>
struct CoeffOps<AlgebraicElement<C>> {
    [[nodiscard]] static AlgebraicElement<C> zero_like(const AlgebraicElement<C>& sample) {
        return AlgebraicElement<C>({}, sample.min_poly());
    }
    [[nodiscard]] static AlgebraicElement<C> one_like(const AlgebraicElement<C>& sample) {
        if (sample.min_poly().empty()) {
            return AlgebraicElement<C>(typename AlgebraicElement<C>::CoeffVec{}, sample.min_poly());
        }
        return AlgebraicElement<C>(
            typename AlgebraicElement<C>::CoeffVec{CoeffOps<C>::one_like(sample.min_poly().front())},
            sample.min_poly());
    }
    [[nodiscard]] static bool is_zero(const AlgebraicElement<C>& c) { return c.is_zero(); }
    [[nodiscard]] static Result<AlgebraicElement<C>> inverse(const AlgebraicElement<C>& c) {
        return c.inverse();
    }
};

// Convenience aliases for the most common towers.
using AlgebraicElementQ      = AlgebraicElement<Rational>;
using AlgebraicTowerTwoLevel = AlgebraicElement<AlgebraicNumber>;

} // namespace algebra
} // namespace cas
