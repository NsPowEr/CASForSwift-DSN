#include "cas/algebraic_number_qx.hpp"
#include "polynomial_internal.hpp"
#include "algebra_internal.hpp"
#include "cas/error.hpp"
#include <algorithm>

namespace cas {
namespace algebra {

struct RatFun {
    RatPoly num;
    RatPoly den;

    RatFun() {
        num = RatPoly({Rational(BigInt(0))});
        den = RatPoly({Rational(BigInt(1))});
    }

    explicit RatFun(Rational c) {
        num = RatPoly({c});
        den = RatPoly({Rational(BigInt(1))});
        normalize();
    }

    explicit RatFun(RatPoly n) : num(std::move(n)) {
        den = RatPoly({Rational(BigInt(1))});
        normalize();
    }

    RatFun(RatPoly n, RatPoly d) : num(std::move(n)), den(std::move(d)) {
        normalize();
    }

    void normalize() {
        normalize_rational_coefficients(num);
        normalize_rational_coefficients(den);
        if (den.is_zero()) {
            den = RatPoly({Rational(BigInt(1))});
            return;
        }
        if (num.is_zero()) {
            den = RatPoly({Rational(BigInt(1))});
            return;
        }
        
        auto [g, s, t] = extended_gcd_rational_poly(num, den);
        if (!g.is_zero() && g.degree() > 0) {
            auto [qn, rn] = div_rem_rational_poly(num, g);
            auto [qd, rd] = div_rem_rational_poly(den, g);
            num = std::move(qn);
            den = std::move(qd);
        } else if (!g.is_zero() && g.degree() == 0) {
            Rational g_const = g.constant_term();
            if (!g_const.numerator().is_zero()) {
                for (auto& c : num.coefficients()) c = c / g_const;
                for (auto& c : den.coefficients()) c = c / g_const;
            }
        }
        
        if (!den.is_zero()) {
            Rational lc = den.leading_coeff();
            if (!lc.numerator().is_zero() && lc != Rational(BigInt(1))) {
                for (auto& c : num.coefficients()) c = c / lc;
                for (auto& c : den.coefficients()) c = c / lc;
            }
        }
        normalize_rational_coefficients(num);
        normalize_rational_coefficients(den);
    }

    [[nodiscard]] bool is_zero() const noexcept {
        return num.is_zero();
    }

    [[nodiscard]] RatFun operator-() const {
        RatPoly n = num;
        for (auto& c : n.coefficients()) c = -c;
        return RatFun(std::move(n), den);
    }

    [[nodiscard]] RatFun operator+(const RatFun& other) const {
        RatPoly n = add_rational_poly(
            mul_rational_poly(num, other.den),
            mul_rational_poly(other.num, den)
        );
        RatPoly d = mul_rational_poly(den, other.den);
        return RatFun(std::move(n), std::move(d));
    }

    [[nodiscard]] RatFun operator-(const RatFun& other) const {
        RatPoly n = sub_rational_poly(
            mul_rational_poly(num, other.den),
            mul_rational_poly(other.num, den)
        );
        RatPoly d = mul_rational_poly(den, other.den);
        return RatFun(std::move(n), std::move(d));
    }

    [[nodiscard]] RatFun operator*(const RatFun& other) const {
        return RatFun(
            mul_rational_poly(num, other.num),
            mul_rational_poly(den, other.den)
        );
    }

    [[nodiscard]] RatFun inverse() const {
        return RatFun(den, num);
    }

    [[nodiscard]] RatFun div(const RatFun& other) const {
        return (*this) * other.inverse();
    }

    [[nodiscard]] bool operator==(const RatFun& other) const {
        return num.coefficients() == other.num.coefficients() &&
               den.coefficients() == other.den.coefficients();
    }
};

struct AlgebraicNumberQx::Impl {
    RatFun a;
    RatFun b;
    RatFun r;

    Impl(RatFun a_val, RatFun b_val, RatFun r_val)
        : a(std::move(a_val)), b(std::move(b_val)), r(std::move(r_val)) {}
};

AlgebraicNumberQx::AlgebraicNumberQx(CoeffVec a_num, CoeffVec a_den,
                                     CoeffVec b_num, CoeffVec b_den,
                                     CoeffVec r_num, CoeffVec r_den) {
    RatFun a_val(RatPoly(std::move(a_num)), RatPoly(std::move(a_den)));
    RatFun b_val(RatPoly(std::move(b_num)), RatPoly(std::move(b_den)));
    RatFun r_val(RatPoly(std::move(r_num)), RatPoly(std::move(r_den)));
    impl_ = std::make_shared<Impl>(std::move(a_val), std::move(b_val), std::move(r_val));
}

AlgebraicNumberQx::AlgebraicNumberQx(CoeffVec a, CoeffVec b, CoeffVec r) {
    RatFun a_val(RatPoly(std::move(a)));
    RatFun b_val(RatPoly(std::move(b)));
    RatFun r_val(RatPoly(std::move(r)));
    impl_ = std::make_shared<Impl>(std::move(a_val), std::move(b_val), std::move(r_val));
}

AlgebraicNumberQx::CoeffVec AlgebraicNumberQx::a_num() const { return impl_->a.num.coefficients(); }
AlgebraicNumberQx::CoeffVec AlgebraicNumberQx::a_den() const { return impl_->a.den.coefficients(); }
AlgebraicNumberQx::CoeffVec AlgebraicNumberQx::b_num() const { return impl_->b.num.coefficients(); }
AlgebraicNumberQx::CoeffVec AlgebraicNumberQx::b_den() const { return impl_->b.den.coefficients(); }
AlgebraicNumberQx::CoeffVec AlgebraicNumberQx::r_num() const { return impl_->r.num.coefficients(); }
AlgebraicNumberQx::CoeffVec AlgebraicNumberQx::r_den() const { return impl_->r.den.coefficients(); }

bool AlgebraicNumberQx::is_zero() const noexcept {
    return impl_->a.is_zero() && impl_->b.is_zero();
}

AlgebraicNumberQx AlgebraicNumberQx::operator-() const {
    return AlgebraicNumberQx(
        (-impl_->a).num.coefficients(), (-impl_->a).den.coefficients(),
        (-impl_->b).num.coefficients(), (-impl_->b).den.coefficients(),
        impl_->r.num.coefficients(), impl_->r.den.coefficients()
    );
}

AlgebraicNumberQx AlgebraicNumberQx::operator+(const AlgebraicNumberQx& other) const {
    RatFun sum_a = impl_->a + other.impl_->a;
    RatFun sum_b = impl_->b + other.impl_->b;
    return AlgebraicNumberQx(
        sum_a.num.coefficients(), sum_a.den.coefficients(),
        sum_b.num.coefficients(), sum_b.den.coefficients(),
        impl_->r.num.coefficients(), impl_->r.den.coefficients()
    );
}

AlgebraicNumberQx AlgebraicNumberQx::operator-(const AlgebraicNumberQx& other) const {
    RatFun diff_a = impl_->a - other.impl_->a;
    RatFun diff_b = impl_->b - other.impl_->b;
    return AlgebraicNumberQx(
        diff_a.num.coefficients(), diff_a.den.coefficients(),
        diff_b.num.coefficients(), diff_b.den.coefficients(),
        impl_->r.num.coefficients(), impl_->r.den.coefficients()
    );
}

AlgebraicNumberQx AlgebraicNumberQx::operator*(const AlgebraicNumberQx& other) const {
    RatFun prod_a = impl_->a * other.impl_->a + impl_->b * other.impl_->b * impl_->r;
    RatFun prod_b = impl_->a * other.impl_->b + impl_->b * other.impl_->a;
    return AlgebraicNumberQx(
        prod_a.num.coefficients(), prod_a.den.coefficients(),
        prod_b.num.coefficients(), prod_b.den.coefficients(),
        impl_->r.num.coefficients(), impl_->r.den.coefficients()
    );
}

bool AlgebraicNumberQx::operator==(const AlgebraicNumberQx& other) const {
    return impl_->a == other.impl_->a && impl_->b == other.impl_->b && impl_->r == other.impl_->r;
}

Result<AlgebraicNumberQx> AlgebraicNumberQx::inverse() const {
    RatFun denom = impl_->a * impl_->a - impl_->b * impl_->b * impl_->r;
    if (denom.is_zero()) {
        return fail<AlgebraicNumberQx>(make_error(CASErrorKind::InvalidArgument, "Division by zero in AlgebraicNumberQx"));
    }
    RatFun inv_a = impl_->a.div(denom);
    RatFun inv_b = (-impl_->b).div(denom);
    return ok(AlgebraicNumberQx(
        inv_a.num.coefficients(), inv_a.den.coefficients(),
        inv_b.num.coefficients(), inv_b.den.coefficients(),
        impl_->r.num.coefficients(), impl_->r.den.coefficients()
    ));
}

Result<AlgebraicNumberQx> AlgebraicNumberQx::div(const AlgebraicNumberQx& other) const {
    auto inv = other.inverse();
    if (inv.is_error()) return fail<AlgebraicNumberQx>(inv.error());
    return ok((*this) * inv.value());
}

Result<AlgebraicNumberQx> AlgebraicNumberQx::pow(std::size_t exponent) const {
    AlgebraicNumberQx result(
        {Rational(1)}, {Rational(1)},
        {Rational(0)}, {Rational(1)},
        impl_->r.num.coefficients(), impl_->r.den.coefficients()
    );
    if (exponent == 0) return ok(result);
    AlgebraicNumberQx base = *this;
    std::size_t power = exponent;
    while (power > 0) {
        if ((power & 1) != 0) {
            result = result * base;
        }
        power >>= 1;
        if (power > 0) {
            base = base * base;
        }
    }
    return ok(result);
}

} // namespace algebra
} // namespace cas
