#include "cas/bigfloat.hpp"
#include "cas/bigint.hpp"
#include "cas/rational.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace cas {

// ── Construction / destruction ────────────────────────────────────────────────

BigFloat::BigFloat(mpfr_prec_t prec) : prec_(std::max(prec, MIN_PREC)) {
    mpfr_init2(val_, prec_);
    mpfr_set_zero(val_, 1);
}

BigFloat::BigFloat(const BigFloat& other) : prec_(other.prec_) {
    mpfr_init2(val_, prec_);
    mpfr_set(val_, other.val_, MPFR_RNDN);
}

BigFloat::BigFloat(BigFloat&& other) noexcept : prec_(other.prec_) {
    // MPFR swap: move internal storage
    mpfr_init2(val_, prec_);
    mpfr_swap(val_, other.val_);
}

BigFloat::~BigFloat() {
    mpfr_clear(val_);
}

BigFloat& BigFloat::operator=(const BigFloat& other) {
    if (this == &other) return *this;
    if (prec_ != other.prec_) {
        mpfr_clear(val_);
        prec_ = other.prec_;
        mpfr_init2(val_, prec_);
    }
    mpfr_set(val_, other.val_, MPFR_RNDN);
    return *this;
}

BigFloat& BigFloat::operator=(BigFloat&& other) noexcept {
    if (this == &other) return *this;
    mpfr_swap(val_, other.val_);
    prec_ = other.prec_;
    return *this;
}

// ── Named constructors ────────────────────────────────────────────────────────

BigFloat BigFloat::from_double(double d, mpfr_prec_t prec) {
    BigFloat result(prec);
    mpfr_set_d(result.val_, d, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::from_integer_string(const std::string& s, mpfr_prec_t prec) {
    BigFloat result(prec);
    mpfr_set_str(result.val_, s.c_str(), 10, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::from_rational_parts(
    const std::string& num_str,
    const std::string& den_str,
    mpfr_prec_t prec) {
    BigFloat num(prec);
    BigFloat den(prec);
    mpfr_set_str(num.val_, num_str.c_str(), 10, MPFR_RNDN);
    mpfr_set_str(den.val_, den_str.c_str(), 10, MPFR_RNDN);
    BigFloat result(prec);
    mpfr_div(result.val_, num.val_, den.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::pi(mpfr_prec_t prec) {
    BigFloat result(prec);
    mpfr_const_pi(result.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::e(mpfr_prec_t prec) {
    BigFloat result(prec);
    BigFloat one = from_double(1.0, prec);
    mpfr_exp(result.val_, one.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::euler_gamma(mpfr_prec_t prec) {
    BigFloat result(prec);
    mpfr_const_euler(result.val_, MPFR_RNDN);
    return result;
}

// ── Helpers ───────────────────────────────────────────────────────────────────

BigFloat BigFloat::make_result(mpfr_prec_t prec) {
    return BigFloat(prec);
}

// ── Arithmetic ────────────────────────────────────────────────────────────────

BigFloat BigFloat::operator+(const BigFloat& rhs) const {
    BigFloat result(std::max(prec_, rhs.prec_));
    mpfr_add(result.val_, val_, rhs.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::operator-(const BigFloat& rhs) const {
    BigFloat result(std::max(prec_, rhs.prec_));
    mpfr_sub(result.val_, val_, rhs.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::operator*(const BigFloat& rhs) const {
    BigFloat result(std::max(prec_, rhs.prec_));
    mpfr_mul(result.val_, val_, rhs.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::operator/(const BigFloat& rhs) const {
    BigFloat result(std::max(prec_, rhs.prec_));
    mpfr_div(result.val_, val_, rhs.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::operator-() const {
    BigFloat result(prec_);
    mpfr_neg(result.val_, val_, MPFR_RNDN);
    return result;
}

// ── Math functions ────────────────────────────────────────────────────────────

BigFloat BigFloat::abs(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_abs(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::sqrt(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_sqrt(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::pow(const BigFloat& base, const BigFloat& exp) {
    BigFloat result(std::max(base.prec_, exp.prec_));
    mpfr_pow(result.val_, base.val_, exp.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::exp(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_exp(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::ln(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_log(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::sin(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_sin(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::cos(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_cos(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::tan(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_tan(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::asin(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_asin(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::acos(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_acos(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::atan(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_atan(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::atan2(const BigFloat& y, const BigFloat& x) {
    BigFloat result(std::max(y.prec_, x.prec_));
    mpfr_atan2(result.val_, y.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::sinh(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_sinh(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::cosh(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_cosh(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::tanh(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_tanh(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::gamma(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_gamma(result.val_, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::lgamma(const BigFloat& x) {
    BigFloat result(x.prec_);
    int sign = 0;
    mpfr_lgamma(result.val_, &sign, x.val_, MPFR_RNDN);
    return result;
}

BigFloat BigFloat::erf(const BigFloat& x) {
    BigFloat result(x.prec_);
    mpfr_erf(result.val_, x.val_, MPFR_RNDN);
    return result;
}

// ── Comparisons ───────────────────────────────────────────────────────────────

bool BigFloat::operator< (const BigFloat& rhs) const noexcept {
    return mpfr_less_p(val_, rhs.val_) != 0;
}
bool BigFloat::operator<=(const BigFloat& rhs) const noexcept {
    return mpfr_lessequal_p(val_, rhs.val_) != 0;
}
bool BigFloat::operator> (const BigFloat& rhs) const noexcept {
    return mpfr_greater_p(val_, rhs.val_) != 0;
}
bool BigFloat::operator>=(const BigFloat& rhs) const noexcept {
    return mpfr_greaterequal_p(val_, rhs.val_) != 0;
}
bool BigFloat::operator==(const BigFloat& rhs) const noexcept {
    return mpfr_equal_p(val_, rhs.val_) != 0;
}
bool BigFloat::operator!=(const BigFloat& rhs) const noexcept {
    return mpfr_equal_p(val_, rhs.val_) == 0;
}

bool BigFloat::is_zero()     const noexcept { return mpfr_zero_p(val_) != 0; }
bool BigFloat::is_nan()      const noexcept { return mpfr_nan_p(val_) != 0; }
bool BigFloat::is_inf()      const noexcept { return mpfr_inf_p(val_) != 0; }
bool BigFloat::is_negative() const noexcept { return mpfr_sgn(val_) < 0; }

// ── Conversion ────────────────────────────────────────────────────────────────

std::string BigFloat::to_string(int decimal_digits) const {
    if (decimal_digits <= 0) decimal_digits = 1;
    // Extra space: sign + digits + decimal point + exponent + null terminator
    const std::size_t buf_size = static_cast<std::size_t>(decimal_digits) + 100U;
    std::vector<char> buf(buf_size);
    mpfr_snprintf(buf.data(), buf_size, "%.*Rg",
        decimal_digits, val_);
    return std::string(buf.data());
}

double BigFloat::to_double() const noexcept {
    return mpfr_get_d(val_, MPFR_RNDN);
}

mpfr_prec_t BigFloat::precision_bits() const noexcept {
    return prec_;
}

} // namespace cas
