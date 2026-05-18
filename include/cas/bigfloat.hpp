#pragma once

#include <mpfr.h>
#include <string>
#include <cstdint>

// Forward declarations
namespace cas { class BigInt; struct Rational; }

namespace cas {

// Precision conversion: decimal significant digits → MPFR bits (with guard bits).
[[nodiscard]] inline mpfr_prec_t decimal_digits_to_bits(unsigned int digits) noexcept {
    return static_cast<mpfr_prec_t>(
        static_cast<double>(digits) * 3.3219280948873626 + 10.0);
}

// BigFloat — arbitrary-precision floating-point backed by MPFR.
// All operations use round-to-nearest (MPFR_RNDN).
// Precision is stored per-object; result precision = max(lhs, rhs) for binary ops.
class BigFloat {
public:
    // Minimum precision: 2 bits (MPFR requirement).
    static constexpr mpfr_prec_t MIN_PREC = 2;
    // Default: 128 bits ≈ 38 decimal digits.
    static constexpr mpfr_prec_t DEFAULT_PREC = 128;

    // ── Constructors / destructor ──────────────────────────────────────────────

    explicit BigFloat(mpfr_prec_t prec = DEFAULT_PREC);
    BigFloat(const BigFloat& other);
    BigFloat(BigFloat&& other) noexcept;
    ~BigFloat();

    BigFloat& operator=(const BigFloat& other);
    BigFloat& operator=(BigFloat&& other) noexcept;

    // ── Named constructors ─────────────────────────────────────────────────────

    [[nodiscard]] static BigFloat from_double(double d,
        mpfr_prec_t prec = DEFAULT_PREC);
    [[nodiscard]] static BigFloat from_integer_string(const std::string& s,
        mpfr_prec_t prec = DEFAULT_PREC);
    [[nodiscard]] static BigFloat from_rational_parts(
        const std::string& numerator_str,
        const std::string& denominator_str,
        mpfr_prec_t prec = DEFAULT_PREC);

    // Math constants at requested precision
    [[nodiscard]] static BigFloat pi(mpfr_prec_t prec = DEFAULT_PREC);
    [[nodiscard]] static BigFloat e(mpfr_prec_t prec = DEFAULT_PREC);
    [[nodiscard]] static BigFloat euler_gamma(mpfr_prec_t prec = DEFAULT_PREC);

    // ── Arithmetic ─────────────────────────────────────────────────────────────

    [[nodiscard]] BigFloat operator+(const BigFloat& rhs) const;
    [[nodiscard]] BigFloat operator-(const BigFloat& rhs) const;
    [[nodiscard]] BigFloat operator*(const BigFloat& rhs) const;
    [[nodiscard]] BigFloat operator/(const BigFloat& rhs) const;
    [[nodiscard]] BigFloat operator-() const;

    // ── Math functions ─────────────────────────────────────────────────────────

    [[nodiscard]] static BigFloat abs(const BigFloat& x);
    [[nodiscard]] static BigFloat sqrt(const BigFloat& x);
    [[nodiscard]] static BigFloat pow(const BigFloat& base, const BigFloat& exp);
    [[nodiscard]] static BigFloat exp(const BigFloat& x);
    [[nodiscard]] static BigFloat ln(const BigFloat& x);

    [[nodiscard]] static BigFloat sin(const BigFloat& x);
    [[nodiscard]] static BigFloat cos(const BigFloat& x);
    [[nodiscard]] static BigFloat tan(const BigFloat& x);
    [[nodiscard]] static BigFloat asin(const BigFloat& x);
    [[nodiscard]] static BigFloat acos(const BigFloat& x);
    [[nodiscard]] static BigFloat atan(const BigFloat& x);
    [[nodiscard]] static BigFloat atan2(const BigFloat& y, const BigFloat& x);

    [[nodiscard]] static BigFloat sinh(const BigFloat& x);
    [[nodiscard]] static BigFloat cosh(const BigFloat& x);
    [[nodiscard]] static BigFloat tanh(const BigFloat& x);

    [[nodiscard]] static BigFloat gamma(const BigFloat& x);
    [[nodiscard]] static BigFloat lgamma(const BigFloat& x);
    [[nodiscard]] static BigFloat erf(const BigFloat& x);

    // ── Comparisons ────────────────────────────────────────────────────────────

    [[nodiscard]] bool operator< (const BigFloat& rhs) const noexcept;
    [[nodiscard]] bool operator<=(const BigFloat& rhs) const noexcept;
    [[nodiscard]] bool operator> (const BigFloat& rhs) const noexcept;
    [[nodiscard]] bool operator>=(const BigFloat& rhs) const noexcept;
    [[nodiscard]] bool operator==(const BigFloat& rhs) const noexcept;
    [[nodiscard]] bool operator!=(const BigFloat& rhs) const noexcept;

    [[nodiscard]] bool is_zero() const noexcept;
    [[nodiscard]] bool is_nan() const noexcept;
    [[nodiscard]] bool is_inf() const noexcept;
    [[nodiscard]] bool is_negative() const noexcept;

    // ── Conversion ─────────────────────────────────────────────────────────────

    // decimal_digits significant decimal digits in g-format (auto fixed/scientific).
    [[nodiscard]] std::string to_string(int decimal_digits) const;
    [[nodiscard]] double to_double() const noexcept;
    [[nodiscard]] mpfr_prec_t precision_bits() const noexcept;

    // Raw access for low-level MPFR interop (e.g. passing to mpfr_* functions).
    [[nodiscard]] const mpfr_t& raw() const noexcept { return val_; }
    [[nodiscard]] mpfr_t& raw() noexcept { return val_; }

private:
    mpfr_t val_;
    mpfr_prec_t prec_;

    // Internal helper — make a result with precision = max(a.prec_, b.prec_)
    [[nodiscard]] static BigFloat make_result(mpfr_prec_t prec);
};

} // namespace cas
