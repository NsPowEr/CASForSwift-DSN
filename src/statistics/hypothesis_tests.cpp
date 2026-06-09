// F7.2-B2 — Hypothesis tests (z, t one/two-sample, chi-squared GoF, F variance).

#include "cas/statistics.hpp"
#include "cas/error_helpers.hpp"

#include <cmath>
#include <numeric>

namespace cas::statistics {

namespace {

CASError make_error(CASErrorKind kind, std::string message) {
    CASError err{};
    err.kind = kind;
    err.message = std::move(message);
    return err;
}

double sample_mean(const std::vector<double>& s) {
    return std::accumulate(s.begin(), s.end(), 0.0) / static_cast<double>(s.size());
}

double sample_var(const std::vector<double>& s, double mean) {
    double acc = 0.0;
    for (double v : s) {
        const double d = v - mean;
        acc += d * d;
    }
    return acc / static_cast<double>(s.size() - 1);
}

}  // namespace

Result<HypothesisTestResult> z_test_one_sample(
    const std::vector<double>& sample, double mu0, double sigma) {
    if (sample.size() < 2U) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument, "z_test_one_sample: need n ≥ 2"));
    }
    if (!(sigma > 0.0)) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument, "z_test_one_sample: sigma must be > 0"));
    }
    const double n = static_cast<double>(sample.size());
    const double mean = sample_mean(sample);
    const double z = (mean - mu0) / (sigma / std::sqrt(n));
    // Two-sided p-value via Normal cdf.
    // p = 2 · (1 − Φ(|z|))
    auto cdf_abs = normal_quantile(0.5, 0.0, 1.0);  // sanity touch (unused)
    (void)cdf_abs;
    // Reuse the standard normal CDF approximation via erf:
    // Φ(z) = 0.5 · (1 + erf(z/√2))
    const double phi = 0.5 * (1.0 + std::erf(std::abs(z) / std::sqrt(2.0)));
    const double p = 2.0 * (1.0 - phi);
    return ok(HypothesisTestResult{
        .statistic = z,
        .p_value = p,
        .df = 0.0,
        .reject_at_alpha_005 = p < 0.05,
    });
}

Result<HypothesisTestResult> t_test_one_sample(
    const std::vector<double>& sample, double mu0) {
    if (sample.size() < 2U) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument, "t_test_one_sample: need n ≥ 2"));
    }
    const double n = static_cast<double>(sample.size());
    const double mean = sample_mean(sample);
    const double var = sample_var(sample, mean);
    if (!(var > 0.0)) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument, "t_test_one_sample: zero variance"));
    }
    const double t = (mean - mu0) / std::sqrt(var / n);
    const double df = n - 1.0;
    auto cdf_res = student_t_cdf(std::abs(t), df);
    if (cdf_res.is_error()) {
        return fail<HypothesisTestResult>(cdf_res.error());
    }
    const double p = 2.0 * (1.0 - cdf_res.value());
    return ok(HypothesisTestResult{
        .statistic = t,
        .p_value = p,
        .df = df,
        .reject_at_alpha_005 = p < 0.05,
    });
}

Result<HypothesisTestResult> t_test_two_sample(
    const std::vector<double>& s1, const std::vector<double>& s2) {
    if (s1.size() < 2U || s2.size() < 2U) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument, "t_test_two_sample: need n ≥ 2 in each"));
    }
    const double n1 = static_cast<double>(s1.size());
    const double n2 = static_cast<double>(s2.size());
    const double m1 = sample_mean(s1);
    const double m2 = sample_mean(s2);
    const double v1 = sample_var(s1, m1);
    const double v2 = sample_var(s2, m2);
    if (!(v1 > 0.0) || !(v2 > 0.0)) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument, "t_test_two_sample: zero variance"));
    }
    const double se = std::sqrt(v1 / n1 + v2 / n2);
    const double t = (m1 - m2) / se;
    // Welch-Satterthwaite df.
    const double a = v1 / n1;
    const double b = v2 / n2;
    const double df = (a + b) * (a + b) /
                      (a * a / (n1 - 1.0) + b * b / (n2 - 1.0));
    auto cdf_res = student_t_cdf(std::abs(t), df);
    if (cdf_res.is_error()) {
        return fail<HypothesisTestResult>(cdf_res.error());
    }
    const double p = 2.0 * (1.0 - cdf_res.value());
    return ok(HypothesisTestResult{
        .statistic = t,
        .p_value = p,
        .df = df,
        .reject_at_alpha_005 = p < 0.05,
    });
}

Result<HypothesisTestResult> chi_squared_goodness_of_fit(
    const std::vector<double>& observed,
    const std::vector<double>& expected) {
    if (observed.empty() || observed.size() != expected.size()) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument,
            "chi_squared_goodness_of_fit: size mismatch or empty"));
    }
    double chi2 = 0.0;
    for (std::size_t i = 0; i < observed.size(); ++i) {
        if (!(expected[i] > 0.0)) {
            return fail<HypothesisTestResult>(make_error(
                CASErrorKind::InvalidArgument,
                "chi_squared_goodness_of_fit: expected count must be > 0"));
        }
        const double d = observed[i] - expected[i];
        chi2 += (d * d) / expected[i];
    }
    const double df = static_cast<double>(observed.size() - 1U);
    auto cdf_res = chi_squared_cdf(chi2, df);
    if (cdf_res.is_error()) {
        return fail<HypothesisTestResult>(cdf_res.error());
    }
    const double p = 1.0 - cdf_res.value();
    return ok(HypothesisTestResult{
        .statistic = chi2,
        .p_value = p,
        .df = df,
        .reject_at_alpha_005 = p < 0.05,
    });
}

Result<HypothesisTestResult> f_test_variance(
    const std::vector<double>& s1, const std::vector<double>& s2) {
    if (s1.size() < 2U || s2.size() < 2U) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument, "f_test_variance: need n ≥ 2 in each"));
    }
    const double m1 = sample_mean(s1);
    const double m2 = sample_mean(s2);
    const double v1 = sample_var(s1, m1);
    const double v2 = sample_var(s2, m2);
    if (!(v2 > 0.0)) {
        return fail<HypothesisTestResult>(make_error(
            CASErrorKind::InvalidArgument, "f_test_variance: denominator variance zero"));
    }
    const double F = v1 / v2;
    const double d1 = static_cast<double>(s1.size() - 1U);
    const double d2 = static_cast<double>(s2.size() - 1U);
    auto cdf_res = f_cdf(F, d1, d2);
    if (cdf_res.is_error()) {
        return fail<HypothesisTestResult>(cdf_res.error());
    }
    // Two-sided p-value for variance ratio.
    const double tail = cdf_res.value();
    const double p = 2.0 * std::min(tail, 1.0 - tail);
    return ok(HypothesisTestResult{
        .statistic = F,
        .p_value = p,
        .df = d1,  // store numerator df; d2 implied by caller's s2
        .reject_at_alpha_005 = p < 0.05,
    });
}

}  // namespace cas::statistics
