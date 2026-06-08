// F7.2-T3 — Discrete distributions: Binomial(n, p), Poisson(λ).
//
// Numerical strategy:
//   - Binomial PMF is computed via a stable log-Gamma reformulation:
//         log P(k) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
//                  + k·log(p) + (n-k)·log(1-p)
//     followed by exp().  This avoids overflow of binomial coefficients
//     for n on the order of 1000 that the naive product form cannot
//     handle.  Edge cases p = 0 or p = 1 short-circuit before the log to
//     avoid log(0).
//   - Poisson PMF uses the same log-Gamma idea:
//         log P(k) = k·log(λ) - λ - lgamma(k+1).
//
// CDFs accumulate PMF terms with a guard against catastrophic cancellation
// at the right tail: once the running sum saturates at 1.0 we stop early.

#include "cas/statistics.hpp"

#include <algorithm>
#include <cmath>

namespace cas::statistics {

namespace {

[[nodiscard]] CASError make_error(CASErrorKind kind, std::string message) {
    CASError err{};
    err.kind = kind;
    err.message = std::move(message);
    return err;
}

}  // namespace

Result<double> binomial_pmf(long long k, long long n, double p) {
    if (n < 0) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "binomial_pmf: n must be non-negative"));
    }
    if (!(p >= 0.0 && p <= 1.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "binomial_pmf: probability p must lie in [0, 1]"));
    }
    if (k < 0 || k > n) return ok(0.0);

    // Boundary probabilities handled explicitly to avoid 0·log(0).
    if (p == 0.0) return ok(k == 0 ? 1.0 : 0.0);
    if (p == 1.0) return ok(k == n ? 1.0 : 0.0);

    const double kd = static_cast<double>(k);
    const double nd = static_cast<double>(n);
    const double log_coeff = std::lgamma(nd + 1.0)
                           - std::lgamma(kd + 1.0)
                           - std::lgamma(nd - kd + 1.0);
    const double log_pmf = log_coeff
                         + kd * std::log(p)
                         + (nd - kd) * std::log(1.0 - p);
    return ok(std::exp(log_pmf));
}

Result<double> binomial_cdf(long long k, long long n, double p) {
    if (n < 0) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "binomial_cdf: n must be non-negative"));
    }
    if (!(p >= 0.0 && p <= 1.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "binomial_cdf: probability p must lie in [0, 1]"));
    }
    if (k < 0) return ok(0.0);
    if (k >= n) return ok(1.0);

    double acc = 0.0;
    for (long long i = 0; i <= k; ++i) {
        auto term = binomial_pmf(i, n, p);
        if (term.is_error()) return term;
        acc += term.value();
        if (acc >= 1.0) {  // numerical saturation guard
            acc = 1.0;
            break;
        }
    }
    return ok(std::min(acc, 1.0));
}

Result<double> poisson_pmf(long long k, double lambda) {
    if (k < 0) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "poisson_pmf: k must be non-negative"));
    }
    if (!(lambda >= 0.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "poisson_pmf: lambda must be non-negative"));
    }
    if (lambda == 0.0) return ok(k == 0 ? 1.0 : 0.0);
    const double kd = static_cast<double>(k);
    const double log_pmf = kd * std::log(lambda) - lambda - std::lgamma(kd + 1.0);
    return ok(std::exp(log_pmf));
}

Result<double> poisson_cdf(long long k, double lambda) {
    if (k < 0) return ok(0.0);
    if (!(lambda >= 0.0)) {
        return fail<double>(make_error(CASErrorKind::InvalidArgument,
            "poisson_cdf: lambda must be non-negative"));
    }
    if (lambda == 0.0) return ok(1.0);
    double acc = 0.0;
    for (long long i = 0; i <= k; ++i) {
        auto term = poisson_pmf(i, lambda);
        if (term.is_error()) return term;
        acc += term.value();
        if (acc >= 1.0) {
            acc = 1.0;
            break;
        }
    }
    return ok(std::min(acc, 1.0));
}

}  // namespace cas::statistics
