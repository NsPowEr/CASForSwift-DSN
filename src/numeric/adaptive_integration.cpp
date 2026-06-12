// adaptive_integration.cpp — F6.D adaptive G7/K15 Gauss-Kronrod integrator.
//
// Implements `integrate_gauss_kronrod_adaptive` with priority-queue panel
// splitting. Convergence and resource caps are driven exclusively by
// CASContext fields (no hardcoded tolerances or depth caps).
//
// References:
//   - R. Piessens, E. de Doncker-Kapenga, C. W. Überhuber, D. K. Kahaner,
//     "QUADPACK: A Subroutine Package for Automatic Integration", Springer,
//     1983, §4 (Algorithm DQAG / GAUSS_KRONROD).
//   - Brent-Zimmermann "Modern Computer Arithmetic" §4.1.3 (adaptive panels).
//
// Algorithm (Adaptive_Numerical_Integration.md §4):
//   1. Compute (I_K, I_G, err) on [a, b] via G7/K15.
//   2. Push panel onto max-heap keyed by err.
//   3. While heap.size() < max_intervals:
//        - if total_error ≤ max(abs_tol, rel_tol·|total_integral|): return success.
//        - pop worst panel; bisect at midpoint; compute G7/K15 on each half;
//          update total_integral and total_error by the difference; push the
//          two new panels.
//   4. On reaching max_intervals without convergence: return success=false
//      with the current best estimate (no silent truncation).
//
// G7/K15 nodes and weights are constexpr (single 18-digit table, QUADPACK §4).
// They are computed analytically as roots of Legendre and Kronrod polynomials.

#include "cas/numeric.hpp"
#include "cas/symbolic.hpp"

#include <cmath>
#include <cstddef>
#include <queue>
#include <vector>

namespace cas::numeric {

namespace {

// QUADPACK §4 Table 1 — Gauss-Kronrod 15-point nodes (symmetric in [-1,1]).
constexpr double kGK15_Nodes[15] = {
    -0.991455371120812639206854697526329,
    -0.949107912342758524526189684047851,
    -0.864864423359769072789712788640926,
    -0.741531185599394439863864773280788,
    -0.586087235467691130294144838258730,
    -0.405845151377397166906606412076961,
    -0.207784955007898467600689403773245,
     0.0,
     0.207784955007898467600689403773245,
     0.405845151377397166906606412076961,
     0.586087235467691130294144838258730,
     0.741531185599394439863864773280788,
     0.864864423359769072789712788640926,
     0.949107912342758524526189684047851,
     0.991455371120812639206854697526329};

// QUADPACK §4 Table 1 — Kronrod weights for the 15-point rule.
constexpr double kGK15_Weights[15] = {
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714,
    0.204432940075298892414161999234649,
    0.190350578064785409913256402421014,
    0.169004726639267902826583426598550,
    0.140653259715525918745189590510238,
    0.104790010322250183839876322541518,
    0.063092092629978553290700663189204,
    0.022935322010529224963732008058970};

// 7-point Gauss-Legendre weights paired with the 7 odd-indexed Kronrod nodes
// (positions 1, 3, 5, 7, 9, 11, 13 inside the 15-point array).
constexpr double kG7_Weights[7] = {
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327,
    0.381830050505118944950369775488975,
    0.279705391489276667901467771423780,
    0.129484966168869693270611432679082};

constexpr int kG7_NodeMap[7] = {1, 3, 5, 7, 9, 11, 13};

struct PanelEstimate {
    double i_kronrod{0.0};
    double i_gauss{0.0};
    double err{0.0};
};

Result<PanelEstimate> gk15_panel(
    const std::function<Result<double>(double)>& f,
    double a, double b) {
    const double half = 0.5 * (b - a);
    const double mid  = 0.5 * (a + b);
    PanelEstimate est;
    double f_vals[15];
    for (int i = 0; i < 15; ++i) {
        const double x = mid + half * kGK15_Nodes[i];
        auto fv = f(x);
        if (fv.is_error()) return fail<PanelEstimate>(fv.error());
        f_vals[i] = fv.value();
        est.i_kronrod += kGK15_Weights[i] * f_vals[i];
    }
    for (int i = 0; i < 7; ++i) {
        est.i_gauss += kG7_Weights[i] * f_vals[kG7_NodeMap[i]];
    }
    est.i_kronrod *= half;
    est.i_gauss   *= half;
    // |G - K| is a conservative bound on the local truncation error.
    est.err = std::abs(est.i_kronrod - est.i_gauss);
    return ok(est);
}

struct Panel {
    double a{0.0};
    double b{0.0};
    double i_kronrod{0.0};
    double err{0.0};
};

struct PanelGreater {
    bool operator()(const Panel& lhs, const Panel& rhs) const noexcept {
        // std::priority_queue is a max-heap by default; sort by err.
        return lhs.err < rhs.err;
    }
};

}  // namespace

Result<GKIntegrationResult> integrate_gauss_kronrod_adaptive(
    const std::function<Result<double>(double)>& f,
    double a,
    double b,
    const symbolic::CASContext& ctx) {

    const double abs_tol = ctx.integration_abs_tol();
    const double rel_tol = ctx.integration_rel_tol();
    const std::size_t max_intervals = ctx.integration_max_intervals();

    GKIntegrationResult out;
    if (!(a < b) && !(a > b)) {
        // Degenerate interval (a == b): integral is zero.
        out.success = true;
        out.interval_count = 1U;
        return ok(out);
    }
    if (a > b) {
        // Integral from b to a is the negative.
        auto flipped = integrate_gauss_kronrod_adaptive(f, b, a, ctx);
        if (flipped.is_error()) return flipped;
        out = flipped.value();
        out.value = -out.value;
        return ok(out);
    }

    auto initial = gk15_panel(f, a, b);
    if (initial.is_error()) return fail<GKIntegrationResult>(initial.error());

    double total_integral = initial.value().i_kronrod;
    double total_error    = initial.value().err;

    std::priority_queue<Panel, std::vector<Panel>, PanelGreater> heap;
    heap.push(Panel{a, b, initial.value().i_kronrod, initial.value().err});

    while (heap.size() < max_intervals) {
        const double abs_int = std::abs(total_integral);
        const double effective_tol = std::max(abs_tol, rel_tol * abs_int);
        if (total_error <= effective_tol) {
            out.value = total_integral;
            out.error_estimate = total_error;
            out.interval_count = heap.size();
            out.success = true;
            return ok(out);
        }

        const Panel worst = heap.top();
        heap.pop();
        const double m = 0.5 * (worst.a + worst.b);

        auto left  = gk15_panel(f, worst.a, m);
        if (left.is_error())  return fail<GKIntegrationResult>(left.error());
        auto right = gk15_panel(f, m, worst.b);
        if (right.is_error()) return fail<GKIntegrationResult>(right.error());

        // Incremental update: replace the popped panel's contribution.
        total_integral += (left.value().i_kronrod + right.value().i_kronrod) - worst.i_kronrod;
        total_error    += (left.value().err + right.value().err) - worst.err;

        heap.push(Panel{worst.a, m, left.value().i_kronrod, left.value().err});
        heap.push(Panel{m, worst.b, right.value().i_kronrod, right.value().err});
    }

    // Resource cap reached without meeting tolerance: report failure but
    // expose the best current estimate so callers can decide.
    out.value = total_integral;
    out.error_estimate = total_error;
    out.interval_count = heap.size();
    out.success = false;
    return ok(out);
}

}  // namespace cas::numeric
