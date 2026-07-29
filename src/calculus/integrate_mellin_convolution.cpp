// A7 §6.7 — Mellin-convolution definite integrator (Meijer_G_Slater.md §6.7).
//
// Closes  ∫_0^∞ K · x^μ · f_1(x) · f_2(x) dx  when each f_i converts (via
// to_meijerg, §5/§3.1) to a Meijer G LINEAR in x, G_i(η_i · x), by the
// Adamchik–Marichev / Meijer (1941) convolution
//
//   ∫_0^∞ G^{m,n}_{p,q}(η x | a;b) · G^{μ',ν'}_{σ,τ}(ω x | c;d) dx
//        = (1/η) · G^{n+μ', m+ν'}_{q+σ, p+τ}(ω/η | ...composition...),
//
// with the parameter composition of §6.7 (verified numerically in the spec,
// mpmath; cross-checked here on the two spec worked examples in the tests).
//
// SOUNDNESS GATE (the unique new risk of §6.7 is applying the exact formula to
// a DIVERGENT integral, which would yield a finite WRONG value):  the formula
// is committed ONLY when convergence at both endpoints is CERTIFIED from the
// G asymptotics (rational parameters ⇒ decidable), else the pattern skips
// (nullopt) and the ∫_0^∞ falls through to the generic Unimplemented — never a
// wrong result (REGOLA ZERO).  The spec's numeric shadow-eval lives in the
// TEST layer (numeric-certify / mpmath): the evaluator is architecturally
// unavailable inside cas_calculus (cas_numeric depends on cas_calculus), so
// the in-engine gate is the convergence certificate; the composition itself
// is exact (every to_meijerg / §6.7 / from_meijerg step is spec-verified).

#include "integrate_definite_patterns.hpp"

#include "cas/extended_real.hpp"
#include "cas/meijerg.hpp"
#include "cas/rational.hpp"

#include <optional>
#include <vector>

namespace cas::calculus {

namespace {

using cas::is_pos_infinity;
using symbolic::CASContext;
using symbolic::MeijerGView;

[[nodiscard]] bool depends_on_var(ExprPtr e, const Symbol& var);

// Real rational value of a literal (Integer/Rational). Anything else (symbol,
// irrational, complex) → nullopt: the convergence certificate can only be
// decided on rational parameters (spec §6.7 "quando i parametri sono
// razionali").
[[nodiscard]] std::optional<Rational> rational_value(ExprPtr e) {
    if (const auto* il = expr_cast<IntegerLit>(e)) return Rational(il->value);
    if (const auto* rl = expr_cast<RationalLit>(e))
        return Rational(rl->numerator, rl->denominator);
    return std::nullopt;
}

[[nodiscard]] bool depends_on_var(ExprPtr e, const Symbol& var) {
    if (e == nullptr) return false;
    if (const auto* s = expr_cast<Symbol>(e)) return s->name == var.name;
    switch (e->kind) {
        case ExprKind::Unary: return depends_on_var(expr_ref<Unary>(e).operand, var);
        case ExprKind::Binary: {
            const auto& b = expr_ref<Binary>(e);
            return depends_on_var(b.left, var) || depends_on_var(b.right, var);
        }
        case ExprKind::Sum:
            for (ExprPtr t : expr_ref<Sum>(e).terms)
                if (depends_on_var(t, var)) return true;
            return false;
        case ExprKind::Product:
            for (ExprPtr f : expr_ref<Product>(e).factors)
                if (depends_on_var(f, var)) return true;
            return false;
        case ExprKind::FuncCall:
            for (ExprPtr a : expr_ref<FuncCall>(e).args)
                if (depends_on_var(a, var)) return true;
            return false;
        default: return false;
    }
}

// z == η · var, with η a var-free expression and the var appearing to the
// first power exactly.  Returns η (nullptr means η = 1).  Conservative: any
// other shape (η·var^r with r≠1, sums, …) → nullopt (not linear ⇒ §6.7 n/a).
[[nodiscard]] std::optional<ExprPtr> match_scaled_var(ExprPtr z, const Symbol& var) {
    if (z == nullptr) return std::nullopt;
    if (const auto* s = expr_cast<Symbol>(z))
        return s->name == var.name ? std::optional<ExprPtr>(nullptr) : std::nullopt;
    // Binary Mul: exactly one side is var, the other var-free.
    if (const auto* bin = expr_cast<Binary>(z); bin != nullptr && bin->op == BinaryOp::Mul) {
        const bool l_is_var = expr_cast<Symbol>(bin->left) &&
            expr_cast<Symbol>(bin->left)->name == var.name;
        const bool r_is_var = expr_cast<Symbol>(bin->right) &&
            expr_cast<Symbol>(bin->right)->name == var.name;
        if (l_is_var && !depends_on_var(bin->right, var)) return std::optional<ExprPtr>(bin->right);
        if (r_is_var && !depends_on_var(bin->left, var)) return std::optional<ExprPtr>(bin->left);
        return std::nullopt;
    }
    // Product: exactly one factor is var (power 1), all others var-free.
    if (const auto* prod = expr_cast<Product>(z)) {
        std::vector<ExprPtr> scale_factors;
        int var_count = 0;
        for (ExprPtr f : prod->factors) {
            if (const auto* s = expr_cast<Symbol>(f); s && s->name == var.name) {
                ++var_count;
            } else if (!depends_on_var(f, var)) {
                scale_factors.push_back(f);
            } else {
                return std::nullopt;  // a factor depends on var but isn't the bare var
            }
        }
        if (var_count != 1) return std::nullopt;
        if (scale_factors.empty()) return std::optional<ExprPtr>(nullptr);
        if (scale_factors.size() == 1U) return std::optional<ExprPtr>(scale_factors.front());
        // Multi-factor η is rare for our linear cores; skip rather than risk a
        // silent η=1 (the caller conservatively bails ⇒ generic Unimplemented).
        return std::nullopt;
    }
    return std::nullopt;
}

struct LinearG {
    MeijerGView view;
    ExprPtr scale;                 // η (var-free), never nullptr (1 materialised)
    std::vector<ExprPtr> prefactor;  // var-free K factors from to_meijerg
};

// core → K · G(η·x) with G linear in x.  nullopt when the shape is outside the
// §6.7-applicable family (caller then skips the pattern).
[[nodiscard]] std::optional<LinearG> extract_linear_g(
    CASContext& ctx, ExprPtr core, const Symbol& var) {
    AstArena& arena = ctx.arena();
    auto g_form = symbolic::to_meijerg(ctx, core);
    if (g_form.is_error()) return std::nullopt;

    const FuncCall* g_node = nullptr;
    std::vector<ExprPtr> prefactor;
    auto absorb = [&](ExprPtr f) -> bool {
        if (const auto* call = expr_cast<FuncCall>(f);
            call != nullptr && call->func_id == BuiltinOp::MeijerG) {
            if (g_node != nullptr) return false;  // two G nodes: not our shape
            g_node = call;
            return true;
        }
        if (depends_on_var(f, var)) return false;  // non-G var-dependent prefactor
        prefactor.push_back(f);
        return true;
    };
    if (const auto* prod = expr_cast<Product>(g_form.value())) {
        for (ExprPtr f : prod->factors)
            if (!absorb(f)) return std::nullopt;
    } else if (!absorb(g_form.value())) {
        return std::nullopt;
    }
    if (g_node == nullptr) return std::nullopt;

    auto view = symbolic::view_meijerg(*g_node);
    if (view.is_error()) return std::nullopt;

    auto scale = match_scaled_var(view.value().z, var);
    if (!scale.has_value()) return std::nullopt;
    ExprPtr eta = scale.value() ? scale.value()
                                : static_cast<ExprPtr>(arena.make<IntegerLit>(BigInt(1)));
    return LinearG{std::move(view.value()), eta, std::move(prefactor)};
}

[[nodiscard]] ExprPtr neg_param(CASContext& ctx, ExprPtr e) {
    AstArena& arena = ctx.arena();
    if (const auto* il = expr_cast<IntegerLit>(e)) return arena.make<IntegerLit>(-il->value);
    if (const auto* rl = expr_cast<RationalLit>(e))
        return arena.make<RationalLit>(-rl->numerator, rl->denominator);
    ExprPtr n = arena.make<Unary>(UnaryOp::Neg, e);
    auto s = ctx.simplify(n);
    return s.is_ok() ? s.value() : n;
}

// 0-endpoint leading exponent of G near x→0⁺:  G ~ x^{min_{k≤m} Re(b_k)}
// (generic b_1..b_m, m≥1).  nullopt when m=0 (different 0-behaviour) or any
// b_k (k≤m) is non-rational ⇒ certificate undecidable here.
[[nodiscard]] std::optional<Rational> zero_exponent(const MeijerGView& v) {
    if (v.m == 0) return std::nullopt;
    std::optional<Rational> lo;
    for (std::size_t k = 0; k < v.m; ++k) {
        auto r = rational_value(v.b[k]);
        if (!r.has_value()) return std::nullopt;
        if (!lo.has_value() || *r < *lo) lo = r;
    }
    return lo;
}

// Sufficient ∞-convergence certificate for the PRODUCT G1·G2 as x→∞.
//
// CRITICAL SOUNDNESS POINT (do not weaken to "p_i < q_i"):  a Meijer G decays
// like a genuine exponential  exp(−η x)  on the positive real axis ONLY when
// q_i − p_i = 1.  For q_i − p_i ≥ 2 the large-x behaviour is
// ALGEBRAIC-OSCILLATORY (e.g. J_0(2√(ωx)) = G^{1,0}_{0,2}(ωx) ~ (ωx)^{−1/4}
// cos(…)), which is NOT absolutely integrable on its own — a product of two
// such factors (J_0·J_0 ~ x^{−1/2}·osc) DIVERGES, and Meijer's exact formula
// would then return a finite WRONG value.  For q_i < p_i the factor GROWS
// exponentially.
//
// So we certify ∞-convergence only when at least one factor is a true Laplace
// kernel (q_i = p_i + 1 ⇒ exp(−η_i x) decay, with η_i > 0 checked separately)
// and the OTHER factor does not grow exponentially (p_j ≤ q_j): exp(−η_i x)
// dominates any algebraic/oscillatory or slower-exponential behaviour, giving
// ABSOLUTE convergence.  The pure-algebraic (p=q) case and the
// oscillatory-only (q−p ≥ 2 with no exp partner) case are left to a follow-up
// (needs the full Meijer δ/c* inequalities) — they return false ⇒ the pattern
// skips, never a wrong result (REGOLA ZERO).
[[nodiscard]] bool infinity_convergent(const MeijerGView& g1, const MeijerGView& g2) {
    const bool decay1 = (g1.q == g1.p + 1);  // exp(−η1 x) kernel
    const bool decay2 = (g2.q == g2.p + 1);
    const bool nogrow1 = (g1.p <= g1.q);
    const bool nogrow2 = (g2.p <= g2.q);
    return (decay1 && nogrow2) || (decay2 && nogrow1);
}

// η decidably a positive rational (branch-safe half-line evaluation, spec §4).
[[nodiscard]] bool positive_rational(ExprPtr eta) {
    auto r = rational_value(eta);
    return r.has_value() && r->numerator().is_positive();
}

}  // namespace

[[nodiscard]] Result<std::optional<ExprPtr>> pattern_mellin_g_convolution(
    const DefiniteContext& dc) {
    // Half-line only:  ∫_0^∞.
    ExprPtr zero = dc.lower;
    const auto* zlit = expr_cast<IntegerLit>(zero);
    if (zlit == nullptr || !zlit->value.is_zero()) return ok(std::optional<ExprPtr>{});
    if (!is_pos_infinity(dc.upper)) return ok(std::optional<ExprPtr>{});

    CASContext& ctx = dc.ctx;
    AstArena& arena = ctx.arena();
    const Symbol& var = dc.var;

    // Split integrand_normalized into: var-free coeff K, global monomial x^μ,
    // and the non-monomial cores (each destined to become a linear G).
    std::vector<ExprPtr> coeff;      // var-free constant factors
    std::vector<ExprPtr> cores;      // non-monomial, var-dependent factors
    Rational mu(BigInt(0));
    bool mu_ok = true;
    auto classify = [&](ExprPtr f) {
        if (!depends_on_var(f, var)) { coeff.push_back(f); return; }
        if (const auto* s = expr_cast<Symbol>(f); s && s->name == var.name) {
            mu = mu + Rational(BigInt(1)); return;
        }
        if (const auto* pw = expr_cast<Binary>(f);
            pw && pw->op == BinaryOp::Pow) {
            if (const auto* s = expr_cast<Symbol>(pw->left); s && s->name == var.name) {
                if (auto e = rational_value(pw->right); e.has_value()) { mu = mu + *e; return; }
                mu_ok = false; return;
            }
        }
        cores.push_back(f);
    };
    if (const auto* prod = expr_cast<Product>(dc.integrand_normalized)) {
        for (ExprPtr f : prod->factors) classify(f);
    } else {
        classify(dc.integrand_normalized);
    }
    if (!mu_ok || cores.size() != 2U) return ok(std::optional<ExprPtr>{});

    // Each core → K_i · G_i(η_i · x).
    auto lg1 = extract_linear_g(ctx, cores[0], var);
    auto lg2 = extract_linear_g(ctx, cores[1], var);
    if (!lg1.has_value() || !lg2.has_value()) return ok(std::optional<ExprPtr>{});

    if (!positive_rational(lg1->scale) || !positive_rational(lg2->scale))
        return ok(std::optional<ExprPtr>{});

    // Fold the global x^μ into G1 via §6.2 power shift:
    //   x^μ · G1(η1 x) = η1^{-μ} · G1(η1 x | a+μ ; b+μ).
    MeijerGView v1 = lg1->view;
    std::vector<ExprPtr> const_factors = std::move(coeff);
    const_factors.insert(const_factors.end(), lg1->prefactor.begin(), lg1->prefactor.end());
    const_factors.insert(const_factors.end(), lg2->prefactor.begin(), lg2->prefactor.end());
    if (!mu.numerator().is_zero()) {
        ExprPtr mu_lit = mu.is_integer()
            ? static_cast<ExprPtr>(arena.make<IntegerLit>(mu.numerator()))
            : static_cast<ExprPtr>(arena.make<RationalLit>(mu.numerator(), mu.denominator()));
        // rebuild the G1 node to feed meijerg_power_shift
        auto g1_node = symbolic::make_meijerg(ctx, v1.m, v1.n, v1.a, v1.b, v1.z);
        if (g1_node.is_error()) return ok(std::optional<ExprPtr>{});
        const auto* g1_call = expr_cast<FuncCall>(g1_node.value());
        if (g1_call == nullptr) return ok(std::optional<ExprPtr>{});
        auto shifted = symbolic::meijerg_power_shift(ctx, *g1_call, mu_lit);
        if (shifted.is_error()) return ok(std::optional<ExprPtr>{});
        const auto* sh_call = expr_cast<FuncCall>(shifted.value());
        if (sh_call == nullptr) return ok(std::optional<ExprPtr>{});
        auto sh_view = symbolic::view_meijerg(*sh_call);
        if (sh_view.is_error()) return ok(std::optional<ExprPtr>{});
        v1 = sh_view.value();
        // constant η1^{-μ}
        ExprPtr neg_mu = arena.make<RationalLit>(
            (Rational(BigInt(0)) - mu).numerator(),
            (Rational(BigInt(0)) - mu).denominator());
        const_factors.push_back(arena.make<Binary>(BinaryOp::Pow, lg1->scale, neg_mu));
    }
    const MeijerGView& v2 = lg2->view;

    // Convergence certificate (both endpoints) — the soundness gate.
    auto z1 = zero_exponent(v1);
    auto z2 = zero_exponent(v2);
    if (!z1.has_value() || !z2.has_value()) return ok(std::optional<ExprPtr>{});
    if (!(*z1 + *z2 + Rational(BigInt(1)) > Rational(BigInt(0))))  // λ0(G1)+λ0(G2) > -1
        return ok(std::optional<ExprPtr>{});
    if (!infinity_convergent(v1, v2)) return ok(std::optional<ExprPtr>{});

    // §6.7 parameter composition.  G1 = G^{m1,n1}_{p1,q1}(η1 x | a;b),
    // G2 = G^{m2,n2}_{p2,q2}(η2 x | c;d).  Result G^{M,N}_{P,Q}(η2/η1):
    //   M = n1+m2,  N = m1+n2,  P = q1+p2,  Q = p1+q2.
    const std::size_t m1 = v1.m, n1 = v1.n, p1 = v1.p, q1 = v1.q;
    const std::size_t m2 = v2.m, n2 = v2.n, p2 = v2.p, q2 = v2.q;
    const std::vector<ExprPtr>& a = v1.a;  // size p1
    const std::vector<ExprPtr>& b = v1.b;  // size q1
    const std::vector<ExprPtr>& c = v2.a;  // size p2
    const std::vector<ExprPtr>& d = v2.b;  // size q2

    // Upper list (size P = q1+p2):
    //   [ -b_1..-b_{m1},  c_1..c_{p2},  -b_{m1+1}..-b_{q1} ]
    std::vector<ExprPtr> upper;
    upper.reserve(q1 + p2);
    for (std::size_t k = 0; k < m1; ++k) upper.push_back(neg_param(ctx, b[k]));
    for (std::size_t k = 0; k < p2; ++k) upper.push_back(c[k]);
    for (std::size_t k = m1; k < q1; ++k) upper.push_back(neg_param(ctx, b[k]));
    // Lower list (size Q = p1+q2):
    //   [ -a_1..-a_{n1},  d_1..d_{q2},  -a_{n1+1}..-a_{p1} ]
    std::vector<ExprPtr> lower;
    lower.reserve(p1 + q2);
    for (std::size_t k = 0; k < n1; ++k) lower.push_back(neg_param(ctx, a[k]));
    for (std::size_t k = 0; k < q2; ++k) lower.push_back(d[k]);
    for (std::size_t k = n1; k < p1; ++k) lower.push_back(neg_param(ctx, a[k]));

    const std::size_t M = n1 + m2;
    const std::size_t N = m1 + n2;

    // Argument ω/η = η2/η1.
    ExprPtr arg = arena.make<Binary>(BinaryOp::Div, lg2->scale, lg1->scale);
    auto arg_s = ctx.simplify(arg);
    if (arg_s.is_ok()) arg = arg_s.value();

    auto g_res = symbolic::make_meijerg(ctx, M, N, std::move(upper), std::move(lower), arg);
    // A PoleOverlap / invalid rebuild means the generic §6.7 form does not
    // apply (confluent case needs the logarithmic variant) → skip, never fail.
    if (g_res.is_error()) return ok(std::optional<ExprPtr>{});

    // Normalise to p ≤ q via the §6.1 inversion identity (DLMF 16.19.1;
    // value-preserving: G^{M,N}_{P,Q}(z) = G^{N,M}_{Q,P}(1/z | 1-b ; 1-a)).
    // The composition (§6.7) can land in the p > q chart depending on which
    // core is labelled G1 (the integral is symmetric, both charts have the
    // SAME value); but slater_expand / the inverse table require p ≤ q, so an
    // un-normalised p > q result would stay an un-folded G even when it equals
    // an elementary closed form.  Inverting once puts it in the foldable chart.
    ExprPtr g_conv = g_res.value();
    if (const auto* gc = expr_cast<FuncCall>(g_conv)) {
        auto gv = symbolic::view_meijerg(*gc);
        if (gv.is_ok() && gv.value().p > gv.value().q) {
            auto inv = symbolic::meijerg_invert_argument(ctx, *gc);
            if (inv.is_ok()) g_conv = inv.value();
        }
    }

    // Assemble  K · (1/η1) · G_conv,  then fold to closed form.
    const_factors.push_back(arena.make<Binary>(BinaryOp::Div,
        arena.make<IntegerLit>(BigInt(1)), lg1->scale));
    const_factors.push_back(g_conv);
    ExprPtr assembled = const_factors.size() == 1U
        ? const_factors.front()
        : arena.make<Product>(std::move(const_factors));

    auto expanded = symbolic::expand_meijerg_nodes(ctx, assembled);
    if (expanded.is_error()) return ok(std::optional<ExprPtr>{});
    auto simplified = ctx.simplify(expanded.value());
    if (simplified.is_error()) return ok(std::optional<ExprPtr>{});

    return ok(std::optional<ExprPtr>(simplified.value()));
}

}  // namespace cas::calculus
