// F7.5.A1 / HC-F75-CYCLOTOMIC-ROOTOF — recognise the canonical roots of
// `(x^n - c^n)/(x - c)` so a RootOf-on-such-a-polynomial compares as
// equal to the corresponding `c·exp(2πi·k/n)` form that Maxima (and
// hand-written corpus answers) emit for cyclotomic equations.
//
// Pattern recognised:
//
//   p(x) = sum_{i=0..d} c^i · x^{d-i}                                     (1)
//        = x^d + c·x^{d-1} + c^2·x^{d-2} + … + c^d
//        = (x^{d+1} - c^{d+1}) / (x - c)
//
// where c is a non-zero rational. The (d+1)-th roots of unity are
// {1, ω, ω^2, …, ω^d} with ω = exp(2πi/(d+1)); the d roots of (1) are
// exactly {c·ω, c·ω^2, …, c·ω^d} (i.e. all (d+1)-th roots of c^{d+1}
// excluding c itself).
//
// The detector is structural and exact (rational arithmetic only); no
// floating point, no closed pattern table beyond the algebraic
// identity. The enumerator emits the d roots as
// `c · exp((2πi·m)/(d+1))` Expr trees with `m = 1..d`, using the
// canonical AST constants MathConstant::Pi and MathConstant::I so the
// downstream `mathematically_equal` comparison can normalise both
// sides without notational hand-tuning.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace cas::algebra {

namespace {

[[nodiscard]] std::optional<Rational> as_rational(ExprPtr e) {
    if (!e) return std::nullopt;
    if (const auto* lit = expr_cast<IntegerLit>(e)) {
        return Rational(lit->value, BigInt(1));
    }
    if (const auto* lit = expr_cast<RationalLit>(e)) {
        return Rational(lit->numerator, lit->denominator);
    }
    return std::nullopt;
}

[[nodiscard]] ExprPtr rational_to_expr(const Rational& r, AstArena& arena) {
    if (r.denominator() == BigInt(1)) {
        return arena.make<IntegerLit>(r.numerator());
    }
    return arena.make<RationalLit>(r.numerator(), r.denominator());
}

[[nodiscard]] int igcd(int a, int b) {
    while (b != 0) { int t = a % b; a = b; b = t; }
    return a < 0 ? -a : a;
}

// Match an ascending integer-coefficient vector against the n-th cyclotomic
// polynomial Φ_n for some n, returning that n. deg(Φ_n) = φ(n), so a degree
// mismatch rules out a candidate immediately. φ(n) ≥ √(n/2) ⇒ n ≤ 2·φ(n)² = 2d²,
// giving a complete, finite search window (capped by ctx.max_cyclotomic_n when
// the user has set a budget). Exact integer comparison — no floating point.
[[nodiscard]] std::optional<int> detect_cyclotomic_index(
    const std::vector<BigInt>& coeffs_asc, symbolic::CASContext& ctx) {
    if (coeffs_asc.size() < 2U) return std::nullopt;        // need degree ≥ 1
    const std::size_t d = coeffs_asc.size() - 1U;
    if (coeffs_asc[d] != BigInt(1)) return std::nullopt;     // Φ_n is monic
    const long long dd = static_cast<long long>(d);
    long long bound = 2 * dd * dd + 2;
    const int cap = ctx.max_cyclotomic_n();
    if (cap > 0 && static_cast<long long>(cap) < bound) bound = cap;
    for (long long n = 1; n <= bound; ++n) {
        algebra::IntPoly phi = algebra::compute_cyclotomic(static_cast<int>(n));
        if (phi.size() != coeffs_asc.size()) continue;       // φ(n) ≠ d
        bool eq = true;
        for (std::size_t k = 0; k < coeffs_asc.size(); ++k) {
            if (phi[k] != coeffs_asc[k]) { eq = false; break; }
        }
        if (eq) return static_cast<int>(n);
    }
    return std::nullopt;
}

}  // namespace

// Recognise p(x) = sum_{i=0..d} c^i · x^{d-i} with c ∈ Q \ {0}.
// Returns the pair {c, n=d+1} on match; nullopt otherwise.
//
// The polynomial argument is given via its univariate coefficient
// vector in ascending power order, exactly what
// `algebra::univariate_coefficients(poly, var, ctx)` produces:
//
//     coeffs[k] = coefficient of x^k          (k = 0..d).
//
// In ascending form the geometric pattern reads coeffs[k] = c^{d-k} so
//   coeffs[d]   = 1                       (monic)
//   coeffs[k+1] = coeffs[k] / c           (ratio 1/c)
[[nodiscard]] std::optional<std::pair<Rational, std::size_t>>
detect_geometric_polynomial_coeffs(const std::vector<ExprPtr>& coeffs_asc) {
    if (coeffs_asc.size() < 3U) return std::nullopt;  // need degree ≥ 2.
    const std::size_t d = coeffs_asc.size() - 1U;

    auto top = as_rational(coeffs_asc[d]);
    if (!top.has_value()) return std::nullopt;
    if (top->numerator() != top->denominator()) return std::nullopt;  // monic.

    auto c_candidate = as_rational(coeffs_asc[d - 1]);
    if (!c_candidate.has_value()) return std::nullopt;
    if (c_candidate->numerator().is_zero()) return std::nullopt;

    const Rational c = *c_candidate;

    // Verify coeffs[k] = c^{d-k} for k = 0..d.
    Rational power = Rational(BigInt(1), BigInt(1));
    for (std::size_t i = 0; i < d; ++i) {
        power = power * c;  // power = c^{i+1}, expected coeffs[d-i-1].
        auto observed = as_rational(coeffs_asc[d - i - 1U]);
        if (!observed.has_value()) return std::nullopt;
        if (*observed != power) return std::nullopt;
    }
    return std::make_pair(c, d + 1U);
}

// Expand a RootOf whose polynomial matches the (x^n - c^n)/(x - c)
// geometric pattern into the d concrete roots `c · exp(2πi·m/n)` for
// m = 1..d. Returns nullopt when the polynomial is not of this form
// or its coefficients are not all rational literals.
[[nodiscard]] std::optional<std::vector<ExprPtr>>
enumerate_geometric_rootof(const RootOf& node, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto coeffs = algebra::univariate_coefficients(
        node.polynomial, node.variable, ctx);
    if (!coeffs.is_ok()) return std::nullopt;
    auto detected = detect_geometric_polynomial_coeffs(coeffs.value());
    if (!detected.has_value()) return std::nullopt;
    const Rational c = detected->first;
    const std::size_t n = detected->second;
    const std::size_t d = n - 1U;

    ExprPtr c_expr = rational_to_expr(c, arena);
    ExprPtr two = arena.make<IntegerLit>(BigInt(2));
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr i = arena.make<Constant>(MathConstant::I);

    // Emit each root twice: once with the natural angle m * 2π/n in
    // [0, 2π), and once with the symmetric angle (m - n) * 2π/n in
    // (-2π, 0) which equals the same root modulo 2πi. Maxima
    // canonicalises angles to (-π, π], so without the negative-angle
    // variant the comparison would miss e.g.
    //     2 * exp(-4πi/5)  (Maxima)  ≡  2 * exp(6πi/5)  (our positive form).
    // mathematically_equal does not yet recognise the
    // `exp(a) == exp(a - 2πi)` identity, so we make both representations
    // available in the candidate set.
    const long long n_signed = static_cast<long long>(n);
    auto emit = [&](long long m_signed) -> ExprPtr {
        ExprPtr m_expr = arena.make<IntegerLit>(BigInt(m_signed));
        ExprPtr n_expr = arena.make<IntegerLit>(BigInt(n_signed));
        ExprPtr num = arena.make<Product>(
            std::vector<ExprPtr>{two, m_expr, pi, i});
        ExprPtr arg = arena.make<Binary>(BinaryOp::Div, num, n_expr);
        ExprPtr exp_term = arena.make<FuncCall>(
            BuiltinOp::Exp, std::vector<ExprPtr>{arg});
        ExprPtr root = arena.make<Product>(
            std::vector<ExprPtr>{c_expr, exp_term});
        // Keep the unsimplified `c * exp(2πim/n)` form so the downstream
        // mathematically_equal pipeline compares cyclotomic forms
        // directly without simplify potentially routing through Ferrari
        // closed-form radicals that diverge structurally from the Maxima
        // exp-canonical output.
        return root;
    };

    std::vector<ExprPtr> roots;
    roots.reserve(2U * d);
    for (std::size_t m = 1U; m <= d; ++m) {
        const long long m_pos = static_cast<long long>(m);
        const long long m_neg = m_pos - n_signed;
        roots.push_back(emit(m_pos));
        roots.push_back(emit(m_neg));
    }
    return roots;
}

// Expand a RootOf whose polynomial is the GENUINE n-th cyclotomic Φ_n into its
// φ(n) primitive roots `exp(2πi·m/n)` with gcd(m,n)=1. Complements the geometric
// enumerator above, which only covers Φ_p for prime p (= (x^p−1)/(x−1)); this
// handles composite n where Φ_n is not geometric (Φ_8 = x⁴+1, Φ_12 = x⁴−x²+1, …).
// Each root is emitted twice (positive angle m·2π/n and its (m−n) negative-angle
// representative) so the comparison aligns with Maxima's (−π, π] canonical range,
// matching the geometric enumerator's convention. nullopt when the polynomial is
// not a cyclotomic Φ_n with integer-literal coefficients.
[[nodiscard]] std::optional<std::vector<ExprPtr>>
enumerate_cyclotomic_rootof(const RootOf& node, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    auto coeffs = algebra::univariate_coefficients(
        node.polynomial, node.variable, ctx);
    if (!coeffs.is_ok()) return std::nullopt;

    std::vector<BigInt> int_coeffs;
    int_coeffs.reserve(coeffs.value().size());
    for (ExprPtr c : coeffs.value()) {
        const auto* il = expr_cast<IntegerLit>(c);
        if (il == nullptr) return std::nullopt;  // Φ_n has integer coefficients
        int_coeffs.push_back(il->value);
    }

    auto n_opt = detect_cyclotomic_index(int_coeffs, ctx);
    if (!n_opt.has_value()) return std::nullopt;
    const int n = *n_opt;

    ExprPtr two = arena.make<IntegerLit>(BigInt(2));
    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr i = arena.make<Constant>(MathConstant::I);
    ExprPtr n_expr = arena.make<IntegerLit>(BigInt(n));
    auto emit = [&](long long m_signed) -> ExprPtr {
        ExprPtr num = arena.make<Product>(std::vector<ExprPtr>{
            two, arena.make<IntegerLit>(BigInt(m_signed)), pi, i});
        ExprPtr arg = arena.make<Binary>(BinaryOp::Div, num, n_expr);
        return arena.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{arg});
    };

    std::vector<ExprPtr> roots;
    for (int m = 1; m < n; ++m) {
        if (igcd(m, n) != 1) continue;  // primitive roots only
        roots.push_back(emit(m));
        roots.push_back(emit(static_cast<long long>(m) - n));
    }
    if (roots.empty()) return std::nullopt;
    return roots;
}

}  // namespace cas::algebra

namespace cas::symbolic {
[[nodiscard]] Result<bool> mathematically_equal(
    ExprPtr lhs, ExprPtr rhs, CASContext& context);
}  // namespace cas::symbolic

namespace cas::algebra {

[[nodiscard]] std::optional<bool> try_rootof_decision(
    ExprPtr lhs, ExprPtr rhs, symbolic::CASContext& ctx) {
    const RootOf* lhs_root = expr_cast<RootOf>(lhs);
    const RootOf* rhs_root = expr_cast<RootOf>(rhs);
    if (lhs_root == nullptr && rhs_root == nullptr) return std::nullopt;

    // Distinct indices on the same minimal polynomial → definitely
    // different roots. Cuts off the polynomial_normal_form path which
    // would otherwise normalise RootOf-RootOf differences to zero by
    // treating both nodes as the same opaque algebraic generator.
    if (lhs_root != nullptr && rhs_root != nullptr &&
        lhs_root->root_index.has_value() &&
        rhs_root->root_index.has_value() &&
        *lhs_root->root_index != *rhs_root->root_index &&
        lhs_root->variable.name == rhs_root->variable.name &&
        structural_equal(lhs_root->polynomial, rhs_root->polynomial)) {
        return false;
    }

    // Exactly-one-side RootOf with a geometric (cyclotomic) minimal
    // polynomial: expand the d closed-form roots and accept iff some
    // enumerated root matches the other side. The "both sides RootOf"
    // case is intentionally NOT handled here — it falls through to the
    // general comparison so set-level callers can still match
    // RootOf(p, k) against another notation for the same root.
    if ((lhs_root == nullptr) != (rhs_root == nullptr)) {
        const RootOf* root = lhs_root != nullptr ? lhs_root : rhs_root;
        ExprPtr other = lhs_root != nullptr ? rhs : lhs;
        auto enumerated = enumerate_geometric_rootof(*root, ctx);
        if (!enumerated.has_value())
            enumerated = enumerate_cyclotomic_rootof(*root, ctx);  // T-025: composite Φ_n
        if (!enumerated.has_value()) return std::nullopt;
        for (ExprPtr candidate : *enumerated) {
            auto eq = symbolic::mathematically_equal(candidate, other, ctx);
            if (eq.is_error()) continue;
            if (eq.value()) return true;
        }
        // Polynomial was geometric but no enumerated root matched — the
        // other side is definitively not a root of this polynomial.
        return false;
    }
    return std::nullopt;
}

}  // namespace cas::algebra
