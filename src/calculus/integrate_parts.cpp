#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <vector>

namespace cas::calculus {

namespace {


[[nodiscard]] int get_ilate_priority(ExprPtr expr) {
    if (!expr) {
        return 100;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        const BuiltinOp func_id = call->func_id;
        // I: Inverse Trig + non-elementary special functions whose derivative
        // collapses the integrand to an elementary form (Erf' = 2/√π·exp(-x²),
        // Erfc' = -2/√π·exp(-x²)). Picking them as u in IBP is the only way
        // to integrate e.g. ∫ x·erf(x)·exp(x²) dx in closed form.
        if (func_id == BuiltinOp::Asin || func_id == BuiltinOp::Acos || func_id == BuiltinOp::Atan ||
            func_id == BuiltinOp::Erf  || func_id == BuiltinOp::Erfc) {
            return 1;
        }
        // F7.5.B1: inverse hyperbolic (parsed as BuiltinOp::Unknown; matched
        // by canonical name).  Same ILATE category as inverse trig — picking
        // them as `u` collapses x·asinh(x)-style integrals to elementary
        // rational forms.
        if (func_id == BuiltinOp::Unknown) {
            const std::string& nm = call->name;
            if (nm == "asinh" || nm == "arcsinh"
                || nm == "acosh" || nm == "arccosh"
                || nm == "atanh" || nm == "arctanh"
                || nm == "acoth" || nm == "arccoth") {
                return 1;
            }
        }
        // L: Logarithmic (Ln and Log are aliases — both natural log).
        if (func_id == BuiltinOp::Ln || func_id == BuiltinOp::Log) {
            return 2;
        }
        // T: Trigonometric
        if (func_id == BuiltinOp::Sin || func_id == BuiltinOp::Cos || func_id == BuiltinOp::Tan ||
            func_id == BuiltinOp::Sec || func_id == BuiltinOp::Csc || func_id == BuiltinOp::Cot) {
            return 4;
        }
        // E: Exponential
        if (func_id == BuiltinOp::Exp) {
            return 5;
        }
    }

    // A: Algebraic (Polynomials, Powers)
    if (expr_is<IntegerLit>(expr) || expr_is<RationalLit>(expr) || expr_is<Symbol>(expr)) {
        return 3;
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) {
            // A power inherits the ILATE class of its base: (log x)^n is still
            // Logarithmic, (asin x)^n Inverse, (sin x)^n Trig, etc., while a
            // purely algebraic base (x^n, (x+1)^n → Symbol/Sum/Product) resolves
            // to Algebraic (3). Without this, log(x)^2 was blanket-classified
            // Algebraic, so ∫ x·(log x)^n dx picked the wrong u and dead-ended
            // (INTEGRATE_NO_STRATEGY) even though it reduces to ∫ x·log x dx.
            return get_ilate_priority(bin->left);
        }
    }
    if (expr_is<Sum>(expr) || expr_is<Product>(expr)) {
        return 3;
    }

    return 6;
}

class IntegrationByPartsGuard {
public:
    IntegrationByPartsGuard(ExprPtr expr, std::size_t max_depth)
        : expr_(expr), max_depth_(max_depth) {}

    [[nodiscard]] Result<void> enter() {
        if (active_stack().size() >= max_depth_) {
            return fail<void>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Integration by parts depth budget exceeded",
                .hint = "serve una strategia diversa dal fallback ricorsivo",
            });
        }

        for (ExprPtr active : active_stack()) {
            if (structural_equal(active, expr_)) {
                return fail<void>(CASError{
                    .kind = CASErrorKind::Unimplemented,
                    .message = "Integration by parts cycle detected",
                    .hint = "il fallback euristico e' rientrato sullo stesso integrando",
                });
            }
        }

        active_stack().push_back(expr_);
        entered_ = true;
        return ok();
    }

    ~IntegrationByPartsGuard() {
        if (!entered_) {
            return;
        }
        active_stack().pop_back();
    }

    IntegrationByPartsGuard(const IntegrationByPartsGuard&) = delete;
    IntegrationByPartsGuard& operator=(const IntegrationByPartsGuard&) = delete;

    // True iff at least one integrate_by_parts is already active on this thread,
    // i.e. the current call is a *nested* recursive IBP rather than the entry.
    [[nodiscard]] static bool nested_active() { return !active_stack().empty(); }

private:
    [[nodiscard]] static std::vector<ExprPtr>& active_stack() {
        thread_local std::vector<ExprPtr> stack;
        return stack;
    }

    ExprPtr expr_;
    std::size_t max_depth_;
    bool entered_ = false;
};

}  // namespace

// Returns true if expr contains a FuncCall(Exp, arg) where arg is not a
// polynomial in var (e.g. exp(1/x), exp(1/x²)).  Used to abort IBP early
// for integrands that the Risch DE solver handles but IBP loops on
// (BUG-HANG-001: growing-denominator infinite IBP recursion).
[[nodiscard]] static bool ibp_has_exp_non_poly_arg(ExprPtr expr, const Symbol& var) {
    if (!expr) return false;
    if (const auto* fc = expr_cast<FuncCall>(expr)) {
        if (fc->func_id == BuiltinOp::Exp && fc->args.size() == 1U) {
            ExprPtr arg = fc->args[0];
            // Walk the arg: if it contains x^n with n<0, Div(a,f(x)), or
            // FuncCall(var-dependent), it is not a polynomial.
            struct Walk {
                const Symbol& v;
                bool found = false;
                void operator()(ExprPtr e) {
                    if (!e || found) return;
                    if (const auto* b = expr_cast<Binary>(e)) {
                        if (b->op == BinaryOp::Div && depends_on(b->right, v)) {
                            found = true; return;
                        }
                        if (b->op == BinaryOp::Pow && depends_on(b->left, v)) {
                            if (const auto* il = expr_cast<IntegerLit>(b->right))
                                if (il->value < BigInt(0)) { found = true; return; }
                            if (expr_is<RationalLit>(b->right)) { found = true; return; }
                            if (depends_on(b->right, v)) { found = true; return; }
                        }
                        (*this)(b->left); (*this)(b->right); return;
                    }
                    if (const auto* u = expr_cast<Unary>(e)) { (*this)(u->operand); return; }
                    if (const auto* s = expr_cast<Sum>(e))
                        for (ExprPtr t : s->terms) (*this)(t);
                    if (const auto* p = expr_cast<Product>(e))
                        for (ExprPtr f : p->factors) (*this)(f);
                    if (const auto* fc2 = expr_cast<FuncCall>(e))
                        for (ExprPtr a : fc2->args) (*this)(a);
                }
            } w{var};
            w(arg);
            if (w.found) return true;
        }
        // Recurse into function args.
        for (ExprPtr a : fc->args) if (ibp_has_exp_non_poly_arg(a, var)) return true;
        return false;
    }
    if (const auto* u = expr_cast<Unary>(expr)) return ibp_has_exp_non_poly_arg(u->operand, var);
    if (const auto* b = expr_cast<Binary>(expr))
        return ibp_has_exp_non_poly_arg(b->left, var) || ibp_has_exp_non_poly_arg(b->right, var);
    if (const auto* s = expr_cast<Sum>(expr))
        { for (ExprPtr t : s->terms) if (ibp_has_exp_non_poly_arg(t, var)) return true; return false; }
    if (const auto* p = expr_cast<Product>(expr))
        { for (ExprPtr f : p->factors) if (ibp_has_exp_non_poly_arg(f, var)) return true; return false; }
    return false;
}

// True if `e` contains a square root (or fractional power) of a var-dependent
// base — i.e. a genuine radical in `var`.
[[nodiscard]] static bool expr_contains_var_radical(ExprPtr e, const Symbol& var) {
    if (!e) return false;
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        if (fc->func_id == BuiltinOp::Sqrt)
            for (ExprPtr a : fc->args) if (depends_on(a, var)) return true;
        for (ExprPtr a : fc->args) if (expr_contains_var_radical(a, var)) return true;
        return false;
    }
    if (const auto* b = expr_cast<Binary>(e)) {
        if (b->op == BinaryOp::Pow && depends_on(b->left, var)
            && expr_is<RationalLit>(b->right)) return true;  // fractional power of var
        return expr_contains_var_radical(b->left, var) || expr_contains_var_radical(b->right, var);
    }
    if (const auto* u = expr_cast<Unary>(e)) return expr_contains_var_radical(u->operand, var);
    if (const auto* s = expr_cast<Sum>(e))
        { for (ExprPtr t : s->terms) if (expr_contains_var_radical(t, var)) return true; return false; }
    if (const auto* p = expr_cast<Product>(e))
        { for (ExprPtr f : p->factors) if (expr_contains_var_radical(f, var)) return true; return false; }
    return false;
}

// True if `base` is a Sum that contains a var-radical, e.g. (x + √(x²+1)).
// Such a denominator needs *conjugate rationalisation* — a bare radical like
// √(1-x²) (FuncCall, not a Sum) does not and is excluded.
[[nodiscard]] static bool is_radical_sum(ExprPtr base, const Symbol& var) {
    return expr_is<Sum>(base) && expr_contains_var_radical(base, var);
}

// True if `e` has a *radical-sum in a denominator* w.r.t. var: a factor
// Pow(base, k<0) (or Div by such) where `base` is a Sum containing a square
// root of var, e.g. (x+√(x²+1))^{-1}.
//
// Anti-hang: IBP cannot reduce these — differentiating the negative power only
// deepens it, so by-parts recursion diverges with a growing negative exponent
// (∫log(x+√(x²+1)) dx degenerates this way once the log factor is consumed).
// They need radical-denominator rationalisation, not IBP.  The clean bare-
// radical case √(1-x²)^{-1} (from ∫asin(x), ∫x·asin(x), …) is NOT matched —
// its base is a FuncCall, not a Sum — and stays on the working sqrt path.
[[nodiscard]] static bool ibp_has_radical_sum_denominator(ExprPtr e, const Symbol& var) {
    if (!e) return false;
    if (const auto* b = expr_cast<Binary>(e)) {
        if (b->op == BinaryOp::Pow) {
            bool neg = false;
            if (const auto* il = expr_cast<IntegerLit>(b->right)) neg = il->value.is_negative();
            else if (const auto* rl = expr_cast<RationalLit>(b->right)) neg = rl->numerator.is_negative();
            if (neg && is_radical_sum(b->left, var)) return true;
        }
        if (b->op == BinaryOp::Div && is_radical_sum(b->right, var)) return true;
        return ibp_has_radical_sum_denominator(b->left, var)
            || ibp_has_radical_sum_denominator(b->right, var);
    }
    if (const auto* u = expr_cast<Unary>(e)) return ibp_has_radical_sum_denominator(u->operand, var);
    if (const auto* s = expr_cast<Sum>(e))
        { for (ExprPtr t : s->terms) if (ibp_has_radical_sum_denominator(t, var)) return true; return false; }
    if (const auto* p = expr_cast<Product>(e))
        { for (ExprPtr f : p->factors) if (ibp_has_radical_sum_denominator(f, var)) return true; return false; }
    if (const auto* fc = expr_cast<FuncCall>(e))
        { for (ExprPtr a : fc->args) if (ibp_has_radical_sum_denominator(a, var)) return true; return false; }
    return false;
}

Result<ExprPtr> integrate_by_parts(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& context) {
    // BUG-HANG-001: IBP on exp(non-polynomial) * rational diverges.
    // The Risch DE solver handles these; skip IBP to prevent infinite loop.
    if (ibp_has_exp_non_poly_arg(expr, var)) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "integrate_by_parts: exp(non-polynomial) factor detected; defer to Risch DE (BUG-HANG-001)",
            .hint = "integrate_risch with Risch DE rational solver",
        });
    }

    // HARDCODE-OF-PASSAGE: HC-IBP-RADSUM-RATIONALIZE
    // Anti-hang: a *nested* IBP whose integrand has a radical-SUM in a
    // denominator (e.g. (x+√(x²+1))^{-1}) cannot terminate — each by-parts step
    // deepens the negative power, so the recursion diverges. This is exactly the
    // ∫log(x+√(x²+1)) dx divergence once the log factor is consumed. Defer (the
    // sub-integrand needs conjugate rationalisation, not IBP). Only nested calls
    // are gated (top-level IBP is exempt, nested_active() false) so the legit
    // first by-parts step that *produces* such a sub-integrand still runs; and a
    // bare radical √(x²+1)^{-1} (FuncCall base, not a Sum) is NOT matched, so the
    // working ∫asin / quadratic-radical path is untouched.
    if (IntegrationByPartsGuard::nested_active()
        && ibp_has_radical_sum_denominator(expr, var)) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "integrate_by_parts: nested IBP on radical-sum-denominator integrand diverges; needs rationalisation (anti-hang)",
            .hint = std::nullopt,
        });
    }

    [[maybe_unused]] const auto hint_guard = context.with_hint(symbolic::SimplifyHints{ .fold_exp_products = true });

    IntegrationByPartsGuard guard(expr, context.max_integrate_by_parts_depth());
    auto guard_result = guard.enter();
    if (guard_result.is_error()) {
        return fail<ExprPtr>(guard_result.error());
    }

    const auto* product = expr_cast<Product>(expr);
    if (!product || product->factors.size() < 2) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::Unimplemented,
            .message = "Integration by parts requires a product of at least two factors",
        });
    }

    std::vector<ExprPtr> factors = product->factors;
    auto it = std::min_element(factors.begin(), factors.end(), [](ExprPtr a, ExprPtr b) {
        return get_ilate_priority(a) < get_ilate_priority(b);
    });

    std::size_t u_idx = std::distance(factors.begin(), it);
    ExprPtr u = factors[u_idx];
    
    std::vector<ExprPtr> dv_factors;
    dv_factors.reserve(factors.size() - 1);
    for (std::size_t i = 0; i < factors.size(); ++i) {
        if (i != u_idx) {
            dv_factors.push_back(factors[i]);
        }
    }

    AstArena& arena = context.arena();
    ExprPtr dv = (dv_factors.size() == 1) ? dv_factors[0] : arena.make<Product>(std::move(dv_factors));

    auto du_res = diff(u, var, 1U, context);
    if (du_res.is_error()) {
        return du_res;
    }
    ExprPtr du = du_res.value();

    auto v_res = integrate(dv, var, context);
    if (v_res.is_error()) {
        return v_res;
    }
    ExprPtr v = v_res.value();

    ExprPtr uv = arena.make<Product>(std::vector<ExprPtr>{u, v});

    ExprPtr vdu = arena.make<Product>(std::vector<ExprPtr>{v, du});

    // HC-F75-B1-IBP-DOUBLE-APPLY fix: simplify the v·du Product before
    // recursive integrate, otherwise patterns like ∫x·log(x)dx produce a
    // sub-integrand (x²/2)·(1/x) which re-routes through integrate_by_parts
    // (both factors algebraic, equal ILATE priority) instead of collapsing
    // to (1/2)·x — leading to recursive by-parts expansion and a 4-term
    // redundant Sum in the final result.
    auto vdu_simp = context.simplify(vdu);
    if (vdu_simp.is_ok()) vdu = vdu_simp.value();

    auto int_vdu_res = integrate(vdu, var, context);
    if (int_vdu_res.is_error()) {
        return int_vdu_res;
    }
    ExprPtr int_vdu = int_vdu_res.value();

    ExprPtr ibp_result = arena.make<Sum>(std::vector<ExprPtr>{
        uv,
        arena.make<Unary>(UnaryOp::Neg, int_vdu)
    });

    // Verification by differentiation: IBP returns uv - ∫v·du; if the
    // recursive ∫v·du did not actually close (e.g. budget exhausted with
    // a partial/zero stub), the algebraic identity D(result) = integrand
    // breaks. Without this guard, callers like ∫ 1/(x·ln(x)) receive
    // the wrong closed form ln(x)^-1·ln|x| instead of falling through
    // to the Risch logarithmic-derivative recognizer.
    auto D_res = diff(ibp_result, var, 1U, context);
    if (D_res.is_ok()) {
        ExprPtr delta = arena.make<Binary>(BinaryOp::Sub, D_res.value(), expr);
        auto delta_tog = algebra::together(delta, context);
        ExprPtr delta_for_simp = delta_tog.is_ok() ? delta_tog.value() : delta;
        auto delta_simp = context.simplify(delta_for_simp);
        bool is_zero = delta_simp.is_ok()
            && expr_is<IntegerLit>(delta_simp.value())
            && expr_ref<IntegerLit>(delta_simp.value()).value.is_zero();
        if (!is_zero) {
            return fail<ExprPtr>(CASError{
                .kind = CASErrorKind::Unimplemented,
                .message = "Integration by parts: recursive sub-integral did not close (D(result) != integrand)",
                .hint = "fallback to algorithmic strategy (Risch/Hermite/RT) required",
            });
        }
    }

    return ok(ibp_result);
}

}  // namespace cas::calculus
