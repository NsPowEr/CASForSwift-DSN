#include "calculus_internal.hpp"
#include "cas/algebra.hpp"
#include "cas/error.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace cas::calculus {

namespace {


[[nodiscard]] int get_ilate_priority(ExprPtr expr) {
    if (!expr) {
        return 100;
    }

    if (const auto* call = expr_cast<FuncCall>(expr)) {
        const BuiltinOp func_id = call->func_id;
        // I: Inverse Trig
        if (func_id == BuiltinOp::Asin || func_id == BuiltinOp::Acos || func_id == BuiltinOp::Atan) {
            return 1;
        }
        // L: Logarithmic
        if (func_id == BuiltinOp::Ln) {
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
            return 3;
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

Result<ExprPtr> integrate_by_parts(
    ExprPtr expr,
    const Symbol& var,
    symbolic::CASContext& context) {
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
