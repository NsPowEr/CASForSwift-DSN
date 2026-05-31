#include "simplify_impl.hpp"
#include "cas/numtheory.hpp"
#include "cas/error_helpers.hpp"

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_funcall_special(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before) {

    if (op == BuiltinOp::Gamma && args.size() == 1U) {
        // Gamma(n) = (n-1)! for positive integer n >= 1
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            if (il->value > BigInt(0)) {
                BigInt n = il->value;
                BigInt result(1);
                for (BigInt k(1); k < n; k += BigInt(1)) result *= k;
                return ok(make_integer(arena_, result));
            }
        }
        // Half-integer: Gamma(1/2) = sqrt(π), recursion via Gamma(z+1) = z*Gamma(z).
        if (const auto* rl = expr_cast<RationalLit>(args.front())) {
            if (rl->denominator == BigInt(2)) {
                BigInt num = rl->numerator;
                Rational z(num, BigInt(2));
                const Rational half(BigInt(1), BigInt(2));
                Rational factor(BigInt(1));
                bool ok_chain = true;
                std::size_t safety = 0;
                const std::size_t safety_max = (context_ != nullptr)
                    ? context_->max_gamma_recursion() : 1024U;
                while (z != half && safety++ < safety_max) {
                    if (z > half) {
                        Rational z_minus_1 = z - Rational(BigInt(1));
                        factor = factor * z_minus_1;
                        z = z_minus_1;
                    } else {
                        if (z.numerator().is_zero()) { ok_chain = false; break; }
                        factor = factor / z;
                        z = z + Rational(BigInt(1));
                    }
                }
                if (ok_chain && safety < safety_max) {
                    ExprPtr sqrt_pi = arena_.make<FuncCall>(
                        BuiltinOp::Sqrt,
                        std::vector<ExprPtr>{arena_.make<Constant>(MathConstant::Pi)});
                    if (factor == Rational(BigInt(1))) return simplify_expr(sqrt_pi);
                    ExprPtr factor_expr = factor.denominator() == BigInt(1)
                        ? static_cast<ExprPtr>(arena_.make<IntegerLit>(factor.numerator()))
                        : static_cast<ExprPtr>(arena_.make<RationalLit>(
                              factor.numerator(), factor.denominator()));
                    return simplify_expr(arena_.make<Product>(
                        std::vector<ExprPtr>{factor_expr, sqrt_pi}));
                }
            }
        }
        // Functional equation: Gamma(z+n) -> z*(z+1)*...*(z+n-1)*Gamma(z)
        if (const auto* sum = expr_cast<Sum>(args.front()); sum && sum->terms.size() >= 2U) {
            int idx = -1;
            BigInt shift(0);
            for (std::size_t i = 0; i < sum->terms.size(); ++i) {
                if (const auto* il = expr_cast<IntegerLit>(sum->terms[i])) {
                    if (idx == -1) { idx = static_cast<int>(i); shift = il->value; }
                    else { idx = -1; break; }
                }
            }
            if (idx >= 0 && !shift.is_zero()) {
                std::vector<ExprPtr> rest;
                rest.reserve(sum->terms.size() - 1U);
                for (std::size_t i = 0; i < sum->terms.size(); ++i) {
                    if (static_cast<int>(i) != idx) rest.push_back(sum->terms[i]);
                }
                ExprPtr z = rest.size() == 1U ? rest.front() : arena_.make<Sum>(std::move(rest));
                std::vector<ExprPtr> factors;
                if (shift > BigInt(0)) {
                    BigInt n = shift;
                    for (BigInt k(0); k < n; k += BigInt(1)) {
                        ExprPtr term = k.is_zero()
                            ? z
                            : static_cast<ExprPtr>(arena_.make<Sum>(std::vector<ExprPtr>{
                                  z, arena_.make<IntegerLit>(k)}));
                        factors.push_back(term);
                    }
                    factors.push_back(arena_.make<FuncCall>(
                        BuiltinOp::Gamma, std::vector<ExprPtr>{z}));
                    return simplify_expr(arena_.make<Product>(std::move(factors)));
                }
                BigInt n = -shift;
                std::vector<ExprPtr> denom_factors;
                for (BigInt k(1); k <= n; k += BigInt(1)) {
                    denom_factors.push_back(arena_.make<Sum>(std::vector<ExprPtr>{
                        z, arena_.make<Unary>(UnaryOp::Neg, arena_.make<IntegerLit>(k))}));
                }
                ExprPtr denom = denom_factors.size() == 1U
                    ? denom_factors.front()
                    : arena_.make<Product>(std::move(denom_factors));
                return simplify_expr(arena_.make<Binary>(
                    BinaryOp::Div,
                    arena_.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{z}),
                    denom));
            }
        }
    }

    // Digamma ψ(z) = Γ'(z)/Γ(z)
    if (op == BuiltinOp::Digamma && args.size() == 1U) {
        auto build_euler_gamma = [&]() -> ExprPtr {
            return arena_.make<Unary>(UnaryOp::Neg,
                arena_.make<Constant>(MathConstant::EulerGamma));
        };
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            if (il->value > BigInt(0)) {
                if (il->value.bit_length() > 16)
                    // F0.8-MIGRATED
                    return make_unimplemented<ExprPtr>(
                        "symbolic", "simplify_funcall_special",
                        "Digamma ψ(n) with positive integer n > 2^16",
                        error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                        "Reduce n or use the asymptotic expansion ψ(n) ≈ ln(n) - 1/(2n) - ...",
                        "F1.x");
                const std::uint64_t n = il->value.to_u64();
                std::vector<ExprPtr> terms;
                terms.push_back(build_euler_gamma());
                for (std::uint64_t k = 1U; k < n; ++k)
                    terms.push_back(arena_.make<RationalLit>(
                        BigInt(1), BigInt(static_cast<std::int64_t>(k))));
                if (terms.size() == 1U) return simplify_expr(terms.front());
                return simplify_expr(arena_.make<Sum>(std::move(terms)));
            }
        }
        if (const auto* rl = expr_cast<RationalLit>(args.front());
            rl != nullptr && rl->denominator == BigInt(2)) {
            BigInt num = rl->numerator;
            if (num > BigInt(0) && (num.to_u64() % 2U) == 1U) {
                const std::uint64_t two_m_plus_1 = num.to_u64();
                const std::uint64_t m = (two_m_plus_1 - 1U) / 2U;
                ExprPtr ln2 = arena_.make<FuncCall>(BuiltinOp::Ln,
                    std::vector<ExprPtr>{make_integer(arena_, BigInt(2))});
                ExprPtr two_ln2 = arena_.make<Product>(std::vector<ExprPtr>{
                    make_integer(arena_, BigInt(2)), ln2});
                std::vector<ExprPtr> terms;
                terms.push_back(build_euler_gamma());
                terms.push_back(arena_.make<Unary>(UnaryOp::Neg, two_ln2));
                for (std::uint64_t k = 1U; k <= m; ++k)
                    terms.push_back(arena_.make<RationalLit>(
                        BigInt(2), BigInt(static_cast<std::int64_t>(2U * k - 1U))));
                return simplify_expr(arena_.make<Sum>(std::move(terms)));
            }
        }
        if (const auto* sum = expr_cast<Sum>(args.front()); sum && sum->terms.size() >= 2U) {
            int idx = -1;
            BigInt shift(0);
            for (std::size_t i = 0; i < sum->terms.size(); ++i) {
                if (const auto* il2 = expr_cast<IntegerLit>(sum->terms[i])) {
                    if (idx == -1) { idx = static_cast<int>(i); shift = il2->value; }
                    else { idx = -1; break; }
                }
            }
            if (idx >= 0 && !shift.is_zero() && !shift.is_negative()
                && shift.bit_length() <= 16) {
                std::vector<ExprPtr> rest;
                rest.reserve(sum->terms.size() - 1U);
                for (std::size_t i = 0; i < sum->terms.size(); ++i) {
                    if (static_cast<int>(i) != idx) rest.push_back(sum->terms[i]);
                }
                ExprPtr z = rest.size() == 1U ? rest.front() : arena_.make<Sum>(std::move(rest));
                const std::uint64_t n = shift.to_u64();
                std::vector<ExprPtr> terms;
                terms.push_back(arena_.make<FuncCall>(BuiltinOp::Digamma,
                    std::vector<ExprPtr>{z}));
                for (std::uint64_t k = 0U; k < n; ++k) {
                    ExprPtr shift_k = (k == 0U)
                        ? z
                        : static_cast<ExprPtr>(arena_.make<Sum>(std::vector<ExprPtr>{
                              z, make_integer(arena_, BigInt(static_cast<std::int64_t>(k)))}));
                    terms.push_back(arena_.make<Binary>(BinaryOp::Div,
                        make_integer(arena_, BigInt(1)), shift_k));
                }
                return simplify_expr(arena_.make<Sum>(std::move(terms)));
            }
        }
    }

    // Polygamma ψ^(n)(z)
    if (op == BuiltinOp::Polygamma && args.size() == 2U) {
        const auto* il_n = expr_cast<IntegerLit>(args[0]);
        if (il_n != nullptr && !il_n->value.is_negative() && il_n->value.bit_length() <= 16) {
            const std::uint64_t n = il_n->value.to_u64();
            if (n == 0U)
                return simplify_expr(arena_.make<FuncCall>(BuiltinOp::Digamma,
                    std::vector<ExprPtr>{args[1]}));
            if (const auto* il_z = expr_cast<IntegerLit>(args[1]);
                il_z != nullptr && il_z->value == BigInt(1)) {
                BigInt fact(1);
                for (std::uint64_t k = 2U; k <= n; ++k)
                    fact *= BigInt(static_cast<std::int64_t>(k));
                ExprPtr fact_expr = make_integer(arena_, fact);
                ExprPtr zeta_expr = arena_.make<FuncCall>(BuiltinOp::Zeta,
                    std::vector<ExprPtr>{make_integer(arena_,
                        BigInt(static_cast<std::int64_t>(n + 1U)))});
                ExprPtr prod = arena_.make<Product>(std::vector<ExprPtr>{fact_expr, zeta_expr});
                if ((n % 2U) == 0U) prod = arena_.make<Unary>(UnaryOp::Neg, prod);
                return simplify_expr(prod);
            }
        }
    }

    // Beta B(x,y) = Gamma(x)*Gamma(y)/Gamma(x+y)
    if (op == BuiltinOp::Beta && args.size() == 2U) {
        ExprPtr x_arg = args[0];
        ExprPtr y_arg = args[1];
        ExprPtr gx  = arena_.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{x_arg});
        ExprPtr gy  = arena_.make<FuncCall>(BuiltinOp::Gamma, std::vector<ExprPtr>{y_arg});
        ExprPtr gxy = arena_.make<FuncCall>(BuiltinOp::Gamma,
            std::vector<ExprPtr>{arena_.make<Sum>(std::vector<ExprPtr>{x_arg, y_arg})});
        ExprPtr num = arena_.make<Product>(std::vector<ExprPtr>{gx, gy});
        return simplify_expr(arena_.make<Binary>(BinaryOp::Div, num, gxy));
    }

    // Pochhammer (rising factorial) (x)_n
    if (op == BuiltinOp::Pochhammer && args.size() == 2U) {
        ExprPtr x_arg = args[0];
        if (const auto* il = expr_cast<IntegerLit>(args[1])) {
            if (il->value.bit_length() > 16)
                // F0.8-MIGRATED
                return make_unimplemented<ExprPtr>(
                    "symbolic", "simplify_funcall_special",
                    "Pochhammer (x)_n with n > 2^16",
                    error::reason_codes::SYMBOLIC_DEGREE_TOO_LARGE,
                    "Cap n via ctx.max_polynomial_degree or use Gamma(x+n)/Gamma(x) identity",
                    "F1.x");
            if (il->value.is_zero()) return ok(make_integer(arena_, BigInt(1)));
            if (!il->value.is_negative()) {
                const std::uint64_t n = il->value.to_u64();
                std::vector<ExprPtr> factors;
                factors.reserve(n);
                for (std::uint64_t k = 0U; k < n; ++k) {
                    ExprPtr term = k == 0U
                        ? x_arg
                        : static_cast<ExprPtr>(arena_.make<Sum>(std::vector<ExprPtr>{
                              x_arg, make_integer(arena_, BigInt(static_cast<std::int64_t>(k)))}));
                    factors.push_back(term);
                }
                ExprPtr prod = factors.size() == 1U
                    ? factors.front()
                    : static_cast<ExprPtr>(arena_.make<Product>(std::move(factors)));
                return simplify_expr(prod);
            }
            // n < 0: (x)_{-m} = 1/((x-1)*(x-2)*...*(x-m))
            BigInt mag = -il->value;
            const std::uint64_t m = mag.to_u64();
            std::vector<ExprPtr> denom_factors;
            denom_factors.reserve(m);
            for (std::uint64_t k = 1U; k <= m; ++k) {
                ExprPtr neg_k = arena_.make<Unary>(UnaryOp::Neg,
                    make_integer(arena_, BigInt(static_cast<std::int64_t>(k))));
                denom_factors.push_back(arena_.make<Sum>(std::vector<ExprPtr>{x_arg, neg_k}));
            }
            ExprPtr denom = denom_factors.size() == 1U
                ? denom_factors.front()
                : static_cast<ExprPtr>(arena_.make<Product>(std::move(denom_factors)));
            return simplify_expr(arena_.make<Binary>(BinaryOp::Div,
                make_integer(arena_, BigInt(1)), denom));
        }
    }

    // Erf — odd function
    if (op == BuiltinOp::Erf && args.size() == 1U) {
        if (is_zero_expr(args.front()))
            return traced_result(RuleId::SimplifyErfZero, target_before, make_integer(arena_, BigInt(0)));
        if (const auto* un = expr_cast<Unary>(args.front()); un && un->op == UnaryOp::Neg) {
            ExprPtr inner_erf = arena_.make<FuncCall>(BuiltinOp::Erf, std::vector<ExprPtr>{un->operand});
            return simplify_expr(arena_.make<Unary>(UnaryOp::Neg, inner_erf));
        }
    }

    // Zeta closed-form via Bernoulli numbers (HC-003 resolved 2026-05-16)
    if (op == BuiltinOp::Zeta && args.size() == 1U) {
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            const BigInt& n = il->value;
            if (n.is_zero())
                return ok(arena_.make<RationalLit>(BigInt(-1), BigInt(2)));
            if (n.is_negative()) {
                if (n.bit_length() > 30)
                    // F0.8-MIGRATED
                    return make_unimplemented<ExprPtr>(
                        "symbolic", "simplify_funcall_special",
                        "Zeta(n) with negative integer n, |n| > 2^30",
                        error::reason_codes::SYMBOLIC_ZETA_OVERFLOW,
                        "Bernoulli number computation requires |n| ≤ 2^30; use ctx.max_bernoulli_index",
                        "F1.x");
                const std::uint64_t magu = (-n).to_u64();
                if (magu % 2U == 0U) return ok(make_integer(arena_, BigInt(0)));
                const unsigned int two_k = static_cast<unsigned int>(magu + 1U);
                Rational b2k = cas::numtheory::bernoulli_number(two_k);
                Rational result = -b2k / Rational(BigInt(static_cast<std::int64_t>(two_k)));
                if (result.denominator() == BigInt(1))
                    return ok(make_integer(arena_, result.numerator()));
                return ok(arena_.make<RationalLit>(result.numerator(), result.denominator()));
            }
            if (n > BigInt(0)) {
                if (n.bit_length() > 30)
                    // F0.8-MIGRATED
                    return make_unimplemented<ExprPtr>(
                        "symbolic", "simplify_funcall_special",
                        "Zeta(n) with positive even integer n > 2^30",
                        error::reason_codes::SYMBOLIC_ZETA_OVERFLOW,
                        "Bernoulli number computation requires n ≤ 2^30; use ctx.max_bernoulli_index",
                        "F1.x");
                const std::uint64_t nu = n.to_u64();
                if (nu % 2U != 0U) return ok(target_before);
                const unsigned int two_k = static_cast<unsigned int>(nu);
                const unsigned int k = two_k / 2U;
                BigInt two_pow(1);
                for (unsigned int i = 0U; i < two_k - 1U; ++i) two_pow *= BigInt(2);
                BigInt fact(1);
                for (unsigned int i = 2U; i <= two_k; ++i)
                    fact *= BigInt(static_cast<std::int64_t>(i));
                Rational b2k = cas::numtheory::bernoulli_number(two_k);
                Rational coeff = b2k * Rational(two_pow, fact);
                if ((k % 2U) == 0U) coeff = -coeff;
                ExprPtr pi_pow = arena_.make<Binary>(
                    BinaryOp::Pow,
                    arena_.make<Constant>(MathConstant::Pi),
                    arena_.make<IntegerLit>(n));
                ExprPtr coeff_expr;
                if (coeff.denominator() == BigInt(1)) {
                    if (coeff.numerator() == BigInt(1))
                        return simplify_expr(pi_pow);
                    coeff_expr = make_integer(arena_, coeff.numerator());
                } else {
                    coeff_expr = arena_.make<RationalLit>(coeff.numerator(), coeff.denominator());
                }
                return simplify_expr(arena_.make<Binary>(BinaryOp::Mul, coeff_expr, pi_pow));
            }
        }
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

} // namespace cas::symbolic::detail
