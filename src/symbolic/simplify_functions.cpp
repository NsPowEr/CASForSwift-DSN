#include "simplify_impl.hpp"
#include "cas/linalg/Matrix.hpp"
#include "cas/numeric.hpp"
#include "cas/algebra.hpp"
#include "../algebra/polynomial_internal.hpp"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>

namespace cas::symbolic::detail {

// ── P2-005: Trig angle reduction ─────────────────────────────────────────────

// Extract rational r such that arg = r * π, or nullopt if arg is not that form.
[[nodiscard]] static std::optional<Rational> try_extract_pi_coefficient(ExprPtr arg) {
    if (!arg) return std::nullopt;

    if (const auto* c = expr_cast<Constant>(arg); c && c->value == MathConstant::Pi)
        return Rational(BigInt(1));

    if (const auto* u = expr_cast<Unary>(arg); u && u->op == UnaryOp::Neg) {
        if (const auto* c = expr_cast<Constant>(u->operand); c && c->value == MathConstant::Pi)
            return Rational(BigInt(-1));
    }

    // Binary(Mul, scalar, Pi) or Binary(Mul, Pi, scalar)
    if (const auto* bin = expr_cast<Binary>(arg); bin && bin->op == BinaryOp::Mul) {
        LiteralRational lr;
        if (const auto* c = expr_cast<Constant>(bin->right); c && c->value == MathConstant::Pi) {
            auto ex = try_get_exact_rational(bin->left, lr);
            if (ex.is_ok() && ex.value()) return lr.value;
        }
        if (const auto* c = expr_cast<Constant>(bin->left); c && c->value == MathConstant::Pi) {
            auto ex = try_get_exact_rational(bin->right, lr);
            if (ex.is_ok() && ex.value()) return lr.value;
        }
    }

    // Product{scalar, Pi} (canonical form after simplification)
    if (const auto* prod = expr_cast<Product>(arg)) {
        ExprPtr scalar_factor = nullptr;
        bool found_pi = false;
        int scalar_count = 0;
        for (ExprPtr f : prod->factors) {
            if (const auto* c = expr_cast<Constant>(f); c && c->value == MathConstant::Pi) {
                found_pi = true;
            } else {
                scalar_factor = f;
                ++scalar_count;
            }
        }
        if (found_pi && scalar_count == 1 && scalar_factor) {
            LiteralRational lr;
            auto ex = try_get_exact_rational(scalar_factor, lr);
            if (ex.is_ok() && ex.value()) return lr.value;
        }
    }

    return std::nullopt;
}

// Reduce rational r to [0, 2) by subtracting multiples of 2.
[[nodiscard]] static Rational reduce_to_period(Rational r) {
    // floor(r/2) = floor(num / (2*den)) with floor toward -inf
    const BigInt two_den = BigInt(2) * r.denominator();
    BigInt q = r.numerator() / two_den;
    BigInt rem = r.numerator() - q * two_den;
    if (rem < BigInt(0)) {
        q = q - BigInt(1);
        rem = rem + two_den;
    }
    return Rational(rem, r.denominator());
}

// Return sin(ref * π) for ref ∈ [0, 1/2], or nullptr if not in table.
// Handles denominators dividing 12: 0, π/6, π/4, π/3, π/2.
[[nodiscard]] static ExprPtr sin_ref_value(Rational ref, AstArena& arena) {
    const Rational zero(BigInt(0));
    const Rational one_sixth(BigInt(1), BigInt(6));
    const Rational one_quarter(BigInt(1), BigInt(4));
    const Rational one_third(BigInt(1), BigInt(3));
    const Rational one_half(BigInt(1), BigInt(2));

    if (ref == zero)       return make_integer(arena, BigInt(0));
    if (ref == one_sixth)  return make_rational(arena, Rational(BigInt(1), BigInt(2)));
    if (ref == one_quarter)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{make_integer(arena, BigInt(2))}),
            make_integer(arena, BigInt(2)));
    if (ref == one_third)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{make_integer(arena, BigInt(3))}),
            make_integer(arena, BigInt(2)));
    if (ref == one_half)   return make_integer(arena, BigInt(1));
    return nullptr;
}

// Return cos(ref * π) for ref ∈ [0, 1/2], or nullptr if not in table.
[[nodiscard]] static ExprPtr cos_ref_value(Rational ref, AstArena& arena) {
    const Rational zero(BigInt(0));
    const Rational one_sixth(BigInt(1), BigInt(6));
    const Rational one_quarter(BigInt(1), BigInt(4));
    const Rational one_third(BigInt(1), BigInt(3));
    const Rational one_half(BigInt(1), BigInt(2));

    if (ref == zero)       return make_integer(arena, BigInt(1));
    if (ref == one_sixth)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{make_integer(arena, BigInt(3))}),
            make_integer(arena, BigInt(2)));
    if (ref == one_quarter)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{make_integer(arena, BigInt(2))}),
            make_integer(arena, BigInt(2)));
    if (ref == one_third)  return make_rational(arena, Rational(BigInt(1), BigInt(2)));
    if (ref == one_half)   return make_integer(arena, BigInt(0));
    return nullptr;
}

// Negate an expression: if it's already 0 return 0, else wrap in Unary(Neg).
[[nodiscard]] static ExprPtr negate_expr(ExprPtr e, AstArena& arena) {
    if (is_zero_expr(e)) return e;
    LiteralRational lr;
    auto ex = try_get_exact_rational(e, lr);
    if (ex.is_ok() && ex.value()) return make_rational(arena, -lr.value);
    return arena.make<Unary>(UnaryOp::Neg, e);
}

// Compute sin(r * π) or cos(r * π) where r is any rational, using quadrant reduction.
// Returns nullptr if reference angle not in table.
[[nodiscard]] static ExprPtr trig_exact_at_pi_multiple(
    BuiltinOp func_id, Rational r, AstArena& arena) {

    const Rational half(BigInt(1), BigInt(2));
    const Rational one(BigInt(1));
    const Rational three_halves(BigInt(3), BigInt(2));

    Rational red = reduce_to_period(r); // red ∈ [0, 2)
    // Quadrant:
    //   Q1 [0, 1/2]:   sin = +sin_ref(red),     cos = +cos_ref(red)
    //   Q2 (1/2, 1]:   sin = +sin_ref(1 - red), cos = -cos_ref(1 - red)
    //   Q3 (1, 3/2]:   sin = -sin_ref(red - 1), cos = -cos_ref(red - 1)
    //   Q4 (3/2, 2):   sin = -sin_ref(2 - red), cos = +cos_ref(2 - red)

    Rational ref;
    int sin_sign = 1;
    int cos_sign = 1;

    if (red <= half) {
        ref = red;
    } else if (red <= one) {
        ref = one - red;
        cos_sign = -1;
    } else if (red <= three_halves) {
        ref = red - one;
        sin_sign = -1;
        cos_sign = -1;
    } else {
        ref = Rational(BigInt(2)) - red;
        sin_sign = -1;
    }

    ExprPtr val = (func_id == BuiltinOp::Sin)
        ? sin_ref_value(ref, arena)
        : cos_ref_value(ref, arena);
    if (!val) return nullptr;

    int sign = (func_id == BuiltinOp::Sin) ? sin_sign : cos_sign;
    if (sign < 0) return negate_expr(val, arena);
    return val;
}

// ─────────────────────────────────────────────────────────────────────────────

// Extract perfect-square factor from n: returns {k, m} where n = k²·m, m squarefree.
[[nodiscard]] static std::pair<BigInt, BigInt> extract_square_factor(BigInt n) {
    BigInt k(1);
    BigInt i(2);
    while (i * i <= n) {
        BigInt i2 = i * i;
        while ((n % i2).is_zero()) {
            k = k * i;
            n = n / i2;
        }
        i = i + BigInt(1);
    }
    return {k, n};
}

// Simplify sqrt(r) for rational r: extract perfect-square factors.
// Returns k*sqrt(m) where k is rational, m is squarefree integer.
// Requires r >= 0.
[[nodiscard]] static Result<ExprPtr> simplify_rational_sqrt(const Rational& r, AstArena& arena) {
    const BigInt& p = r.numerator();
    const BigInt& q = r.denominator();
    if (p.is_zero()) return ok(arena.make<IntegerLit>(BigInt(0)));
    // r must be normalized and positive (caller ensures this)
    auto [p_out, p_rem] = extract_square_factor(p);
    auto [q_out, q_rem] = extract_square_factor(q);
    // sqrt(r) = (p_out/q_out) * sqrt(p_rem/q_rem)
    //         = (p_out/(q_out*q_rem)) * sqrt(p_rem*q_rem)
    // Since r is normalized, gcd(p,q)=1 => gcd(p_rem,q_rem)=1 => p_rem*q_rem is squarefree.
    BigInt final_radicand = p_rem * q_rem;
    Rational coeff(p_out, q_out * q_rem);

    ExprPtr coeff_expr;
    if (coeff.denominator() == BigInt(1)) {
        coeff_expr = arena.make<IntegerLit>(coeff.numerator());
    } else {
        coeff_expr = arena.make<RationalLit>(coeff.numerator(), coeff.denominator());
    }

    if (final_radicand == BigInt(1)) return ok(coeff_expr);

    ExprPtr sqrt_expr = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{arena.make<IntegerLit>(final_radicand)});
    if (coeff.numerator() == BigInt(1) && coeff.denominator() == BigInt(1))
        return ok(sqrt_expr);
    return ok(arena.make<Binary>(BinaryOp::Mul, coeff_expr, sqrt_expr));
}

[[nodiscard]] static BigInt integer_sqrt(const BigInt& n) {
    if (n.is_zero()) return BigInt(0);
    static const BigInt one(1);
    if (n == one) return one;
    
    // Initial guess: bit length / 2
    BigInt x = one.shift_left_bits((n.bit_length() + 1) / 2);
    while (true) {
        BigInt y = (x + n / x) / BigInt(2);
        if (y >= x) return x;
        x = std::move(y);
    }
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const FuncCall& node) {
    std::vector<ExprPtr> args;
    args.reserve(node.args.size());
    for (std::size_t i = 0; i < node.args.size(); ++i) {
        auto simplify_arg = [&]() -> Result<ExprPtr> {
            if (trace_enabled_) {
                std::vector<ExprPtr> current_args = node.args;
                for (std::size_t j = 0; j < args.size(); ++j) current_args[j] = args[j];
                ScopedFrame frame(*this, [this, current_args = std::move(current_args), i, name = node.name](ExprPtr value) mutable {
                    current_args[i] = value;
                    return arena_.make<FuncCall>(name, std::move(current_args));
                });
                return simplify_expr(node.args[i]);
            }
            return simplify_expr(node.args[i]);
        };
        auto res = simplify_arg();
        if (res.is_error()) return res;
        args.push_back(res.value());
    }

    const ExprPtr target_before = expr_ptr_sequence_identical(args, node.args) ? original : (trace_enabled_ ? arena_.make<FuncCall>(node.func_id, args) : ExprPtr{});

    if (rewrite_provider_ != nullptr && may_rewrite_function_call(node.func_id, args)) {
        ExprPtr rewrite_target = expr_ptr_sequence_identical(args, node.args) ? original : arena_.make<FuncCall>(node.func_id, args);
        auto rewritten = rewrite_provider_->try_rewrite(rewrite_target, arena_, assumptions_, context_);
        if (rewritten.is_ok() && rewritten.value() != rewrite_target) {
            append_trace(RuleId::RewriteProviderApplied, rewrite_target, rewritten.value());
            return simplify_expr(rewritten.value());
        }
    }

    if (node.func_id == BuiltinOp::Sin && args.size() == 1U) {
        if (is_zero_expr(args.front())) return traced_result(RuleId::SimplifySinZero, target_before, make_integer(arena_, BigInt(0)));
        if (is_constant_expr(args.front(), MathConstant::Pi)) return traced_result(RuleId::SimplifySinPi, target_before, make_integer(arena_, BigInt(0)));
        
        // sin(asin(x)) -> x
        if (const auto* call = expr_cast<FuncCall>(args.front()); call && call->func_id == BuiltinOp::Asin) {
            return ok(call->args[0]);
        }
        
        // P2-005: sin(r*π) exact for r = k/12
        if (auto coeff = try_extract_pi_coefficient(args.front())) {
            if (ExprPtr val = trig_exact_at_pi_multiple(BuiltinOp::Sin, *coeff, arena_)) {
                auto simplified = simplify_expr(val);
                if (simplified.is_ok()) return traced_result(RuleId::Unknown, target_before, simplified.value());
            }
        }
    }
    if (node.func_id == BuiltinOp::Cos && args.size() == 1U) {
        if (is_zero_expr(args.front())) return traced_result(RuleId::SimplifyCosZero, target_before, make_integer(arena_, BigInt(1)));
        if (is_constant_expr(args.front(), MathConstant::Pi)) return traced_result(RuleId::SimplifyCosPi, target_before, make_integer(arena_, BigInt(-1)));

        // cos(acos(x)) -> x
        if (const auto* call = expr_cast<FuncCall>(args.front()); call && call->func_id == BuiltinOp::Acos) {
            return ok(call->args[0]);
        }

        // P2-005: cos(r*π) exact for r = k/12
        if (auto coeff = try_extract_pi_coefficient(args.front())) {
            if (ExprPtr val = trig_exact_at_pi_multiple(BuiltinOp::Cos, *coeff, arena_)) {
                auto simplified = simplify_expr(val);
                if (simplified.is_ok()) return traced_result(RuleId::Unknown, target_before, simplified.value());
            }
        }
    }
    if (node.func_id == BuiltinOp::Tan && args.size() == 1U) {
        // tan(atan(x)) -> x
        if (const auto* call = expr_cast<FuncCall>(args.front()); call && call->func_id == BuiltinOp::Atan) {
            return ok(call->args[0]);
        }
    }

    if (node.func_id == BuiltinOp::Exp && args.size() == 1U) {
        if (is_zero_expr(args.front())) return traced_result(RuleId::SimplifyExpZero, target_before, make_integer(arena_, BigInt(1)));
        if (is_one_expr(args.front())) return traced_result(RuleId::SimplifyExpOne, target_before, make_constant(arena_, MathConstant::E));
        
        // exp(ln(x)) -> x
        if (const auto* ln_call = expr_cast<FuncCall>(args.front()); ln_call && ln_call->func_id == BuiltinOp::Ln && ln_call->args.size() == 1U) {
            return ok(ln_call->args[0]);
        }

        if (is_constant_expr(args.front(), MathConstant::Infinity)) return traced_result(RuleId::Unknown, target_before, make_constant(arena_, MathConstant::Infinity));
        if (expr_is<Unary>(args.front()) && expr_ref<Unary>(args.front()).op == UnaryOp::Neg && is_constant_expr(expr_ref<Unary>(args.front()).operand, MathConstant::Infinity)) {
            return traced_result(RuleId::Unknown, target_before, make_integer(arena_, BigInt(0)));
        }
        if (const auto* sum = expr_cast<Sum>(args.front())) {
            std::vector<ExprPtr> factors;
            for (ExprPtr term : sum->terms) factors.push_back(arena_.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{term}));
            auto rewritten = simplify_product_factors(factors, arena_.make<Product>(factors));
            if (rewritten.is_ok()) { append_trace(RuleId::SimplifyExpSum, target_before, rewritten.value()); return rewritten; }
        }
    }

    if (node.func_id == BuiltinOp::Ln && args.size() == 1U) {
        if (is_zero_expr(args.front())) return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "ln(0) is undefined"));
        if (is_one_expr(args.front())) return traced_result(RuleId::SimplifyLnOne, target_before, make_integer(arena_, BigInt(0)));
        if (is_constant_expr(args.front(), MathConstant::E)) return traced_result(RuleId::SimplifyLnE, target_before, make_integer(arena_, BigInt(1)));
        if (is_constant_expr(args.front(), MathConstant::Infinity)) return traced_result(RuleId::Unknown, target_before, make_constant(arena_, MathConstant::Infinity));
        
        // ln(e^x) -> x
        if (const auto* power = expr_cast<Binary>(args.front()); power != nullptr && power->op == BinaryOp::Pow && is_constant_expr(power->left, MathConstant::E)) {
            return traced_result(RuleId::SimplifyLnExp, target_before, power->right);
        }
        // ln(exp(x)) -> x
        if (const auto* exp_call = expr_cast<FuncCall>(args.front()); exp_call && exp_call->func_id == BuiltinOp::Exp && exp_call->args.size() == 1U) {
            return traced_result(RuleId::SimplifyLnExp, target_before, exp_call->args[0]);
        }

        // ln(a*b) -> ln(a) + ln(b) for a,b > 0
        if (const auto* prod = expr_cast<Product>(args.front())) {
            bool all_pos = true;
            for (auto f : prod->factors) if (!is_known_positive(f)) { all_pos = false; break; }
            if (all_pos) {
                std::vector<ExprPtr> ln_factors;
                for (auto f : prod->factors) {
                    auto res = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{f}));
                    if (res.is_error()) return res;
                    ln_factors.push_back(res.value());
                }
                return simplify_expr(arena_.make<Sum>(std::move(ln_factors)));
            }
        }
        // ln(sqrt(x)) = (1/2)*ln(x)  — valido per x > 0, identità esatta
        if (const auto* sqrt_call = expr_cast<FuncCall>(args.front());
            sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U) {
            ExprPtr half = make_rational(arena_, Rational(BigInt(1), BigInt(2)));
            ExprPtr ln_inner = arena_.make<FuncCall>(BuiltinOp::Ln, sqrt_call->args);
            return simplify_expr(arena_.make<Binary>(BinaryOp::Mul, half, ln_inner));
        }
        // Branch cut principal value: ln(-x) = ln(x) + I*pi for x > 0
        if (const auto* neg = expr_cast<Unary>(args.front()); neg != nullptr && neg->op == UnaryOp::Neg) {
            ExprPtr inner = neg->operand;
            bool inner_pos = is_known_positive(inner) || is_constant_expr(inner, MathConstant::E);
            if (!inner_pos) {
                LiteralRational rat;
                auto ex = try_get_exact_rational(inner, rat);
                if (ex.is_ok() && ex.value() && !rat.value.numerator().is_negative() && !rat.value.numerator().is_zero()) {
                    inner_pos = true;
                }
            }
            if (inner_pos) {
                ExprPtr ln_inner;
                if (is_one_expr(inner)) {
                    ln_inner = make_integer(arena_, BigInt(0));
                } else {
                    auto r = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{inner}));
                    if (r.is_error()) return r;
                    ln_inner = r.value();
                }
                ExprPtr i_pi = arena_.make<Binary>(BinaryOp::Mul,
                    make_constant(arena_, MathConstant::I),
                    make_constant(arena_, MathConstant::Pi));
                return simplify_expr(arena_.make<Binary>(BinaryOp::Add, ln_inner, i_pi));
            }
        }
    }
    if (node.func_id == BuiltinOp::Asin && args.size() == 1U) {
        // asin(sin(x)) -> x if -pi/2 <= x <= pi/2
        if (const auto* call = expr_cast<FuncCall>(args.front()); call && call->func_id == BuiltinOp::Sin) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi_2 = arena_.make<Binary>(BinaryOp::Div, arena_.make<Constant>(MathConstant::Pi), make_integer(arena_, 2));
                ExprPtr neg_pi_2 = arena_.make<Unary>(UnaryOp::Neg, pi_2);
                if (assumptions_->is_greater_equal(x, neg_pi_2) && assumptions_->is_greater_equal(pi_2, x)) {
                    return ok(x);
                }
            }
        }
    }
    if (node.func_id == BuiltinOp::Acos && args.size() == 1U) {
        // acos(cos(x)) -> x if 0 <= x <= pi
        if (const auto* call = expr_cast<FuncCall>(args.front()); call && call->func_id == BuiltinOp::Cos) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi = arena_.make<Constant>(MathConstant::Pi);
                ExprPtr zero = make_integer(arena_, 0);
                if (assumptions_->is_greater_equal(x, zero) && assumptions_->is_greater_equal(pi, x)) {
                    return ok(x);
                }
            }
        }
    }
    if (node.func_id == BuiltinOp::Atan && args.size() == 1U) {
        // atan(tan(x)) -> x if -pi/2 < x < pi/2
        if (const auto* call = expr_cast<FuncCall>(args.front()); call && call->func_id == BuiltinOp::Tan) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi_2 = arena_.make<Binary>(BinaryOp::Div, arena_.make<Constant>(MathConstant::Pi), make_integer(arena_, 2));
                ExprPtr neg_pi_2 = arena_.make<Unary>(UnaryOp::Neg, pi_2);
                if (assumptions_->is_greater(x, neg_pi_2) && assumptions_->is_greater(pi_2, x)) {
                    return ok(x);
                }
            }
        }
        // Special values: atan(0)=0, atan(1)=pi/4, atan(-1)=-pi/4,
        //                 atan(sqrt(3))=pi/3, atan(1/sqrt(3))=pi/6
        if (is_zero_expr(args.front())) {
            return ok(make_integer(arena_, BigInt(0)));
        }
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            if (il->value == BigInt(1)) {
                return simplify_expr(arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi),
                    make_integer(arena_, BigInt(4))));
            }
            if (il->value == BigInt(-1)) {
                return simplify_expr(arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Binary>(BinaryOp::Div,
                        arena_.make<Constant>(MathConstant::Pi),
                        make_integer(arena_, BigInt(4)))));
            }
        }
        // atan is odd: atan(-x) = -atan(x)
        if (const auto* un = expr_cast<Unary>(args.front()); un && un->op == UnaryOp::Neg) {
            ExprPtr inner = arena_.make<FuncCall>(
                BuiltinOp::Atan, std::vector<ExprPtr>{un->operand});
            return simplify_expr(arena_.make<Unary>(UnaryOp::Neg, inner));
        }
    }
    if (node.func_id == BuiltinOp::Sqrt && args.size() == 1U) {
        // Denesting sqrt(a + b*sqrt(c))
        if (const auto* sum = expr_cast<Sum>(args.front()); sum && sum->terms.size() == 2) {
            LiteralRational rat_a, rat_b, rat_c;
            ExprPtr a_ptr = nullptr, b_ptr = nullptr, c_ptr = nullptr;
            
            // Try to find 'a' and 'b*sqrt(c)'
            for (auto term : sum->terms) {
                if (auto ex = try_get_exact_rational(term, rat_a); ex.is_ok() && ex.value()) {
                    a_ptr = term;
                } else if (const auto* prod = expr_cast<Product>(term)) {
                    // Look for rational coefficient and a single sqrt factor
                    Rational b_coeff(1);
                    ExprPtr c_val = nullptr;
                    bool found_sqrt = false;
                    for (ExprPtr f : prod->factors) {
                        LiteralRational lr;
                        if (auto ex = try_get_exact_rational(f, lr); ex.is_ok() && ex.value()) {
                            b_coeff *= lr.value;
                        } else if (const auto* sqrt_c = expr_cast<FuncCall>(f); sqrt_c && sqrt_c->func_id == BuiltinOp::Sqrt && !found_sqrt) {
                            c_val = sqrt_c->args[0];
                            found_sqrt = true;
                        } else {
                            // Other factors found, not a simple b*sqrt(c)
                            found_sqrt = false;
                            break;
                        }
                    }
                    if (found_sqrt && c_val) {
                        LiteralRational lr_c;
                        if (auto ex_c = try_get_exact_rational(c_val, lr_c); ex_c.is_ok() && ex_c.value()) {
                            rat_b.value = b_coeff;
                            rat_c.value = lr_c.value;
                            b_ptr = make_rational(arena_, b_coeff);
                            c_ptr = c_val;
                        }
                    }
                }
            }

            if (a_ptr && b_ptr && c_ptr) {
                // sqrt(a + b*sqrt(c)) = sqrt((a + sqrt(a^2 - b^2*c))/2) + sign(b)*sqrt((a - sqrt(a^2 - b^2*c))/2)
                Rational a = rat_a.value;
                Rational b = rat_b.value;
                Rational c = rat_c.value;
                Rational discriminant = a*a - b*b*c;
                if (discriminant >= Rational(0)) {
                    BigInt d_num = discriminant.numerator();
                    BigInt d_den = discriminant.denominator();
                    BigInt s_num = integer_sqrt(d_num);
                    BigInt s_den = integer_sqrt(d_den);
                    if (s_num * s_num == d_num && s_den * s_den == d_den) {
                        Rational s(s_num, s_den);
                        Rational x = (a + s) / Rational(2);
                        Rational y = (a - s) / Rational(2);
                        
                        auto sqrt_x = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{make_rational(arena_, x)}));
                        auto sqrt_y = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{make_rational(arena_, y)}));
                        
                        if (sqrt_x.is_ok() && sqrt_y.is_ok()) {
                            ExprPtr res;
                            if (b >= Rational(0)) {
                                res = arena_.make<Sum>(std::vector<ExprPtr>{sqrt_x.value(), sqrt_y.value()});
                            } else {
                                res = arena_.make<Binary>(BinaryOp::Sub, sqrt_x.value(), sqrt_y.value());
                            }
                            return simplify_expr(res);
                        }
                    }
                }
            }
        }

        LiteralRational rat;
        auto exact = try_get_exact_rational(args.front(), rat);
        if (exact.is_error()) return fail<ExprPtr>(exact.error());
        if (exact.is_ok() && exact.value()) {
            if (rat.value.numerator().is_zero()) {
                return traced_result(RuleId::Unknown, target_before, make_integer(arena_, BigInt(0)));
            }
            if (rat.value == Rational(BigInt(1))) {
                return traced_result(RuleId::Unknown, target_before, make_integer(arena_, BigInt(1)));
            }
            
            if (rat.value.numerator().is_negative()) {
                auto pos_rat = -rat.value;
                auto sqrt_pos = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{make_rational(arena_, pos_rat)}));
                if (sqrt_pos.is_ok()) {
                    auto product = simplify_expr(arena_.make<Binary>(
                        BinaryOp::Mul,
                        arena_.make<Constant>(MathConstant::I),
                        sqrt_pos.value()));
                    if (product.is_error()) return product;
                    return traced_result(RuleId::Unknown, target_before, product.value());
                }
            }

            // Valutazione quadrati perfetti
            auto num_sqrt = integer_sqrt(rat.value.numerator());
            auto den_sqrt = integer_sqrt(rat.value.denominator());
            if (num_sqrt * num_sqrt == rat.value.numerator() && den_sqrt * den_sqrt == rat.value.denominator()) {
                return traced_result(RuleId::Unknown, target_before, make_rational(arena_, Rational(num_sqrt, den_sqrt)));
            }
            // Estrazione fattori quadrati parziali: sqrt(k²·m) → k·sqrt(m)
            auto denested = simplify_rational_sqrt(rat.value, arena_);
            if (denested.is_ok()) {
                return traced_result(RuleId::Unknown, target_before, denested.value());
            }
        }

        if (is_known_negative(args.front())) {
            auto negated_arg = simplify_expr(arena_.make<Unary>(UnaryOp::Neg, args.front()));
            if (negated_arg.is_ok()) {
                auto sqrt_pos = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{negated_arg.value()}));
                if (sqrt_pos.is_ok()) {
                    ExprPtr res = arena_.make<Binary>(BinaryOp::Mul, arena_.make<Constant>(MathConstant::I), sqrt_pos.value());
                    return traced_result(RuleId::Unknown, target_before, simplify_expr(res).value());
                }
            }
        }
        if (const auto* power = expr_cast<Binary>(args.front()); power != nullptr && power->op == BinaryOp::Pow) {
            if (auto exp = try_get_integer_exponent(power->right); exp.has_value() && *exp == BigInt(2)) {
                if (is_known_nonnegative(power->left)) { append_assumption(target_before); return traced_result(RuleId::SimplifySqrtSquare, target_before, power->left); }
                return traced_result(RuleId::SimplifySqrtSquare, target_before, arena_.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{power->left}));
            }
        }
        // sqrt(sqrt(x)) -> x^(1/4)
        if (const auto* inner = expr_cast<FuncCall>(args.front()); inner && inner->func_id == BuiltinOp::Sqrt) {
            return simplify_expr(arena_.make<Binary>(BinaryOp::Pow, inner->args[0], make_rational(arena_, Rational(BigInt(1), BigInt(4)))));
        }
    }

    if (node.func_id == BuiltinOp::Sin || node.func_id == BuiltinOp::Cos) {
    }


    if (node.func_id == BuiltinOp::Gamma && args.size() == 1U) {
        // Gamma(n) = (n-1)! for positive integers n >= 1
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            if (il->value > BigInt(0)) {
                BigInt n = il->value;
                BigInt result(1);
                for (BigInt k(1); k < n; k += BigInt(1)) {
                    result *= k;
                }
                return ok(make_integer(arena_, result));
            }
        }
        // Half-integer Gamma values:  Gamma(1/2) = sqrt(pi),
        // Gamma(n + 1/2) = (2n)! / (4^n * n!) * sqrt(pi)   for non-negative integer n.
        // Gamma(-n + 1/2) computed via recursion Gamma(z) = Gamma(z+1)/z down to Gamma(1/2).
        if (const auto* rl = expr_cast<RationalLit>(args.front())) {
            if (rl->denominator == BigInt(2)) {
                // arg = num/2 where num is odd (otherwise it would simplify to an integer).
                BigInt num = rl->numerator;
                // Reduce to Gamma(1/2) by recursion:
                //   Gamma(z+1) = z * Gamma(z)  =>  Gamma(z) = Gamma(z+1)/z, for z != 0.
                // Walk from arg back / forward to 1/2, accumulating multiplicative factors.
                Rational z(num, BigInt(2));
                const Rational half(BigInt(1), BigInt(2));
                // Build accumulator as a rational factor f, then multiply by sqrt(pi).
                Rational factor(BigInt(1));
                bool ok_chain = true;
                int safety = 0;
                // HARDCODE-OF-PASSAGE: HC-002 — Gamma half-integer recursion bound.
                // Fix: esporre ctx.max_gamma_recursion(). Vedi HARDCODE_LEDGER.md.
                const int safety_max = 1024;
                while (z != half && safety++ < safety_max) {
                    if (z > half) {
                        // Gamma(z) = (z-1) * Gamma(z-1)
                        Rational z_minus_1 = z - Rational(BigInt(1));
                        factor = factor * z_minus_1;
                        z = z_minus_1;
                    } else {
                        // z < 1/2 (still half-integer):  Gamma(z) = Gamma(z+1) / z, z != 0.
                        if (z.numerator().is_zero()) { ok_chain = false; break; }
                        factor = factor / z;
                        z = z + Rational(BigInt(1));
                    }
                }
                if (ok_chain && safety < safety_max) {
                    ExprPtr sqrt_pi = arena_.make<FuncCall>(
                        BuiltinOp::Sqrt,
                        std::vector<ExprPtr>{arena_.make<Constant>(MathConstant::Pi)});
                    if (factor == Rational(BigInt(1))) {
                        return simplify_expr(sqrt_pi);
                    }
                    ExprPtr factor_expr = factor.denominator() == BigInt(1)
                        ? static_cast<ExprPtr>(arena_.make<IntegerLit>(factor.numerator()))
                        : static_cast<ExprPtr>(arena_.make<RationalLit>(factor.numerator(), factor.denominator()));
                    return simplify_expr(
                        arena_.make<Product>(std::vector<ExprPtr>{factor_expr, sqrt_pi}));
                }
            }
        }
        // Functional equation:  Gamma(z + n)  -> (z)(z+1)...(z+n-1) * Gamma(z)
        // when the argument is a Sum and the integer summand can be stripped.
        if (const auto* sum = expr_cast<Sum>(args.front()); sum && sum->terms.size() >= 2U) {
            // Find a single IntegerLit term we can strip.
            int idx = -1;
            BigInt shift(0);
            for (std::size_t i = 0; i < sum->terms.size(); ++i) {
                if (const auto* il = expr_cast<IntegerLit>(sum->terms[i])) {
                    if (idx == -1) { idx = static_cast<int>(i); shift = il->value; }
                    else { idx = -1; break; }  // multiple integer terms: leave alone
                }
            }
            if (idx >= 0 && !shift.is_zero()) {
                // Reconstruct z = sum without that term.
                std::vector<ExprPtr> rest;
                rest.reserve(sum->terms.size() - 1U);
                for (std::size_t i = 0; i < sum->terms.size(); ++i) {
                    if (static_cast<int>(i) != idx) rest.push_back(sum->terms[i]);
                }
                ExprPtr z = rest.size() == 1U ? rest.front()
                                              : arena_.make<Sum>(std::move(rest));
                std::vector<ExprPtr> factors;
                if (shift > BigInt(0)) {
                    // Gamma(z + n) = z*(z+1)*...*(z+n-1) * Gamma(z), n positive integer.
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
                // shift < 0: Gamma(z - n) = Gamma(z) / [(z-1)(z-2)...(z-n)]
                BigInt n = -shift;
                std::vector<ExprPtr> denom_factors;
                for (BigInt k(1); k <= n; k += BigInt(1)) {
                    denom_factors.push_back(
                        arena_.make<Sum>(std::vector<ExprPtr>{
                            z,
                            arena_.make<Unary>(UnaryOp::Neg, arena_.make<IntegerLit>(k))}));
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

    // P3-005: Re, Im, Conj over a + b*I form
    // Returns {real_part, imag_part} from Sum/Product containing Constant(I)
    auto extract_complex = [&](ExprPtr expr) -> std::optional<std::pair<ExprPtr, ExprPtr>> {
        ExprPtr real_part = nullptr;
        ExprPtr imag_part = nullptr;
        ExprPtr zero = make_integer(arena_, BigInt(0));

        auto is_imaginary_unit = [](ExprPtr e) {
            const auto* c = expr_cast<Constant>(e);
            return c != nullptr && c->value == MathConstant::I;
        };

        auto extract_imag_factor = [&](ExprPtr e) -> ExprPtr {
            // e = I  → imag factor = 1
            if (is_imaginary_unit(e)) return make_integer(arena_, BigInt(1));
            // e = -I → imag factor = -1
            if (const auto* u = expr_cast<Unary>(e); u && u->op == UnaryOp::Neg && is_imaginary_unit(u->operand))
                return make_integer(arena_, BigInt(-1));
            // e = Product{..., I} → imag factor = Product{remaining}
            if (const auto* prod = expr_cast<Product>(e)) {
                std::vector<ExprPtr> non_i;
                bool found_i = false;
                for (ExprPtr f : prod->factors) {
                    if (!found_i && is_imaginary_unit(f)) { found_i = true; continue; }
                    non_i.push_back(f);
                }
                if (found_i) {
                    if (non_i.size() == 1U) return non_i[0];
                    if (non_i.empty()) return make_integer(arena_, BigInt(1));
                    return arena_.make<Product>(std::move(non_i));
                }
            }
            return nullptr;
        };

        if (const auto* sum = expr_cast<Sum>(expr)) {
            for (ExprPtr term : sum->terms) {
                ExprPtr maybe_b = extract_imag_factor(term);
                if (maybe_b) {
                    imag_part = maybe_b;
                } else {
                    real_part = term;
                }
            }
            if (!real_part) real_part = zero;
            if (!imag_part) return std::nullopt;
            return std::make_pair(real_part, imag_part);
        }
        // bare I
        ExprPtr maybe_b = extract_imag_factor(expr);
        if (maybe_b) return std::make_pair(zero, maybe_b);
        return std::nullopt;
    };

    if (args.size() == 1U) {
        auto parts = extract_complex(args.front());

        if (node.func_id == BuiltinOp::Re) {
            if (parts) return ok(parts->first);
            // Re(real) = real if no imaginary component detected
        }

        if (node.func_id == BuiltinOp::Im) {
            if (parts) return ok(parts->second);
        }

        if (node.func_id == BuiltinOp::Conj && parts) {
            // conj(a + b*I) = a - b*I
            ExprPtr neg_b_i = arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Unary>(UnaryOp::Neg, parts->second),
                arena_.make<Constant>(MathConstant::I)
            });
            auto conj_expr = arena_.make<Sum>(std::vector<ExprPtr>{parts->first, neg_b_i});
            auto simp = simplify_expr(conj_expr);
            if (simp.is_ok()) return simp;
        }
    }

    if (node.func_id == BuiltinOp::Erf && args.size() == 1U) {
        if (is_zero_expr(args.front())) return traced_result(RuleId::SimplifyErfZero, target_before, make_integer(arena_, BigInt(0)));
        // erf is odd:  erf(-x) = -erf(x)
        if (const auto* un = expr_cast<Unary>(args.front()); un && un->op == UnaryOp::Neg) {
            ExprPtr inner_erf = arena_.make<FuncCall>(
                BuiltinOp::Erf, std::vector<ExprPtr>{un->operand});
            auto neg = arena_.make<Unary>(UnaryOp::Neg, inner_erf);
            return simplify_expr(neg);
        }
    }
    // Zeta special values:
    //   zeta(0) = -1/2,  zeta(-2n) = 0 for n >= 1 (trivial zeros),
    //   zeta(-1) = -1/12, zeta(-3) = 1/120, zeta(-5) = -1/252  (via Bernoulli),
    //   zeta(2) = pi^2/6, zeta(4) = pi^4/90, zeta(6) = pi^6/945,
    //   zeta(8) = pi^8/9450, zeta(10) = pi^10/93555.
    // (Closed-form: zeta(2n) = (-1)^(n+1) * (2pi)^(2n) * B_{2n} / (2 (2n)!))
    if (node.func_id == BuiltinOp::Zeta && args.size() == 1U) {
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            const BigInt& n = il->value;
            if (n.is_zero()) {
                return ok(arena_.make<RationalLit>(BigInt(-1), BigInt(2)));
            }
            // Negative even integers: trivial zeros.
            if (n.is_negative()) {
                BigInt mag = -n;
                bool is_even = (mag.to_u64() % 2U) == 0U;
                if (is_even) {
                    return ok(make_integer(arena_, BigInt(0)));
                }
                // Negative odd: zeta(-(2k-1)) = -B_{2k} / (2k).  We only do small cases.
                std::uint64_t magu = mag.to_u64();
                if (magu == 1U) return ok(arena_.make<RationalLit>(BigInt(-1), BigInt(12)));
                if (magu == 3U) return ok(arena_.make<RationalLit>(BigInt(1), BigInt(120)));
                if (magu == 5U) return ok(arena_.make<RationalLit>(BigInt(-1), BigInt(252)));
                if (magu == 7U) return ok(arena_.make<RationalLit>(BigInt(1), BigInt(240)));
            }
            // Positive even integers: closed form via pi-power.
            // HARDCODE-OF-PASSAGE: HC-003 — Zeta lookup table per nu in {2,4,6,8,10,12}.
            // Viola REGOLA ZERO: scelta deliberata di lookup invece di Bernoulli generale.
            // Fix: esporre bernoulli_numbers() (oggi static in summation.cpp:63) come API
            // pubblica e usare formula closed-form ζ(2n) = (-1)^(n+1)·(2π)^(2n)·B_{2n}/(2(2n)!)
            // per nu arbitrario. Vedi HARDCODE_LEDGER.md.
            if (n > BigInt(0)) {
                std::uint64_t nu = n.to_u64();
                if (nu % 2U == 0U && nu <= 12U) {
                    // Small lookup of even-zeta denominators (numerator = 1).
                    BigInt denom;
                    switch (nu) {
                        case 2U:  denom = BigInt(6); break;
                        case 4U:  denom = BigInt(90); break;
                        case 6U:  denom = BigInt(945); break;
                        case 8U:  denom = BigInt(9450); break;
                        case 10U: denom = BigInt(93555); break;
                        case 12U: denom = BigInt::parse("638512875").value(); break;
                        default: denom = BigInt(0);
                    }
                    if (!denom.is_zero()) {
                        ExprPtr pi_pow = arena_.make<Binary>(
                            BinaryOp::Pow,
                            arena_.make<Constant>(MathConstant::Pi),
                            arena_.make<IntegerLit>(n));
                        // For zeta(12) the numerator differs; recheck:
                        //   zeta(2) = pi^2/6        num=1
                        //   zeta(4) = pi^4/90       num=1
                        //   zeta(6) = pi^6/945      num=1
                        //   zeta(8) = pi^8/9450     num=1
                        //   zeta(10) = pi^10/93555  num=1
                        //   zeta(12) = 691*pi^12/638512875   num=691
                        if (nu == 12U) {
                            ExprPtr num = arena_.make<IntegerLit>(BigInt(691));
                            ExprPtr prod = arena_.make<Product>(std::vector<ExprPtr>{num, pi_pow});
                            return simplify_expr(arena_.make<Binary>(
                                BinaryOp::Div, prod, arena_.make<IntegerLit>(denom)));
                        }
                        return simplify_expr(arena_.make<Binary>(
                            BinaryOp::Div, pi_pow, arena_.make<IntegerLit>(denom)));
                    }
                }
            }
        }
    }
    if (node.func_id == BuiltinOp::Abs && args.size() == 1U) {
        if (is_known_nonnegative(args.front())) {
            return traced_result(RuleId::SimplifyAbsPositive, target_before, args.front());
        }
        if (is_known_negative(args.front())) {
            auto neg = arena_.make<Unary>(UnaryOp::Neg, args.front());
            auto simplified = simplify_expr(neg);
            if (simplified.is_ok()) return traced_result(RuleId::SimplifyAbsNegative, target_before, simplified.value());
        }
        // abs(abs(x)) -> abs(x)
        if (const auto* inner = expr_cast<FuncCall>(args.front()); inner && inner->func_id == BuiltinOp::Abs) {
            return traced_result(RuleId::SimplifyAbsAbs, target_before, args.front());
        }
        // abs(-x) -> abs(x)
        if (const auto* unary = expr_cast<Unary>(args.front()); unary && unary->op == UnaryOp::Neg) {
            return simplify_expr(arena_.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{unary->operand}));
        }
        // abs(I) = 1
        if (const auto* c = expr_cast<Constant>(args.front()); c && c->value == MathConstant::I) {
            return ok(make_integer(arena_, BigInt(1)));
        }
        // abs(a + b*I) = sqrt(a^2 + b^2)
        {
            auto parts = extract_complex(args.front());
            if (parts) {
                ExprPtr a_sq = arena_.make<Binary>(BinaryOp::Pow, parts->first, make_integer(arena_, BigInt(2)));
                ExprPtr b_sq = arena_.make<Binary>(BinaryOp::Pow, parts->second, make_integer(arena_, BigInt(2)));
                ExprPtr sum = arena_.make<Sum>(std::vector<ExprPtr>{a_sq, b_sq});
                ExprPtr rooted = arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{sum});
                return simplify_expr(rooted);
            }
        }
    }
    // arg(z) — principal argument.  Only handles the cases we can fully decide:
    //   arg(0) = 0  (convention)
    //   arg(positive real) = 0
    //   arg(negative real) = pi
    //   arg(i)  = pi/2,   arg(-i) = -pi/2
    //   arg(a + b*i) numeric, a > 0:  atan(b/a)
    if (node.func_id == BuiltinOp::Arg && args.size() == 1U) {
        if (is_zero_expr(args.front())) return ok(make_integer(arena_, BigInt(0)));
        if (is_known_positive(args.front())) return ok(make_integer(arena_, BigInt(0)));
        if (is_known_negative(args.front())) {
            return ok(arena_.make<Constant>(MathConstant::Pi));
        }
        // arg(i) = pi/2
        if (const auto* c = expr_cast<Constant>(args.front()); c && c->value == MathConstant::I) {
            return simplify_expr(arena_.make<Binary>(
                BinaryOp::Div,
                arena_.make<Constant>(MathConstant::Pi),
                make_integer(arena_, BigInt(2))));
        }
        // arg(-i) = -pi/2
        if (const auto* un = expr_cast<Unary>(args.front()); un && un->op == UnaryOp::Neg) {
            if (const auto* c = expr_cast<Constant>(un->operand); c && c->value == MathConstant::I) {
                return simplify_expr(arena_.make<Unary>(
                    UnaryOp::Neg,
                    arena_.make<Binary>(BinaryOp::Div,
                        arena_.make<Constant>(MathConstant::Pi),
                        make_integer(arena_, BigInt(2)))));
            }
        }
        // a + b*i form
        auto parts = extract_complex(args.front());
        if (parts) {
            ExprPtr a = parts->first;
            ExprPtr b = parts->second;
            // a > 0 → arctan(b/a)
            if (is_known_positive(a)) {
                ExprPtr ratio = arena_.make<Binary>(BinaryOp::Div, b, a);
                ExprPtr atan = arena_.make<FuncCall>(BuiltinOp::Atan, std::vector<ExprPtr>{ratio});
                return simplify_expr(atan);
            }
            // a = 0
            if (is_zero_expr(a)) {
                if (is_known_positive(b)) {
                    return simplify_expr(arena_.make<Binary>(
                        BinaryOp::Div,
                        arena_.make<Constant>(MathConstant::Pi),
                        make_integer(arena_, BigInt(2))));
                }
                if (is_known_negative(b)) {
                    return simplify_expr(arena_.make<Unary>(
                        UnaryOp::Neg,
                        arena_.make<Binary>(BinaryOp::Div,
                            arena_.make<Constant>(MathConstant::Pi),
                            make_integer(arena_, BigInt(2)))));
                }
            }
            // a < 0:  arg = atan(b/a) + sign(b)*pi
            if (is_known_negative(a)) {
                ExprPtr ratio = arena_.make<Binary>(BinaryOp::Div, b, a);
                ExprPtr atan = arena_.make<FuncCall>(BuiltinOp::Atan, std::vector<ExprPtr>{ratio});
                ExprPtr pi_expr = arena_.make<Constant>(MathConstant::Pi);
                if (is_known_nonnegative(b)) {
                    return simplify_expr(arena_.make<Sum>(std::vector<ExprPtr>{atan, pi_expr}));
                }
                if (is_known_negative(b)) {
                    return simplify_expr(arena_.make<Sum>(std::vector<ExprPtr>{
                        atan, arena_.make<Unary>(UnaryOp::Neg, pi_expr)}));
                }
            }
        }
    }
    // ln(I) = I * pi / 2
    if (node.func_id == BuiltinOp::Ln && args.size() == 1U) {
        if (const auto* c = expr_cast<Constant>(args.front()); c && c->value == MathConstant::I) {
            ExprPtr i_pi = arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Constant>(MathConstant::I),
                arena_.make<Constant>(MathConstant::Pi)});
            return simplify_expr(
                arena_.make<Binary>(BinaryOp::Div, i_pi, make_integer(arena_, BigInt(2))));
        }
        // ln(-1) = I*pi  (principal branch).  Detect the literal -1.
        if (const auto* il = expr_cast<IntegerLit>(args.front()); il && il->value == BigInt(-1)) {
            return simplify_expr(arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Constant>(MathConstant::I),
                arena_.make<Constant>(MathConstant::Pi)}));
        }
    }
    if (node.func_id == BuiltinOp::Sign && args.size() == 1U) {
        if (is_known_positive(args.front())) {
            return traced_result(RuleId::SimplifySignPositive, target_before, make_integer(arena_, BigInt(1)));
        }
        if (is_known_negative(args.front())) {
            return traced_result(RuleId::SimplifySignNegative, target_before, make_integer(arena_, BigInt(-1)));
        }
        if (is_zero_expr(args.front())) {
            return traced_result(RuleId::SimplifySignZero, target_before, make_integer(arena_, BigInt(0)));
        }
    }

    if (context_ != nullptr && args.size() == 1U && expr_is<Matrix>(args.front())) {
        const auto& m_node = expr_ref<Matrix>(args.front());
        cas::linalg::MatrixExpr m_expr(m_node.rows, m_node.cols, m_node.elements);

        if (node.func_id == BuiltinOp::Det) {
            auto res = cas::linalg::determinant(m_expr, *context_);
            if (res.is_ok()) return ok(res.value());
        }
        if (node.func_id == BuiltinOp::Rank) {
            auto res = cas::linalg::rank(m_expr, *context_);
            if (res.is_ok()) return ok(make_integer(arena_, static_cast<long long>(res.value())));
        }
        if (node.func_id == BuiltinOp::Trace) {
            auto res = cas::linalg::trace(m_expr, *context_);
            if (res.is_ok()) return ok(res.value());
        }
        if (node.func_id == BuiltinOp::Inv) {
            auto res = cas::linalg::inverse(m_expr, *context_);
            if (res.is_ok()) return ok(arena_.make<Matrix>(res.value().rows(), res.value().cols(), res.value().elements()));
        }
        if (node.func_id == BuiltinOp::Transpose) {
            auto res = cas::linalg::transpose(m_expr);
            if (res.is_ok()) return ok(arena_.make<Matrix>(res.value().rows(), res.value().cols(), res.value().elements()));
        }
    }

    if (node.func_id == BuiltinOp::N && args.size() == 1U) {
        if (const auto* matrix = expr_cast<Matrix>(args.front())) {
            std::vector<ExprPtr> numeric_elements;
            for (auto elem : matrix->elements) {
                auto val = cas::numeric::eval(elem);
                if (val.is_ok()) {
                    numeric_elements.push_back(arena_.make<DecimalLit>(val.value()));
                } else {
                    numeric_elements.push_back(elem);
                }
            }
            return ok(arena_.make<Matrix>(matrix->rows, matrix->cols, std::move(numeric_elements)));
        }

        auto val = cas::numeric::eval(args.front());
        if (val.is_ok()) {
            return ok(arena_.make<DecimalLit>(val.value()));
        }
    }

    if (expr_ptr_sequence_identical(args, node.args)) return ok(original);
    return ok(arena_.make<FuncCall>(node.func_id, std::move(args)));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Integral& node) { return simplify_passthrough(original, node); }
Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Derivative& node) { return simplify_passthrough(original, node); }
Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Limit& node) { return simplify_passthrough(original, node); }

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const RootOf& node) {
    auto simplified_poly = simplify_expr(node.polynomial);
    if (simplified_poly.is_error()) return simplified_poly;

    if (context_) {
        // Risolviamo esplicitamente solo fino a max_rootof_explicit_degree (default 2) per evitare swell
        auto poly_res = cas::algebra::parse_polynomial(simplified_poly.value(), node.variable, *context_);
        if (poly_res.is_ok() && cas::algebra::poly_degree(poly_res.value()) <= context_->max_rootof_explicit_degree()) {
            auto roots_res = cas::algebra::solve_polynomial(simplified_poly.value(), node.variable, *context_);
            if (roots_res.is_ok()) {
                const auto& roots = roots_res.value();
                if (roots.size() == 1) {
                    if (roots[0]->kind != ExprKind::RootOf) return ok(roots[0]);
                } else if (node.root_index.has_value() && *node.root_index < roots.size()) {
                    if (roots[*node.root_index]->kind != ExprKind::RootOf) return ok(roots[*node.root_index]);
                }
            }
        }
    }

    if (simplified_poly.value() == node.polynomial) return ok(original);
    return ok(arena_.make<RootOf>(simplified_poly.value(), node.variable, node.root_index));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Matrix& node) {
    std::vector<ExprPtr> elements;
    elements.reserve(node.elements.size());
    for (std::size_t i = 0; i < node.elements.size(); ++i) {
        auto simplify_el = [&]() -> Result<ExprPtr> {
            if (trace_enabled_) {
                std::vector<ExprPtr> current = node.elements;
                for (std::size_t j = 0; j < elements.size(); ++j) current[j] = elements[j];
                ScopedFrame frame(*this, [this, current = std::move(current), i, r = node.rows, c = node.cols](ExprPtr value) mutable {
                    current[i] = value; return arena_.make<Matrix>(r, c, std::move(current));
                });
                return simplify_expr(node.elements[i]);
            }
            return simplify_expr(node.elements[i]);
        };
        auto res = simplify_el();
        if (res.is_error()) return res;
        elements.push_back(res.value());
    }
    if (expr_ptr_sequence_identical(elements, node.elements)) return ok(original);
    return ok(arena_.make<Matrix>(node.rows, node.cols, std::move(elements)));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const SeriesExp& node) {
    std::vector<std::pair<long long, ExprPtr>> simplified_terms;
    simplified_terms.reserve(node.terms.size());
    bool changed = false;

    // Simplify point
    auto simplified_point = [&]() -> Result<ExprPtr> {
        if (trace_enabled_) {
            ScopedFrame frame(*this, [this, &node](ExprPtr value) {
                return arena_.make<SeriesExp>(node.var, value, node.terms, node.order);
            });
            return simplify_expr(node.point);
        }
        return simplify_expr(node.point);
    }();
    if (simplified_point.is_error()) return simplified_point;
    if (simplified_point.value() != node.point) changed = true;

    // Simplify terms
    for (std::size_t i = 0; i < node.terms.size(); ++i) {
        auto simplified_coeff = [&]() -> Result<ExprPtr> {
            if (trace_enabled_) {
                ScopedFrame frame(*this, [this, &node, i](ExprPtr value) {
                    auto terms = node.terms;
                    terms[i].second = value;
                    return arena_.make<SeriesExp>(node.var, node.point, std::move(terms), node.order);
                });
                return simplify_expr(node.terms[i].second);
            }
            return simplify_expr(node.terms[i].second);
        }();
        if (simplified_coeff.is_error()) return simplified_coeff;
        if (simplified_coeff.value() != node.terms[i].second) changed = true;
        simplified_terms.push_back({node.terms[i].first, simplified_coeff.value()});
    }

    if (!changed) return ok(original);
    return ok(arena_.make<SeriesExp>(node.var, simplified_point.value(), std::move(simplified_terms), node.order));
}

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Quantity& node) {
    auto res = simplify_expr(node.value);
    if (res.is_error()) return res;
    if (res.value() == node.value) return ok(original);
    return ok(arena_.make<Quantity>(res.value(), node.dimensions));
}

template <typename Node>
Result<ExprPtr> Simplifier::simplify_passthrough(ExprPtr original, const Node&) { return ok(original); }

bool Simplifier::may_rewrite_function_call(BuiltinOp op, const std::vector<ExprPtr>& args) const {
    if (args.size() != 1U) return false;
    if (op == BuiltinOp::Tan) return true;
    if (is_parity_rewrite_function(op)) {
        if (is_zero_expr(args.front())) return op == BuiltinOp::Sin || op == BuiltinOp::Cos || op == BuiltinOp::Tan || op == BuiltinOp::Sinh || op == BuiltinOp::Tanh;
        const auto* unary = expr_cast<Unary>(args.front());
        return unary != nullptr && unary->op == UnaryOp::Neg;
    }
    if (op == BuiltinOp::Exp) return is_zero_expr(args.front()) || is_one_expr(args.front());
    if (op == BuiltinOp::Ln) {
        if (is_one_expr(args.front()) || is_constant_expr(args.front(), MathConstant::E)) return true;
        if (const auto* quot = expr_cast<Binary>(args.front()); quot != nullptr && quot->op == BinaryOp::Div) return is_known_positive(quot->left) && is_known_positive(quot->right);
        if (const auto* pow = expr_cast<Binary>(args.front()); pow != nullptr && pow->op == BinaryOp::Pow) return is_constant_expr(pow->left, MathConstant::E) || is_known_positive(pow->left);
        if (const auto* prod = expr_cast<Product>(args.front())) return std::all_of(prod->factors.begin(), prod->factors.end(), [this](ExprPtr f) { return is_known_positive(f); });
        const auto* sqrt = expr_cast<FuncCall>(args.front());
        return sqrt != nullptr && sqrt->func_id == BuiltinOp::Sqrt && sqrt->args.size() == 1U && is_known_positive(sqrt->args.front());
    }
    if (op != BuiltinOp::Sqrt) {
        if (op == BuiltinOp::Det || op == BuiltinOp::Rank || op == BuiltinOp::Trace || op == BuiltinOp::Inv || op == BuiltinOp::Transpose || op == BuiltinOp::N) {
            return !args.empty() && (expr_is<Matrix>(args.front()) || op == BuiltinOp::N);
        }
        return false;
    }
    if (const auto* quot = expr_cast<Binary>(args.front()); quot != nullptr && quot->op == BinaryOp::Div) return is_known_nonnegative(quot->left) && is_known_positive(quot->right);
    const auto* prod = expr_cast<Product>(args.front());
    return prod != nullptr && std::all_of(prod->factors.begin(), prod->factors.end(), [this](ExprPtr f) { return is_known_nonnegative(f); });
}

bool Simplifier::may_rewrite_sum_terms(const std::vector<ExprPtr>& terms) const {
    ExprPtr s_arg, c_arg;
    for (ExprPtr term : terms) {
        const auto* pow = expr_cast<Binary>(term);
        if (pow == nullptr || pow->op != BinaryOp::Pow) continue;
        if (auto exp = try_get_integer_exponent(pow->right); !exp.has_value() || *exp != BigInt(2)) continue;
        const auto* call = expr_cast<FuncCall>(pow->left);
        if (call == nullptr || call->args.size() != 1U) continue;
        if (call->func_id == BuiltinOp::Sin) {
            if (c_arg && call->args.front() == c_arg) return true;
            if (!s_arg) s_arg = call->args.front();
        } else if (call->func_id == BuiltinOp::Cos) {
            if (s_arg && call->args.front() == s_arg) return true;
            if (!c_arg) c_arg = call->args.front();
        }
    }
    return false;
}

bool Simplifier::may_rewrite_power(ExprPtr base, ExprPtr exponent) const {
    if (!is_constant_expr(base, MathConstant::E)) return false;
    const auto* ln = expr_cast<FuncCall>(exponent);
    return ln != nullptr && ln->func_id == BuiltinOp::Ln && ln->args.size() == 1U && is_known_positive(ln->args.front());
}

bool Simplifier::is_known_positive(ExprPtr expr) const {
    if (!expr) return false;
    
    // 1. Check Assumptions (Primary Choice)
    if (assumptions_ != nullptr && assumptions_->is_positive(expr)) return true;
    
    // 2. Fallback for literals and constants
    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value()) {
        return !rat.value.numerator().is_zero() && !rat.value.numerator().is_negative();
    }
    if (const auto* constant = expr_cast<Constant>(expr)) return is_known_positive_constant(constant->value);
    
    // 3. Structural fallback (if assumptions_ is null or missed it)
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Pow) return is_known_positive(bin->left);
        if (bin->op == BinaryOp::Div) return is_known_positive(bin->left) && is_known_positive(bin->right);
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        if (prod->factors.empty()) return false;
        for (ExprPtr f : prod->factors) if (!is_known_positive(f)) return false;
        return true;
    }
    
    return false;
}

bool Simplifier::is_known_nonnegative(ExprPtr expr) const {
    if (!expr) return false;
    
    // 1. Check Assumptions (Primary Choice)
    if (assumptions_ != nullptr && assumptions_->is_nonnegative(expr)) return true;
    
    // 2. Check if it's strictly positive
    if (is_known_positive(expr)) return true;

    // 3. Fallback for literals and constants
    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value()) {
        return !rat.value.numerator().is_negative();
    }
    if (const auto* constant = expr_cast<Constant>(expr)) return is_known_nonnegative_constant(constant->value);
    
    // 4. Structural fallback
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Div) return is_known_nonnegative(bin->left) && is_known_positive(bin->right);
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        if (prod->factors.empty()) return false;
        for (ExprPtr f : prod->factors) if (!is_known_nonnegative(f)) return false;
        return true;
    }
    
    return false;
}

bool Simplifier::is_known_negative(ExprPtr expr) const {
    if (!expr) return false;
    if (assumptions_ != nullptr && assumptions_->is_negative(expr)) return true;

    LiteralRational rat;
    if (auto exact = try_get_exact_rational(expr, rat); exact.is_ok() && exact.value()) return rat.value.numerator().is_negative();
    if (const auto* bin = expr_cast<Binary>(expr)) {
        if (bin->op == BinaryOp::Div) return (is_known_negative(bin->left) && is_known_positive(bin->right)) ||
                                            (is_known_positive(bin->left) && is_known_negative(bin->right));
    }
    if (const auto* prod = expr_cast<Product>(expr)) {
        int negative_count = 0;
        for (ExprPtr f : prod->factors) {
            if (is_known_negative(f)) negative_count++;
            else if (!is_known_positive(f)) return false;
        }
        return (negative_count % 2 != 0);
    }
    if (const auto* unary = expr_cast<Unary>(expr); unary != nullptr && unary->op == UnaryOp::Neg) {
        return is_known_positive(unary->operand);
    }
    return false;
}

} // namespace cas::symbolic::detail
