#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

// ── P2-005 / L2-10: Trig angle reduction helpers ────────────────────────────

[[nodiscard]] static std::optional<Rational> try_extract_pi_coefficient(ExprPtr arg) {
    if (!arg) return std::nullopt;

    if (const auto* c = expr_cast<Constant>(arg); c && c->value == MathConstant::Pi)
        return Rational(BigInt(1));

    if (const auto* u = expr_cast<Unary>(arg); u && u->op == UnaryOp::Neg) {
        if (const auto* c = expr_cast<Constant>(u->operand); c && c->value == MathConstant::Pi)
            return Rational(BigInt(-1));
    }

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

[[nodiscard]] static Rational reduce_to_period(Rational r) {
    const BigInt two_den = BigInt(2) * r.denominator();
    BigInt q = r.numerator() / two_den;
    BigInt rem = r.numerator() - q * two_den;
    if (rem < BigInt(0)) {
        q = q - BigInt(1);
        rem = rem + two_den;
    }
    return Rational(rem, r.denominator());
}

// Base values for constructible angles.  These are Gauss-construction
// primitives for Fermat primes p=3,5 — NOT an exhaustive lookup table.
// Arbitrary angles reach via half-angle recursion in cos/sin_ref_value.

[[nodiscard]] static ExprPtr cos_pi_over_5(AstArena& arena) {
    ExprPtr sqrt5 = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_integer(arena, BigInt(5))});
    ExprPtr sum = arena.make<Sum>(std::vector<ExprPtr>{make_integer(arena, BigInt(1)), sqrt5});
    return arena.make<Binary>(BinaryOp::Div, sum, make_integer(arena, BigInt(4)));
}

[[nodiscard]] static ExprPtr cos_2pi_over_5(AstArena& arena) {
    ExprPtr sqrt5 = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_integer(arena, BigInt(5))});
    ExprPtr neg_one = arena.make<Unary>(UnaryOp::Neg, make_integer(arena, BigInt(1)));
    ExprPtr sum = arena.make<Sum>(std::vector<ExprPtr>{sqrt5, neg_one});
    return arena.make<Binary>(BinaryOp::Div, sum, make_integer(arena, BigInt(4)));
}

[[nodiscard]] static ExprPtr sin_pi_over_5(AstArena& arena) {
    ExprPtr sqrt5 = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_integer(arena, BigInt(5))});
    ExprPtr two_sqrt5 = arena.make<Product>(std::vector<ExprPtr>{
        make_integer(arena, BigInt(2)), sqrt5});
    ExprPtr neg_two_sqrt5 = arena.make<Unary>(UnaryOp::Neg, two_sqrt5);
    ExprPtr inner_sum = arena.make<Sum>(std::vector<ExprPtr>{
        make_integer(arena, BigInt(10)), neg_two_sqrt5});
    ExprPtr outer = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{inner_sum});
    return arena.make<Binary>(BinaryOp::Div, outer, make_integer(arena, BigInt(4)));
}

[[nodiscard]] static ExprPtr sin_2pi_over_5(AstArena& arena) {
    ExprPtr sqrt5 = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_integer(arena, BigInt(5))});
    ExprPtr two_sqrt5 = arena.make<Product>(std::vector<ExprPtr>{
        make_integer(arena, BigInt(2)), sqrt5});
    ExprPtr inner_sum = arena.make<Sum>(std::vector<ExprPtr>{
        make_integer(arena, BigInt(10)), two_sqrt5});
    ExprPtr outer = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{inner_sum});
    return arena.make<Binary>(BinaryOp::Div, outer, make_integer(arena, BigInt(4)));
}

[[nodiscard]] static ExprPtr sin_ref_value_table(Rational ref, AstArena& arena) {
    const Rational zero(BigInt(0));
    const Rational one_tenth(BigInt(1), BigInt(10));
    const Rational one_sixth(BigInt(1), BigInt(6));
    const Rational one_fifth(BigInt(1), BigInt(5));
    const Rational one_quarter(BigInt(1), BigInt(4));
    const Rational three_tenths(BigInt(3), BigInt(10));
    const Rational one_third(BigInt(1), BigInt(3));
    const Rational two_fifths(BigInt(2), BigInt(5));
    const Rational one_half(BigInt(1), BigInt(2));

    if (ref == zero)       return make_integer(arena, BigInt(0));
    if (ref == one_tenth)  return cos_2pi_over_5(arena);
    if (ref == one_sixth)  return make_rational(arena, Rational(BigInt(1), BigInt(2)));
    if (ref == one_fifth)  return sin_pi_over_5(arena);
    if (ref == one_quarter)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{make_integer(arena, BigInt(2))}),
            make_integer(arena, BigInt(2)));
    if (ref == three_tenths) return cos_pi_over_5(arena);
    if (ref == one_third)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{make_integer(arena, BigInt(3))}),
            make_integer(arena, BigInt(2)));
    if (ref == two_fifths) return sin_2pi_over_5(arena);
    if (ref == one_half)   return make_integer(arena, BigInt(1));
    return nullptr;
}

[[nodiscard]] static ExprPtr cos_ref_value_table(Rational ref, AstArena& arena) {
    const Rational zero(BigInt(0));
    const Rational one_tenth(BigInt(1), BigInt(10));
    const Rational one_sixth(BigInt(1), BigInt(6));
    const Rational one_fifth(BigInt(1), BigInt(5));
    const Rational one_quarter(BigInt(1), BigInt(4));
    const Rational three_tenths(BigInt(3), BigInt(10));
    const Rational one_third(BigInt(1), BigInt(3));
    const Rational two_fifths(BigInt(2), BigInt(5));
    const Rational one_half(BigInt(1), BigInt(2));

    if (ref == zero)       return make_integer(arena, BigInt(1));
    if (ref == one_tenth)  return sin_2pi_over_5(arena);
    if (ref == one_sixth)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{make_integer(arena, BigInt(3))}),
            make_integer(arena, BigInt(2)));
    if (ref == one_fifth)  return cos_pi_over_5(arena);
    if (ref == one_quarter)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{make_integer(arena, BigInt(2))}),
            make_integer(arena, BigInt(2)));
    if (ref == three_tenths) return sin_pi_over_5(arena);
    if (ref == one_third)  return make_rational(arena, Rational(BigInt(1), BigInt(2)));
    if (ref == two_fifths) return cos_2pi_over_5(arena);
    if (ref == one_half)   return make_integer(arena, BigInt(0));
    return nullptr;
}

// Forward declarations for mutual recursion.
[[nodiscard]] static ExprPtr cos_ref_value(Rational ref, AstArena& arena);
[[nodiscard]] static ExprPtr sin_ref_value(Rational ref, AstArena& arena);
[[nodiscard]] static ExprPtr try_angle_combination(Rational ref, BuiltinOp func, AstArena& arena);
[[nodiscard]] static ExprPtr trig_exact_at_pi_multiple(BuiltinOp func_id, Rational r, AstArena& arena);

// L2-10: half-angle recursion — cos(α/2) = sqrt((1+cos α)/2).
// Arbitrary depth via doubling; no fixed power-of-two table.
[[nodiscard]] static ExprPtr cos_ref_value_half_angle(Rational ref, AstArena& arena) {
    const Rational zero(BigInt(0));
    const Rational one_half(BigInt(1), BigInt(2));
    if (ref < zero || ref > one_half) return nullptr;
    if (ExprPtr table_hit = cos_ref_value_table(ref, arena)) return table_hit;
    std::vector<Rational> path;
    Rational current = ref;
    ExprPtr base = nullptr;
    constexpr int kMaxHalfAngleDepth = 32;
    for (int depth = 0; depth < kMaxHalfAngleDepth; ++depth) {
        Rational doubled = current * Rational(BigInt(2));
        if (doubled > one_half) {
            Rational complement = Rational(BigInt(1)) - doubled;
            if (complement < zero || complement > one_half) break;
            ExprPtr cos_compl = cos_ref_value_table(complement, arena);
            if (cos_compl != nullptr) {
                base = arena.make<Unary>(UnaryOp::Neg, cos_compl);
                break;
            }
            break;
        }
        path.push_back(current);
        current = doubled;
        if ((base = cos_ref_value_table(current, arena)) != nullptr) break;
    }
    if (base == nullptr) return nullptr;
    ExprPtr cos_val = base;
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
        ExprPtr one_plus = arena.make<Sum>(std::vector<ExprPtr>{
            make_integer(arena, BigInt(1)), cos_val});
        ExprPtr half_arg = arena.make<Binary>(BinaryOp::Div, one_plus,
            make_integer(arena, BigInt(2)));
        cos_val = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{half_arg});
    }
    return cos_val;
}

// L2-10: Chebyshev T_p applied to cos(π/q) gives cos(p·π/q).
[[nodiscard]] static ExprPtr cos_ref_value(Rational ref, AstArena& arena) {
    if (ExprPtr direct = cos_ref_value_half_angle(ref, arena)) return direct;
    const BigInt& p = ref.numerator();
    const BigInt& q = ref.denominator();
    if (p <= BigInt(1) || q <= BigInt(1)) {
        // p=1 case not reachable by half-angle: try angle combination (e.g. cos(π/15)).
        if (p == BigInt(1) && q > BigInt(1))
            return try_angle_combination(ref, BuiltinOp::Cos, arena);
        return nullptr;
    }
    if (p.bit_length() > 16) return nullptr;
    Rational one_over_q(BigInt(1), q);
    ExprPtr cos_q = cos_ref_value_half_angle(one_over_q, arena);
    if (cos_q == nullptr) {
        // Base angle not in half-angle reachable set; try combination path.
        if (ExprPtr combo = try_angle_combination(ref, BuiltinOp::Cos, arena))
            return combo;
        return nullptr;
    }
    const std::uint64_t p_u = p.to_u64();
    ExprPtr T_prev = make_integer(arena, BigInt(1));
    ExprPtr T_curr = cos_q;
    for (std::uint64_t k = 1U; k < p_u; ++k) {
        ExprPtr two_x_Tk = arena.make<Product>(std::vector<ExprPtr>{
            make_integer(arena, BigInt(2)), cos_q, T_curr});
        ExprPtr neg_prev = arena.make<Unary>(UnaryOp::Neg, T_prev);
        ExprPtr next = arena.make<Sum>(std::vector<ExprPtr>{two_x_Tk, neg_prev});
        T_prev = T_curr;
        T_curr = next;
    }
    return T_curr;
}

[[nodiscard]] static ExprPtr sin_ref_value(Rational ref, AstArena& arena) {
    if (ExprPtr table_hit = sin_ref_value_table(ref, arena)) return table_hit;
    const Rational zero(BigInt(0));
    const Rational one_half(BigInt(1), BigInt(2));
    if (ref < zero || ref > one_half) return nullptr;
    Rational doubled = ref * Rational(BigInt(2));
    ExprPtr cos_2ref = nullptr;
    if (doubled > one_half) {
        Rational complement = Rational(BigInt(1)) - doubled;
        ExprPtr cos_c = cos_ref_value(complement, arena);
        if (cos_c != nullptr) cos_2ref = arena.make<Unary>(UnaryOp::Neg, cos_c);
    } else {
        cos_2ref = cos_ref_value(doubled, arena);
    }
    if (cos_2ref == nullptr) {
        // cos(2·ref·π) not reachable; try angle combination for sin.
        return try_angle_combination(ref, BuiltinOp::Sin, arena);
    }
    ExprPtr neg = arena.make<Unary>(UnaryOp::Neg, cos_2ref);
    ExprPtr one_minus = arena.make<Sum>(std::vector<ExprPtr>{make_integer(arena, BigInt(1)), neg});
    ExprPtr half_arg = arena.make<Binary>(BinaryOp::Div, one_minus, make_integer(arena, BigInt(2)));
    return arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{half_arg});
}

// L2-10: angle subtraction formula for constructible angles not reachable by halving/Chebyshev.
// Handles denominators like 15 = lcm(3,5) by expressing ref = r1 - r2 where r1 is a base
// table angle (den ≤ 10) and r2 has strictly smaller denominator than ref.
// cos(A-B) = cos(A)cos(B) + sin(A)sin(B), sin(A-B) = sin(A)cos(B) - cos(A)sin(B).
[[nodiscard]] static ExprPtr try_angle_combination(Rational ref, BuiltinOp func, AstArena& arena) {
    const Rational zero(BigInt(0));
    if (ref.denominator() < BigInt(7)) return nullptr;

    // Enumerates rational multiples of π with small denominators as r1 candidates.
    static const std::pair<int,int> kBaseAngles[] = {
        {1,2},{1,3},{2,3},{1,4},{3,4},{1,5},{2,5},{3,5},{4,5},
        {1,6},{5,6},{1,10},{3,10},{7,10},{9,10}
    };
    for (auto [k1_i, d1_i] : kBaseAngles) {
        auto r1 = Rational{BigInt{k1_i}, BigInt{d1_i}};
        auto r2 = r1 - ref;
        if (r2 <= zero) continue;
        if (r2.denominator() >= ref.denominator()) continue;

        ExprPtr c1 = trig_exact_at_pi_multiple(BuiltinOp::Cos, r1, arena);
        ExprPtr s1 = trig_exact_at_pi_multiple(BuiltinOp::Sin, r1, arena);
        ExprPtr c2 = trig_exact_at_pi_multiple(BuiltinOp::Cos, r2, arena);
        ExprPtr s2 = trig_exact_at_pi_multiple(BuiltinOp::Sin, r2, arena);
        if (!c1 || !s1 || !c2 || !s2) continue;

        if (func == BuiltinOp::Cos) {
            // cos(r1-r2)·π = cos(r1·π)cos(r2·π) + sin(r1·π)sin(r2·π)
            return arena.make<Sum>(std::vector<ExprPtr>{
                arena.make<Product>(std::vector<ExprPtr>{c1, c2}),
                arena.make<Product>(std::vector<ExprPtr>{s1, s2})});
        } else {
            // sin(r1-r2)·π = sin(r1·π)cos(r2·π) - cos(r1·π)sin(r2·π)
            return arena.make<Sum>(std::vector<ExprPtr>{
                arena.make<Product>(std::vector<ExprPtr>{s1, c2}),
                arena.make<Unary>(UnaryOp::Neg,
                    arena.make<Product>(std::vector<ExprPtr>{c1, s2}))});
        }
    }
    return nullptr;
}

[[nodiscard]] static ExprPtr negate_expr(ExprPtr e, AstArena& arena) {
    if (is_zero_expr(e)) return e;
    LiteralRational lr;
    auto ex = try_get_exact_rational(e, lr);
    if (ex.is_ok() && ex.value()) return make_rational(arena, -lr.value);
    return arena.make<Unary>(UnaryOp::Neg, e);
}

// Quadrant reduction: sin/cos(r·π) for any rational r.
[[nodiscard]] static ExprPtr trig_exact_at_pi_multiple(
    BuiltinOp func_id, Rational r, AstArena& arena) {

    const Rational half(BigInt(1), BigInt(2));
    const Rational one(BigInt(1));
    const Rational three_halves(BigInt(3), BigInt(2));

    Rational red = reduce_to_period(r);
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

// ── Simplifier::simplify_funcall_trig ────────────────────────────────────────

Result<ExprPtr> Simplifier::simplify_funcall_trig(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before) {

    if (op == BuiltinOp::Sin && args.size() == 1U) {
        if (is_zero_expr(args.front()))
            return traced_result(RuleId::SimplifySinZero, target_before, make_integer(arena_, BigInt(0)));
        if (is_constant_expr(args.front(), MathConstant::Pi))
            return traced_result(RuleId::SimplifySinPi, target_before, make_integer(arena_, BigInt(0)));
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Asin)
            return ok(call->args[0]);
        if (auto coeff = try_extract_pi_coefficient(args.front())) {
            if (ExprPtr val = trig_exact_at_pi_multiple(BuiltinOp::Sin, *coeff, arena_)) {
                auto simplified = simplify_expr(val);
                if (simplified.is_ok())
                    return traced_result(RuleId::Unknown, target_before, simplified.value());
            }
        }
    }
    if (op == BuiltinOp::Cos && args.size() == 1U) {
        if (is_zero_expr(args.front()))
            return traced_result(RuleId::SimplifyCosZero, target_before, make_integer(arena_, BigInt(1)));
        if (is_constant_expr(args.front(), MathConstant::Pi))
            return traced_result(RuleId::SimplifyCosPi, target_before, make_integer(arena_, BigInt(-1)));
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Acos)
            return ok(call->args[0]);
        if (auto coeff = try_extract_pi_coefficient(args.front())) {
            if (ExprPtr val = trig_exact_at_pi_multiple(BuiltinOp::Cos, *coeff, arena_)) {
                auto simplified = simplify_expr(val);
                if (simplified.is_ok())
                    return traced_result(RuleId::Unknown, target_before, simplified.value());
            }
        }
    }
    if (op == BuiltinOp::Tan && args.size() == 1U) {
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Atan)
            return ok(call->args[0]);
    }

    if (op == BuiltinOp::Asin && args.size() == 1U) {
        // asin(sin(x)) -> x if -π/2 ≤ x ≤ π/2
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Sin) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi_2 = arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi), make_integer(arena_, BigInt(2)));
                ExprPtr neg_pi_2 = arena_.make<Unary>(UnaryOp::Neg, pi_2);
                if (assumptions_->is_greater_equal(x, neg_pi_2)
                    && assumptions_->is_greater_equal(pi_2, x))
                    return ok(x);
            }
        }
    }
    if (op == BuiltinOp::Acos && args.size() == 1U) {
        // acos(cos(x)) -> x if 0 ≤ x ≤ π
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Cos) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi = arena_.make<Constant>(MathConstant::Pi);
                ExprPtr zero = make_integer(arena_, BigInt(0));
                if (assumptions_->is_greater_equal(x, zero)
                    && assumptions_->is_greater_equal(pi, x))
                    return ok(x);
            }
        }
    }
    if (op == BuiltinOp::Atan && args.size() == 1U) {
        // atan(tan(x)) -> x if -π/2 < x < π/2
        if (const auto* call = expr_cast<FuncCall>(args.front());
            call && call->func_id == BuiltinOp::Tan) {
            ExprPtr x = call->args[0];
            if (assumptions_) {
                ExprPtr pi_2 = arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi), make_integer(arena_, BigInt(2)));
                ExprPtr neg_pi_2 = arena_.make<Unary>(UnaryOp::Neg, pi_2);
                if (assumptions_->is_greater(x, neg_pi_2)
                    && assumptions_->is_greater(pi_2, x))
                    return ok(x);
            }
        }
        if (is_zero_expr(args.front())) return ok(make_integer(arena_, BigInt(0)));
        if (const auto* il = expr_cast<IntegerLit>(args.front())) {
            if (il->value == BigInt(1))
                return simplify_expr(arena_.make<Binary>(BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi),
                    make_integer(arena_, BigInt(4))));
            if (il->value == BigInt(-1))
                return simplify_expr(arena_.make<Unary>(UnaryOp::Neg,
                    arena_.make<Binary>(BinaryOp::Div,
                        arena_.make<Constant>(MathConstant::Pi),
                        make_integer(arena_, BigInt(4)))));
        }
        // atan odd: atan(-x) = -atan(x)
        if (const auto* un = expr_cast<Unary>(args.front()); un && un->op == UnaryOp::Neg) {
            ExprPtr inner_atan = arena_.make<FuncCall>(
                BuiltinOp::Atan, std::vector<ExprPtr>{un->operand});
            return simplify_expr(arena_.make<Unary>(UnaryOp::Neg, inner_atan));
        }
    }

    // L2-07: Addition formula — sin/cos(x + kπ/n) where kπ/n has an exact value.
    // Only triggered when the argument is a Sum containing at least one extractable π-term.
    // The non-π residual becomes the "x" argument; cos/sin of the π-part evaluate to exact
    // algebraic numbers, so the result is always a simplification, never an expansion spiral.
    if ((op == BuiltinOp::Sin || op == BuiltinOp::Cos) && args.size() == 1U) {
        if (const auto* sum_arg = expr_cast<Sum>(args.front())) {
            Rational pi_total(BigInt(0));
            std::vector<ExprPtr> non_pi_terms;
            for (ExprPtr term : sum_arg->terms) {
                if (auto c = try_extract_pi_coefficient(term))
                    pi_total = pi_total + *c;
                else
                    non_pi_terms.push_back(term);
            }
            if (!non_pi_terms.empty() && !(pi_total == Rational(BigInt(0)))) {
                ExprPtr cv = trig_exact_at_pi_multiple(BuiltinOp::Cos, pi_total, arena_);
                ExprPtr sv = trig_exact_at_pi_multiple(BuiltinOp::Sin, pi_total, arena_);
                if (cv && sv) {
                    ExprPtr x = (non_pi_terms.size() == 1U)
                        ? non_pi_terms[0]
                        : arena_.make<Sum>(std::move(non_pi_terms));
                    auto cv_s = simplify_expr(cv);
                    auto sv_s = simplify_expr(sv);
                    if (cv_s.is_ok() && sv_s.is_ok()) {
                        ExprPtr sin_x = arena_.make<FuncCall>(BuiltinOp::Sin, std::vector<ExprPtr>{x});
                        ExprPtr cos_x = arena_.make<FuncCall>(BuiltinOp::Cos, std::vector<ExprPtr>{x});
                        ExprPtr t1, t2, expanded;
                        if (op == BuiltinOp::Sin) {
                            // sin(x + kπ) = sin(x)cos(kπ) + cos(x)sin(kπ)
                            t1 = arena_.make<Product>(std::vector<ExprPtr>{cv_s.value(), sin_x});
                            t2 = arena_.make<Product>(std::vector<ExprPtr>{sv_s.value(), cos_x});
                        } else {
                            // cos(x + kπ) = cos(x)cos(kπ) - sin(x)sin(kπ)
                            t1 = arena_.make<Product>(std::vector<ExprPtr>{cv_s.value(), cos_x});
                            ExprPtr neg_sv = arena_.make<Unary>(UnaryOp::Neg, sv_s.value());
                            t2 = arena_.make<Product>(std::vector<ExprPtr>{neg_sv, sin_x});
                        }
                        expanded = arena_.make<Sum>(std::vector<ExprPtr>{t1, t2});
                        return simplify_expr(expanded);
                    }
                }
            }
        }
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

} // namespace cas::symbolic::detail
