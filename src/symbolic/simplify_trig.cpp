#include "simplify_impl.hpp"
#include "simplify_trig_chebyshev_impl.hpp"
#include "simplify_trig_tables_impl.hpp"
#include "cas/formatter.hpp"
#include <cstdlib>
#include <iostream>

namespace cas::symbolic::detail {
static bool is_debug_chebyshev() { static const char* v = std::getenv("CAS_DEBUG_CHEBYSHEV"); return v && std::string(v) == "1"; }
}



// F1.4c — Stack depth guard for try_angle_combination.
// Uses a thread-local counter to prevent mutual recursion:
//   try_angle_combination → trig_exact_at_pi_multiple
//     → cos_ref_value / sin_ref_value → try_angle_combination
// Maximum recursion depth: kTrigCombinationMaxDepth (= 3).
// At depth ≥ limit, try_angle_combination returns nullptr so the
// caller falls through to the RootOf generator.
//
// Rationale for depth=3:
//   Depth 0: outer call for target angle (e.g. cos(π/7))
//   Depth 1: recursive call for sub-angle r2 = r1 - ref within the
//            combination loop (r2 has smaller denominator than ref)
//   Depth 2: second-level sub-angle (rare but occurs for composite q)
//   Depth 3: guaranteed termination — denominators strictly decrease
//            at each depth level so combinatorial explosion is bounded.
// Reference: CLAUDE.md §REGOLA ZERO; HARDCODE_LEDGER HPP-014-c.

namespace cas::symbolic::detail {

// Thread-local depth guard for try_angle_combination.
// Increment on entry, decrement on exit (RAII guard below).
// Zero when no combination is in progress.
static thread_local int s_trig_combination_depth = 0;

// RAII guard: increments depth on construction, decrements on destruction.
struct TrigCombinationDepthGuard {
    TrigCombinationDepthGuard()  noexcept { ++s_trig_combination_depth; }
    ~TrigCombinationDepthGuard() noexcept { --s_trig_combination_depth; }
    TrigCombinationDepthGuard(const TrigCombinationDepthGuard&) = delete;
    TrigCombinationDepthGuard& operator=(const TrigCombinationDepthGuard&) = delete;
};

// Maximum recursion depth for try_angle_combination.
// Mathematical justification: denominators strictly decrease at each
// recursion level (r2.denom < ref.denom ≤ kBaseAngleMaxDenom=60).
// With 3 levels we cover all constructible angle subtractions reachable
// without triggering exponential fan-out on non-constructible denominators
// (q=7, q=17, q=11, ...).
// Configurable here; see also HARDCODE_LEDGER HPP-014c.
constexpr int kTrigCombinationMaxDepth = 3;

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
    const BigInt& p = ref.numerator(); const BigInt& q = ref.denominator();
    if (is_debug_chebyshev()) std::cerr << "[CHEBYSHEV_DEBUG] cos_ref_value: ref = " << p.decimal() << "/" << q.decimal() << "\n";
    if (ExprPtr direct = cos_ref_value_half_angle(ref, arena)) {
        if (is_debug_chebyshev()) std::cerr << "  [DIRECT HALVING SUCCESS] ref = " << p.decimal() << "/" << q.decimal() << "\n";
        return direct;
    }
    // Co-function reduction: cos(ref·π) = sin((1/2 − ref)·π).  For ref = p/q
    // with p > q/3, the complement (q−2p)/(2q) reduces to a STRICTLY smaller
    // numerator (e.g. 7/16 → 1/16), reachable by the cheap sin half-angle path —
    // avoiding the O(p) Chebyshev T_p recurrence below, which otherwise builds an
    // enormous expression that simplifies right at the depth/timeout budget and
    // bails non-deterministically under suite-wide load (root cause of T-004,
    // the flaky cos(7π/16)).  The strict-decrease guard bounds the mutual
    // cos/sin recursion and never makes an already-cheap case worse.  Both ref
    // and its complement lie in [0,1/2], so the identity holds with no sign flip.
    {
        const Rational cofunc = Rational(BigInt(1), BigInt(2)) - ref;  // auto-reduced
        if (cofunc.numerator() > BigInt(0) && cofunc.numerator() < p) {
            if (ExprPtr s = sin_ref_value(cofunc, arena)) {
                if (is_debug_chebyshev()) std::cerr << "  [COFUNCTION cos→sin] ref = " << p.decimal() << "/" << q.decimal() << " → sin " << cofunc.numerator().decimal() << "/" << cofunc.denominator().decimal() << "\n";
                return s;
            }
        }
    }
    if (p <= BigInt(1) || q <= BigInt(1)) {
        if (p == BigInt(1) && q > BigInt(1)) {
            if (q.bit_length() <= 16U) {
                const auto q_u64 = q.to_u64();
                if (q_u64 == 17U || q_u64 == 257U || q_u64 == 65537U) {
                    if (q_u64 <= static_cast<std::uint64_t>(kCosPolyMaxQ)) {
                        if (is_debug_chebyshev()) std::cerr << "  [FERMAT PRIME ROOTOF] q = " << q_u64 << "\n";
                        return build_rootof_cos_pi_q(static_cast<int>(q_u64), arena);
                    }
                }
            }
            if (ExprPtr combo = try_angle_combination(ref, BuiltinOp::Cos, arena)) {
                if (is_debug_chebyshev()) std::cerr << "  [COMBO SUCCESS] ref = " << p.decimal() << "/" << q.decimal() << "\n";
                return combo;
            }
            if (q.bit_length() <= 16U) {
                const auto q_u64 = q.to_u64();
                if (q_u64 >= 2U && q_u64 <= static_cast<std::uint64_t>(kCosPolyMaxQ)) {
                    if (is_debug_chebyshev()) std::cerr << "  [NON-CONSTRUCTIBLE ROOTOF] q = " << q_u64 << "\n";
                    return build_rootof_cos_pi_q(static_cast<int>(q_u64), arena);
                }
            }
            return nullptr;
        }
        return nullptr;
    }
    if (p.bit_length() > 16) return nullptr;
    Rational one_over_q(BigInt(1), q);
    ExprPtr cos_q = cos_ref_value_half_angle(one_over_q, arena);
    if (cos_q == nullptr) {
        if (ExprPtr combo = try_angle_combination(ref, BuiltinOp::Cos, arena)) {
            if (is_debug_chebyshev()) std::cerr << "  [BASE COMBO SUCCESS] ref = " << p.decimal() << "/" << q.decimal() << "\n";
            return combo;
        }
        {
            const auto q_u64 = q.to_u64();
            if (q_u64 >= 2U && q_u64 <= static_cast<std::uint64_t>(kCosPolyMaxQ)) {
                if (is_debug_chebyshev()) std::cerr << "  [FALLBACK BASE ROOTOF] q = " << q_u64 << "\n";
                cos_q = build_rootof_cos_pi_q(static_cast<int>(q_u64), arena);
            }
        }
        if (cos_q == nullptr) return nullptr;
    }
    const std::uint64_t p_u = p.to_u64();
    if (is_debug_chebyshev()) std::cerr << "[CHEBYSHEV_DEBUG] Chebyshev loop: p = " << p_u << ", cos_q = " << formatter::TextFormatter{}.format(cos_q) << "\n";
    ExprPtr T_prev = make_integer(arena, BigInt(1));
    ExprPtr T_curr = cos_q;
    for (std::uint64_t k = 1U; k < p_u; ++k) {
        ExprPtr two_x_Tk = arena.make<Product>(std::vector<ExprPtr>{make_integer(arena, BigInt(2)), cos_q, T_curr});
        ExprPtr neg_prev = arena.make<Unary>(UnaryOp::Neg, T_prev);
        ExprPtr next = arena.make<Sum>(std::vector<ExprPtr>{two_x_Tk, neg_prev});
        T_prev = T_curr; T_curr = next;
        if (is_debug_chebyshev()) std::cerr << "  k = " << k << ", T_curr = " << formatter::TextFormatter{}.format(T_curr) << "\n";
    }
    if (is_debug_chebyshev()) std::cerr << "[CHEBYSHEV_DEBUG] cos_ref_value result: " << formatter::TextFormatter{}.format(T_curr) << "\n";
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

// L2-10 / HPP-014 resolved: angle subtraction formula with algorithmically-generated
// base angle set.  Replaces the closed kBaseAngles[] table (HPP-014) with all
// constructible p/q ∈ (0,1] (gcd(p,q)=1) reachable by cos_ref_value_half_angle
// for denominators q ≤ kBaseAngleMaxDenom.
//
// Constructibility criterion (Gauss, Disquisitiones Arithmeticae §VII):
//   cos(π/q) is expressible in nested radicals over Q ⟺
//   q = 2^a · ∏ distinct Fermat primes from {3,5,17,257,65537}.
//   The half-angle recursion in cos_ref_value_half_angle is sound for exactly
//   these denominators (it bottoms out at the base table values for q ∈ {1,2,3,4,5,6,10}).
//   All other p/q are reachable as combinations of constructible angles via
//   the subtraction formula, or fall through to the RootOf generator.
//
// Max denominator: kBaseAngleMaxDenom (configurable at compile time; default 60
// covers all constructible denominators ≤ 60 = lcm(3,4,5) and handles test
// cases cos(π/15), cos(π/17), cos(2π/7), sin(π/30) and similar).
// Runtime cost: O(kBaseAngleMaxDenom²) angle checks per call; each check is O(1)
// table lookup or O(log q) half-angle recursion — acceptable for q ≤ 60.
// For larger q, the RootOf fallback in cos_ref_value handles the remainder.
//
// Tracked in HARDCODE_LEDGER.md: HPP-014 → Resolved (this commit).
// Remaining open item: wire CASContext::max_trig_exact_denom() through the
// static call chain so the limit is user-configurable at runtime (see
// HARDCODE_LEDGER.md entry HPP-014-RT).
constexpr int kBaseAngleMaxDenom = 120;

[[nodiscard]] static ExprPtr try_angle_combination(Rational ref, BuiltinOp func, AstArena& arena) {
    if (is_debug_chebyshev()) std::cerr << "[CHEBYSHEV_DEBUG] try_angle_combination: ref = " << ref.numerator().decimal() << "/" << ref.denominator().decimal() << ", func = " << (func == BuiltinOp::Cos ? "Cos" : "Sin") << ", depth = " << s_trig_combination_depth << "\n";
    if (s_trig_combination_depth >= kTrigCombinationMaxDepth) {
        if (is_debug_chebyshev()) std::cerr << "  [DEPTH LIMIT HIT] depth = " << s_trig_combination_depth << "\n";
        return nullptr;
    }
    const TrigCombinationDepthGuard depth_guard;
    const Rational zero(BigInt(0));
    if (ref.denominator() < BigInt(3)) return nullptr;
    auto igcd = [](int a, int b) -> int { while (b) { a %= b; std::swap(a, b); } return a; };
    for (int q1 = 2; q1 <= kBaseAngleMaxDenom; ++q1) {
        for (int p1 = 1; p1 < q1; ++p1) {
            if (igcd(p1, q1) != 1) continue;
            Rational r1{BigInt(p1), BigInt(q1)};
            Rational r2 = r1 - ref;
            if (r2 <= zero) continue;
            if (r2.denominator() >= ref.denominator()) continue;
            ExprPtr c1 = trig_exact_at_pi_multiple(BuiltinOp::Cos, r1, arena);
            ExprPtr s1 = trig_exact_at_pi_multiple(BuiltinOp::Sin, r1, arena);
            ExprPtr c2 = trig_exact_at_pi_multiple(BuiltinOp::Cos, r2, arena);
            ExprPtr s2 = trig_exact_at_pi_multiple(BuiltinOp::Sin, r2, arena);
            if (!c1 || !s1 || !c2 || !s2) continue;
            if (func == BuiltinOp::Cos) {
                ExprPtr res = arena.make<Sum>(std::vector<ExprPtr>{arena.make<Product>(std::vector<ExprPtr>{c1, c2}), arena.make<Product>(std::vector<ExprPtr>{s1, s2})});
                if (is_debug_chebyshev()) std::cerr << "  [SUCCESS] found Cos combo: " << p1 << "/" << q1 << " - " << r2.numerator().decimal() << "/" << r2.denominator().decimal() << "\n";
                return res;
            } else {
                ExprPtr res = arena.make<Sum>(std::vector<ExprPtr>{arena.make<Product>(std::vector<ExprPtr>{s1, c2}), arena.make<Unary>(UnaryOp::Neg, arena.make<Product>(std::vector<ExprPtr>{c1, s2}))});
                if (is_debug_chebyshev()) std::cerr << "  [SUCCESS] found Sin combo: " << p1 << "/" << q1 << " - " << r2.numerator().decimal() << "/" << r2.denominator().decimal() << "\n";
                return res;
            }
        }
    }
    if (is_debug_chebyshev()) std::cerr << "  [FAILED] no combo found for ref = " << ref.numerator().decimal() << "/" << ref.denominator().decimal() << "\n";
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

    // Arc-trig: asin, acos, atan inverse-function rules delegated to
    // simplify_trig_inverse.cpp (split from simplify_trig.cpp for ≤500 LOC).
    if (op == BuiltinOp::Asin || op == BuiltinOp::Acos || op == BuiltinOp::Atan)
        return simplify_funcall_arc_trig(original, op, std::move(args), target_before);

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
