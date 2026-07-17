#include "simplify_impl.hpp"

namespace cas::symbolic::detail {

// Extract {real, imag} from a + b*I form expressions.
// Returns nullopt if no imaginary component detected.
[[nodiscard]] static std::optional<std::pair<ExprPtr, ExprPtr>>
extract_complex_parts(ExprPtr expr, AstArena& arena) {
    ExprPtr zero = arena.make<IntegerLit>(BigInt(0));

    if (const auto* complex = expr_cast<ComplexLit>(expr)) {
        return std::make_pair(
            make_rational(arena, Rational(complex->re_num, complex->re_den)),
            make_rational(arena, Rational(complex->im_num, complex->im_den)));
    }

    auto is_imaginary_unit = [](ExprPtr e) -> bool {
        const auto* c = expr_cast<Constant>(e);
        return c != nullptr && c->value == MathConstant::I;
    };

    auto extract_imag_factor = [&](ExprPtr e) -> ExprPtr {
        if (const auto* complex = expr_cast<ComplexLit>(e)) {
             // (a+bi) is not purely imaginary unless a=0.
             if (complex->re_num.is_zero())
                 return make_rational(arena, Rational(complex->im_num, complex->im_den));
             return nullptr;
        }
        if (is_imaginary_unit(e))
            return arena.make<IntegerLit>(BigInt(1));
        if (const auto* u = expr_cast<Unary>(e);
            u && u->op == UnaryOp::Neg && is_imaginary_unit(u->operand))
            return arena.make<IntegerLit>(BigInt(-1));
        if (const auto* prod = expr_cast<Product>(e)) {
            std::vector<ExprPtr> non_i;
            bool found_i = false;
            for (ExprPtr f : prod->factors) {
                if (!found_i && is_imaginary_unit(f)) { found_i = true; continue; }
                non_i.push_back(f);
            }
            if (found_i) {
                if (non_i.size() == 1U) return non_i[0];
                if (non_i.empty()) return arena.make<IntegerLit>(BigInt(1));
                return arena.make<Product>(std::move(non_i));
            }
        }
        return nullptr;
    };

    if (const auto* sum = expr_cast<Sum>(expr)) {
        ExprPtr real_part = nullptr;
        ExprPtr imag_part = nullptr;
        for (ExprPtr term : sum->terms) {
            ExprPtr maybe_b = extract_imag_factor(term);
            if (maybe_b) imag_part = maybe_b;
            else         real_part = term;
        }
        if (!real_part) real_part = zero;
        if (!imag_part) return std::nullopt;
        return std::make_pair(real_part, imag_part);
    }

    ExprPtr maybe_b = extract_imag_factor(expr);
    if (maybe_b) return std::make_pair(zero, maybe_b);
    
    LiteralRational lr;
    if (try_get_exact_rational(expr, lr).is_ok() && lr.exact) {
        return std::make_pair(expr, zero);
    }

    return std::nullopt;
}

// ── Simplifier::simplify_funcall_complex ──────────────────────────────────────

Result<ExprPtr> Simplifier::simplify_funcall_complex(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before) {

    if (args.size() == 1U) {
        if (op == BuiltinOp::Re) {
            LiteralComplex comp;
            if (try_get_exact_complex(args.front(), comp).is_ok() && comp.exact) {
                return ok(make_rational(arena_, comp.value.real()));
            }
        }
        if (op == BuiltinOp::Im) {
            LiteralComplex comp;
            if (try_get_exact_complex(args.front(), comp).is_ok() && comp.exact) {
                return ok(make_rational(arena_, comp.value.imag()));
            }
        }
        if (op == BuiltinOp::Conj) {
            LiteralComplex comp;
            if (try_get_exact_complex(args.front(), comp).is_ok() && comp.exact) {
                return ok(make_complex(arena_, comp.value.conjugate()));
            }
        }

        auto parts = extract_complex_parts(args.front(), arena_);

        if (op == BuiltinOp::Re) {
            if (parts) return ok(parts->first);
        }

        if (op == BuiltinOp::Im) {
            if (parts) return ok(parts->second);
        }

        if (op == BuiltinOp::Conj && parts) {
            ExprPtr neg_b_i = arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Unary>(UnaryOp::Neg, parts->second),
                arena_.make<Constant>(MathConstant::I)
            });
            auto conj_expr = arena_.make<Sum>(std::vector<ExprPtr>{parts->first, neg_b_i});
            auto simp = simplify_expr(conj_expr);
            if (simp.is_ok()) return simp;
        }

        if (op == BuiltinOp::Abs) {
            if (is_known_nonnegative(args.front()))
                return traced_result(RuleId::SimplifyAbsPositive, target_before, args.front());
            if (is_known_negative(args.front())) {
                auto neg = arena_.make<Unary>(UnaryOp::Neg, args.front());
                auto simplified = simplify_expr(neg);
                if (simplified.is_ok())
                    return traced_result(RuleId::SimplifyAbsNegative, target_before, simplified.value());
            }
            if (const auto* inner = expr_cast<FuncCall>(args.front());
                inner && inner->func_id == BuiltinOp::Abs)
                return traced_result(RuleId::SimplifyAbsAbs, target_before, args.front());
            if (const auto* unary = expr_cast<Unary>(args.front());
                unary && unary->op == UnaryOp::Neg)
                return simplify_expr(arena_.make<FuncCall>(
                    BuiltinOp::Abs, std::vector<ExprPtr>{unary->operand}));
            if (const auto* c = expr_cast<Constant>(args.front());
                c && c->value == MathConstant::I)
                return ok(make_integer(arena_, BigInt(1)));
            // A31 fase 2 (Domain_Conditions_Propagation.md §10.3.R3):
            // abs(b^(2k)) -> b^(2k) for literal positive even exponent 2k,
            // exact when b is real (refused today because e.g. b=i gives
            // b^2 = -1). Opt-in: rewrite and register Real(b).
            if (const auto* power = expr_cast<Binary>(args.front());
                power != nullptr && power->op == BinaryOp::Pow
                && context_ != nullptr && context_->conditional_domain_rules()) {
                const auto* exp_lit = expr_cast<IntegerLit>(power->right);
                if (exp_lit != nullptr && !exp_lit->value.is_negative()
                    && !exp_lit->value.is_zero()
                    && (exp_lit->value % BigInt(2)).is_zero()) {
                    auto cond = context_->emit_side_condition(
                        DomainConditionKind::Real, power->left);
                    if (cond.is_error()) return fail<ExprPtr>(cond.error());
                    return traced_result(RuleId::SimplifyAbsPositive,
                        target_before, args.front());
                }
            }
            // abs(a + b*I) = sqrt(a^2 + b^2)
            if (parts) {
                ExprPtr a_sq = arena_.make<Binary>(BinaryOp::Pow,
                    parts->first, make_integer(arena_, BigInt(2)));
                ExprPtr b_sq = arena_.make<Binary>(BinaryOp::Pow,
                    parts->second, make_integer(arena_, BigInt(2)));
                ExprPtr sum_expr = arena_.make<Sum>(std::vector<ExprPtr>{a_sq, b_sq});
                ExprPtr rooted = arena_.make<FuncCall>(
                    BuiltinOp::Sqrt, std::vector<ExprPtr>{sum_expr});
                return simplify_expr(rooted);
            }
        }

        if (op == BuiltinOp::Arg) {
            if (is_zero_expr(args.front()))
                return ok(make_integer(arena_, BigInt(0)));
            if (is_known_positive(args.front()))
                return ok(make_integer(arena_, BigInt(0)));
            if (is_known_negative(args.front()))
                return ok(arena_.make<Constant>(MathConstant::Pi));
            if (const auto* c = expr_cast<Constant>(args.front());
                c && c->value == MathConstant::I) {
                return simplify_expr(arena_.make<Binary>(
                    BinaryOp::Div,
                    arena_.make<Constant>(MathConstant::Pi),
                    make_integer(arena_, BigInt(2))));
            }
            if (const auto* un = expr_cast<Unary>(args.front());
                un && un->op == UnaryOp::Neg) {
                if (const auto* c = expr_cast<Constant>(un->operand);
                    c && c->value == MathConstant::I) {
                    return simplify_expr(arena_.make<Unary>(
                        UnaryOp::Neg,
                        arena_.make<Binary>(BinaryOp::Div,
                            arena_.make<Constant>(MathConstant::Pi),
                            make_integer(arena_, BigInt(2)))));
                }
            }
            if (parts) {
                ExprPtr a = parts->first;
                ExprPtr b = parts->second;
                if (is_known_positive(a)) {
                    ExprPtr ratio = arena_.make<Binary>(BinaryOp::Div, b, a);
                    return simplify_expr(arena_.make<FuncCall>(
                        BuiltinOp::Atan, std::vector<ExprPtr>{ratio}));
                }
                if (is_zero_expr(a)) {
                    if (is_known_positive(b))
                        return simplify_expr(arena_.make<Binary>(
                            BinaryOp::Div,
                            arena_.make<Constant>(MathConstant::Pi),
                            make_integer(arena_, BigInt(2))));
                    if (is_known_negative(b))
                        return simplify_expr(arena_.make<Unary>(
                            UnaryOp::Neg,
                            arena_.make<Binary>(BinaryOp::Div,
                                arena_.make<Constant>(MathConstant::Pi),
                                make_integer(arena_, BigInt(2)))));
                }
                if (is_known_negative(a)) {
                    ExprPtr ratio = arena_.make<Binary>(BinaryOp::Div, b, a);
                    ExprPtr atan_expr = arena_.make<FuncCall>(
                        BuiltinOp::Atan, std::vector<ExprPtr>{ratio});
                    ExprPtr pi_expr = arena_.make<Constant>(MathConstant::Pi);
                    if (is_known_nonnegative(b))
                        return simplify_expr(arena_.make<Sum>(
                            std::vector<ExprPtr>{atan_expr, pi_expr}));
                    if (is_known_negative(b))
                        return simplify_expr(arena_.make<Sum>(std::vector<ExprPtr>{
                            atan_expr, arena_.make<Unary>(UnaryOp::Neg, pi_expr)}));
                }
            }
        }

        if (op == BuiltinOp::Sign) {
            if (is_known_positive(args.front()))
                return traced_result(RuleId::SimplifySignPositive, target_before,
                    make_integer(arena_, BigInt(1)));
            if (is_known_negative(args.front()))
                return traced_result(RuleId::SimplifySignNegative, target_before,
                    make_integer(arena_, BigInt(-1)));
            if (is_zero_expr(args.front()))
                return traced_result(RuleId::SimplifySignZero, target_before,
                    make_integer(arena_, BigInt(0)));
        }
    }

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

} // namespace cas::symbolic::detail
