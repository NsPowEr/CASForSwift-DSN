#include "simplify_impl.hpp"
#include "simplify_branch_cut.hpp"

namespace cas::symbolic::detail {

// ── Simplifier::simplify_funcall_exp_log_sqrt ─────────────────────────────────
//
// W9.3 split: the Sqrt branch (and its rational-sqrt / denesting helpers) was
// extracted to simplify_sqrt.cpp + simplify_sqrt_helpers.cpp to keep this
// translation unit under the 500-line anti-monolith limit. This function now
// handles Exp and Ln/Log directly and delegates Sqrt to simplify_funcall_sqrt.

Result<ExprPtr> Simplifier::simplify_funcall_exp_log_sqrt(
    ExprPtr original, BuiltinOp op, std::vector<ExprPtr> args, ExprPtr target_before) {

    if (op == BuiltinOp::Exp && args.size() == 1U) {
        if (is_zero_expr(args.front()))
            return traced_result(RuleId::SimplifyExpZero, target_before, make_integer(arena_, BigInt(1)));
        if (is_one_expr(args.front()))
            return traced_result(RuleId::SimplifyExpOne, target_before, make_constant(arena_, MathConstant::E));
        // exp(ln(x)) = x is valid ONLY for x > 0 (principal branch of
        // ln has a cut along the negative real axis). Pre-fix this
        // rule was applied unconditionally — wrong for symbolic x.
        // Reference: Bronstein "Symbolic Integration" §3.3.
        if (const auto* ln_call = expr_cast<FuncCall>(args.front());
            ln_call && ln_call->func_id == BuiltinOp::Ln && ln_call->args.size() == 1U) {
            ExprPtr ln_arg = ln_call->args[0];
            if (is_known_positive(ln_arg)) {
                return ok(ln_arg);
            }
            // Otherwise keep symbolic exp(ln(arg)).
        }
        // exp(c · ln(x)) → x^c when x > 0 (principal branch).
        // L2-19 (positivity gating) extended to scaled-log form. The Product
        // argument must contain exactly one FuncCall(Ln, x) with x > 0; the
        // remaining factors form the exponent c. Verified by re-exponentiation
        // identity exp(c·ln(x)) = x^c on the principal branch (Bronstein §3.3).
        if (const auto* prod_arg = expr_cast<Product>(args.front())) {
            ExprPtr ln_arg_inner = nullptr;
            std::vector<ExprPtr> c_factors;
            int ln_count = 0;
            for (ExprPtr f : prod_arg->factors) {
                if (const auto* lc = expr_cast<FuncCall>(f);
                    lc && lc->func_id == BuiltinOp::Ln && lc->args.size() == 1U) {
                    ln_arg_inner = lc->args[0];
                    ++ln_count;
                } else {
                    c_factors.push_back(f);
                }
            }
            if (ln_count == 1 && ln_arg_inner && is_known_positive(ln_arg_inner)) {
                ExprPtr c_expr;
                if (c_factors.empty()) c_expr = make_integer(arena_, BigInt(1));
                else if (c_factors.size() == 1U) c_expr = c_factors[0];
                else c_expr = arena_.make<Product>(std::move(c_factors));
                return simplify_expr(arena_.make<Binary>(BinaryOp::Pow, ln_arg_inner, c_expr));
            }
        }
        // exp(I · θ) → cos(θ) + I·sin(θ) (Euler's formula, principal branch).
        // Detected when argument is a Product containing exactly one factor
        // equal to MathConstant::I; the remaining factors compose θ. Safe on
        // the principal branch for any θ ∈ R; if θ has a complex part it
        // further splits via exp(α+iβ) = exp(α)(cos β + i sin β).
        if (const auto* prod_arg = expr_cast<Product>(args.front())) {
            int i_count = 0;
            bool neg_unit = false;
            std::vector<ExprPtr> theta_factors;
            for (ExprPtr f : prod_arg->factors) {
                if (const auto* cc = expr_cast<Constant>(f);
                    cc && cc->value == MathConstant::I) {
                    ++i_count;
                    continue;
                }
                // A purely imaginary coefficient captures any scaled factor of I.
                if (const auto* cl = expr_cast<ComplexLit>(f);
                    cl != nullptr && cl->re_num.is_zero() && !cl->im_num.is_zero()) {
                    ++i_count;
                    if (cl->im_num == BigInt(-1) && cl->im_den == BigInt(1)) {
                        neg_unit = !neg_unit;
                    } else if (!(cl->im_num == BigInt(1) && cl->im_den == BigInt(1))) {
                        if (cl->im_den == BigInt(1)) {
                            theta_factors.push_back(make_integer(arena_, cl->im_num));
                        } else {
                            theta_factors.push_back(arena_.make<RationalLit>(cl->im_num, cl->im_den));
                        }
                    }
                    continue;
                }
                theta_factors.push_back(f);
            }
            if (i_count == 1) {
                if (neg_unit) {
                    theta_factors.push_back(make_integer(arena_, BigInt(-1)));
                }
                ExprPtr theta;
                if (theta_factors.empty()) theta = make_integer(arena_, BigInt(1));
                else if (theta_factors.size() == 1U) theta = theta_factors[0];
                else theta = arena_.make<Product>(std::move(theta_factors));
                ExprPtr cos_theta = arena_.make<FuncCall>(BuiltinOp::Cos,
                    std::vector<ExprPtr>{theta});
                ExprPtr sin_theta = arena_.make<FuncCall>(BuiltinOp::Sin,
                    std::vector<ExprPtr>{theta});
                ExprPtr i_sin = arena_.make<Product>(std::vector<ExprPtr>{
                    arena_.make<Constant>(MathConstant::I), sin_theta});
                return simplify_expr(arena_.make<Sum>(std::vector<ExprPtr>{
                    cos_theta, i_sin}));
            }
        }
        if (is_constant_expr(args.front(), MathConstant::Infinity))
            return traced_result(RuleId::Unknown, target_before, make_constant(arena_, MathConstant::Infinity));
        if (expr_is<Unary>(args.front())
            && expr_ref<Unary>(args.front()).op == UnaryOp::Neg
            && is_constant_expr(expr_ref<Unary>(args.front()).operand, MathConstant::Infinity))
            return traced_result(RuleId::Unknown, target_before, make_integer(arena_, BigInt(0)));
        if (const auto* sum = expr_cast<Sum>(args.front())) {
            std::vector<ExprPtr> factors;
            for (ExprPtr term : sum->terms)
                factors.push_back(arena_.make<FuncCall>(BuiltinOp::Exp,
                    std::vector<ExprPtr>{term}));
            auto rewritten = simplify_product_factors(factors, arena_.make<Product>(factors));
            if (rewritten.is_ok()) {
                append_trace(RuleId::SimplifyExpSum, target_before, rewritten.value());
                return rewritten;
            }
        }
    }

    if ((op == BuiltinOp::Ln || op == BuiltinOp::Log) && args.size() == 1U) {
        if (is_zero_expr(args.front()))
            return fail<ExprPtr>(make_error(CASErrorKind::Undefined, "ln(0) is undefined"));
        if (is_one_expr(args.front()))
            return traced_result(RuleId::SimplifyLnOne, target_before, make_integer(arena_, BigInt(0)));
        if (is_constant_expr(args.front(), MathConstant::E))
            return traced_result(RuleId::SimplifyLnE, target_before, make_integer(arena_, BigInt(1)));
        if (is_constant_expr(args.front(), MathConstant::Infinity))
            return traced_result(RuleId::Unknown, target_before, make_constant(arena_, MathConstant::Infinity));

        // F8.0-6.2: branch-cut bookkeeping for ln(exp(z)).
        // Default (legacy):  ln(exp(z)) → z  unconditionally.
        // When ctx.branch_cut_aware_logexp() == true AND z is not known real:
        //   ln(exp(z)) → z + 2πi·K(z)
        // K(z) is preserved verbatim by the simplifier (Task 6.1).
        auto build_branch_aware_logexp = [&](ExprPtr z) -> ExprPtr {
            const bool aware = context_ != nullptr
                            && context_->branch_cut_aware_logexp();
            if (!aware) return z;
            const bool z_is_real = context_ != nullptr
                                && context_->assumptions().is_real(z);
            if (z_is_real) return z;
            ExprPtr K_z = arena_.make<FuncCall>(BuiltinOp::UnwindingNumber,
                std::vector<ExprPtr>{z});
            ExprPtr two_pi_i = arena_.make<Product>(std::vector<ExprPtr>{
                make_integer(arena_, BigInt(2)),
                make_constant(arena_, MathConstant::Pi),
                make_constant(arena_, MathConstant::I)});
            ExprPtr offset = arena_.make<Binary>(BinaryOp::Mul, two_pi_i, K_z);
            return arena_.make<Binary>(BinaryOp::Add, z, offset);
        };
        // ln(e^x) -> x   (or x + 2πi·K(x) when branch-aware)
        if (const auto* power = expr_cast<Binary>(args.front());
            power != nullptr && power->op == BinaryOp::Pow
            && is_constant_expr(power->left, MathConstant::E))
            return traced_result(RuleId::SimplifyLnExp, target_before,
                build_branch_aware_logexp(power->right));
        // ln(exp(x)) -> x  (or x + 2πi·K(x) when branch-aware)
        if (const auto* exp_call = expr_cast<FuncCall>(args.front());
            exp_call && exp_call->func_id == BuiltinOp::Exp && exp_call->args.size() == 1U)
            return traced_result(RuleId::SimplifyLnExp, target_before,
                build_branch_aware_logexp(exp_call->args[0]));

        // ln(a*b) -> ln(a) + ln(b) for a,b > 0
        // F8.0-6.2 / Task 20 BC-3 (Branch_Cut_Propagation.md §2 rule 4):
        //   ln(z1·z2) = ln(z1) + ln(z2) - 2πi · K(ln(z1) + ln(z2))
        // When ctx.strict_branch_cuts() is on and at least one factor is not
        // provably positive, emit the full unwinding-corrected form instead
        // of the unsafe positive-only reduction.
        if (const auto* prod = expr_cast<Product>(args.front())) {
            bool all_pos = true;
            for (auto f : prod->factors) if (!is_known_positive(f)) { all_pos = false; break; }
            if (all_pos) {
                std::vector<ExprPtr> ln_factors;
                for (auto f : prod->factors) {
                    auto res = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln,
                        std::vector<ExprPtr>{f}));
                    if (res.is_error()) return res;
                    ln_factors.push_back(res.value());
                }
                return simplify_expr(arena_.make<Sum>(std::move(ln_factors)));
            }
            const bool strict = (context_ != nullptr) && context_->strict_branch_cuts();
            if (strict && prod->factors.size() >= 2U) {
                // Build  Σ ln(fᵢ)  +  Σ_{i<j}  −2πi · K(ln(fᵢ) + ln(fⱼ))
                // The pairwise K(·) terms encode the multi-factor unwinding by
                // associativity:  K(a+b+c) telescopes through partial sums.
                std::vector<ExprPtr> terms;
                terms.reserve(prod->factors.size() + (prod->factors.size() - 1U));
                for (auto f : prod->factors) {
                    terms.push_back(arena_.make<FuncCall>(BuiltinOp::Ln,
                        std::vector<ExprPtr>{f}));
                }
                // Pairwise correction Σᵢ −2πi · K(ln(fᵢ) + ln(fᵢ₊₁))
                for (std::size_t i = 0; i + 1U < prod->factors.size(); ++i) {
                    terms.push_back(symbolic::branch_cut::make_log_product_correction(
                        prod->factors[i], prod->factors[i + 1U], arena_));
                }
                return ok(arena_.make<Sum>(std::move(terms)));
            }
        }
        // ln(z1 / z2) -> ln(z1) - ln(z2) (or with corrections)
        // F8.0-6.2 / Task 20 BC-3b (Branch_Cut_Propagation.md §2 rule 5):
        //   ln(z1 / z2) = ln(z1) - ln(z2) - 2πi · K(ln(z1) - ln(z2))
        if (const auto* div = expr_cast<Binary>(args.front());
            div != nullptr && div->op == BinaryOp::Div) {
            bool all_pos = is_known_positive(div->left) && is_known_positive(div->right);
            if (all_pos) {
                auto ln_a = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{div->left}));
                if (ln_a.is_error()) return ln_a;
                auto ln_b = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{div->right}));
                if (ln_b.is_error()) return ln_b;
                return simplify_expr(arena_.make<Binary>(BinaryOp::Sub, ln_a.value(), ln_b.value()));
            }
            const bool strict = (context_ != nullptr) && context_->strict_branch_cuts();
            if (strict) {
                ExprPtr ln_a = arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{div->left});
                ExprPtr ln_b = arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{div->right});
                ExprPtr diff = arena_.make<Binary>(BinaryOp::Sub, ln_a, ln_b);
                ExprPtr correction = symbolic::branch_cut::make_log_quotient_correction(div->left, div->right, arena_);
                return simplify_expr(arena_.make<Binary>(BinaryOp::Add, diff, correction));
            }
            // Legacy / non-strict default: ln(z1) - ln(z2)
            auto ln_a = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{div->left}));
            if (ln_a.is_error()) return ln_a;
            auto ln_b = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{div->right}));
            if (ln_b.is_error()) return ln_b;
            return simplify_expr(arena_.make<Binary>(BinaryOp::Sub, ln_a.value(), ln_b.value()));
        }

        // ln(sqrt(x)) = (1/2)*ln(x) — identità esatta
        if (const auto* sqrt_call = expr_cast<FuncCall>(args.front());
            sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt
            && sqrt_call->args.size() == 1U) {
            ExprPtr half = make_rational(arena_, Rational(BigInt(1), BigInt(2)));
            ExprPtr ln_inner = arena_.make<FuncCall>(BuiltinOp::Ln, sqrt_call->args);
            return simplify_expr(arena_.make<Binary>(BinaryOp::Mul, half, ln_inner));
        }
        // Branch cut principal value: ln(-x) = ln(x) + I·π for x > 0
        if (const auto* neg = expr_cast<Unary>(args.front());
            neg != nullptr && neg->op == UnaryOp::Neg) {
            ExprPtr inner = neg->operand;
            bool inner_pos = is_known_positive(inner) || is_constant_expr(inner, MathConstant::E);
            if (!inner_pos) {
                LiteralRational rat;
                auto ex = try_get_exact_rational(inner, rat);
                if (ex.is_ok() && ex.value()
                    && !rat.value.numerator().is_negative()
                    && !rat.value.numerator().is_zero())
                    inner_pos = true;
            }
            if (inner_pos) {
                ExprPtr ln_inner;
                if (is_one_expr(inner)) {
                    ln_inner = make_integer(arena_, BigInt(0));
                } else {
                    auto r = simplify_expr(arena_.make<FuncCall>(BuiltinOp::Ln,
                        std::vector<ExprPtr>{inner}));
                    if (r.is_error()) return r;
                    ln_inner = r.value();
                }
                ExprPtr i_pi = arena_.make<Binary>(BinaryOp::Mul,
                    make_constant(arena_, MathConstant::I),
                    make_constant(arena_, MathConstant::Pi));
                return simplify_expr(arena_.make<Binary>(BinaryOp::Add, ln_inner, i_pi));
            }
        }
        // Complex principal branch: ln(I) = I·π/2
        if (const auto* c = expr_cast<Constant>(args.front()); c && c->value == MathConstant::I) {
            ExprPtr i_pi = arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Constant>(MathConstant::I),
                arena_.make<Constant>(MathConstant::Pi)});
            return simplify_expr(arena_.make<Binary>(
                BinaryOp::Div, i_pi, make_integer(arena_, BigInt(2))));
        }
        // Complex principal branch: ln(-I) = -I·π/2.
        // Handled separately from ln(-x) for x>0 because i is neither
        // positive nor negative real; this is the unique branch-cut
        // edge case at arg(-I) = -π/2.
        if (const auto* un = expr_cast<Unary>(args.front());
            un && un->op == UnaryOp::Neg) {
            const auto* c = expr_cast<Constant>(un->operand);
            if (c && c->value == MathConstant::I) {
                ExprPtr neg_i_pi = arena_.make<Product>(std::vector<ExprPtr>{
                    arena_.make<IntegerLit>(BigInt(-1)),
                    arena_.make<Constant>(MathConstant::I),
                    arena_.make<Constant>(MathConstant::Pi)});
                return simplify_expr(arena_.make<Binary>(
                    BinaryOp::Div, neg_i_pi, make_integer(arena_, BigInt(2))));
            }
        }
        // Post-F1.6: arguments may arrive as ComplexLit after canonicalization.
        // HC-F16-LN-COMPLEX-FULL: General formula ln(a + b·I) = ln|z| + I·arg(z).
        // Naturally handles all complex literals, delegating to robust Abs and Arg.
        if (expr_is<ComplexLit>(args.front())) {
            ExprPtr z = args.front();
            ExprPtr abs_z = arena_.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{z});
            ExprPtr arg_z = arena_.make<FuncCall>(BuiltinOp::Arg, std::vector<ExprPtr>{z});
            ExprPtr ln_abs = arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{abs_z});
            ExprPtr i_arg = arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Constant>(MathConstant::I), arg_z});
            return simplify_expr(arena_.make<Sum>(std::vector<ExprPtr>{ln_abs, i_arg}));
        }
        // ln(-1) = I·π  (also covered by ln(-x) above, but explicit for clarity)
        if (const auto* il = expr_cast<IntegerLit>(args.front()); il && il->value == BigInt(-1)) {
            return simplify_expr(arena_.make<Product>(std::vector<ExprPtr>{
                arena_.make<Constant>(MathConstant::I),
                arena_.make<Constant>(MathConstant::Pi)}));
        }
        // L2-08: ln(a + b·I) = ln|z| + I·arg(z)  (principal branch).
        // Triggered only when the argument is a Sum containing an imaginary term b·I.
        // abs(z) and arg(z) dispatch back to simplify_funcall_complex which already
        // handles the quadrant logic — no duplication.
        if (const auto* sum_z = expr_cast<Sum>(args.front())) {
            auto is_i_unit = [](ExprPtr e) -> bool {
                const auto* c = expr_cast<Constant>(e);
                return c != nullptr && c->value == MathConstant::I;
            };
            auto extract_imag_coeff = [&](ExprPtr term) -> ExprPtr {
                // Post-F1.6 a standalone `i` canonicalizes to ComplexLit(0,1),
                // so a Sum like √3 + i carries a pure-imaginary ComplexLit term.
                // Recognize it (re=0) and return its imaginary coefficient, matching
                // extract_complex_parts in simplify_complex.cpp. (A ComplexLit with
                // nonzero real part is fully handled by the ln(ComplexLit) branch
                // above and never reaches this Sum path.)
                if (const auto* cl = expr_cast<ComplexLit>(term)) {
                    if (cl->re_num.is_zero() && !cl->im_num.is_zero())
                        return make_rational(arena_, Rational(cl->im_num, cl->im_den));
                    return nullptr;
                }
                if (is_i_unit(term)) return make_integer(arena_, BigInt(1));
                if (const auto* prod = expr_cast<Product>(term)) {
                    bool found_i = false;
                    std::vector<ExprPtr> others;
                    for (ExprPtr f : prod->factors) {
                        if (!found_i && is_i_unit(f)) { found_i = true; continue; }
                        others.push_back(f);
                    }
                    if (found_i) {
                        if (others.empty()) return make_integer(arena_, BigInt(1));
                        if (others.size() == 1U) return others[0];
                        return arena_.make<Product>(std::move(others));
                    }
                }
                return nullptr;
            };
            ExprPtr real_part = nullptr, imag_part = nullptr;
            for (ExprPtr term : sum_z->terms) {
                ExprPtr b = extract_imag_coeff(term);
                if (b && !imag_part) imag_part = b;
                else if (!b && !real_part) real_part = term;
            }
            if (real_part && imag_part) {
                // ln(a + b·I) = ln(abs(a+b·I)) + I·arg(a+b·I)
                ExprPtr z = args.front();
                ExprPtr abs_z = arena_.make<FuncCall>(BuiltinOp::Abs, std::vector<ExprPtr>{z});
                ExprPtr arg_z = arena_.make<FuncCall>(BuiltinOp::Arg, std::vector<ExprPtr>{z});
                ExprPtr ln_abs = arena_.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{abs_z});
                ExprPtr i_arg = arena_.make<Product>(std::vector<ExprPtr>{
                    arena_.make<Constant>(MathConstant::I), arg_z});
                ExprPtr result = arena_.make<Sum>(std::vector<ExprPtr>{ln_abs, i_arg});
                return simplify_expr(result);
            }
        }
    }

    if (op == BuiltinOp::Sqrt && args.size() == 1U)
        return simplify_funcall_sqrt(original, std::move(args), target_before);

    const auto& orig_args = expr_ref<FuncCall>(original).args;
    if (expr_ptr_sequence_identical(args, orig_args)) return ok(original);
    return ok(arena_.make<FuncCall>(op, std::move(args)));
}

} // namespace cas::symbolic::detail
