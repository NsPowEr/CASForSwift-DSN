// HC-F8-MONOLITH-WAIVER tier-1 split (2026-06-11): estratto da
// simplify_arithmetic.cpp.  Dispatch Pow nodes; dependencies via simplify_impl.hpp.

#include "simplify_impl.hpp"
#include "simplify_branch_cut.hpp"
#include "../algebra/polynomial_internal.hpp"

namespace cas::symbolic::detail {

// ── Trig power linearization (Chebyshev / DeMoivre) ──────────────────────

// Returns linearized form of sin^n or cos^n using multiple-angle identities.
// Odd n=2m+1:  trig^n = (1/4^m) * sum_{j=0}^{m} c_j * trig((n-2j)*arg)
// Even n=2m:   trig^n = (1/4^m) * [C(n,m) + 2*sum_{j=0}^{m-1} c_j * cos((n-2j)*arg)]

Result<ExprPtr> Simplifier::simplify_power(ExprPtr base, ExprPtr exponent, ExprPtr target_before) {
    if (!target_before && trace_enabled_) {
        target_before = arena_.make<Binary>(BinaryOp::Pow, base, exponent);
    }

    // F7.5.F1 Phase 2 — extended-real arithmetic propagation for Pow.
    // Covers ±∞^n, ComplexInf^n, 0^0, 1^∞, Indeterminate^x, x^Indeterminate.
    // Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Extended_Real_AST.md
    if (auto ext = try_simplify_pow_extended_real(base, exponent, arena_); ext) {
        return traced_result(RuleId::SimplifyPowerOne, target_before, *ext);
    }

    if (rewrite_provider_ != nullptr && may_rewrite_power(base, exponent)) {
        ExprPtr rewrite_target = target_before ? target_before : arena_.make<Binary>(BinaryOp::Pow, base, exponent);
        auto rewritten = rewrite_provider_->try_rewrite(rewrite_target, arena_, assumptions_, context_);
        if (rewritten.is_error()) return rewritten;
        if (rewritten.value() != rewrite_target) {
            append_trace(RuleId::RewriteProviderApplied, rewrite_target, rewritten.value());
            return simplify_expr(rewritten.value());
        }
    }

    LiteralRational exp_rat;
    LiteralComplex base_comp;
    auto b_c_exact = try_get_exact_complex(base, base_comp);
    if (b_c_exact.is_error()) return fail<ExprPtr>(b_c_exact.error());
    auto e_exact = try_get_exact_rational(exponent, exp_rat);
    if (e_exact.is_error()) return fail<ExprPtr>(e_exact.error());

    if (b_c_exact.value() && e_exact.value()) {
        if (base_comp.value.is_zero() && exp_rat.value.numerator().is_zero()) {
            return ok(arena_.make<Binary>(BinaryOp::Pow, base, exponent));
        }
        if (exp_rat.value.is_integer()) {
            const BigInt power = exp_rat.value.numerator();
            if (power.is_zero()) return ok(make_integer(arena_, BigInt(1)));

            ComplexRational res = ComplexRational::one();
            ComplexRational b = base_comp.value;
            BigInt p = power.is_negative() ? -power : power;

            while (!p.is_zero()) {
                if ((p % BigInt(2)) == BigInt(1)) res = res * b;
                p /= BigInt(2);
                if (!p.is_zero()) b = b * b;
            }

            if (power.is_negative()) {
                auto inv_res = ComplexRational::one().divide(res);
                if (inv_res.is_error()) return fail<ExprPtr>(inv_res.error());
                return ok(make_complex(arena_, inv_res.value()));
            }
            return ok(make_complex(arena_, res));
        }
    }

    if (is_constant_expr(base, MathConstant::Infinity) && e_exact.value() && exp_rat.value.is_integer()) {
        const BigInt power = exp_rat.value.numerator();
        if (power.is_negative()) {
            return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(0)));
        }
        if (!power.is_zero()) {
            return traced_result(RuleId::SimplifyPowerOne, target_before, base);
        }
    }

    // x^0 = 1 is valid only when x ≠ 0. The limit 0^0 is mathematically
    // indeterminate (lim_{x→0+} x^x = 1, but 0^0 has no agreed value).
    // We keep 0^0 symbolic — the rational fast-path above already does
    // this for the literal case; here we guard the structural case
    // where base is structurally zero but not necessarily literal.
    if (is_zero_expr(exponent)) {
        if (is_zero_expr(base)) {
            // Keep Pow(0, 0) symbolic; caller (e.g. limit engine via
            // L'Hôpital or Gruntz) decides the indeterminate form.
            return ok(arena_.make<Binary>(BinaryOp::Pow, base, exponent));
        }
        return traced_result(RuleId::SimplifyPowerZero, target_before, make_integer(arena_, BigInt(1)));
    }
    if (is_one_expr(exponent)) return traced_result(RuleId::SimplifyPowerOne, target_before, base);
    if (is_one_expr(base)) return ok(base);

    // (−x)^n for integer n: even → x^n, odd → −(x^n). Letting the sign escape the
    // power lets downstream consumers collapse e.g. (−√3)² → 3 (needed by abs(z),
    // arg(z) and the trig special-value tables). Exact for any integer n; the I /
    // −I fast-paths below still apply to x^n via the recursive call.
    if (const auto* neg_base = expr_cast<Unary>(base);
        neg_base != nullptr && neg_base->op == UnaryOp::Neg) {
        if (auto n = try_get_integer_exponent(exponent); n.has_value()) {
            auto inner = simplify_power(neg_base->operand, exponent);
            if (inner.is_error()) return inner;
            if ((*n % BigInt(2)) == BigInt(0))
                return traced_result(RuleId::SimplifyPowerOne, target_before, inner.value());
            return traced_result(RuleId::SimplifyPowerOne, target_before,
                arena_.make<Unary>(UnaryOp::Neg, inner.value()));
        }
    }

    if (e_exact.value() && exp_rat.value.is_integer()) {
        // Chebyshev/DeMoivre linearization: sin/cos^n → multiple-angle sum (n ≥ 2).
        // Self-contained; extracted to simplify_arithmetic_power_trig.cpp.
        if (auto trig = try_linearize_trig_power(base, exp_rat.value.numerator())) {
            return *trig;
        }
    }

    // RootOf algebraic reduction: RootOf(P, x, k)^n -> reduce using P(R) = 0
    if (const auto* root = expr_cast<RootOf>(base)) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value() && context_) {
            const BigInt n_val = *maybe_n;
            if (n_val > BigInt(0)) {
                auto poly_res = cas::algebra::parse_polynomial(root->polynomial, root->variable, *context_);
                if (poly_res.is_ok()) {
                    auto poly = poly_res.value();
                    std::size_t d = cas::algebra::poly_degree(poly);
                    if (d > 0 && n_val >= BigInt(static_cast<long long>(d))) {
                        // Create monomial x^n
                        cas::algebra::PolyExpr x_n = cas::algebra::poly_make_monomial(make_integer(arena_, BigInt(1)), static_cast<std::size_t>(n_val.to_u64()));
                        auto div_res = cas::algebra::divide_poly_with_remainder(x_n, poly, *context_);
                        if (div_res.is_ok()) {
                            auto reduced_poly_expr = cas::algebra::polynomial_to_expr(div_res.value().remainder, root->variable, *context_);
                            if (reduced_poly_expr.is_ok()) {
                                // Sostituisci Symbol(var) con RootOf
                                auto result = substitute(reduced_poly_expr.value(), root->variable, base, *context_);
                                if (result.is_ok()) {
                                    append_trace(RuleId::Unknown, target_before, result.value());
                                    return simplify_expr(result.value());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Quantity power: Quantity(v, d)^n -> Quantity(v^n, d*n) for integer n
    if (const auto* q = expr_cast<Quantity>(base)) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            long long n_val = n.to_u64();
            if (n.is_negative()) n_val = -static_cast<long long>((-n).to_u64());

            SIDimensions dims = q->dimensions;
            dims.m = static_cast<int16_t>(dims.m * n_val);
            dims.kg = static_cast<int16_t>(dims.kg * n_val);
            dims.s = static_cast<int16_t>(dims.s * n_val);
            dims.A = static_cast<int16_t>(dims.A * n_val);
            dims.K = static_cast<int16_t>(dims.K * n_val);
            dims.mol = static_cast<int16_t>(dims.mol * n_val);
            dims.cd = static_cast<int16_t>(dims.cd * n_val);

            auto inner_pow = simplify_power(q->value, exponent);
            if (inner_pow.is_error()) return inner_pow;

            if (dims.is_dimensionless()) return inner_pow;
            return ok(arena_.make<Quantity>(inner_pow.value(), dims));
        }
    }

    // x^(1/2) → sqrt(x)  (canonical form)
    {
        LiteralRational exp_half;
        auto e_exact_half = try_get_exact_rational(exponent, exp_half);
        if (e_exact_half.is_ok() && e_exact_half.value() &&
            exp_half.value.numerator() == BigInt(1) && exp_half.value.denominator() == BigInt(2)) {
            ExprPtr sqrt_expr = arena_.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{base});
            return simplify_expr(sqrt_expr);
        }
    }
    // HC-KV-04: x^(p/2), p odd → sqrt(x)^p  (p<0: 1/sqrt(x)^|p|).
    {
        LiteralRational eq;
        auto eqr = try_get_exact_rational(exponent, eq);
        if (eqr.is_ok() && eqr.value()
            && eq.value.denominator() == BigInt(2) && !eq.value.is_integer()) {
            const BigInt p = eq.value.numerator();
            ExprPtr sb = arena_.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{base});
            ExprPtr out = p.is_negative()
                ? arena_.make<Binary>(BinaryOp::Div, make_integer(arena_, BigInt(1)),
                    arena_.make<Binary>(BinaryOp::Pow, sb,
                        make_integer(arena_, -p)))
                : arena_.make<Binary>(BinaryOp::Pow, sb, make_integer(arena_, p));
            return simplify_expr(out);
        }
    }
    if (is_zero_expr(base)) {
        // 0^x = 0 whenever x is PROVABLY positive (literal or assumption): exact
        // and unconditional (no A31 side-condition). Covers 0^(1/2)=0; 0^0 and
        // 0^(negative) stay symbolic (dedicated branches / conditional path below).
        if (is_known_positive(exponent)) {
            return traced_result(RuleId::SimplifyZeroPowerPositive, target_before,
                make_integer(arena_, BigInt(0)));
        }
        LiteralRational exp_rat_check;
        auto exp_check = try_get_exact_rational(exponent, exp_rat_check);
        if (exp_check.is_ok() && exp_check.value() && exp_rat_check.value.is_integer() && !exp_rat_check.value.numerator().is_negative() && !exp_rat_check.value.numerator().is_zero()) {
            return traced_result(RuleId::SimplifyZeroPowerPositive, target_before, make_integer(arena_, BigInt(0)));
        }
        // A31 fase 2 (Domain_Conditions_Propagation.md §10.3.R4): 0^e -> 0
        // for a SYMBOLIC exponent, exact when Re(e) > 0 (§1). The vocabulary
        // has no real-part predicate (§3.2), so the stronger Positive(e) is
        // registered — sound, possibly over-restrictive. Literal exponents
        // stay with the exact branch above; a provably negative or zero
        // exponent keeps the refusal (contradiction guard; 0^0 and 0^-k are
        // handled by the dedicated branches elsewhere in this function).
        if (context_ != nullptr && context_->conditional_domain_rules()
            && !(exp_check.is_ok() && exp_check.value())
            && !is_zero_expr(exponent) && !is_known_negative(exponent)) {
            auto cond = context_->emit_side_condition(
                DomainConditionKind::Positive, exponent);
            if (cond.is_error()) return fail<ExprPtr>(cond.error());
            return traced_result(RuleId::SimplifyZeroPowerPositive,
                target_before, make_integer(arena_, BigInt(0)));
        }
    }

    if (is_constant_expr(base, MathConstant::E)) {
        // E^(a + b + ...) = E^a * E^b * ...  (linearizzazione esponenziale)
        if (const auto* sum_exp = expr_cast<Sum>(exponent); sum_exp != nullptr && !sum_exp->terms.empty()) {
            std::vector<ExprPtr> factors;
            factors.reserve(sum_exp->terms.size());
            for (ExprPtr term : sum_exp->terms) {
                factors.push_back(arena_.make<Binary>(BinaryOp::Pow, base, term));
            }
            return simplify_expr(arena_.make<Product>(std::move(factors)));
        }

        // E^(ln(x)) = x: valido SOLO se x è dimostrabile positivo.
        // Per x simbolico senza assumption, exp(ln(x)) NON è uguale a x
        // nel campo complesso (branch cut su ramo principale del logaritmo).
        // Pre-fix (bug): logica invertita `if (!is_known_nonpositive)`
        // applicava la cancellazione anche a simboli ignoti.
        // Riferimento math: Bronstein "Symbolic Integration" §3.3.
        const auto* call = expr_cast<FuncCall>(exponent);
        if (call != nullptr
            && (call->func_id == BuiltinOp::Ln || call->func_id == BuiltinOp::Log)
            && call->args.size() == 1U) {
            ExprPtr arg = call->args.front();
            if (is_known_positive(arg)) {
                return traced_result(RuleId::SimplifyExpLnPositive, target_before, arg);
            }
            // A31 fase 2 (Domain_Conditions_Propagation.md §10.3.R1): on the
            // principal branch E^(ln z) = z is exact for EVERY z != 0 (ln is
            // a right inverse of exp, cut included). Opt-in: rewrite and
            // register NonZero(z). Mirrors the FuncCall(Exp) site in
            // simplify_exp_log.cpp.
            if (context_ != nullptr && context_->conditional_domain_rules()
                && !is_zero_expr(arg)) {
                auto cond = context_->emit_side_condition(
                    DomainConditionKind::NonZero, arg);
                if (cond.is_error()) return fail<ExprPtr>(cond.error());
                return traced_result(RuleId::SimplifyExpLnPositive,
                    target_before, arg);
            }
            // Altrimenti: mantieni forma simbolica E^(ln(arg)).
        }
    }

    // I^n simplification
    if (const auto* constant = expr_cast<Constant>(base); constant != nullptr && constant->value == MathConstant::I) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            BigInt rem = n % BigInt(4);
            if (rem.is_negative()) rem += BigInt(4);
            const long long r = rem.to_u64();
            if (r == 0) return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(1)));
            if (r == 1) return traced_result(RuleId::SimplifyPowerOne, target_before, base);
            if (r == 2) return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(-1)));
            if (r == 3) return traced_result(RuleId::SimplifyPowerOne, target_before, arena_.make<Unary>(UnaryOp::Neg, base));
        }
    }
    if (const auto* unary = expr_cast<Unary>(base);
        unary != nullptr && unary->op == UnaryOp::Neg && is_constant_expr(unary->operand, MathConstant::I)) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            BigInt rem = n % BigInt(4);
            if (rem.is_negative()) rem += BigInt(4);
            const long long r = rem.to_u64();
            if (r == 0) return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(1)));
            if (r == 1) return traced_result(RuleId::SimplifyPowerOne, target_before, base);
            if (r == 2) return traced_result(RuleId::SimplifyPowerOne, target_before, make_integer(arena_, BigInt(-1)));
            if (r == 3) return traced_result(RuleId::SimplifyPowerOne, target_before, unary->operand);
        }
    }

    // sqrt(A)^n simplification
    if (const auto* sqrt_call = expr_cast<FuncCall>(base);
        sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            if (!n.is_zero()) {
                ExprPtr arg = sqrt_call->args.front();
                BigInt k = n / BigInt(2);
                BigInt rem = n % BigInt(2);

                if (rem.is_zero()) {
                    // (sqrt(A))^(2k) = A^k  (works for both positive and negative k)
                    return simplify_power(arg, make_integer(arena_, k));
                } else if (n != BigInt(-1) && !k.is_zero()) {
                    // (sqrt(A))^(2k+1) = A^k * sqrt(A)
                    // Under Euclidean division for negative n:
                    // e.g. -3 / 2 = -2, -3 % 2 = 1.
                    // (sqrt(A))^-3 = A^-2 * sqrt(A). This is correct and contains no negative sqrt.
                    // We skip n == -1 to avoid infinite loop.
                    auto ak = simplify_power(arg, make_integer(arena_, k));
                    if (ak.is_error()) return ak;
                    return simplify_product_factors({ak.value(), base}, ExprPtr{}, false);
                }
            }
        }
    }

    // (a · b · c · ...)^n  →  a^n · b^n · c^n · ...   for any integer n ≠ 0.
    // Distribution over Product is always valid for integer n: no branch-cut
    // issues, and the number of factors is preserved (no monomial blow-up as
    // with Pow of a Sum).  Negative exponents are essential for canonical
    // cancellation: without distribution, `(d_0·d_1·d_2)^-1` stays as
    // Pow(Product, -1) and downstream Product-flattening cannot combine it
    // with positive Pows on the same symbols (Bareiss-Edmonds inverse
    // extraction depends on this).  Reference: HC-F4-INV-SYMBOLIC-CANONICAL.
    //
    // Guards:
    //   * |n| ≤ 20 to prevent BigInt explosion on numeric factors.
    //   * For n < 0, require all factors known-nonzero: otherwise
    //     simplify_power on a zero factor returns "rational division by
    //     zero" (legitimate for 0^-1) and aborts the caller — limit/Gruntz
    //     relies on the unsimplified Pow(Product, -1) form to manipulate
    //     intermediates that vanish at the limit point.
    if (const auto* product = expr_cast<Product>(base); product != nullptr) {
        if (auto maybe_n = try_get_integer_exponent(exponent); maybe_n.has_value()) {
            const BigInt n = *maybe_n;
            if (!n.is_zero()) {
                BigInt abs_n = n.is_negative() ? -n : n;
                if (abs_n <= BigInt(20)) {
                    bool safe = !n.is_negative();
                    if (!safe) {
                        // Negative exponent: require every factor to be
                        // structurally non-zero (literal nonzero, Symbol,
                        // Constant) or assumed nonzero.  Symbolic factors
                        // that could vanish (e.g. Sum) stay un-distributed
                        // so limit/Gruntz can manipulate the Pow(Product, n)
                        // form without triggering 0^(-n).
                        safe = true;
                        for (ExprPtr factor : product->factors) {
                            bool nz = false;
                            if (const auto* il = expr_cast<IntegerLit>(factor)) nz = !il->value.is_zero();
                            else if (const auto* rl = expr_cast<RationalLit>(factor)) nz = !rl->numerator.is_zero();
                            else if (expr_is<Symbol>(factor) || expr_is<Constant>(factor)) nz = true;
                            else if (is_assumed_nonzero(factor)) nz = true;
                            if (!nz) { safe = false; break; }
                        }
                    }
                    if (safe) {
                        std::vector<ExprPtr> factors;
                        factors.reserve(product->factors.size());
                        for (ExprPtr factor : product->factors) {
                            auto p = simplify_power(factor, exponent);
                            if (p.is_error()) return p;
                            factors.push_back(p.value());
                        }
                        return simplify_product_factors(factors);
                    }
                    // Full distribution unsafe (a symbolic factor may vanish):
                    // keep Pow(Product, n) intact and fall through. NOTE: pulling
                    // out just the numeric coefficient here — (c·rest)^n →
                    // c^n·Pow(rest,n) — is mathematically exact but BREAKS Gruntz:
                    // its mrv machinery relies on the single un-split Pow node for
                    // vanishing intermediates (regression: AcidTest Gruntz limit →
                    // "ComplexRational: division by zero"). The numeric-factor
                    // normalization that asin/atan round-trips need belongs at the
                    // Sum/like-term layer, not in this universal power path. → T-054.
                }
            }
        }
    }

    if (const auto* outer = expr_cast<Binary>(base); outer != nullptr && outer->op == BinaryOp::Pow) {
        bool safe = false;
        auto c_int = try_get_integer_exponent(exponent);
        if (c_int.has_value()) {
            safe = true;
        } else if (is_known_positive(outer->left)) {
            safe = true;
        }

        if (safe) {
            auto new_exp = simplify_expr(arena_.make<Binary>(BinaryOp::Mul, outer->right, exponent));
            if (new_exp.is_ok()) {
                auto rewritten = simplify_power(outer->left, new_exp.value());
                if (rewritten.is_ok()) {
                    append_trace(RuleId::SimplifyFlattenNestedPowers, target_before, rewritten.value());
                    return rewritten;
                }
            }
        } else {
            // F8.0-6.2 / Task 20 BC-2 (Branch_Cut_Propagation.md §2 rule 3):
            //   (z^a)^b = z^(a·b) · e^(2πi · b · K(a · ln(z)))
            // When ctx.strict_branch_cuts() is active and the base is not
            // provably positive, the legacy `(a^b)^c → a^(b·c)` flatten loses
            // the unwinding correction. Emit the corrected product instead so
            // the identity stays algebraically exact in the complex plane.
            const bool strict = (context_ != nullptr) && context_->strict_branch_cuts();
            if (strict) {
                auto new_exp = simplify_expr(arena_.make<Binary>(BinaryOp::Mul, outer->right, exponent));
                if (new_exp.is_ok()) {
                    auto flattened = simplify_power(outer->left, new_exp.value());
                    if (flattened.is_ok()) {
                        ExprPtr correction = branch_cut::make_pow_of_pow_correction(
                            outer->left, outer->right, exponent, arena_);
                        ExprPtr corrected = arena_.make<Binary>(BinaryOp::Mul,
                            flattened.value(), correction);
                        append_trace(RuleId::SimplifyFlattenNestedPowers, target_before, corrected);
                        return ok(corrected);
                    }
                }
                // Strict mode: never silently drop the correction. If we can't
                // build the flattened form, preserve the original (z^a)^b.
                append_trace(RuleId::SimplifyFlattenNestedPowers, target_before, base);
                return ok(arena_.make<Binary>(BinaryOp::Pow, base, exponent));
            }
            // Legacy fallback for literal rationals (already covers some cases)
            LiteralRational l_exp, r_exp;
            auto l_exact = try_get_exact_rational(outer->right, l_exp);
            auto r_exact = try_get_exact_rational(exponent, r_exp);
            if (l_exact.is_ok() && l_exact.value() && r_exact.is_ok() && r_exact.value()) {
                // (a^b)^c where b, c are rational. Legacy behaviour — strict
                // mode handles the branch-cut-aware path above.
                auto rewritten = simplify_power(outer->left, make_rational(arena_, l_exp.value * r_exp.value));
                if (rewritten.is_ok()) {
                    append_trace(RuleId::SimplifyFlattenNestedPowers, target_before, rewritten.value());
                    return rewritten;
                }
            }
        }
    }

    // sqrt(A)^2 = A for any real A (principal branch: sqrt(-k) = i*sqrt(k), (i*sqrt(k))^2 = -k = A).
    // Only apply when A is a known real constant (nonneg or known-negative rational).
    if (const auto* sqrt_call = expr_cast<FuncCall>(base);
        sqrt_call != nullptr && sqrt_call->func_id == BuiltinOp::Sqrt && sqrt_call->args.size() == 1U) {
        LiteralRational exp_r;
        auto e_ex = try_get_exact_rational(exponent, exp_r);
        if (e_ex.is_ok() && e_ex.value() && exp_r.value.is_integer() && exp_r.value.numerator() == BigInt(2)) {
            // (sqrt(x))^2 = x is always true for principal complex branch.
            return traced_result(RuleId::SimplifyFlattenNestedPowers, target_before, sqrt_call->args.front());
        }
    }

    // GAP #6: exp(a)^b -> exp(a*b)
    ExprPtr exp_base = base;
    bool neg_sign = false;
    if (const auto* neg = expr_cast<Unary>(base); neg != nullptr && neg->op == UnaryOp::Neg) {
        exp_base = neg->operand;
        neg_sign = true;
    }

    if (const auto* exp_call = expr_cast<FuncCall>(exp_base);
        exp_call != nullptr && exp_call->func_id == BuiltinOp::Exp && exp_call->args.size() == 1U) {
        auto new_arg = simplify_expr(arena_.make<Binary>(BinaryOp::Mul, exp_call->args.front(), exponent));
        if (new_arg.is_ok()) {
            ExprPtr res = arena_.make<FuncCall>(BuiltinOp::Exp, std::vector<ExprPtr>{new_arg.value()});
            if (neg_sign) {
                if (auto exp_val = try_get_integer_exponent(exponent); exp_val.has_value()) {
                    if ((*exp_val % BigInt(2)) != BigInt(0)) res = arena_.make<Unary>(UnaryOp::Neg, res);
                    // if even, remains positive
                } else {
                    // complex case, skip for now
                    return ok(arena_.make<Binary>(BinaryOp::Pow, base, exponent));
                }
            }
            return simplify_expr(res);
        }
    }

    if (expr_is<Binary>(target_before)) {
        const auto& before = expr_ref<Binary>(target_before);
        if (before.op == BinaryOp::Pow && before.left == base && before.right == exponent) return ok(target_before);
    }
    return ok(arena_.make<Binary>(BinaryOp::Pow, base, exponent));
}

}  // namespace cas::symbolic::detail
