#include "simplify_arithmetic_chain_impl.hpp"
#include <algorithm>
#include <optional>

// simplify_node(Product) + simplify_product_factors implementation.
// Sum handling is in simplify_arithmetic_chain_sum.cpp.
// F1.4 like-term helpers are in simplify_arithmetic_chain_liketerm.cpp.

namespace cas::symbolic::detail {

Result<ExprPtr> Simplifier::simplify_node(ExprPtr original, const Product& node) {
    return simplify_product_factors(node.factors, original);
}

Result<ExprPtr> Simplifier::simplify_product_factors(
    const std::vector<ExprPtr>& factors,
    ExprPtr target_before,
    bool inputs_are_simplified)
{
    if (!target_before && trace_enabled_)
        target_before = make_product_target(factors);

    // Step 1: recursively simplify each factor; detect zero early.
    // EXCEPTION: Gamma(z+n) with integer n is kept UNsimplified so that Step 4
    // can pair it directly with another Gamma for reflection.  Without this guard
    // simplify_expr expands Gamma(z+n) → z*(z+1)*…*Gamma(z) which Step 8 then
    // distributes into a Sum, hiding the Gamma from the reflection detector.
    // Returns true when f is Gamma(z+n) or Gamma(z-n) with non-zero integer n.
    // Recognizes both Sum and Binary(Add/Sub) forms that arise from parsing.
    auto is_deferred_gamma = [](ExprPtr f) -> bool {
        const auto* fc = expr_cast<FuncCall>(f);
        if (!fc || fc->func_id != BuiltinOp::Gamma || fc->args.size() != 1U)
            return false;
        ExprPtr arg = fc->args[0];
        // Sum form: Sum([..., IntegerLit(n), ...])
        if (const auto* sum = expr_cast<Sum>(arg)) {
            for (ExprPtr t : sum->terms) {
                if (const auto* il = expr_cast<IntegerLit>(t))
                    if (!il->value.is_zero()) return true;
            }
            return false;
        }
        // Binary Add/Sub with one integer operand: z+n, n+z, z-n, n-z
        // Also handles -n forms like Unary(Neg, IntegerLit(k)).
        auto is_nonzero_int = [](ExprPtr e) -> bool {
            if (const auto* il = expr_cast<IntegerLit>(e))
                return !il->value.is_zero();
            if (const auto* u = expr_cast<Unary>(e))
                if (u->op == UnaryOp::Neg)
                    if (const auto* il = expr_cast<IntegerLit>(u->operand))
                        return !il->value.is_zero();
            return false;
        };
        if (const auto* bin = expr_cast<Binary>(arg)) {
            if (bin->op == BinaryOp::Add || bin->op == BinaryOp::Sub) {
                if (is_nonzero_int(bin->left) || is_nonzero_int(bin->right))
                    return true;
            }
        }
        // Plain integer literal Gamma(n) — not shifted, no need to defer
        return false;
    };

    std::vector<ExprPtr> initial_factors;
    bool has_zero = false;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        auto simplify_current = [&]() -> Result<ExprPtr> {
            if (!inputs_are_simplified && is_deferred_gamma(factors[i]))
                return ok(factors[i]);  // keep Gamma(z+n) unsimplified for Step 4
            return inputs_are_simplified ? ok(factors[i]) : simplify_expr(factors[i]);
        };
        Result<ExprPtr> simplified =
            fail<ExprPtr>(make_error(CASErrorKind::InternalError, "Unreachable"));
        if (trace_enabled_) {
            std::vector<ExprPtr> current_factors = factors;
            for (std::size_t j = 0;
                 j < initial_factors.size() && j < current_factors.size(); ++j)
                current_factors[j] = initial_factors[j];
            ScopedFrame frame(*this,
                [this, current_factors = std::move(current_factors), i](ExprPtr value) mutable {
                    current_factors[i] = value;
                    return make_product_target(current_factors);
                });
            simplified = simplify_current();
        } else {
            simplified = simplify_current();
        }
        if (simplified.is_error()) return simplified;
        if (is_zero_expr(simplified.value())) has_zero = true;
        initial_factors.push_back(simplified.value());
    }

    if (has_zero) {
        for (ExprPtr factor : initial_factors) {
            if (is_constant_expr(factor, MathConstant::Infinity))
                return fail<ExprPtr>(
                    make_error(CASErrorKind::Undefined, "0 * Infinity is undefined"));
            if (const auto* u = expr_cast<Unary>(factor);
                u != nullptr && u->op == UnaryOp::Neg) {
                if (is_constant_expr(u->operand, MathConstant::Infinity))
                    return fail<ExprPtr>(
                        make_error(CASErrorKind::Undefined, "0 * Infinity is undefined"));
            }
        }
        return traced_result(RuleId::SimplifyMultiplyByZero,
            target_before, make_integer(arena_, BigInt(0)));
    }

    // Step 2: flatten nested Products and Div nodes.
    std::vector<ExprPtr> flat_factors;
    std::function<Result<void>(ExprPtr, bool)> flatten =
        [&](ExprPtr f, bool invert) -> Result<void> {
        if (const auto* prod = expr_cast<Product>(f)) {
            for (ExprPtr factor : prod->factors) {
                auto res = flatten(factor, invert);
                if (res.is_error()) return res;
            }
        } else if (const auto* div = expr_cast<Binary>(f);
                   div != nullptr && div->op == BinaryOp::Div) {
            if (auto res_l = flatten(div->left, invert);   res_l.is_error()) return res_l;
            if (auto res_r = flatten(div->right, !invert); res_r.is_error()) return res_r;
        } else if (invert) {
            auto inv = simplify_power(f, make_integer(arena_, BigInt(-1)));
            if (inv.is_error()) return fail<void>(inv.error());
            flat_factors.push_back(inv.value());
        } else {
            flat_factors.push_back(f);
        }
        return ok();
    };
    for (ExprPtr f : initial_factors)
        if (auto res = flatten(f, false); res.is_error())
            return fail<ExprPtr>(res.error());

    // Step 3: separate numeric coefficient, imaginary unit accumulation,
    //         infinity tracking, and symbolic base/exponent pairs.
    ComplexRational coefficient = ComplexRational::one();
    std::vector<std::pair<ExprPtr, BigInt>> symbolic;
    bool has_infinity = false;
    int infinity_sign = 1;

    for (ExprPtr f : flat_factors) {
        if (auto timeout = check_timeout(); timeout.is_error())
            return fail<ExprPtr>(timeout.error());

        int sign = 1;
        ExprPtr core = f;
        if (const auto* u = expr_cast<Unary>(f);
            u != nullptr && u->op == UnaryOp::Neg) {
            sign = -1;
            core = u->operand;
        }
        if (is_constant_expr(core, MathConstant::Infinity)) {
            has_infinity = true;
            infinity_sign *= sign;
            continue;
        }

        if (sign == -1) {
            coefficient = coefficient * ComplexRational(Rational(BigInt(-1)));
        }

        LiteralComplex comp;
        auto exact = try_get_exact_complex(core, comp);
        if (exact.is_ok() && exact.value()) {
            coefficient = coefficient * comp.value;
            continue;
        }

        if (is_constant_expr(core, MathConstant::I)) {
            coefficient = coefficient * ComplexRational::imag_unit();
            continue;
        }

        if (const auto* binary = expr_cast<Binary>(core);
            binary != nullptr && binary->op == BinaryOp::Pow) {
            if (is_constant_expr(binary->left, MathConstant::I)) {
                if (auto exponent = try_get_integer_exponent(binary->right);
                    exponent.has_value()) {
                    // I^n = I % 4
                    BigInt rem = *exponent % BigInt(4);
                    if (rem.is_negative()) rem += BigInt(4);
                    long long r = rem.to_u64();
                    if (r == 1) coefficient = coefficient * ComplexRational::imag_unit();
                    else if (r == 2) coefficient = coefficient * ComplexRational(Rational(BigInt(-1)));
                    else if (r == 3) coefficient = coefficient * ComplexRational(Rational(BigInt(0)), Rational(BigInt(-1)));
                    continue;
                }
            }
            if (auto exponent = try_get_integer_exponent(binary->right);
                exponent.has_value()) {
                symbolic.push_back({binary->left, *exponent});
                continue;
            }
        }
        symbolic.push_back({core, BigInt(1)});
    }

    if (has_infinity) {
        if (coefficient.real().numerator().is_negative() || coefficient.imag().numerator().is_negative()) 
            infinity_sign *= -1; // Simplified infinity sign check
        ExprPtr inf = arena_.make<Constant>(MathConstant::Infinity);
        if (infinity_sign < 0)
            return ok(arena_.make<Unary>(UnaryOp::Neg, inf));
        return ok(inf);
    }

    if (coefficient.is_zero())
        return traced_result(RuleId::SimplifyMultiplyByZero,
            target_before, make_integer(arena_, BigInt(0)));
    merge_symbolic_factors(symbolic);

    // Step 4: L3-04 Gamma reflection identity.
    //   Γ(z)·Γ(1-z) = π/sin(πz)    [m=1 case]
    //   Γ(z)·Γ(-z)  = -π/(z·sin(πz)) [m=0 case]
    // Fires only when two Gamma factors sum to an integer m ∈ {0, 1}.
    {
        auto gamma_arg = [](const std::pair<ExprPtr, BigInt>& p) -> ExprPtr {
            if (p.second != BigInt(1)) return nullptr;
            const auto* fc = expr_cast<FuncCall>(p.first);
            if (!fc || fc->func_id != BuiltinOp::Gamma || fc->args.size() != 1U)
                return nullptr;
            return fc->args[0];
        };
        bool reflected_any = false;
        bool keep_scanning = true;
        while (keep_scanning) {
            keep_scanning = false;
            std::vector<std::size_t> gamma_idx;
            for (std::size_t i = 0; i < symbolic.size(); ++i)
                if (gamma_arg(symbolic[i]) != nullptr) gamma_idx.push_back(i);
            for (std::size_t a = 0; a < gamma_idx.size() && !keep_scanning; ++a) {
                for (std::size_t b = a + 1; b < gamma_idx.size() && !keep_scanning; ++b) {
                    ExprPtr za = gamma_arg(symbolic[gamma_idx[a]]);
                    ExprPtr zb = gamma_arg(symbolic[gamma_idx[b]]);
                    if (!za || !zb) continue;
                    auto sum_simp = simplify_expr(
                        arena_.make<Sum>(std::vector<ExprPtr>{za, zb}));
                    if (sum_simp.is_error()) continue;
                    const auto* m_lit = expr_cast<IntegerLit>(sum_simp.value());
                    if (!m_lit) continue;
                    const bool m_is_one  = (m_lit->value == BigInt(1));
                    const bool m_is_zero = m_lit->value.is_zero();
                    if (!m_is_one && !m_is_zero) continue;
                    ExprPtr pi_const = arena_.make<Constant>(MathConstant::Pi);
                    auto pi_z_simp = simplify_expr(
                        arena_.make<Product>(std::vector<ExprPtr>{pi_const, za}));
                    if (pi_z_simp.is_error()) continue;
                    auto sin_simp = simplify_expr(
                        arena_.make<FuncCall>(BuiltinOp::Sin,
                            std::vector<ExprPtr>{pi_z_simp.value()}));
                    if (sin_simp.is_error()) continue;
                    std::size_t ia = gamma_idx[a], ib = gamma_idx[b];
                    if (ia > ib) std::swap(ia, ib);
                    symbolic.erase(symbolic.begin() + ib);
                    symbolic.erase(symbolic.begin() + ia);
                    symbolic.push_back({pi_const,          BigInt(1)});
                    symbolic.push_back({sin_simp.value(),  BigInt(-1)});
                    if (m_is_zero) {
                        coefficient = coefficient * ComplexRational(Rational(BigInt(-1)));
                        symbolic.push_back({za, BigInt(-1)});
                    }
                    reflected_any = true;
                    keep_scanning = true;
                }
            }
        }
        if (reflected_any) merge_symbolic_factors(symbolic);
    }

    // Step 4.5: Rational cancellation of negated-Sum pairs.
    // Detects S^ka · (-S)^kb → (-1)^kb · S^(ka+kb), removing both when ka+kb=0.
    // This handles (x+1)·(-x-1)^{-1} = -1 without Step 8 expanding first.
    // General: any two symbolic factors whose bases satisfy simplify(Sa+Sb)==0.
    {
        bool cancelled_any = true;
        while (cancelled_any) {
            cancelled_any = false;
            const std::size_t sz = symbolic.size();
            for (std::size_t ii = 0; ii < sz && !cancelled_any; ++ii) {
                ExprPtr sa = symbolic[ii].first;
                BigInt   ka = symbolic[ii].second;
                for (std::size_t jj = ii + 1; jj < sz; ++jj) {
                    ExprPtr sb = symbolic[jj].first;
                    BigInt   kb = symbolic[jj].second;
                    // Quick pre-filter: skip trivially-unrelated types.
                    // Both must be non-Symbol (Symbols handled by monomial collector).
                    if (expr_is<Symbol>(sa) && expr_is<Symbol>(sb)) continue;
                    // Test Sa + Sb == 0.
                    auto probe = simplify_expr(
                        arena_.make<Sum>(std::vector<ExprPtr>{sa, sb}));
                    if (probe.is_error()) continue;
                    if (!is_zero_expr(probe.value())) continue;
                    // Sb = -Sa.  Absorption: Sa^ka · (-Sa)^kb = (-1)^kb · Sa^(ka+kb).
                    // (-1)^kb flips sign iff |kb| is odd.
                    BigInt abs_kb = kb.is_negative() ? BigInt(0) - kb : kb;
                    BigInt two(2);
                    if (!((abs_kb % two).is_zero()))
                        coefficient = coefficient * ComplexRational(Rational(BigInt(-1)));
                    BigInt new_exp = ka + kb;
                    // Remove jj first (higher index).
                    symbolic.erase(symbolic.begin() + static_cast<std::ptrdiff_t>(jj));
                    if (new_exp.is_zero()) {
                        symbolic.erase(symbolic.begin() + static_cast<std::ptrdiff_t>(ii));
                    } else {
                        symbolic[ii].second = new_exp;
                    }
                    cancelled_any = true;
                    break;
                }
            }
        }
    }

    // Step 5: L2-07 double-angle compaction: sin(x)·cos(x) → (1/2)·sin(2x).
    // After merge_symbolic_factors, canonical order places Cos before Sin (Cos
    // priority 77 < Sin priority 78 in term_order.cpp), so the pair may appear
    // as (cos,sin) rather than (sin,cos).  We search all unordered pairs.
    {
        std::optional<std::pair<std::size_t, std::size_t>> sc_pair;
        ExprPtr sc_arg = nullptr;
        for (std::size_t ii = 0; ii < symbolic.size() && !sc_pair; ++ii) {
            if (symbolic[ii].second != BigInt(1)) continue;
            const auto* fi = expr_cast<FuncCall>(symbolic[ii].first);
            if (!fi || fi->args.size() != 1U) continue;
            const bool fi_sin = (fi->func_id == BuiltinOp::Sin);
            const bool fi_cos = (fi->func_id == BuiltinOp::Cos);
            if (!fi_sin && !fi_cos) continue;
            for (std::size_t jj = ii + 1; jj < symbolic.size(); ++jj) {
                if (symbolic[jj].second != BigInt(1)) continue;
                const auto* fj = expr_cast<FuncCall>(symbolic[jj].first);
                if (!fj || fj->args.size() != 1U) continue;
                const bool fj_sin = (fj->func_id == BuiltinOp::Sin);
                const bool fj_cos = (fj->func_id == BuiltinOp::Cos);
                if (!fj_sin && !fj_cos) continue;
                // Need one Sin and one Cos with the same argument.
                if (!((fi_sin && fj_cos) || (fi_cos && fj_sin))) continue;
                if (structural_equal(fi->args[0], fj->args[0])) {
                    sc_pair = {ii, jj};
                    sc_arg  = fi->args[0];
                    break;
                }
            }
        }
        if (sc_pair) {
            const auto [si, ci] = *sc_pair;
            coefficient = coefficient * ComplexRational(Rational(BigInt(1), BigInt(2)));
            ExprPtr two_arg = arena_.make<Binary>(BinaryOp::Mul,
                make_integer(arena_, BigInt(2)), sc_arg);
            ExprPtr sin2x = arena_.make<FuncCall>(BuiltinOp::Sin,
                std::vector<ExprPtr>{two_arg});
            symbolic.erase(symbolic.begin() + ci);
            symbolic.erase(symbolic.begin() + si);
            auto sin2x_s = simplify_expr(sin2x);
            symbolic.push_back({sin2x_s.is_ok() ? sin2x_s.value() : sin2x, BigInt(1)});
        }
    }

    // Step 6: L3-08 Quantity multiplication — combine Quantity^1 factors.
    {
        std::vector<std::size_t> qty_idx;
        for (std::size_t i = 0; i < symbolic.size(); ++i)
            if (symbolic[i].second == BigInt(1)
                && expr_is<Quantity>(symbolic[i].first))
                qty_idx.push_back(i);
        if (qty_idx.size() >= 2) {
            SIDimensions combined_dim{};
            std::vector<ExprPtr> values;
            for (std::size_t k : qty_idx) {
                const auto& q = expr_ref<Quantity>(symbolic[k].first);
                values.push_back(q.value);
                combined_dim.m   += q.dimensions.m;
                combined_dim.kg  += q.dimensions.kg;
                combined_dim.s   += q.dimensions.s;
                combined_dim.A   += q.dimensions.A;
                combined_dim.K   += q.dimensions.K;
                combined_dim.mol += q.dimensions.mol;
                combined_dim.cd  += q.dimensions.cd;
            }
            for (auto it = qty_idx.rbegin(); it != qty_idx.rend(); ++it)
                symbolic.erase(symbolic.begin() + *it);
            ExprPtr combined_value;
            if (values.size() == 1U) {
                combined_value = values[0];
            } else {
                auto val_prod = simplify_product_factors(values, ExprPtr{}, false);
                combined_value = val_prod.is_ok()
                    ? val_prod.value()
                    : arena_.make<Product>(std::move(values));
            }
            symbolic.push_back({arena_.make<Quantity>(combined_value, combined_dim),
                                 BigInt(1)});
        }
    }

    // Step 6.5: HC-F4-QR-SYMBOLIC-TIMEOUT fix — sqrt(a)·sqrt(a) → a quando a è
    // strutturalmente non-negativo (somma/prodotto di quadrati, Pow esponente
    // pari, FuncCall sempre non-negativo, costante non-negativa) oppure
    // dichiarato tale via assumptions. Evita di delegare a Step 7 che produce
    // sqrt(a²) e poi tenta perfect-square detection via factorize — costo
    // esponenziale su rationali grandi (8×8 QR random).
    {
        auto is_struct_nonneg = [](ExprPtr e, auto&& self) -> bool {
            if (!e) return false;
            if (const auto* il = expr_cast<IntegerLit>(e)) return !il->value.is_negative();
            if (const auto* rl = expr_cast<RationalLit>(e)) return !rl->numerator.is_negative();
            if (const auto* b = expr_cast<Binary>(e); b && b->op == BinaryOp::Pow) {
                if (const auto* il = expr_cast<IntegerLit>(b->right)) {
                    if (!il->value.is_negative() && (il->value % BigInt(2)).is_zero()) return true;
                }
            }
            if (const auto* sum = expr_cast<Sum>(e)) {
                for (auto t : sum->terms) if (!self(t, self)) return false;
                return true;
            }
            if (const auto* prod = expr_cast<Product>(e)) {
                for (auto f : prod->factors) if (!self(f, self)) return false;
                return true;
            }
            if (const auto* fc = expr_cast<FuncCall>(e)) {
                if (fc->func_id == BuiltinOp::Sqrt
                    || fc->func_id == BuiltinOp::Abs
                    || fc->func_id == BuiltinOp::Exp
                    || fc->func_id == BuiltinOp::Cosh) return true;
            }
            return false;
        };
        auto known_nonneg = [&](ExprPtr e) -> bool {
            if (is_struct_nonneg(e, is_struct_nonneg)) return true;
            if (assumptions_ && assumptions_->is_nonnegative(e)) return true;
            return false;
        };
        bool merged_any = true;
        while (merged_any) {
            merged_any = false;
            for (std::size_t i = 0; i < symbolic.size() && !merged_any; ++i) {
                if (symbolic[i].second != BigInt(1)) continue;
                const auto* fa = expr_cast<FuncCall>(symbolic[i].first);
                if (!fa || fa->func_id != BuiltinOp::Sqrt || fa->args.size() != 1U) continue;
                for (std::size_t j = i + 1; j < symbolic.size(); ++j) {
                    if (symbolic[j].second != BigInt(1)) continue;
                    const auto* fb = expr_cast<FuncCall>(symbolic[j].first);
                    if (!fb || fb->func_id != BuiltinOp::Sqrt || fb->args.size() != 1U) continue;
                    if (fa->args[0] != fb->args[0]) continue;
                    if (!known_nonneg(fa->args[0])) continue;
                    ExprPtr arg = fa->args[0];
                    symbolic.erase(symbolic.begin() + j);
                    symbolic.erase(symbolic.begin() + i);
                    symbolic.push_back({arg, BigInt(1)});
                    merged_any = true;
                    break;
                }
            }
        }
    }

    // Step 7: L1-12 sqrt(a)·sqrt(b) → sqrt(a·b) for non-negative rational args.
    {
        auto get_rat_pos = [](ExprPtr e) -> std::optional<Rational> {
            const auto* fc = expr_cast<FuncCall>(e);
            if (!fc || fc->func_id != BuiltinOp::Sqrt || fc->args.size() != 1U)
                return std::nullopt;
            if (const auto* il = expr_cast<IntegerLit>(fc->args[0])) {
                if (il->value.is_negative()) return std::nullopt;
                return Rational(il->value, BigInt(1));
            }
            if (const auto* rl = expr_cast<RationalLit>(fc->args[0])) {
                if (rl->numerator.is_negative()) return std::nullopt;
                return Rational(rl->numerator, rl->denominator);
            }
            return std::nullopt;
        };
        bool merged_any = true;
        while (merged_any) {
            merged_any = false;
            for (std::size_t i = 0; i < symbolic.size() && !merged_any; ++i) {
                if (symbolic[i].second != BigInt(1)) continue;
                auto ra = get_rat_pos(symbolic[i].first);
                if (!ra) continue;
                for (std::size_t j = i + 1; j < symbolic.size(); ++j) {
                    if (symbolic[j].second != BigInt(1)) continue;
                    auto rb = get_rat_pos(symbolic[j].first);
                    if (!rb) continue;
                    Rational prod = (*ra) * (*rb);
                    ExprPtr arg = (prod.denominator() == BigInt(1))
                        ? static_cast<ExprPtr>(arena_.make<IntegerLit>(prod.numerator()))
                        : static_cast<ExprPtr>(arena_.make<RationalLit>(
                            prod.numerator(), prod.denominator()));
                    auto new_sqrt = simplify_expr(
                        arena_.make<FuncCall>(BuiltinOp::Sqrt,
                            std::vector<ExprPtr>{arg}));
                    ExprPtr replacement = new_sqrt.is_ok() ? new_sqrt.value()
                        : arena_.make<FuncCall>(BuiltinOp::Sqrt,
                            std::vector<ExprPtr>{arg});
                    symbolic.erase(symbolic.begin() + j);
                    symbolic.erase(symbolic.begin() + i);
                    LiteralRational rep_rat;
                    auto rep_check = try_get_exact_rational(replacement, rep_rat);
                    if (rep_check.is_ok() && rep_check.value()) {
                        coefficient = coefficient * ComplexRational(std::move(rep_rat.value));
                    } else {
                        symbolic.push_back({replacement, BigInt(1)});
                    }
                    merged_any = true;
                    break;
                }
            }
        }
    }

    // Step 8: distribute coefficient over a Sum factor (if present).
    {
        ExprPtr sum_factor = nullptr;
        std::size_t sum_idx = 0;
        bool found_sum = false;
        for (std::size_t i = 0; i < symbolic.size(); ++i) {
            if (symbolic[i].second == BigInt(1)
                && expr_is<Sum>(symbolic[i].first)) {
                sum_factor = symbolic[i].first;
                sum_idx    = i;
                found_sum  = true;
                break;
            }
        }
        // Second scan needed when the first loop exited early on coefficient check.
        if (!found_sum
            && !(coefficient == ComplexRational::one())
            && std::any_of(symbolic.begin(), symbolic.end(),
                [](const auto& p) {
                    return expr_is<Sum>(p.first) && p.second == BigInt(1);
                }))
        {
            for (std::size_t i = 0; i < symbolic.size(); ++i) {
                if (symbolic[i].second == BigInt(1)
                    && expr_is<Sum>(symbolic[i].first)) {
                    sum_factor = symbolic[i].first;
                    sum_idx    = i;
                    found_sum  = true;
                    break;
                }
            }
        }

        if (found_sum) {
            const auto* sum = expr_cast<Sum>(sum_factor);
            std::vector<std::pair<ExprPtr, BigInt>> other_sym = symbolic;
            other_sym.erase(other_sym.begin() + sum_idx);

            std::vector<ExprPtr> distributed_terms;
            for (ExprPtr term : sum->terms) {
                std::vector<ExprPtr> factors_for_term;
                if (!(coefficient == ComplexRational::one()))
                    factors_for_term.push_back(make_complex(arena_, coefficient));
                factors_for_term.push_back(term);
                for (const auto& [base, exp] : other_sym) {
                    factors_for_term.push_back(exp == BigInt(1)
                        ? base
                        : arena_.make<Binary>(BinaryOp::Pow,
                            base, make_integer(arena_, exp)));
                }
                auto prod = simplify_product_factors(factors_for_term, ExprPtr{}, false);
                if (prod.is_error()) return prod;
                distributed_terms.push_back(prod.value());
            }
            return simplify_sum_terms(distributed_terms, target_before, true);
        }
    }

    // Step 9: assemble final normalized result.
    std::vector<ExprPtr> normalized;
    bool is_neg = (coefficient == ComplexRational(Rational(BigInt(-1))));
    if (!is_neg && (!(coefficient == ComplexRational::one()) || symbolic.empty()))
        normalized.push_back(make_complex(arena_, coefficient));
    for (const auto& [base, exp] : symbolic) {
        if (exp.is_zero()) continue;
        normalized.push_back(exp == BigInt(1)
            ? base
            : arena_.make<Binary>(BinaryOp::Pow, base, make_integer(arena_, exp)));
    }
    normalized.erase(
        std::remove_if(normalized.begin(), normalized.end(),
            [](ExprPtr e) { return is_one_expr(e); }),
        normalized.end());

    if (normalized.empty())
        return traced_result(RuleId::SimplifyMultiplyByOne,
            target_before, make_complex(arena_, coefficient));

    ExprPtr result = normalized.size() == 1U
        ? normalized.front()
        : arena_.make<Product>(std::move(normalized));

    if (is_neg)
        return traced_result(RuleId::SimplifyMultiplyByOne,
            target_before, arena_.make<Unary>(UnaryOp::Neg, result));
    return ok(result);
}

} // namespace cas::symbolic::detail
