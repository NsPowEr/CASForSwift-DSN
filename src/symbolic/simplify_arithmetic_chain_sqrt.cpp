// Sqrt-pair collapses applied during Product simplification.
//
// Two phases run back-to-back:
//
//   1. sqrt(a) · sqrt(a) → a   (when `a` is provably non-negative).
//      `is_struct_nonneg` recognises sums/products of structurally
//      non-negative terms, even-power Pow, Sqrt / Abs / Exp / Cosh
//      results, provably-positive real constants (pi, e, gamma), and
//      non-negative literals; the assumption store can certify additional
//      symbols.
//
//   2. sqrt(a) · sqrt(b) → sqrt(a·b)  for a, b PROVABLY non-negative (A50:
//      generalized beyond bare literal rationals — see below).
//      The product is re-simplified so perfect squares collapse into
//      the numeric coefficient.
//
// Both passes iterate to a fixed point and mutate `symbolic` /
// `coefficient` in place.  The first pass is gated by HC-F4-QR-SYMBOLIC-
// TIMEOUT to avoid an exponential blow-up when later steps would
// otherwise rebuild `sqrt(a^2)` and attempt perfect-square factoring.

#include "simplify_arithmetic_chain_impl.hpp"

namespace cas::symbolic::detail {

Result<void> Simplifier::collapse_sqrt_pairs(
    std::vector<std::pair<ExprPtr, BigInt>>& symbolic,
    ComplexRational& coefficient)
{
    // Shared by both phases (A50): proves an expression is non-negative
    // compositionally, so Phase 2 below can reuse the exact same soundness
    // argument as Phase 1 instead of a narrower literal-only check.
    // `Constant` case: pi/e/gamma are the only MathConstant values that are
    // real AND positive (CLAUDE.md's own exception for exact math constants);
    // `i`/the infinities/NaN/Indeterminate are deliberately excluded — they
    // are not ordinary non-negative reals.
    auto is_struct_nonneg = [](ExprPtr e, auto&& self) -> bool {
        if (!e) return false;
        if (const auto* il = expr_cast<IntegerLit>(e)) return !il->value.is_negative();
        if (const auto* rl = expr_cast<RationalLit>(e)) return !rl->numerator.is_negative();
        if (const auto* c = expr_cast<Constant>(e)) {
            return c->value == MathConstant::Pi
                || c->value == MathConstant::E
                || c->value == MathConstant::EulerGamma;
        }
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

    // Phase 1: sqrt(a) * sqrt(a) → a  for structurally / assumption-derived
    // non-negative `a`.
    {
        bool changed = true;
        while (changed) {
            changed = false;
            std::vector<std::pair<ExprPtr, BigInt>> next_symbolic;
            for (const auto& entry : symbolic) {
                const auto* fa = expr_cast<FuncCall>(entry.first);
                if (fa && fa->func_id == BuiltinOp::Sqrt && fa->args.size() == 1U && known_nonneg(fa->args[0])) {
                    ExprPtr arg = fa->args[0];
                    BigInt k = entry.second;
                    BigInt q(0);
                    BigInt r(0);
                    bool should_decompose = false;

                    if ((k % BigInt(2)).is_zero()) {
                        q = k / BigInt(2);
                        r = BigInt(0);
                        should_decompose = true;
                    } else if (k >= BigInt(3)) {
                        q = k / BigInt(2);
                        r = BigInt(1);
                        should_decompose = true;
                    } else if (k <= BigInt(-3)) {
                        q = k / BigInt(2);
                        r = BigInt(-1);
                        should_decompose = true;
                    }

                    if (should_decompose) {
                        if (!q.is_zero()) {
                            auto res_power = simplify_power(arg, make_integer(arena_, q));
                            if (res_power.is_error()) return fail<void>(res_power.error());
                            LiteralRational rat;
                            bool is_rat = false;
                            auto rat_check = try_get_exact_rational(res_power.value(), rat);
                            if (rat_check.is_ok() && rat_check.value()) {
                                is_rat = true;
                            }
                            if (is_rat) {
                                coefficient = coefficient * ComplexRational(std::move(rat.value));
                            } else {
                                next_symbolic.push_back({arg, q});
                            }
                            changed = true;
                        }
                        if (!r.is_zero()) {
                            next_symbolic.push_back({entry.first, r});
                        }
                    } else {
                        next_symbolic.push_back(entry);
                    }
                } else {
                    next_symbolic.push_back(entry);
                }
            }
            if (changed) {
                symbolic = std::move(next_symbolic);
                merge_symbolic_factors(symbolic);
            }
        }
    }

    // Phase 2: sqrt(a) * sqrt(b) → sqrt(a*b) for a, b PROVABLY non-negative.
    //
    // A50 — generalizes the old rational-only rule via `extract_monomial`
    // (the same rational-coefficient / symbolic-factor decomposition Step 4
    // uses for like-term collection, `simplify_arithmetic.cpp`): `a`/`b` no
    // longer need to be bare literal rationals, just provably non-negative —
    // e.g. `sqrt(pi/2)*sqrt(2) -> sqrt(pi)`, needed to let `D(F)=f` close on
    // `integrate(exp(2*x^2))` (the general gaussian, A > 0, non-perfect-
    // square case). Soundness: the extracted rational coefficient must be
    // >= 0 AND every remaining symbolic factor (at its own extracted
    // exponent) must be provably non-negative via `known_nonneg` — reusing
    // the exact machinery Phase 1 already relies on, not a new ad-hoc rule.
    // Old behaviour (bare literal a, b) is the special case with an empty
    // symbolic remainder, so nothing regresses.
    {
        auto nonneg_monomial = [&](ExprPtr e) -> std::optional<MonomialTerm> {
            const auto* fc = expr_cast<FuncCall>(e);
            if (!fc || fc->func_id != BuiltinOp::Sqrt || fc->args.size() != 1U)
                return std::nullopt;
            auto term = extract_monomial(fc->args[0]);
            if (term.is_error() || !term.value().has_value()) return std::nullopt;
            const MonomialTerm& mt = *term.value();
            if (mt.coefficient.numerator().is_negative()) return std::nullopt;
            for (const auto& factor_exp : mt.key.factors) {
                if ((factor_exp.second % BigInt(2)).is_zero()) continue;
                if (!known_nonneg(factor_exp.first)) return std::nullopt;
            }
            return mt;
        };
        bool merged_any = true;
        while (merged_any) {
            merged_any = false;
            for (std::size_t i = 0; i < symbolic.size() && !merged_any; ++i) {
                if (symbolic[i].second != BigInt(1)) continue;
                auto ma = nonneg_monomial(symbolic[i].first);
                if (!ma) continue;
                for (std::size_t j = i + 1; j < symbolic.size(); ++j) {
                    if (symbolic[j].second != BigInt(1)) continue;
                    auto mb = nonneg_monomial(symbolic[j].first);
                    if (!mb) continue;

                    // A50 — la fusione deve produrre una VERA riduzione
                    // (coefficiente non banale, o un fattore condiviso le cui
                    // potenze si consolidano), non la semplice concatenazione
                    // di due simboli scorrelati sotto un unico radicale. Senza
                    // questa guardia, `sqrt(a)*sqrt(b)` per simboli a,b
                    // positivi-per-assunzione ma senza altra relazione veniva
                    // ricombinato in `sqrt(a*b)` — disfacendo l'espansione
                    // opposta e deliberata di `builtin_rewrite_algebraic.cpp`
                    // (canonicalizza `sqrt(prodotto non-negativo)` nella forma
                    // SEPARATA), che serve a valle (es. cancellazioni QR,
                    // `test_sqrt_fold.cpp`). Le due direzioni oscillerebbero.
                    const bool nontrivial_coefficient =
                        !(ma->coefficient == Rational(BigInt(1)))
                        || !(mb->coefficient == Rational(BigInt(1)));
                    bool shares_base = false;
                    for (const auto& fa : ma->key.factors) {
                        for (const auto& fb : mb->key.factors) {
                            if (structural_equal(fa.first, fb.first)) { shares_base = true; break; }
                        }
                        if (shares_base) break;
                    }
                    if (!nontrivial_coefficient && !shares_base) continue;

                    Rational new_coeff = ma->coefficient * mb->coefficient;
                    std::vector<std::pair<ExprPtr, BigInt>> merged_factors = ma->key.factors;
                    merged_factors.insert(merged_factors.end(),
                        mb->key.factors.begin(), mb->key.factors.end());
                    merge_symbolic_factors(merged_factors);

                    std::vector<ExprPtr> arg_factors;
                    if (!(new_coeff == Rational(BigInt(1)))) {
                        arg_factors.push_back((new_coeff.denominator() == BigInt(1))
                            ? static_cast<ExprPtr>(arena_.make<IntegerLit>(new_coeff.numerator()))
                            : static_cast<ExprPtr>(arena_.make<RationalLit>(
                                new_coeff.numerator(), new_coeff.denominator())));
                    }
                    for (const auto& factor_exp : merged_factors) {
                        arg_factors.push_back((factor_exp.second == BigInt(1))
                            ? factor_exp.first
                            : static_cast<ExprPtr>(arena_.make<Binary>(BinaryOp::Pow,
                                factor_exp.first, make_integer(arena_, factor_exp.second))));
                    }
                    ExprPtr new_arg = arg_factors.empty()
                        ? static_cast<ExprPtr>(arena_.make<IntegerLit>(new_coeff.numerator()))
                        : (arg_factors.size() == 1U
                            ? arg_factors.front()
                            : static_cast<ExprPtr>(arena_.make<Product>(arg_factors)));

                    auto new_sqrt = simplify_expr(
                        arena_.make<FuncCall>(BuiltinOp::Sqrt,
                            std::vector<ExprPtr>{new_arg}));
                    ExprPtr replacement = new_sqrt.is_ok() ? new_sqrt.value()
                        : arena_.make<FuncCall>(BuiltinOp::Sqrt,
                            std::vector<ExprPtr>{new_arg});
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
    return ok();
}

}  // namespace cas::symbolic::detail
