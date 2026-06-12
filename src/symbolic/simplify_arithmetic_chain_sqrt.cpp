// Sqrt-pair collapses applied during Product simplification.
//
// Two phases run back-to-back:
//
//   1. sqrt(a) · sqrt(a) → a   (when `a` is provably non-negative).
//      `is_struct_nonneg` recognises sums/products of structurally
//      non-negative terms, even-power Pow, Sqrt / Abs / Exp / Cosh
//      results, and non-negative literals; the assumption store can
//      certify additional symbols.
//
//   2. sqrt(a) · sqrt(b) → sqrt(a·b)  for non-negative rational a, b.
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
    // Phase 1: sqrt(a) * sqrt(a) → a  for structurally / assumption-derived
    // non-negative `a`.
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

    // Phase 2: sqrt(a) * sqrt(b) → sqrt(a*b) for non-negative rational a, b.
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
    return ok();
}

}  // namespace cas::symbolic::detail
