// F7.5.A4 — Hyperbolic + circular reciprocal/quotient normalization,
//           plus log/ln synonym unification.
//
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Sech_Csch_Identity.md.
//
// Bottom-up rewrite that replaces sech/csch/coth/tanh function calls
// (parsed as FuncCall("sech", …) with BuiltinOp::Custom, since these
// operators are not enumerated in BuiltinOp — extending the enum would
// touch 76 switch statements under -Wswitch -Werror, deferred to
// Fase 8) with their canonical form in terms of cosh/sinh:
//
//   sech(u) → 1/cosh(u)
//   csch(u) → 1/sinh(u)
//   coth(u) → cosh(u)/sinh(u)
//   tanh(u) → sinh(u)/cosh(u)
//
// Applied to both sides of `mathematically_equal` before structural and
// algebraic comparison. Structural sharing preserved when no rewrite
// fires in a subtree.

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "algebra_internal.hpp"

#include <string>
#include <utility>
#include <vector>

namespace cas::algebra {

ExprPtr hyperbolic_normalize(ExprPtr expr, AstArena& arena) {
    if (!expr) return expr;
    switch (expr->kind) {
        case ExprKind::FuncCall: {
            const auto* fc = expr_cast<FuncCall>(expr);
            if (fc == nullptr) return expr;
            if (fc->args.size() != 1U) {
                std::vector<ExprPtr> new_args;
                new_args.reserve(fc->args.size());
                bool changed = false;
                for (auto& a : fc->args) {
                    auto na = hyperbolic_normalize(a, arena);
                    if (na.get() != a.get()) changed = true;
                    new_args.push_back(na);
                }
                if (!changed) return expr;
                return arena.make<FuncCall>(fc->name, std::move(new_args));
            }
            ExprPtr u = hyperbolic_normalize(fc->args[0], arena);
            const std::string& nm = fc->name;
            auto make_call = [&](const char* name) {
                std::vector<ExprPtr> a; a.push_back(u);
                return arena.make<FuncCall>(std::string(name), std::move(a));
            };
            auto one = arena.make<IntegerLit>(BigInt(1));
            if (nm == "sech") {
                return arena.make<Binary>(BinaryOp::Div, one, make_call("cosh"));
            }
            if (nm == "csch") {
                return arena.make<Binary>(BinaryOp::Div, one, make_call("sinh"));
            }
            if (nm == "coth") {
                return arena.make<Binary>(BinaryOp::Div,
                                         make_call("cosh"), make_call("sinh"));
            }
            if (nm == "tanh") {
                return arena.make<Binary>(BinaryOp::Div,
                                         make_call("sinh"), make_call("cosh"));
            }
            // Circular reciprocals (Maxima-Giac notation parity).
            if (nm == "sec") {
                return arena.make<Binary>(BinaryOp::Div, one, make_call("cos"));
            }
            if (nm == "csc") {
                return arena.make<Binary>(BinaryOp::Div, one, make_call("sin"));
            }
            if (nm == "cot") {
                return arena.make<Binary>(BinaryOp::Div,
                                         make_call("cos"), make_call("sin"));
            }
            // log/ln synonym: parser maps `log` → BuiltinOp::Log and `ln`
            // → BuiltinOp::Ln; rest of the engine treats them as the same
            // natural log, but FuncCall::name still differs. Unify on "ln".
            if (nm == "log") {
                return make_call("ln");
            }
            if (u.get() == fc->args[0].get()) return expr;
            std::vector<ExprPtr> new_args; new_args.push_back(u);
            return arena.make<FuncCall>(fc->name, std::move(new_args));
        }
        case ExprKind::Unary: {
            const auto* un = expr_cast<Unary>(expr);
            if (un == nullptr) return expr;
            auto inner = hyperbolic_normalize(un->operand, arena);
            if (inner.get() == un->operand.get()) return expr;
            return arena.make<Unary>(un->op, inner);
        }
        case ExprKind::Binary: {
            const auto* bi = expr_cast<Binary>(expr);
            if (bi == nullptr) return expr;
            auto l = hyperbolic_normalize(bi->left, arena);
            auto r = hyperbolic_normalize(bi->right, arena);
            if (l.get() == bi->left.get() && r.get() == bi->right.get()) return expr;
            return arena.make<Binary>(bi->op, l, r);
        }
        case ExprKind::Sum: {
            const auto* su = expr_cast<Sum>(expr);
            if (su == nullptr) return expr;
            std::vector<ExprPtr> nt;
            nt.reserve(su->terms.size());
            bool changed = false;
            for (auto& t : su->terms) {
                auto nx = hyperbolic_normalize(t, arena);
                if (nx.get() != t.get()) changed = true;
                nt.push_back(nx);
            }
            if (!changed) return expr;
            return arena.make<Sum>(std::move(nt));
        }
        case ExprKind::Product: {
            const auto* pr = expr_cast<Product>(expr);
            if (pr == nullptr) return expr;
            std::vector<ExprPtr> nf;
            nf.reserve(pr->factors.size());
            bool changed = false;
            for (auto& f : pr->factors) {
                auto nx = hyperbolic_normalize(f, arena);
                if (nx.get() != f.get()) changed = true;
                nf.push_back(nx);
            }
            if (!changed) return expr;
            return arena.make<Product>(std::move(nf));
        }
        default:
            return expr;
    }
}

}  // namespace cas::algebra
