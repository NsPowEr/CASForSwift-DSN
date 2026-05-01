#include <cas/ast.hpp>
#include <string>
#include <vector>

namespace cas::symbolic {

static int get_precedence(ExprKind kind, ExprPtr e = ExprPtr{}) {
    switch (kind) {
        case ExprKind::Binary: {
            if (const auto* b = expr_cast<Binary>(e)) {
                if (b->op == BinaryOp::Pow) return 60;
                if (b->op == BinaryOp::Mul) return 40;
                if (b->op == BinaryOp::Div) return 40;
                if (b->op == BinaryOp::Add) return 30;
                if (b->op == BinaryOp::Sub) return 30;
            }
            return 50;
        }
        case ExprKind::Product: return 40;
        case ExprKind::Sum: return 30;
        case ExprKind::FuncCall: {
            if (const auto* f = expr_cast<FuncCall>(e)) {
                if (f->name == "exp") return 80;
                if (f->name == "ln" || f->name == "log") return 75;
                if (f->name == "sin" || f->name == "cos") return 70;
            }
            return 25;
        }
        case ExprKind::Symbol: return 10;
        case ExprKind::IntegerLit: return 5;
        case ExprKind::RationalLit: return 5;
        default: return 20;
    }
}

static std::vector<ExprPtr> get_children(ExprPtr e) {
    if (!e) return {};
    if (const auto* s = expr_cast<Sum>(e)) return s->terms;
    if (const auto* p = expr_cast<Product>(e)) return p->factors;
    if (const auto* b = expr_cast<Binary>(e)) return {b->left, b->right};
    if (const auto* u = expr_cast<Unary>(e)) return {u->operand};
    if (const auto* f = expr_cast<FuncCall>(e)) return f->args;
    return {};
}

bool is_lpo_less(ExprPtr s, ExprPtr t);

bool is_lpo_less(ExprPtr s, ExprPtr t) {
    if (s == t) return false;

    auto s_kind = expr_kind(s);
    auto t_kind = expr_kind(t);
    auto s_args = get_children(s);
    auto t_args = get_children(t);

    // 1. Esiste ti tale che s <= ti
    for (auto ti : t_args) {
        if (s == ti || is_lpo_less(s, ti)) return true;
    }

    // 2. f < g e per ogni si, si < t
    if (get_precedence(s_kind, s) < get_precedence(t_kind, t)) {
        for (auto si : s_args) {
            if (!is_lpo_less(si, t)) return false;
        }
        return true;
    }

    // 3. f = g e {s1...sn} <lex {t1...tm} e per ogni si, si < t
    if (s_kind == t_kind) {
        if (s_args.size() != t_args.size()) return s_args.size() < t_args.size();
        
        bool lex_less = false;
        for (size_t i = 0; i < s_args.size(); ++i) {
            if (is_lpo_less(s_args[i], t_args[i])) {
                lex_less = true;
                break;
            }
            if (is_lpo_less(t_args[i], s_args[i])) return false;
        }
        
        if (lex_less) {
            for (auto si : s_args) {
                if (!is_lpo_less(si, t)) return false;
            }
            return true;
        }
    }

    return false;
}

} // namespace cas::symbolic
