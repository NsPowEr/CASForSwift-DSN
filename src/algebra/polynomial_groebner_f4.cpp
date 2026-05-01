#include "polynomial_groebner_f4.hpp"
#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "cas/ast_debug.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <vector>
#include <numeric>
#include <iostream>

namespace cas::algebra {

// --- 1. Comparatori Monomiali ---

struct MonomialLexComparator {
    bool operator()(const Monomial& a, const Monomial& b) const {
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] > b[i]) return true;
            if (a[i] < b[i]) return false;
        }
        return false;
    }
};

struct MonomialGRevLexComparator {
    bool operator()(const Monomial& a, const Monomial& b) const {
        unsigned int da = 0; for (unsigned int x : a) da += x;
        unsigned int db = 0; for (unsigned int x : b) db += x;
        if (da != db) return da > db;
        for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
            if (a[i] != b[i]) return a[i] < b[i];
        }
        return false;
    }
};

using MonomialComparator = std::function<bool(const Monomial&, const Monomial&)>;

static MonomialComparator get_comparator(MonomialOrder order) {
    if (order == MonomialOrder::Lex) return MonomialLexComparator();
    return MonomialGRevLexComparator();
}

Monomial PolyF4::leading_monomial(MonomialOrder order) const {
    if (terms.empty()) return {};
    auto comp = get_comparator(order);
    auto it = terms.begin();
    auto best_it = it;
    for (++it; it != terms.end(); ++it) {
        if (comp(it->first, best_it->first)) best_it = it;
    }
    return best_it->first;
}

Rational PolyF4::leading_coefficient(MonomialOrder order) const {
    if (terms.empty()) return Rational(0);
    auto comp = get_comparator(order);
    auto it = terms.begin();
    auto best_it = it;
    for (++it; it != terms.end(); ++it) {
        if (comp(it->first, best_it->first)) best_it = it;
    }
    return best_it->second;
}

void PolyF4::make_monic(MonomialOrder order) {
    if (terms.empty()) return;
    Rational lc = leading_coefficient(order);
    if (lc.numerator().is_zero()) return;
    for (auto& [mon, coeff] : terms) coeff = coeff / lc;
}

// --- 2. Matrice di Macaulay ---

class MacaulayMatrix {
public:
    explicit MacaulayMatrix(MonomialComparator comp) : comp_(comp), monomial_to_col(comp) {}

    std::vector<Monomial> col_monomials;
    std::map<Monomial, size_t, MonomialComparator> monomial_to_col;
    std::vector<PolyF4> poly_rows; // Store polynomials directly, build matrix at the end

    void add_polynomial(const PolyF4& poly) {
        poly_rows.push_back(poly);
        for (const auto& [mon, coeff] : poly.terms) {
            register_monomial(mon);
        }
    }

    void register_monomial(const Monomial& mon) {
        if (monomial_to_col.find(mon) == monomial_to_col.end()) {
            col_monomials.push_back(mon);
            std::sort(col_monomials.begin(), col_monomials.end(), comp_);
            monomial_to_col.clear();
            for (size_t i = 0; i < col_monomials.size(); ++i) {
                monomial_to_col[col_monomials[i]] = i;
            }
        }
    }

    void gaussian_elimination() {
        if (poly_rows.empty()) return;
        
        // Build rows now that col_monomials is stable
        std::vector<std::map<size_t, Rational>> rows;
        for (const auto& poly : poly_rows) {
            std::map<size_t, Rational> row;
            for (const auto& [mon, coeff] : poly.terms) {
                row[monomial_to_col[mon]] = coeff;
            }
            rows.push_back(std::move(row));
        }

        size_t num_cols = col_monomials.size();
        size_t pivot_row = 0;
        for (size_t j = 0; j < num_cols && pivot_row < rows.size(); ++j) {
            size_t sel = rows.size();
            for (size_t i = pivot_row; i < rows.size(); ++i) {
                if (rows[i].count(j) && !rows[i][j].numerator().is_zero()) { sel = i; break; }
            }
            if (sel == rows.size()) continue;
            std::swap(rows[pivot_row], rows[sel]);
            
            Rational pivot_val = rows[pivot_row][j];
            for (auto& [col, val] : rows[pivot_row]) val = val / pivot_val;

            for (size_t i = 0; i < rows.size(); ++i) {
                if (i != pivot_row && rows[i].count(j) && !rows[i][j].numerator().is_zero()) {
                    Rational factor = rows[i][j];
                    for (const auto& [c, v] : rows[pivot_row]) {
                        rows[i][c] = rows[i][c] - factor * v;
                        if (rows[i][c].numerator().is_zero()) rows[i].erase(c);
                    }
                }
            }
            pivot_row++;
        }

        // Store back reduced rows
        poly_rows.clear();
        for (const auto& row : rows) {
            if (row.empty()) continue;
            PolyF4 p;
            for (const auto& [col_idx, coeff] : row) {
                p.terms[col_monomials[col_idx]] = coeff;
            }
            poly_rows.push_back(std::move(p));
        }
    }

    std::vector<PolyF4> get_reduced_polynomials() const {
        return poly_rows;
    }

private:
    MonomialComparator comp_;
};

// --- 3. Algoritmo F4 ---

static Monomial lcm(const Monomial& a, const Monomial& b) {
    Monomial res(a.size());
    for (size_t i = 0; i < a.size(); ++i) res[i] = std::max(a[i], b[i]);
    return res;
}

static bool divides(const Monomial& a, const Monomial& b) {
    for (size_t i = 0; i < a.size(); ++i) if (a[i] > b[i]) return false;
    return true;
}

struct Pair {
    size_t i, j;
    Monomial deg_lcm;
};

std::vector<PolyF4> f4_groebner(std::vector<PolyF4> G, MonomialOrder order) {
    size_t n_vars = 0;
    if (!G.empty()) n_vars = G[0].leading_monomial(order).size();
    
    for (auto& g : G) g.make_monic(order);

    MonomialComparator comp = get_comparator(order);
    std::vector<Pair> P;
    for (size_t i = 0; i < G.size(); ++i) {
        for (size_t j = i + 1; j < G.size(); ++j) {
            Monomial lm1 = G[i].leading_monomial(order);
            Monomial lm2 = G[j].leading_monomial(order);
            if (lm1.empty() || lm2.empty()) continue;
            
            bool coprime = true;
            for (size_t k = 0; k < n_vars; k++) {
                if (lm1[k] > 0 && lm2[k] > 0) { coprime = false; break; }
            }
            if (coprime) continue;

            P.push_back({i, j, lcm(lm1, lm2)});
        }
    }

    while (!P.empty()) {
        unsigned int min_deg = 0xFFFFFFFF;
        for (const auto& p : P) {
            unsigned int d = 0; for (unsigned int x : p.deg_lcm) d += x;
            if (d < min_deg) min_deg = d;
        }

        std::vector<Pair> L;
        std::vector<Pair> next_P;
        for (const auto& p : P) {
            unsigned int d = 0; for (unsigned int x : p.deg_lcm) d += x;
            if (d == min_deg) L.push_back(p);
            else next_P.push_back(p);
        }
        P = std::move(next_P);

        std::set<Monomial, MonomialComparator> done(comp);
        std::set<Monomial, MonomialComparator> todo(comp);
        MacaulayMatrix matrix(comp);

        for (const auto& p : L) {
            Monomial lm1 = G[p.i].leading_monomial(order);
            Monomial lm2 = G[p.j].leading_monomial(order);
            if (lm1.size() != n_vars || lm2.size() != n_vars) {
                throw std::runtime_error("Monomial size mismatch");
            }
            Monomial t1(n_vars), t2(n_vars);
            for (size_t k = 0; k < n_vars; ++k) {
                t1[k] = p.deg_lcm[k] - lm1[k];
                t2[k] = p.deg_lcm[k] - lm2[k];
            }
            
            auto add_mult = [&](const PolyF4& poly, const Monomial& t) {
                PolyF4 tp;
                for (const auto& [mon, coeff] : poly.terms) {
                    Monomial nm(n_vars);
                    for (size_t k = 0; k < n_vars; ++k) nm[k] = mon[k] + t[k];
                    tp.terms[nm] = coeff;
                    if (done.find(nm) == done.end()) todo.insert(nm);
                }
                matrix.add_polynomial(tp);
            };
            
            add_mult(G[p.i], t1);
            add_mult(G[p.j], t2);
        }

        while (!todo.empty()) {
            Monomial m = *todo.begin();
            todo.erase(todo.begin());
            done.insert(m);

            for (const auto& g : G) {
                Monomial lm = g.leading_monomial(order);
                if (!lm.empty() && divides(lm, m)) {
                    Monomial t(n_vars);
                    for (size_t k = 0; k < n_vars; ++k) t[k] = m[k] - lm[k];
                    
                    PolyF4 tg;
                    for (const auto& [mon, coeff] : g.terms) {
                        Monomial nm(n_vars);
                        for (size_t k = 0; k < n_vars; ++k) nm[k] = mon[k] + t[k];
                        tg.terms[nm] = coeff;
                        if (done.find(nm) == done.end()) todo.insert(nm);
                    }
                    matrix.add_polynomial(tg);
                    break;
                }
            }
        }

        matrix.gaussian_elimination();
        std::vector<PolyF4> reduced = matrix.get_reduced_polynomials();
        for (auto& h : reduced) {
            Monomial lm_h = h.leading_monomial(order);
            if (lm_h.empty()) continue;
            
            bool is_new = true;
            for (const auto& g : G) {
                if (g.leading_monomial(order) == lm_h) { is_new = false; break; }
            }
            if (is_new) {
                h.make_monic(order);
                G.push_back(h);
                size_t new_idx = G.size() - 1;
                for (size_t i = 0; i < new_idx; ++i) {
                    Monomial lm1 = G[i].leading_monomial(order);
                    bool coprime = true;
                    for (size_t k = 0; k < n_vars; k++) {
                        if (lm1[k] > 0 && lm_h[k] > 0) { coprime = false; break; }
                    }
                    if (coprime) continue;
                    P.push_back({i, new_idx, lcm(lm1, lm_h)});
                }
            }
        }
    }
    return G;
}

// --- 4. Inter-riduzione ---

void inter_reduce(std::vector<PolyF4>& G, MonomialOrder order) {
    if (G.empty()) return;
    size_t n_vars = G[0].leading_monomial(order).size();
    auto comp = get_comparator(order);

    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < G.size(); ++i) {
            PolyF4 f = G[i];
            // Reduce f mod G \ {G[i]}
            bool f_changed = false;
            auto it = f.terms.begin();
            while (it != f.terms.end()) {
                Monomial m = it->first;
                bool reduced = false;
                for (size_t j = 0; j < G.size(); ++j) {
                    if (i == j) continue;
                    Monomial lmj = G[j].leading_monomial(order);
                    if (divides(lmj, m)) {
                        Rational factor = it->second / G[j].leading_coefficient(order);
                        Monomial t(n_vars);
                        for (size_t k = 0; k < n_vars; ++k) t[k] = m[k] - lmj[k];
                        
                        // f = f - factor * G[j] * t
                        for (const auto& [mon_j, coeff_j] : G[j].terms) {
                            Monomial nm(n_vars);
                            for (size_t k = 0; k < n_vars; ++k) nm[k] = mon_j[k] + t[k];
                            f.terms[nm] = f.terms[nm] - factor * coeff_j;
                            if (f.terms[nm].numerator().is_zero()) f.terms.erase(nm);
                        }
                        f_changed = true;
                        reduced = true;
                        it = f.terms.begin(); // restart scan
                        break;
                    }
                }
                if (!reduced) ++it;
            }

            if (f_changed) {
                if (f.is_zero()) {
                    G.erase(G.begin() + i);
                    changed = true;
                    break;
                }
                f.make_monic();
                G[i] = std::move(f);
                changed = true;
            }
        }
    }
}

static Rational decimal_to_rational(const std::string& decimal) {
    size_t dot_pos = decimal.find('.');
    if (dot_pos == std::string::npos) {
        auto res = BigInt::parse(decimal);
        return Rational(res.is_ok() ? res.value() : BigInt(0));
    }
    std::string integral = decimal.substr(0, dot_pos);
    std::string fractional = decimal.substr(dot_pos + 1);
    auto res_num = BigInt::parse(integral + fractional);
    BigInt num = res_num.is_ok() ? res_num.value() : BigInt(0);
    BigInt den = BigInt(1);
    for (size_t i = 0; i < fractional.length(); ++i) den = den * BigInt(10);
    return Rational(num, den);
}

// --- 5. Conversioni ---

Result<PolyF4> expr_to_f4(ExprPtr expr, const std::vector<Symbol>& vars, symbolic::CASContext& ctx) {
    if (!expr) return fail<PolyF4>(make_error(CASErrorKind::InvalidArgument, "Null expression in expr_to_f4"));
    
    size_t n = vars.size();
    auto get_var_idx = [&](const std::string& name) -> std::optional<size_t> {
        for (size_t i = 0; i < n; ++i) if (vars[i].name == name) return i;
        return std::nullopt;
    };

    if (const auto* il = expr_cast<IntegerLit>(expr)) {
        PolyF4 p; p.terms[Monomial(n, 0)] = Rational(il->value); return ok(p);
    }
    if (const auto* rl = expr_cast<RationalLit>(expr)) {
        PolyF4 p; p.terms[Monomial(n, 0)] = Rational(rl->numerator, rl->denominator); return ok(p);
    }
    if (const auto* dl = expr_cast<DecimalLit>(expr)) {
        auto rat = decimal_to_rational(dl->text);
        PolyF4 p; p.terms[Monomial(n, 0)] = rat; return ok(p);
    }
    if (const auto* sym = expr_cast<Symbol>(expr)) {
        auto idx = get_var_idx(sym->name);
        if (idx) {
            PolyF4 p; Monomial mon(n, 0); mon[*idx] = 1;
            p.terms[mon] = Rational(1); return ok(p);
        }
        return fail<PolyF4>(make_error(CASErrorKind::Unimplemented, "Parametric constants not supported in expr_to_f4: " + sym->name));
    }
    if (const auto* un = expr_cast<Unary>(expr)) {
        if (un->op == UnaryOp::Neg) {
            auto r = expr_to_f4(un->operand, vars, ctx);
            if (r.is_error()) return r;
            for (auto& [mon, coeff] : r.value().terms) coeff = -coeff;
            return r;
        }
    }
    if (const auto* bin = expr_cast<Binary>(expr)) {
        auto lhs = expr_to_f4(bin->left, vars, ctx);
        auto rhs = expr_to_f4(bin->right, vars, ctx);
        if (lhs.is_error()) return lhs;
        if (rhs.is_error()) return rhs;

        if (bin->op == BinaryOp::Add) {
            for (const auto& [mon, coeff] : rhs.value().terms) {
                lhs.value().terms[mon] = lhs.value().terms[mon] + coeff;
                if (lhs.value().terms[mon].numerator().is_zero()) lhs.value().terms.erase(mon);
            }
            return lhs;
        }
        if (bin->op == BinaryOp::Sub) {
            for (const auto& [mon, coeff] : rhs.value().terms) {
                lhs.value().terms[mon] = lhs.value().terms[mon] - coeff;
                if (lhs.value().terms[mon].numerator().is_zero()) lhs.value().terms.erase(mon);
            }
            return lhs;
        }
        if (bin->op == BinaryOp::Mul) {
            PolyF4 res;
            for (const auto& [mon_l, coeff_l] : lhs.value().terms) {
                for (const auto& [mon_r, coeff_r] : rhs.value().terms) {
                    Monomial mon(n);
                    for (size_t k = 0; k < n; ++k) mon[k] = mon_l[k] + mon_r[k];
                    res.terms[mon] = res.terms[mon] + coeff_l * coeff_r;
                    if (res.terms[mon].numerator().is_zero()) res.terms.erase(mon);
                }
            }
            return ok(res);
        }
        if (bin->op == BinaryOp::Pow) {
            const auto* il = expr_cast<IntegerLit>(bin->right);
            if (il && !il->value.is_negative()) {
                long long exp_val = std::stoll(il->value.decimal());
                PolyF4 res; res.terms[Monomial(n, 0)] = Rational(1);
                PolyF4 base = lhs.value();
                while (exp_val > 0) {
                    if (exp_val % 2 == 1) {
                        PolyF4 next_res;
                        for (const auto& [ml, cl] : res.terms) {
                            for (const auto& [mb, cb] : base.terms) {
                                Monomial m(n); for (size_t k = 0; k < n; ++k) m[k] = ml[k] + mb[k];
                                next_res.terms[m] = next_res.terms[m] + cl * cb;
                                if (next_res.terms[m].numerator().is_zero()) next_res.terms.erase(m);
                            }
                        }
                        res = std::move(next_res);
                    }
                    PolyF4 next_base;
                    for (const auto& [m1, c1] : base.terms) {
                        for (const auto& [m2, c2] : base.terms) {
                            Monomial m(n); for (size_t k = 0; k < n; ++k) m[k] = m1[k] + m2[k];
                            next_base.terms[m] = next_base.terms[m] + c1 * c2;
                            if (next_base.terms[m].numerator().is_zero()) next_base.terms.erase(m);
                        }
                    }
                    base = std::move(next_base);
                    exp_val /= 2;
                }
                return ok(res);
            }
        }
    }
    if (const auto* s = expr_cast<Sum>(expr)) {
        PolyF4 res;
        for (auto term : s->terms) {
            auto r = expr_to_f4(term, vars, ctx);
            if (r.is_error()) return r;
            for (const auto& [mon, coeff] : r.value().terms) {
                res.terms[mon] = res.terms[mon] + coeff;
                if (res.terms[mon].numerator().is_zero()) res.terms.erase(mon);
            }
        }
        return ok(res);
    }
    if (const auto* p = expr_cast<Product>(expr)) {
        PolyF4 res; res.terms[Monomial(n, 0)] = Rational(1);
        for (auto factor : p->factors) {
            auto r = expr_to_f4(factor, vars, ctx);
            if (r.is_error()) return r;
            PolyF4 next_res;
            for (const auto& [ml, cl] : res.terms) {
                for (const auto& [mr, cr] : r.value().terms) {
                    Monomial m(n); for (size_t k = 0; k < n; ++k) m[k] = ml[k] + mr[k];
                    next_res.terms[m] = next_res.terms[m] + cl * cr;
                    if (next_res.terms[m].numerator().is_zero()) next_res.terms.erase(m);
                }
            }
            res = std::move(next_res);
        }
        return ok(res);
    }

    return fail<PolyF4>(make_error(CASErrorKind::Unimplemented, "expr_to_f4 unsupported: " + debug_print(expr)));
}

Result<ExprPtr> f4_to_expr(const PolyF4& p, const std::vector<Symbol>& vars, symbolic::CASContext& ctx) {
    std::vector<ExprPtr> sum_terms;
    AstArena& arena = ctx.arena();
    for (const auto& [mon, coeff] : p.terms) {
        ExprPtr term_expr = make_rational_expr(arena, coeff);
        for (size_t i = 0; i < vars.size(); ++i) {
            if (mon[i] > 0) {
                ExprPtr var_expr = arena.make<Symbol>(vars[i].name);
                if (mon[i] > 1) {
                    var_expr = arena.make<Binary>(BinaryOp::Pow, var_expr, arena.make<IntegerLit>(BigInt(mon[i])));
                }
                term_expr = arena.make<Binary>(BinaryOp::Mul, term_expr, var_expr);
            }
        }
        sum_terms.push_back(term_expr);
    }
    if (sum_terms.empty()) return ok(arena.make<IntegerLit>(BigInt(0)));
    if (sum_terms.size() == 1) return ok(sum_terms[0]);
    return ok(arena.make<Sum>(sum_terms));
}

// --- 6. Solver via Shape Lemma ---

Result<std::vector<std::vector<ExprPtr>>> solve_nonlinear_system_f4(
    const std::vector<ExprPtr>& equations,
    const std::vector<Symbol>& variables,
    symbolic::CASContext& ctx) {
    
    try {
        size_t n = variables.size();
        MonomialOrder order = MonomialOrder::Lex;
        std::vector<PolyF4> F;
        for (ExprPtr eq : equations) {
            auto p = expr_to_f4(eq, variables, ctx);
            if (p.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(p.error());
            F.push_back(p.value());
        }

        // Lex order GB for Shape Lemma
        auto G = f4_groebner(std::move(F), order);

        inter_reduce(G, order);

        if (G.empty()) return ok(std::vector<std::vector<ExprPtr>>{});

        // Shape Lemma: try to find a polynomial in only the last variable
        std::optional<PolyF4> pure_last;
        for (const auto& g : G) {
            bool only_last = true;
            for (const auto& [mon, coeff] : g.terms) {
                for (size_t i = 0; i < n - 1; ++i) {
                    if (mon[i] > 0) { only_last = false; break; }
                }
                if (!only_last) break;
            }
            if (only_last) { pure_last = g; break; }
        }

        if (!pure_last) {
            return fail<std::vector<std::vector<ExprPtr>>>(make_error(CASErrorKind::Unimplemented, "Shape lemma solver: pure polynomial not found"));
        }

        // Solve for the last variable
        auto last_expr_res = f4_to_expr(*pure_last, variables, ctx);
        if (last_expr_res.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(last_expr_res.error());

        auto last_roots_res = solve_polynomial(last_expr_res.value(), variables.back(), ctx);
        if (last_roots_res.is_error()) return fail<std::vector<std::vector<ExprPtr>>>(last_roots_res.error());

        std::vector<std::vector<ExprPtr>> all_solutions;
        for (ExprPtr root : last_roots_res.value()) {
            std::vector<ExprPtr> sol(n);
            sol[n - 1] = root;
            all_solutions.push_back(std::move(sol));
        }

        for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
            std::vector<std::vector<ExprPtr>> next_solutions;
            for (const auto& sol : all_solutions) {
                std::vector<ExprPtr> roots_for_i;
                for (const auto& g : G) {
                    Monomial lm = g.leading_monomial(order);
                    if (lm.size() > static_cast<size_t>(i) && lm[i] > 0) {
                        bool other_vars_ok = true;
                        for (size_t k = 0; k < n; ++k) {
                            if (k < static_cast<size_t>(i) && lm[k] > 0) { other_vars_ok = false; break; }
                        }
                        if (other_vars_ok) {
                            auto g_expr_res = f4_to_expr(g, variables, ctx);
                            if (g_expr_res.is_error()) continue;
                            ExprPtr g_expr = g_expr_res.value();
                            
                            for (size_t k = i + 1; k < n; ++k) {
                                auto sub = symbolic::substitute(g_expr, variables[k], sol[k], ctx);
                                if (sub.is_ok()) g_expr = sub.value();
                            }
                            auto solved = solve_polynomial(g_expr, variables[i], ctx);
                            if (solved.is_ok() && !solved.value().empty()) {
                                roots_for_i = solved.value();
                                break;
                            }
                        }
                    }
                }
                for (ExprPtr r : roots_for_i) {
                    std::vector<ExprPtr> new_sol = sol;
                    new_sol[i] = r;
                    next_solutions.push_back(std::move(new_sol));
                }
            }
            all_solutions = std::move(next_solutions);
        }

        return ok(all_solutions);
    } catch (const std::exception& e) {
        std::cerr << "EXCEPTION in solve_nonlinear_system_f4: " << e.what() << std::endl;
        throw;
    }
}

} // namespace cas::algebra
