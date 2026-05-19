#include "polynomial_groebner_f4.hpp"
#include "polynomial_groebner_f4_internal.hpp"
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
#include <unordered_map>
#include <vector>

namespace cas::algebra {

// --- 1. Comparatori e Hash Monomiali ---

struct MonomialHash {
    std::size_t operator()(const Monomial& v) const {
        std::size_t seed = v.size();
        for (auto x : v) {
            seed ^= std::hash<unsigned int>{}(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

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
    explicit MacaulayMatrix(MonomialComparator comp) : comp_(comp) {}

    std::vector<Monomial> col_monomials;
    std::unordered_map<Monomial, size_t, MonomialHash> monomial_to_col;
    std::vector<PolyF4> poly_rows; // Store polynomials directly, build matrix at the end

    void add_polynomial(const PolyF4& poly) {
        poly_rows.push_back(poly);
        for (const auto& [mon, coeff] : poly.terms) {
            register_monomial(mon);
        }
    }

    void register_monomial(const Monomial& mon) {
        // O(1) average insertion
        monomial_to_col.emplace(mon, SIZE_MAX);
    }

    void gaussian_elimination() {
        if (poly_rows.empty()) return;

        // Assign column indices once, in sorted order.
        // We must sort the collected monomials using comp_ (largest first).
        col_monomials.clear();
        col_monomials.reserve(monomial_to_col.size());
        for (const auto& [mon, _] : monomial_to_col) {
            col_monomials.push_back(mon);
        }
        std::sort(col_monomials.begin(), col_monomials.end(), comp_);

        for (size_t i = 0; i < col_monomials.size(); ++i) {
            monomial_to_col[col_monomials[i]] = i;
        }

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

struct Pair {
    size_t i, j;
    Monomial deg_lcm;
};

static Monomial lcm(const Monomial& a, const Monomial& b) {
    Monomial res(a.size());
    for (size_t i = 0; i < a.size(); ++i) res[i] = std::max(a[i], b[i]);
    return res;
}

static bool divides(const Monomial& a, const Monomial& b) {
    for (size_t i = 0; i < a.size(); ++i) if (a[i] > b[i]) return false;
    return true;
}

Result<std::vector<PolyF4>> f4_groebner(std::vector<PolyF4> G, MonomialOrder order) {
    // Termination follows from the Hilbert basis theorem (every ideal in
    // K[x_1, ..., x_n] is finitely generated) and Buchberger termination
    // theorem: under sugar-ordered Gebauer-Moeller pair pruning the pair
    // queue empties in finite time. No artificial batch counter is
    // therefore load-bearing for correctness; the previous
    //
    //     if (++batches > kMaxF4Batches=2048)
    //         return buchberger_groebner(buchberger_seed, order);
    //
    // bail has been removed so F4 runs to completion on its own pair
    // queue. Memory-side guards on the Macaulay matrix (rows / monomial
    // count) remain below as a hardware safety net; on overflow we still
    // fall back to Buchberger, which (post-Sugar) handles the seed at
    // theoretical minimum basis cardinality.
    constexpr std::size_t kMaxMacaulayRows = 512;
    constexpr std::size_t kMaxMacaulayMonomials = 512;
    constexpr std::size_t kMaxPendingMonomials = 1024;
    size_t n_vars = 0;
    if (!G.empty()) n_vars = G[0].leading_monomial(order).size();
    
    for (auto& g : G) g.make_monic(order);
    const std::vector<PolyF4> buchberger_seed = G;

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
                return fail<std::vector<PolyF4>>(make_error(
                    CASErrorKind::InvalidArgument,
                    "F4 input contains inconsistent monomial arity"));
            }
            Monomial t1(n_vars), t2(n_vars);
            for (size_t k = 0; k < n_vars; ++k) {
                t1[k] = p.deg_lcm[k] - lm1[k];
                t2[k] = p.deg_lcm[k] - lm2[k];
            }
            
            auto add_mult = [&](const PolyF4& poly, const Monomial& t) -> Result<void> {
                if (matrix.poly_rows.size() >= kMaxMacaulayRows ||
                    matrix.monomial_to_col.size() >= kMaxMacaulayMonomials ||
                    todo.size() + done.size() >= kMaxPendingMonomials) {
                    return fail<void>(make_error(
                        CASErrorKind::Timeout,
                        "F4 Macaulay matrix construction exceeded resource guard"));
                }
                PolyF4 tp;
                for (const auto& [mon, coeff] : poly.terms) {
                    Monomial nm(n_vars);
                    for (size_t k = 0; k < n_vars; ++k) nm[k] = mon[k] + t[k];
                    tp.terms[nm] = coeff;
                    if (done.find(nm) == done.end()) todo.insert(nm);
                }
                matrix.add_polynomial(tp);
                return ok();
            };
            
            auto add_left = add_mult(G[p.i], t1);
            if (add_left.is_error()) return detail::buchberger_groebner(buchberger_seed, order);
            auto add_right = add_mult(G[p.j], t2);
            if (add_right.is_error()) return detail::buchberger_groebner(buchberger_seed, order);
        }

        while (!todo.empty()) {
            if (matrix.poly_rows.size() >= kMaxMacaulayRows ||
                matrix.monomial_to_col.size() >= kMaxMacaulayMonomials ||
                todo.size() + done.size() >= kMaxPendingMonomials) {
                return detail::buchberger_groebner(buchberger_seed, order);
            }
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
    return ok(std::move(G));
}

} // namespace cas::algebra
