#include "simplify_arithmetic_chain_impl.hpp"
#include <algorithm>

// F1.4: symbolic like-term collection helpers.
//
// Collects a·x + b·x → (a+b)·x for SYMBOLIC coefficients a, b (not just
// rational literals).
//
// Algorithm (coefficient / monomial decomposition):
//   Given a Sum {t₁, t₂, …, tₙ}, for each pair (ti, tj) compute the
//   intersection of their symbolic factor sets.  If the intersection is
//   non-empty AND no residual overlaps the shared monomial (canonicality
//   guard), replace with (qi + qj)·shared_monomial.
//
// The residual-overlap guard is critical for idempotence: without it,
//   x^3 + x → (x^2+1)*x which re-expands to x^3+x on the next simplify.
// The guard rejects any merge where a residual base appears in the shared
// monomial, preventing non-canonical factored forms for polynomial sums.

namespace cas::symbolic::detail {

// Internal struct (definition-local; the public API is the three functions
// declared in simplify_arithmetic_chain_impl.hpp).
struct SymbolicMonomial {
    std::vector<std::pair<ExprPtr, BigInt>> factors;
    ExprPtr coefficient;  // nullptr means coefficient is 1
};

bool decompose_term(
    ExprPtr expr,
    Rational& coeff_out,
    std::vector<std::pair<ExprPtr, BigInt>>& factors_out)
{
    coeff_out = Rational(BigInt(1));
    factors_out.clear();

    // Pure rational scalar
    {
        LiteralRational lr;
        auto ex = try_get_exact_rational(expr, lr);
        if (ex.is_ok() && ex.value()) {
            coeff_out = lr.value;
            return true;
        }
    }

    auto push_factor = [&](ExprPtr f) {
        if (const auto* bin = expr_cast<Binary>(f);
            bin && bin->op == BinaryOp::Pow) {
            if (auto exp = try_get_integer_exponent(bin->right); exp.has_value()) {
                factors_out.push_back({bin->left, *exp});
                return;
            }
        }
        factors_out.push_back({f, BigInt(1)});
    };

    // Single non-numeric factor
    if (expr_is<Symbol>(expr) || expr_is<FuncCall>(expr) || expr_is<Constant>(expr)) {
        push_factor(expr);
        return true;
    }

    // NOTA (A48): questo ramo NON ricorre di proposito. `-(2·x)` resta
    // coeff=-1 / fattori=[Product(2,x)] invece di coeff=-2 / fattori=[x].
    // L'asimmetria e' voluta: far ricorrere il ramo qui cambia le DECISIONI DI
    // MERGE di `try_merge_symbolic_like_terms`, che comincerebbe a fondere
    // coppie prima ignorate. Misurato: `-pi - pi*x` veniva fattorizzato in
    // `pi*(-1-x)` e la riduzione trigonometrica — che cerca la forma
    // distribuita — smetteva di vedere l'argomento, rompendo la riflessione
    // di Gamma (`gamma(x+2)*gamma(-1-x)` non arrivava piu' a `pi/sin(pi*x)`).
    // Il fabbisogno di A48 e' un'altra cosa — raccogliere COEFFICIENTI opposti
    // — ed e' risolto localmente in `collect_coefficient_terms`, che scarta il
    // segno prima di decomporre senza toccare questa euristica.
    if (const auto* u = expr_cast<Unary>(expr); u && u->op == UnaryOp::Neg) {
        coeff_out = Rational(BigInt(-1));
        expr = u->operand;
        LiteralRational lr;
        auto ex = try_get_exact_rational(expr, lr);
        if (ex.is_ok() && ex.value()) { coeff_out *= lr.value; return true; }
        push_factor(expr);
        return true;
    }

    if (const auto* prod = expr_cast<Product>(expr)) {
        for (ExprPtr f : prod->factors) {
            LiteralRational lr;
            auto ex = try_get_exact_rational(f, lr);
            if (ex.is_ok() && ex.value()) {
                coeff_out *= lr.value;
            } else {
                push_factor(f);
            }
        }
        std::sort(factors_out.begin(), factors_out.end(),
            [](const auto& a, const auto& b) {
                int c = canonical_compare(a.first, b.first);
                return c != 0 ? c < 0 : a.second < b.second;
            });
        return true;
    }

    if (const auto* bin = expr_cast<Binary>(expr);
        bin && bin->op == BinaryOp::Pow) {
        if (auto exp = try_get_integer_exponent(bin->right); exp.has_value()) {
            factors_out.push_back({bin->left, *exp});
            return true;
        }
    }

    // Cannot decompose further — leave as single opaque factor with exponent 1.
    factors_out.push_back({expr, BigInt(1)});
    return true;
}

ExprPtr build_coeff_monomial(
    ExprPtr coeff_expr,
    Rational numeric_coeff,
    const std::vector<std::pair<ExprPtr, BigInt>>& mono_factors,
    AstArena& arena)
{
    std::vector<ExprPtr> parts;
    bool neg = numeric_coeff.numerator().is_negative();
    Rational abs_coeff = neg ? -numeric_coeff : numeric_coeff;

    if (coeff_expr) {
        if (!(abs_coeff == Rational(BigInt(1)))) {
            parts.push_back(make_rational(arena, abs_coeff));
        }
        parts.push_back(coeff_expr);
    } else {
        if (!(abs_coeff == Rational(BigInt(1))) || mono_factors.empty()) {
            parts.push_back(make_rational(arena, abs_coeff));
        }
    }

    for (const auto& [base, exp] : mono_factors) {
        if (exp == BigInt(1)) {
            parts.push_back(base);
        } else {
            parts.push_back(arena.make<Binary>(BinaryOp::Pow, base,
                make_integer(arena, exp)));
        }
    }

    ExprPtr result;
    if (parts.empty()) {
        result = make_integer(arena, BigInt(1));
    } else if (parts.size() == 1U) {
        result = parts[0];
    } else {
        result = arena.make<Product>(std::move(parts));
    }
    return neg ? arena.make<Unary>(UnaryOp::Neg, result) : result;
}

namespace {

// A48 — raccoglie i termini simili DENTRO il Sum dei coefficienti.
//
// Il sito chiamante (simplify_arithmetic_chain_sum.cpp, Step 5) vieta — a
// ragione — di ri-simplificare il PRODOTTO fuso: rientrerebbe in
// `simplify_product_factors` e riaprirebbe il ciclo distribuzione/
// fattorizzazione. Quel divieto pero' NON riguarda il Sum dei coefficienti,
// che invece deve essere in forma raccolta come qualunque altro Sum.
//
// Senza questa raccolta il coefficiente veniva costruito grezzo (splice +
// sort, nessuna cancellazione) e `A - A` non si annullava: con
// A = (2x-2)·e^{2x} il risultato era `(-2 - 2x + 2x + 2)·e^{2x}` invece di 0,
// per QUALUNQUE A = polinomio·trascendente (repro misurato in A48).
//
// L'algoritmo e' quello dello Step 4 — chiave = lista di fattori canonica,
// valore = coefficiente razionale accumulato — applicato qui in locale. Non
// costruisce alcun Product, quindi non puo' innescare il ciclo che il divieto
// del chiamante protegge. `decompose_term` restituisce i fattori gia' ordinati
// canonicamente, quindi il confronto elemento per elemento e' lecito.
//
// Ritorna i termini superstiti: vuoto significa che il coefficiente si annulla
// per intero, cioe' che il termine fuso e' 0.
[[nodiscard]] std::vector<ExprPtr> collect_coefficient_terms(
    const std::vector<ExprPtr>& coeff_terms, AstArena& arena)
{
    struct Bucket {
        std::vector<std::pair<ExprPtr, BigInt>> factors;
        Rational coeff;
    };
    std::vector<Bucket> buckets;
    buckets.reserve(coeff_terms.size());

    for (ExprPtr term : coeff_terms) {
        // Scarta i `Neg` PRIMA di decomporre. `decompose_term` si ferma al
        // primo livello per non alterare le decisioni di merge (vedi la nota
        // sul suo ramo Neg), ma qui serve la chiave profonda: senza,
        // `-(2·x)` darebbe fattori=[Product(2,x)] e `2·x` fattori=[x], due
        // chiavi diverse per termini opposti, e `-(2x) + 2x` non cancellerebbe.
        // Lo scarto e' locale a questa raccolta, quindi non ha effetti
        // sull'euristica F1.4.
        Rational sign(BigInt(1));
        ExprPtr stripped = term;
        while (const auto* neg = expr_cast<Unary>(stripped)) {
            if (neg->op != UnaryOp::Neg) break;
            sign = -sign;
            stripped = neg->operand;
        }
        Rational c;
        std::vector<std::pair<ExprPtr, BigInt>> f;
        if (!decompose_term(stripped, c, f)) {
            buckets.push_back(Bucket{{{term, BigInt(1)}}, Rational(BigInt(1))});
            continue;
        }
        c *= sign;
        bool merged = false;
        for (Bucket& b : buckets) {
            if (b.factors.size() != f.size()) continue;
            bool same = true;
            for (std::size_t k = 0; k < f.size(); ++k) {
                if (!structural_equal(b.factors[k].first, f[k].first)
                    || b.factors[k].second != f[k].second) {
                    same = false;
                    break;
                }
            }
            if (same) {
                b.coeff += c;
                merged = true;
                break;
            }
        }
        if (!merged) buckets.push_back(Bucket{std::move(f), c});
    }

    std::vector<ExprPtr> collected;
    collected.reserve(buckets.size());
    for (const Bucket& b : buckets) {
        if (b.coeff.numerator().is_zero()) continue;
        collected.push_back(build_coeff_monomial(nullptr, b.coeff, b.factors, arena));
    }
    return collected;
}

}  // namespace

bool try_merge_symbolic_like_terms(
    std::vector<ExprPtr>& terms,
    AstArena& arena)
{
    const std::size_t n = terms.size();
    for (std::size_t i = 0; i < n; ++i) {
        Rational ci;
        std::vector<std::pair<ExprPtr, BigInt>> fi;
        if (!decompose_term(terms[i], ci, fi)) continue;
        if (fi.empty()) continue;  // pure scalar, skip

        for (std::size_t j = i + 1; j < n; ++j) {
            Rational cj;
            std::vector<std::pair<ExprPtr, BigInt>> fj;
            if (!decompose_term(terms[j], cj, fj)) continue;
            if (fj.empty()) continue;

            // Compute intersection of factor sets (structural equality for
            // interned nodes is pointer comparison — O(1)).
            std::vector<std::pair<ExprPtr, BigInt>> shared;
            for (const auto& [base_i, exp_i] : fi) {
                for (const auto& [base_j, exp_j] : fj) {
                    if (structural_equal(base_i, base_j)) {
                        BigInt min_exp = exp_i < exp_j ? exp_i : exp_j;
                        if (!min_exp.is_zero() && !min_exp.is_negative()) {
                            shared.push_back({base_i, min_exp});
                        }
                        break;
                    }
                }
            }
            if (shared.empty()) continue;

            // Polynomial-context guard: if ALL shared bases are plain Symbols
            // or powers of Symbols, skip.  In that context the monomial
            // collector in Step 4 handles like-term collection; running F1.4
            // on top would produce factored forms (e.g. (b+z)*x) that the
            // distribute step (Step 8) immediately expands back, causing
            // non-confluence with the distributive property tests.
            // F1.4 is only useful when at least one shared base is a
            // transcendental / radical expression (FuncCall / non-Symbol Power).
            {
                bool any_nonpoly_shared = false;
                for (const auto& [sb, _se] : shared) {
                    if (!expr_is<Symbol>(sb)) {
                        any_nonpoly_shared = true;
                        break;
                    }
                }
                if (!any_nonpoly_shared) continue;
            }

            // Compute residual factors after removing the shared part.
            auto subtract_shared = [&](const std::vector<std::pair<ExprPtr, BigInt>>& src)
                -> std::vector<std::pair<ExprPtr, BigInt>>
            {
                std::vector<std::pair<ExprPtr, BigInt>> res;
                for (const auto& [base, exp] : src) {
                    BigInt remaining = exp;
                    for (const auto& [sb, se] : shared) {
                        if (structural_equal(base, sb)) {
                            remaining = remaining - se;
                            break;
                        }
                    }
                    if (!remaining.is_zero()) res.push_back({base, remaining});
                }
                return res;
            };

            auto ri = subtract_shared(fi);
            auto rj = subtract_shared(fj);

            // Canonicality guard: reject merge if any residual base appears in
            // the shared monomial.  Prevents x^3 + x → (x^2+1)*x (would break
            // idempotence and the polynomial coefficient-extraction invariant).
            // E.g. x^3+x: shared={x:1}, ri={x:2}. x is in both → reject.
            auto residual_overlaps_shared =
                [&](const std::vector<std::pair<ExprPtr, BigInt>>& r) -> bool {
                for (const auto& [rb, _r] : r) {
                    for (const auto& [sb, _s] : shared) {
                        if (structural_equal(rb, sb)) return true;
                    }
                }
                return false;
            };
            if (residual_overlaps_shared(ri) || residual_overlaps_shared(rj)) continue;

            // Build quotients qi = ci*∏ri, qj = cj*∏rj.
            auto build_quotient = [&](Rational c,
                const std::vector<std::pair<ExprPtr, BigInt>>& r) -> ExprPtr
            {
                return build_coeff_monomial(nullptr, c, r, arena);
            };

            ExprPtr qi = build_quotient(ci, ri);
            ExprPtr qj = build_quotient(cj, rj);

            // Form (qi + qj) * shared_monomial.
            // Sort the coefficient sum canonically so the result is in the
            // same normal form on repeated simplify passes (idempotency).
            // A39: a quotient can itself be a Sum — e.g. collecting
            // (1-2x^2)*e - 1*e over the shared factor e gives qi = 1-2x^2.
            // Nesting it whole would produce Sum([Sum([...]), -1]), which the
            // strict-canonicity checker rejects and which no later pass
            // flattens (the Sum flatten in simplify_sum_terms runs before this
            // collection). Splice such terms in instead, before the sort, so
            // the canonical order is computed on the flat term list.
            std::vector<ExprPtr> coeff_terms;
            coeff_terms.reserve(2U);
            for (ExprPtr q : {qi, qj}) {
                if (const auto* nested = expr_cast<Sum>(q)) {
                    coeff_terms.insert(coeff_terms.end(),
                        nested->terms.begin(), nested->terms.end());
                } else {
                    coeff_terms.push_back(q);
                }
            }
            // A48: raccogli PRIMA di ordinare — altrimenti il coefficiente
            // resta grezzo e `A - A` non si annulla (vedi
            // collect_coefficient_terms).
            coeff_terms = collect_coefficient_terms(coeff_terms, arena);
            if (coeff_terms.empty()) {
                // Il coefficiente si annulla per intero: il termine fuso e' 0,
                // quindi i due addendi spariscono entrambi. Lasciare uno zero
                // letterale in lista lo porterebbe fino al Sum finale, perche'
                // dopo lo Step 5 non c'e' un'altra raccolta.
                terms.erase(terms.begin() + static_cast<std::ptrdiff_t>(j));
                terms.erase(terms.begin() + static_cast<std::ptrdiff_t>(i));
                return true;
            }
            std::sort(coeff_terms.begin(), coeff_terms.end(),
                [](ExprPtr a, ExprPtr b) {
                    int d = polynomial_degree(a) - polynomial_degree(b);
                    return d != 0 ? d > 0 : canonical_compare(a, b) < 0;
                });
            ExprPtr coeff_sum = coeff_terms.size() == 1U
                ? coeff_terms.front()
                : arena.make<Sum>(std::move(coeff_terms));

            // Build the list of factor expressions for the merged product.
            // Each entry is (base_for_sort, full_factor_expr).
            // base_for_sort is what merge_symbolic_factors uses as sort key:
            // for a plain expr it is the expr itself; for base^exp it is base.
            // Sorting here mirrors merge_symbolic_factors so that the Product
            // node created below already has canonical factor order.  Without
            // this sort pass-1 produces Product{coeff_sum, base^exp} while
            // pass-2 (via simplify_product_factors → merge_symbolic_factors)
            // produces Product{base^exp, coeff_sum}, violating idempotency.
            struct FactorEntry {
                ExprPtr sort_key;   // base used by canonical_compare
                ExprPtr factor_expr;
            };
            std::vector<FactorEntry> factor_entries;
            // A48 — prima della raccolta `coeff_sum` era SEMPRE un Sum (>= 2
            // termini spliciati), quindi poteva entrare nel Product come unico
            // fattore senza problemi. Ora puo' collassare a un singolo termine,
            // e quel termine puo' essere:
            //   * il letterale 1        -> non va aggiunto affatto (un fattore
            //                              `1` non e' canonico);
            //   * un Product (es. 2·a)  -> va SPLICIATO, non annidato: un
            //                              Product dentro un Product viola la
            //                              canonicita' e rompe l'idempotenza
            //                              (misurato su (a+b)·sin x + (a−b)·sin x,
            //                              che produceva Product{Product{2,a}, sin x}
            //                              e faceva scattare il canary A34).
            auto push_factor_entry = [&](ExprPtr f) {
                ExprPtr key = f;
                if (const auto* pw = expr_cast<Binary>(f);
                    pw != nullptr && pw->op == BinaryOp::Pow) {
                    key = pw->left;
                }
                factor_entries.push_back({key, f});
            };
            const auto* coeff_literal = expr_cast<IntegerLit>(coeff_sum);
            const bool coeff_is_one =
                coeff_literal != nullptr && coeff_literal->value == BigInt(1);
            if (!coeff_is_one) {
                if (const auto* coeff_product = expr_cast<Product>(coeff_sum)) {
                    for (ExprPtr f : coeff_product->factors) push_factor_entry(f);
                } else {
                    push_factor_entry(coeff_sum);
                }
            }
            for (const auto& [base, exp] : shared) {
                ExprPtr full = (exp == BigInt(1))
                    ? base
                    : arena.make<Binary>(BinaryOp::Pow, base, make_integer(arena, exp));
                factor_entries.push_back({base, full});
            }
            // Stable sort preserves relative order among equal-key entries,
            // matching the stable deterministic output of merge_symbolic_factors.
            std::stable_sort(factor_entries.begin(), factor_entries.end(),
                [](const FactorEntry& lhs, const FactorEntry& rhs) {
                    return canonical_compare(lhs.sort_key, rhs.sort_key) < 0;
                });
            std::vector<ExprPtr> merged_parts;
            merged_parts.reserve(factor_entries.size());
            for (const auto& fe : factor_entries)
                merged_parts.push_back(fe.factor_expr);
            ExprPtr merged = merged_parts.size() == 1U
                ? merged_parts[0]
                : arena.make<Product>(std::move(merged_parts));

            terms[i] = merged;
            terms.erase(terms.begin() + static_cast<std::ptrdiff_t>(j));
            return true;
        }
    }
    return false;
}

} // namespace cas::symbolic::detail
