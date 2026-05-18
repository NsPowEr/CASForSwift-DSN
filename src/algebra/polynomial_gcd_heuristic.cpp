#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "cas/numtheory.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace cas::algebra {

namespace {

using Monomial = std::vector<unsigned int>;
using SparsePoly = std::map<Monomial, BigInt>;

// Kronecker substitution: x_i = B^{(D+1)^{i-1}}
// Mapping multivariate polynomial to BigInt
Result<BigInt> kronecker_evaluate(const MultivariatePolynomial& poly, const std::vector<Symbol>& vars, const BigInt& B, std::size_t D_plus_1) {
    BigInt total(0);
    for (const auto& term : poly.terms()) {
        BigInt term_val = term.coefficient;
        [[maybe_unused]] std::size_t exponent_sum = 0;
        std::size_t weight = 1;
        
        // Calcoliamo l'esponente combinato k = sum e_i * (D+1)^{i-1}
        std::size_t k = 0;
        for (std::size_t i = 0; i < vars.size(); ++i) {
            unsigned int e_i = 0;
            for (const auto& factor : term.factors) {
                if (factor.first.name == vars[i].name) {
                    e_i = factor.second;
                    break;
                }
            }
            k += e_i * weight;
            weight *= D_plus_1;
        }
        
        // term_val *= B^k
        if (k > 0) {
            term_val *= bigint_pow_nonnegative(B, k);
        }
        total += term_val;
    }
    return ok(std::move(total));
}

// Reconstruction from BigInt in base B
Result<MultivariatePolynomial> kronecker_reconstruct(BigInt g, const std::vector<Symbol>& vars, const BigInt& B, std::size_t D_plus_1) {
    std::vector<MultivariateTerm> terms;
    
    // We need to write g in base B with signed digits in [-B/2, B/2]
    BigInt half_B = B / BigInt(2);
    
    std::size_t k = 0;
    while (!g.is_zero()) {
        BigInt rem = g % B;
        g /= B;
        
        if (rem > half_B) {
            rem -= B;
            g += BigInt(1);
        } else if (rem < -half_B) {
            rem += B;
            g -= BigInt(1);
        }
        
        if (!rem.is_zero()) {
            // Decodifica k in (e_1, ..., e_n) base D+1
            std::vector<std::pair<Symbol, unsigned int>> factors;
            std::size_t temp_k = k;
            for (std::size_t i = 0; i < vars.size(); ++i) {
                unsigned int e_i = static_cast<unsigned int>(temp_k % D_plus_1);
                temp_k /= D_plus_1;
                if (e_i > 0) {
                    factors.push_back({vars[i], e_i});
                }
            }
            
            terms.push_back({
                .coefficient = rem,
                .factors = std::move(factors)
            });
        }
        k++;
    }
    
    return ok(MultivariatePolynomial(std::move(terms)));
}

[[nodiscard]] std::vector<Symbol> collect_variables(
    const MultivariatePolynomial& lhs,
    const MultivariatePolynomial& rhs) {
    std::vector<Symbol> vars = lhs.variables();
    for (const auto& var : rhs.variables()) {
        if (std::find_if(vars.begin(), vars.end(), [&](const Symbol& existing) {
                return existing.name == var.name;
            }) == vars.end()) {
            vars.push_back(var);
        }
    }
    std::sort(vars.begin(), vars.end(), [](const Symbol& a, const Symbol& b) {
        return a.name < b.name;
    });
    return vars;
}

[[nodiscard]] SparsePoly to_sparse(const MultivariatePolynomial& poly, const std::vector<Symbol>& vars) {
    SparsePoly sparse;
    for (const auto& term : poly.terms()) {
        Monomial monomial(vars.size(), 0U);
        for (const auto& [symbol, exponent] : term.factors) {
            auto it = std::find_if(vars.begin(), vars.end(), [&](const Symbol& var) {
                return var.name == symbol.name;
            });
            if (it != vars.end()) {
                monomial[static_cast<std::size_t>(std::distance(vars.begin(), it))] += exponent;
            }
        }
        sparse[monomial] += term.coefficient;
        if (sparse[monomial].is_zero()) {
            sparse.erase(monomial);
        }
    }
    return sparse;
}

[[nodiscard]] bool monomial_divides(const Monomial& divisor, const Monomial& dividend) {
    for (std::size_t i = 0U; i < divisor.size(); ++i) {
        if (divisor[i] > dividend[i]) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] Monomial monomial_quotient(const Monomial& dividend, const Monomial& divisor) {
    Monomial quotient(dividend.size(), 0U);
    for (std::size_t i = 0U; i < dividend.size(); ++i) {
        quotient[i] = dividend[i] - divisor[i];
    }
    return quotient;
}

void add_sparse_term(SparsePoly& poly, const Monomial& monomial, const BigInt& coefficient) {
    if (coefficient.is_zero()) {
        return;
    }
    poly[monomial] += coefficient;
    if (poly[monomial].is_zero()) {
        poly.erase(monomial);
    }
}

void subtract_sparse_term_multiple(
    SparsePoly& remainder,
    const SparsePoly& divisor,
    const Monomial& monomial,
    const BigInt& coefficient) {
    for (const auto& [base_monomial, base_coefficient] : divisor) {
        Monomial product(base_monomial.size(), 0U);
        for (std::size_t i = 0U; i < base_monomial.size(); ++i) {
            product[i] = base_monomial[i] + monomial[i];
        }
        add_sparse_term(remainder, product, -(base_coefficient * coefficient));
    }
}

[[nodiscard]] bool divides_exactly(
    const MultivariatePolynomial& dividend,
    const MultivariatePolynomial& divisor,
    const std::vector<Symbol>& vars) {
    if (divisor.is_zero()) {
        return false;
    }

    SparsePoly remainder = to_sparse(dividend, vars);
    const SparsePoly divisor_sparse = to_sparse(divisor, vars);
    if (divisor_sparse.empty()) {
        return false;
    }

    const auto [divisor_lm, divisor_lc] = *std::prev(divisor_sparse.end());
    const std::size_t max_steps = (remainder.size() + 1U) * (divisor_sparse.size() + 1U) * (vars.size() + 1U);
    std::size_t steps = 0U;

    while (!remainder.empty()) {
        if (++steps > max_steps) {
            return false;
        }

        const auto [remainder_lm, remainder_lc] = *std::prev(remainder.end());
        if (!monomial_divides(divisor_lm, remainder_lm)) {
            return false;
        }
        if ((remainder_lc % divisor_lc) != BigInt(0)) {
            return false;
        }

        const BigInt term_coefficient = remainder_lc / divisor_lc;
        const Monomial term_monomial = monomial_quotient(remainder_lm, divisor_lm);
        subtract_sparse_term_multiple(remainder, divisor_sparse, term_monomial, term_coefficient);
    }

    return true;
}

[[nodiscard]] bool verify_gcd_candidate(
    const MultivariatePolynomial& P,
    const MultivariatePolynomial& Q,
    const MultivariatePolynomial& G) {
    std::vector<Symbol> all_vars = collect_variables(P, Q);
    for (const auto& var : G.variables()) {
        if (std::find_if(all_vars.begin(), all_vars.end(), [&](const Symbol& existing) {
                return existing.name == var.name;
            }) == all_vars.end()) {
            all_vars.push_back(var);
        }
    }
    std::sort(all_vars.begin(), all_vars.end(), [](const Symbol& a, const Symbol& b) {
        return a.name < b.name;
    });
    return divides_exactly(P, G, all_vars) && divides_exactly(Q, G, all_vars);
}

} // namespace

Result<MultivariatePolynomial> gcd_heuristic(const MultivariatePolynomial& P, const MultivariatePolynomial& Q) {
    if (P.is_zero()) return ok(Q);
    if (Q.is_zero()) return ok(P);

    auto vars_p = P.variables();
    auto vars_q = Q.variables();
    std::vector<Symbol> all_vars = vars_p;
    for (const auto& v : vars_q) {
        if (std::find_if(all_vars.begin(), all_vars.end(), [&](const Symbol& s){ return s.name == v.name; }) == all_vars.end()) {
            all_vars.push_back(v);
        }
    }
    std::sort(all_vars.begin(), all_vars.end(), [](const Symbol& a, const Symbol& b){ return a.name < b.name; });

    // Determine D = max degree in any variable
    std::size_t D = 0;
    for (const auto& v : all_vars) {
        std::size_t deg_p = 0;
        for (const auto& t : P.terms()) {
            for (const auto& f : t.factors) if (f.first.name == v.name) deg_p = std::max(deg_p, (std::size_t)f.second);
        }
        std::size_t deg_q = 0;
        for (const auto& t : Q.terms()) {
            for (const auto& f : t.factors) if (f.first.name == v.name) deg_q = std::max(deg_q, (std::size_t)f.second);
        }
        D = std::max({D, deg_p, deg_q});
    }

    // Determine max coefficient magnitude
    BigInt max_coeff(0);
    for (const auto& t : P.terms()) max_coeff = std::max(max_coeff, t.coefficient.abs());
    for (const auto& t : Q.terms()) max_coeff = std::max(max_coeff, t.coefficient.abs());

    // Mignotte bound (Knuth TAOCP vol2 §4.6.2): GCD coefficients bounded by
    // sqrt(D+1) * 2^D * max_coeff. We use the strict envelope
    //   B = 2 * (max_coeff + 1) * 2^(D+1) + 1
    // which dominates Mignotte for all (D, max_coeff) and guarantees no
    // collision in the Kronecker substitution image. The +1 keeps signed-digit
    // reconstruction unambiguous (digit > B/2 ⇒ negative branch).
    // Minimum structural floor B ≥ 3 is required so that signed digits
    // {-1, 0, +1} are representable; for D=0, max_coeff=0 this is the only
    // case where the formula would otherwise produce B<3.
    BigInt B = max_coeff + BigInt(1);
    B = B * BigInt(2);
    for (std::size_t i = 0; i <= D; ++i) B = B * BigInt(2);
    B = B + BigInt(1);
    if (B < BigInt(3)) B = BigInt(3);

    std::size_t D_plus_1 = D + 1;
    
    auto p_val = kronecker_evaluate(P, all_vars, B, D_plus_1);
    if (p_val.is_error()) return fail<MultivariatePolynomial>(p_val.error());
    
    auto q_val = kronecker_evaluate(Q, all_vars, B, D_plus_1);
    if (q_val.is_error()) return fail<MultivariatePolynomial>(q_val.error());
    
    BigInt g_val = gcd(p_val.value().abs(), q_val.value().abs());
    
    auto G = kronecker_reconstruct(g_val, all_vars, B, D_plus_1);
    if (G.is_error()) return G;
    
    // GCDHEU is heuristic, so the reconstructed candidate must be certified by
    // exact multivariate division before it can leave the algorithm boundary.
    if (verify_gcd_candidate(P, Q, G.value())) {
        return G;
    }
    
    return fail<MultivariatePolynomial>(make_error(CASErrorKind::InternalError, "GCDHEU failed to reconstruct valid GCD"));
}

} // namespace cas::algebra
