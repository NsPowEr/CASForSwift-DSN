// F3.2 — Sparse multivariate polynomial over Z (MPoly): exact arithmetic.
// See factor_multivariate_internal.hpp for the contract.

#include "factor_multivariate_internal.hpp"

#include <algorithm>
#include <utility>

namespace cas::algebra {

namespace {

[[nodiscard]] Monomial mono_add(const Monomial& a, const Monomial& b) {
    Monomial r(a.size(), 0U);
    for (std::size_t i = 0; i < a.size(); ++i) {
        r[i] = a[i] + b[i];
    }
    return r;
}

[[nodiscard]] std::optional<Monomial> mono_sub(const Monomial& a, const Monomial& b) {
    Monomial r(a.size(), 0U);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] < b[i]) {
            return std::nullopt;
        }
        r[i] = a[i] - b[i];
    }
    return r;
}

}  // namespace

void MPoly::prune() {
    for (auto it = terms.begin(); it != terms.end();) {
        if (it->second.is_zero()) {
            it = terms.erase(it);
        } else {
            ++it;
        }
    }
}

MPoly mpoly_zero() { return MPoly{}; }

MPoly mpoly_constant(const BigInt& c, std::size_t nvars) {
    MPoly p;
    if (!c.is_zero()) {
        p.terms.emplace(Monomial(nvars, 0U), c);
    }
    return p;
}

MPoly mpoly_add(const MPoly& a, const MPoly& b) {
    MPoly r = a;
    for (const auto& [mono, coeff] : b.terms) {
        BigInt& slot = r.terms[mono];
        slot += coeff;
        if (slot.is_zero()) {
            r.terms.erase(mono);
        }
    }
    return r;
}

MPoly mpoly_sub(const MPoly& a, const MPoly& b) {
    MPoly r = a;
    for (const auto& [mono, coeff] : b.terms) {
        BigInt& slot = r.terms[mono];
        slot -= coeff;
        if (slot.is_zero()) {
            r.terms.erase(mono);
        }
    }
    return r;
}

MPoly mpoly_scale(const MPoly& a, const BigInt& s) {
    MPoly r;
    if (s.is_zero()) {
        return r;
    }
    for (const auto& [mono, coeff] : a.terms) {
        r.terms.emplace(mono, coeff * s);
    }
    return r;
}

MPoly mpoly_mul(const MPoly& a, const MPoly& b) {
    MPoly r;
    for (const auto& [ma, ca] : a.terms) {
        for (const auto& [mb, cb] : b.terms) {
            Monomial mm = mono_add(ma, mb);
            BigInt prod = ca * cb;
            BigInt& slot = r.terms[mm];
            slot += prod;
            if (slot.is_zero()) {
                r.terms.erase(mm);
            }
        }
    }
    return r;
}

bool mpoly_equal(const MPoly& a, const MPoly& b) {
    return a.terms == b.terms;
}

unsigned int mpoly_degree_in(const MPoly& p, std::size_t var_index) {
    unsigned int d = 0U;
    for (const auto& [mono, coeff] : p.terms) {
        (void)coeff;
        d = std::max(d, mono[var_index]);
    }
    return d;
}

BigInt mpoly_integer_content(const MPoly& p) {
    BigInt g(0);
    for (const auto& [mono, coeff] : p.terms) {
        (void)mono;
        g = gcd(g, coeff.abs());
    }
    return g.is_zero() ? BigInt(1) : g;
}

MPoly mpoly_coeff_in(const MPoly& p, std::size_t var_index, unsigned int deg) {
    MPoly r;
    for (const auto& [mono, coeff] : p.terms) {
        if (mono[var_index] == deg) {
            Monomial m = mono;
            m[var_index] = 0U;
            r.terms[m] += coeff;
        }
    }
    r.prune();
    return r;
}

MPoly mpoly_leading_coeff_in(const MPoly& p, std::size_t var_index) {
    if (p.is_zero()) {
        return MPoly{};
    }
    return mpoly_coeff_in(p, var_index, mpoly_degree_in(p, var_index));
}

MPoly mpoly_eval_var(const MPoly& p, std::size_t var_index, const BigInt& value) {
    MPoly r;
    for (const auto& [mono, coeff] : p.terms) {
        BigInt c = coeff;
        unsigned int e = mono[var_index];
        if (e > 0U) {
            c *= bigint_pow_nonnegative(value, static_cast<std::size_t>(e));
        }
        Monomial m = mono;
        m[var_index] = 0U;
        BigInt& slot = r.terms[m];
        slot += c;
        if (slot.is_zero()) {
            r.terms.erase(m);
        }
    }
    return r;
}

std::optional<IntPoly> mpoly_to_intpoly(const MPoly& p, std::size_t var_index) {
    IntPoly result;
    for (const auto& [mono, coeff] : p.terms) {
        for (std::size_t i = 0; i < mono.size(); ++i) {
            if (i != var_index && mono[i] != 0U) {
                return std::nullopt;  // genuinely multivariate
            }
        }
        std::size_t d = static_cast<std::size_t>(mono[var_index]);
        if (result.size() <= d) {
            result.coefficients().resize(d + 1U, BigInt(0));
        }
        result[d] += coeff;
    }
    normalize_integer_poly(result);
    return result;
}

MPoly mpoly_from_intpoly(const IntPoly& f, std::size_t var_index, std::size_t nvars) {
    MPoly r;
    for (std::size_t d = 0; d < f.size(); ++d) {
        if (f[d].is_zero()) {
            continue;
        }
        Monomial m(nvars, 0U);
        m[var_index] = static_cast<unsigned int>(d);
        r.terms.emplace(std::move(m), f[d]);
    }
    return r;
}

// Exact division over Z via leading-monomial reduction (lex on the dense
// exponent vector, which std::map already provides through Monomial ordering).
std::optional<MPoly> mpoly_exact_div(const MPoly& a, const MPoly& b) {
    if (b.is_zero()) {
        return std::nullopt;
    }
    if (a.is_zero()) {
        return MPoly{};
    }
    // Leading term of b = largest monomial under lex (map is sorted ascending).
    const auto b_lead = std::prev(b.terms.end());
    const Monomial& b_lm = b_lead->first;
    const BigInt& b_lc = b_lead->second;

    MPoly rem = a;
    MPoly quot;
    while (!rem.is_zero()) {
        const auto r_lead = std::prev(rem.terms.end());
        const Monomial& r_lm = r_lead->first;
        const BigInt& r_lc = r_lead->second;

        auto delta = mono_sub(r_lm, b_lm);
        if (!delta.has_value()) {
            return std::nullopt;  // not divisible
        }
        BigInt r = r_lc % b_lc;
        if (!r.is_zero()) {
            return std::nullopt;  // coefficient not divisible
        }
        BigInt q = r_lc / b_lc;
        // quot += q * x^delta
        BigInt& qslot = quot.terms[*delta];
        qslot += q;
        if (qslot.is_zero()) {
            quot.terms.erase(*delta);
        }
        // rem -= (q * x^delta) * b
        MPoly term;
        term.terms.emplace(*delta, q);
        rem = mpoly_sub(rem, mpoly_mul(term, b));
    }
    return quot;
}

}  // namespace cas::algebra
