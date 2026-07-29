#include "polynomial_internal.hpp"
#include "cas/ast.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>
#include <vector>

namespace cas::algebra {

namespace {

// Trial-division primality test.
[[nodiscard]] bool is_prime(int n) {
    if (n < 2) return false;
    if (n < 4) return true;          // 2, 3
    if (n % 2 == 0) return false;
    for (int i = 3; static_cast<long long>(i) * i <= n; i += 2) {
        if (n % i == 0) return false;
    }
    return true;
}

// Distinct prime factors of n (ascending), via trial division.
[[nodiscard]] std::vector<int> prime_factors_distinct(int n) {
    std::vector<int> ps;
    for (int p = 2; static_cast<long long>(p) * p <= n; ++p) {
        if (n % p == 0) {
            ps.push_back(p);
            while (n % p == 0) n /= p;
        }
    }
    if (n > 1) ps.push_back(n);
    return ps;
}

[[nodiscard]] std::vector<int> divisors_of(int n) {
    std::vector<int> d;
    for (int i = 1; static_cast<long long>(i) * i <= n; ++i) {
        if (n % i == 0) {
            d.push_back(i);
            if (i != n / i) d.push_back(n / i);
        }
    }
    std::sort(d.begin(), d.end());
    return d;
}

// Exact polynomial division a / b in Z[x] (assumes b monic and b | a, which
// holds for the cyclotomic identity Φ_k(x^p) = Φ_k(x) · Φ_{kp}(x) when p ∤ k).
[[nodiscard]] IntPoly poly_div_exact(const IntPoly& a, const IntPoly& b) {
    if (b.empty() || a.size() < b.size()) return IntPoly{};
    std::vector<BigInt> r = a.coefficients();
    const std::size_t da = r.size() - 1U;
    const std::size_t db = b.degree();
    std::vector<BigInt> q(da - db + 1U, BigInt(0));
    for (int i = static_cast<int>(da - db); i >= 0; --i) {
        q[i] = r[i + db];  // b is monic (leading coeff = 1)
        for (std::size_t j = 0; j <= db; ++j) {
            r[i + j] -= q[i] * b[j];
        }
    }
    IntPoly out(std::move(q));
    out.normalize([](const BigInt& v) { return v.is_zero(); });
    return out;
}

// Substitute x → x^p:  Φ(x) ↦ Φ(x^p).  Coefficient at degree i moves to i·p.
[[nodiscard]] IntPoly subst_x_pow(const IntPoly& f, int p) {
    if (f.empty() || p <= 1) return f;
    std::vector<BigInt> c(f.degree() * static_cast<std::size_t>(p) + 1U, BigInt(0));
    for (std::size_t i = 0; i < f.size(); ++i) {
        c[i * static_cast<std::size_t>(p)] = f[i];
    }
    return IntPoly(std::move(c));  // leading coeff preserved (nonzero)
}

// Φ_p(x) = 1 + x + … + x^(p-1) for prime p.
[[nodiscard]] IntPoly cyclotomic_prime(int p) {
    return IntPoly(std::vector<BigInt>(static_cast<std::size_t>(p), BigInt(1)));
}

bool poly_equal(const IntPoly& a, const IntPoly& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

// Inverse totient: all n with φ(n) = d (ascending, deduplicated).
//
// A prime p can divide n only if (p-1) | φ(n) = d.  Building n in strictly
// increasing prime order makes each canonical factorisation appear once.
// For a prime power p^j, the φ-contribution is (p-1)·p^(j-1).  Every n with
// φ(n) = d satisfies n ≤ max(6, 2d²) (Rosser–Schoenfeld), which bounds the
// search and prevents int overflow.
void inverse_totient_rec(long long d_rem, int min_prime, long long cur_n,
                         long long n_cap, std::vector<long long>& out) {
    if (d_rem == 1) {
        out.push_back(cur_n);
        // Do NOT return: when 2 has not yet been used (min_prime < 2, i.e. the
        // top frame for d = 1), the prime 2 still has φ-contribution 1 and must
        // be appended so that φ(2n) = φ(n) solutions (e.g. n = 2) are found.
        // For d ≥ 2 every frame reaching d_rem == 1 has min_prime ≥ 2, so the
        // p > min_prime guard below makes the loop a no-op (no duplicates).
    }
    for (int s : divisors_of(static_cast<int>(d_rem))) {
        const int p = s + 1;                 // p-1 = s must divide d_rem
        if (p <= min_prime || !is_prime(p)) continue;
        long long rem = d_rem / s;            // consume (p-1) for p^1
        long long pk = p;                     // p^j
        if (cur_n * pk <= n_cap) {
            inverse_totient_rec(rem, p, cur_n * pk, n_cap, out);
        }
        // Higher prime powers p^j (j ≥ 2): extra φ-factor p each time.
        while (rem % p == 0) {
            rem /= p;
            pk *= p;
            if (cur_n * pk > n_cap) break;
            inverse_totient_rec(rem, p, cur_n * pk, n_cap, out);
        }
    }
}

[[nodiscard]] std::vector<int> inverse_totient(int d) {
    const long long n_cap = std::max<long long>(6, 2LL * d * d);
    std::vector<long long> raw;
    inverse_totient_rec(d, 1, 1, n_cap, raw);
    std::sort(raw.begin(), raw.end());
    raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
    std::vector<int> out;
    out.reserve(raw.size());
    for (long long n : raw) out.push_back(static_cast<int>(n));
    return out;
}

} // namespace

// Public: Φ_n(x), the n-th cyclotomic polynomial (monic, degree φ(n)).
//
// Squarefree reduction keeps every intermediate at degree O(φ(n)), so there is
// no order-n materialisation and no OOM cap (this lifts the former deg ≤ 724
// limitation, ledger A5-LARGECYCLO):
//
//   1. Φ_rad(n)(x) via the fold  Φ_{k·p}(x) = Φ_k(x^p) / Φ_k(x)  over the
//      distinct primes p of n (each step's intermediate has degree ≤ 2·φ(n)).
//   2. Φ_n(x) = Φ_rad(n)(x^(n/rad(n)))  (Apostol, "Introduction to Analytic
//      Number Theory", Thm 8.13).
//
// Returns the empty IntPoly only for n ≤ 0 (invalid input).
IntPoly compute_cyclotomic(int n) {
    if (n <= 0) return IntPoly{};
    if (n == 1) return IntPoly({BigInt(-1), BigInt(1)});  // x - 1

    const std::vector<int> primes = prime_factors_distinct(n);

    IntPoly phi;
    int rad = 1;
    for (int p : primes) {
        if (rad == 1) {
            phi = cyclotomic_prime(p);                 // Φ_p
        } else {
            // p ∤ rad here, so Φ_{rad·p}(x) = Φ_rad(x^p) / Φ_rad(x).
            phi = poly_div_exact(subst_x_pow(phi, p), phi);
        }
        rad *= p;
    }

    const int m = n / rad;                              // n / rad(n) = ∏ p_i^(a_i-1)
    if (m > 1) phi = subst_x_pow(phi, m);               // Φ_n(x) = Φ_rad(x^m)
    return phi;
}

// Public: if `poly` is the n-th cyclotomic polynomial, return n (smallest such),
// else nullopt.  Only the candidates {n : φ(n) = deg(poly)} are tested
// (inverse totient), instead of scanning every n up to 2d².
//
// `max_n ≥ 0` restricts the search to n ≤ max_n (from ctx.max_cyclotomic_n());
// max_n < 0 means "no extra cap" — the inverse-totient set is already finite.
std::optional<int> is_cyclotomic(const IntPoly& poly, int max_n) {
    if (poly.empty()) return std::nullopt;
    if (poly.leading_coeff() != BigInt(1)) return std::nullopt;

    const std::size_t deg = poly.degree();
    if (deg == 0) return std::nullopt;

    // deg fits comfortably in int for practical CAS inputs (2·d² must not
    // overflow the inverse-totient search; d < ~32000 keeps 2d² < INT_MAX).
    if (deg > static_cast<std::size_t>(30000)) return std::nullopt;
    const int d = static_cast<int>(deg);

    std::vector<int> candidates = inverse_totient(d);   // ascending → smallest n first
    for (int n : candidates) {
        if (max_n >= 0 && n > max_n) continue;
        IntPoly phi = compute_cyclotomic(n);
        if (phi.degree() == deg && poly_equal(poly, phi)) return n;
    }
    return std::nullopt;
}

std::vector<ExprPtr> cyclotomic_roots(int n, const Symbol& /*var*/, AstArena& arena) {
    std::vector<ExprPtr> roots;
    // Primitive n-th roots of unity: exp(2*pi*i*k/n) where gcd(k, n) == 1
    auto gcd = [](int a, int b) {
        while (b) { a %= b; std::swap(a, b); }
        return a;
    };

    ExprPtr pi = arena.make<Constant>(MathConstant::Pi);
    ExprPtr i_const = arena.make<Constant>(MathConstant::I);
    ExprPtr two = arena.make<IntegerLit>(BigInt(2));
    ExprPtr n_lit = arena.make<IntegerLit>(BigInt(n));

    // Handle n=1 (phi(1)=1, root is 1)
    if (n == 1) {
        roots.push_back(arena.make<IntegerLit>(BigInt(1)));
        return roots;
    }
    // Handle n=2 (phi(2)=1, root is -1)
    if (n == 2) {
        roots.push_back(arena.make<IntegerLit>(BigInt(-1)));
        return roots;
    }

    for (int k = 1; k < n; ++k) {
        if (gcd(k, n) == 1) {
            // exp(2 * pi * i * k / n)
            ExprPtr k_lit = arena.make<IntegerLit>(BigInt(k));
            ExprPtr exponent = arena.make<Binary>(BinaryOp::Mul,
                arena.make<Binary>(BinaryOp::Mul, two, pi),
                arena.make<Binary>(BinaryOp::Mul, i_const,
                    arena.make<Binary>(BinaryOp::Div, k_lit, n_lit)));

            // We use exp function
            roots.push_back(arena.make<FuncCall>("exp", std::vector<ExprPtr>{exponent}));
        }
    }

    return roots;
}

} // namespace cas::algebra
