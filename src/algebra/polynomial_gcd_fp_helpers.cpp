#include "polynomial_gcd_fp_internal.hpp"
#include "cas/numtheory.hpp"

#include <algorithm>
#include <set>

namespace cas::algebra::fp_helpers {

BigInt pos_mod(const BigInt& a, const BigInt& m) {
    BigInt r = a % m;
    if (r.is_negative()) r += m.abs();
    return r;
}

BigInt bigint_gcd(BigInt a, BigInt b) {
    a = a.abs(); b = b.abs();
    while (!b.is_zero()) { BigInt t = a % b; a = std::move(b); b = std::move(t); }
    return a;
}

BigInt centered_repr(const BigInt& r, const BigInt& M) {
    BigInt two_r = r + r;
    return (two_r > M) ? (r - M) : r;
}

BSparsePoly reduce_sparse_mod_p(const BSparsePoly& sp, const BigInt& p) {
    BSparsePoly r;
    for (const auto& [mono, c] : sp) {
        BigInt cr = pos_mod(c, p);
        if (!cr.is_zero()) r[mono] = cr;
    }
    return r;
}

BigInt sparse_inf_norm(const BSparsePoly& sp) {
    BigInt best(0);
    for (const auto& [_, c] : sp) {
        BigInt ac = c.abs();
        if (ac > best) best = ac;
    }
    return best;
}

std::size_t deg_in_var(const BSparsePoly& sp, std::size_t vi) {
    std::size_t d = 0U;
    for (const auto& [m, _] : sp)
        if (vi < m.size()) d = std::max<std::size_t>(d, m[vi]);
    return d;
}

BSparsePoly eval_var_mod_p(const BSparsePoly& sp, std::size_t var_idx,
                            const BigInt& val, const BigInt& p) {
    BSparsePoly out;
    for (const auto& [mono, c] : sp) {
        BigInt new_c = c;
        if (var_idx < mono.size() && mono[var_idx] > 0U) {
            BigInt acc(1);
            for (unsigned int i = 0; i < mono[var_idx]; ++i)
                acc = pos_mod(acc * val, p);
            new_c = pos_mod(new_c * acc, p);
        }
        if (new_c.is_zero()) continue;
        BMonomial m2 = mono;
        if (var_idx < m2.size()) m2[var_idx] = 0U;
        out[m2] += new_c;
        out[m2] = pos_mod(out[m2], p);
        if (out[m2].is_zero()) out.erase(m2);
    }
    return out;
}

std::optional<BSparsePoly> univariate_sparse_gcd_fp(
    const BSparsePoly& A, const BSparsePoly& B,
    std::size_t var_idx, const BigInt& p) {
    auto to_dense = [&](const BSparsePoly& sp) -> std::vector<BigInt> {
        std::size_t d = deg_in_var(sp, var_idx);
        std::vector<BigInt> v(d + 1U, BigInt(0));
        for (const auto& [mono, c] : sp) {
            std::size_t k = (var_idx < mono.size()) ? mono[var_idx] : 0U;
            v[k] = pos_mod(v[k] + c, p);
        }
        while (!v.empty() && v.back().is_zero()) v.pop_back();
        return v;
    };
    auto a = to_dense(A);
    auto b = to_dense(B);
    while (!b.empty()) {
        BigInt inv_lc = pos_mod(b.back(), p);
        auto inv = numtheory::modular_inverse(inv_lc, p);
        if (inv.is_error()) return std::nullopt;
        BigInt ilc = inv.value();
        std::vector<BigInt> rem = a;
        while (!rem.empty() && rem.size() >= b.size()) {
            std::size_t dd = rem.size() - b.size();
            BigInt factor = pos_mod(rem.back() * ilc, p);
            for (std::size_t i = 0; i < b.size(); ++i)
                rem[i + dd] = pos_mod(rem[i + dd] - factor * b[i], p);
            while (!rem.empty() && rem.back().is_zero()) rem.pop_back();
        }
        a = std::move(b);
        b = std::move(rem);
    }
    if (a.empty()) return BSparsePoly{};
    auto inv = numtheory::modular_inverse(pos_mod(a.back(), p), p);
    if (inv.is_error()) return std::nullopt;
    BigInt ilc = inv.value();
    BSparsePoly out;
    for (std::size_t k = 0; k < a.size(); ++k) {
        BigInt c = pos_mod(a[k] * ilc, p);
        if (c.is_zero()) continue;
        BMonomial m;
        m.assign(var_idx + 1U, 0U);
        m[var_idx] = static_cast<unsigned int>(k);
        out[m] = c;
    }
    return out;
}

std::optional<BSparsePoly> lagrange_interp_fp(
    const std::vector<BigInt>& bs, const std::vector<BSparsePoly>& vals,
    std::size_t var_idx, std::size_t target_deg, const BigInt& p) {
    const std::size_t n = bs.size();
    if (n != target_deg + 1U || n == 0U) return std::nullopt;
    std::set<BMonomial> all_monos;
    for (const auto& v : vals)
        for (const auto& [m, _] : v) all_monos.insert(m);
    BSparsePoly result;
    std::vector<std::vector<BigInt>> Li_polys(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<BigInt> num = { BigInt(1) };
        BigInt den(1);
        for (std::size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            std::vector<BigInt> nn(num.size() + 1U, BigInt(0));
            for (std::size_t k = 0; k < num.size(); ++k) {
                nn[k + 1U] = pos_mod(nn[k + 1U] + num[k], p);
                nn[k]      = pos_mod(nn[k]     - num[k] * bs[j], p);
            }
            num = std::move(nn);
            den = pos_mod(den * pos_mod(bs[i] - bs[j], p), p);
        }
        auto inv_d = numtheory::modular_inverse(den, p);
        if (inv_d.is_error()) return std::nullopt;
        for (auto& c : num) c = pos_mod(c * inv_d.value(), p);
        Li_polys[i] = std::move(num);
    }
    for (const auto& m : all_monos) {
        for (std::size_t i = 0; i < n; ++i) {
            auto it = vals[i].find(m);
            if (it == vals[i].end() || it->second.is_zero()) continue;
            const BigInt& vi = it->second;
            for (std::size_t k = 0; k < Li_polys[i].size(); ++k) {
                if (Li_polys[i][k].is_zero()) continue;
                BMonomial m2 = m;
                if (m2.size() <= var_idx) m2.resize(var_idx + 1U, 0U);
                m2[var_idx] = static_cast<unsigned int>(k);
                result[m2] = pos_mod(result[m2] + vi * Li_polys[i][k], p);
                if (result[m2].is_zero()) result.erase(m2);
            }
        }
    }
    return result;
}

std::map<std::size_t, BSparsePoly> layers_by_var(
    const BSparsePoly& sp, std::size_t var_idx) {
    std::map<std::size_t, BSparsePoly> out;
    for (const auto& [m, c] : sp) {
        std::size_t k = (var_idx < m.size()) ? m[var_idx] : 0U;
        BMonomial m2 = m;
        if (var_idx < m2.size()) m2[var_idx] = 0U;
        out[k][m2] = c;
    }
    return out;
}

BSparsePoly reassemble_layers(
    const std::map<std::size_t, BSparsePoly>& layers,
    std::size_t var_idx, const BigInt& p) {
    BSparsePoly out;
    for (const auto& [k, layer] : layers) {
        for (const auto& [m, c] : layer) {
            if (c.is_zero()) continue;
            BMonomial nm = m;
            if (nm.size() <= var_idx) nm.resize(var_idx + 1U, 0U);
            nm[var_idx] = static_cast<unsigned int>(k);
            out[nm] = pos_mod(c, p);
            if (out[nm].is_zero()) out.erase(nm);
        }
    }
    return out;
}

BSparsePoly mul_mod_p(const BSparsePoly& A, const BSparsePoly& B,
                      const BigInt& p) {
    BSparsePoly out;
    if (A.empty() || B.empty()) return out;
    std::size_t n = 0;
    for (const auto& [m, _] : A) n = std::max(n, m.size());
    for (const auto& [m, _] : B) n = std::max(n, m.size());
    for (const auto& [ma, ca] : A) {
        for (const auto& [mb, cb] : B) {
            BMonomial nm(n, 0U);
            for (std::size_t i = 0; i < n; ++i) {
                unsigned int ea = (i < ma.size()) ? ma[i] : 0U;
                unsigned int eb = (i < mb.size()) ? mb[i] : 0U;
                nm[i] = ea + eb;
            }
            BigInt v = pos_mod(out[nm] + ca * cb, p);
            if (v.is_zero()) out.erase(nm);
            else out[nm] = v;
        }
    }
    return out;
}

std::optional<BSparsePoly> exact_div_fp(
    const BSparsePoly& A, const BSparsePoly& B, const BigInt& p) {
    if (B.empty()) return std::nullopt;
    if (A.empty()) return BSparsePoly{};
    std::size_t n = 0;
    for (const auto& [m, _] : A) n = std::max(n, m.size());
    for (const auto& [m, _] : B) n = std::max(n, m.size());
    BSparsePoly rem = A;
    BSparsePoly quo;
    auto [dlm, dlc] = *std::prev(B.end());
    auto inv = numtheory::modular_inverse(pos_mod(dlc, p), p);
    if (inv.is_error()) return std::nullopt;
    const BigInt ilc = inv.value();
    const std::size_t budget = (rem.size() + 1U) * (B.size() + 1U) + 16U;
    std::size_t steps = 0;
    while (!rem.empty()) {
        if (++steps > budget) return std::nullopt;
        auto [rlm, rlc] = *std::prev(rem.end());
        BMonomial qm(n, 0U);
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t ri = (i < rlm.size()) ? rlm[i] : 0U;
            std::size_t di = (i < dlm.size()) ? dlm[i] : 0U;
            if (ri < di) return std::nullopt;
            qm[i] = static_cast<unsigned int>(ri - di);
        }
        BigInt qc = pos_mod(rlc * ilc, p);
        quo[qm] = qc;
        for (const auto& [bm, bc] : B) {
            BMonomial nm(n, 0U);
            for (std::size_t i = 0; i < n; ++i)
                nm[i] = qm[i] + ((i < bm.size()) ? bm[i] : 0U);
            BigInt v = pos_mod(rem[nm] - pos_mod(qc * bc, p), p);
            if (v.is_zero()) rem.erase(nm);
            else rem[nm] = v;
        }
    }
    return quo;
}

}  // namespace cas::algebra::fp_helpers
