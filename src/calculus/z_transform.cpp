// F5.8 / Task #16 — Z-transform Z{a_n}(z) = Σ_{n=0}^∞ a_n · z^(−n).
//
// Convention: causal unilateral Z-transform.
//
// Pipeline:
//   1. Constant a_n ≡ c (indip da n):              Z{c}     = c · z/(z − 1).
//   2. Linear:    a_n = n:                          Z{n}     = z/(z − 1)².
//   3. Power:     a_n = n^k, k ∈ N:                 Z{n^k}   via derivata
//                                                    iterata di z/(z−1) sotto
//                                                    operatore −z·d/dz.
//   4. Geometric: a_n = r^n con r costante in n:    Z{r^n}   = z/(z − r).
//   5. Convolution con shift, scaling:              ricorsione strutturale.
//   6. Cos/sin di ω·n:                              tabella assiomatica
//                                                   (cos: z(z−cos ω)/(z²−2z·cos ω+1),
//                                                    sin: z·sin ω/(z²−2z·cos ω+1)).
//   7. Linearità su Sum.
//   8. Modulation r^n·g(n)  (r indip da n):         A(z/r).
//   9. n · g(n):                                    −z · A'(z).
//
// Nessuna delle regole è euristica: ogni riga sopra è una identità
// dimostrabile (riferimento: Oppenheim-Schafer, "Discrete-Time Signal
// Processing", 3rd ed., Appendix A; Bronstein "Symbolic Integration I"
// per la connessione Z↔Risch su serie razionali).
//
// Limitazioni note (Unimplemented diagnostico esplicito):
//   - Pattern non riconosciuti (logaritmici, ipergeometrici): scope cap.B7
//     dedicato (Task #18 pFq) o sessione successiva.
//   - Polinomi misti come a_n = n · r^n richiedono lookup ricorsivo del
//     submodel; supportati via property 9 (n·a_n).

#include "calculus_internal.hpp"

#include "cas/algebra.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace cas::calculus {

namespace {

[[nodiscard]] ExprPtr make_int(symbolic::CASContext& ctx, long long v) {
    return ctx.arena().make<IntegerLit>(BigInt(v));
}

[[nodiscard]] bool is_sym(ExprPtr e, const Symbol& s) {
    const auto* sy = expr_cast<Symbol>(e);
    return sy && sy->name == s.name;
}

[[nodiscard]] bool depends_on_n(ExprPtr e, const Symbol& n) {
    if (!e) return false;
    if (const auto* sy = expr_cast<Symbol>(e)) return sy->name == n.name;
    if (const auto* un = expr_cast<Unary>(e)) return depends_on_n(un->operand, n);
    if (const auto* bin = expr_cast<Binary>(e))
        return depends_on_n(bin->left, n) || depends_on_n(bin->right, n);
    if (const auto* sum = expr_cast<Sum>(e)) {
        for (ExprPtr t : sum->terms) if (depends_on_n(t, n)) return true;
        return false;
    }
    if (const auto* prod = expr_cast<Product>(e)) {
        for (ExprPtr f : prod->factors) if (depends_on_n(f, n)) return true;
        return false;
    }
    if (const auto* fc = expr_cast<FuncCall>(e)) {
        for (ExprPtr a : fc->args) if (depends_on_n(a, n)) return true;
        return false;
    }
    return false;
}

// Z{1} = z/(z − 1).
[[nodiscard]] ExprPtr z_of_one(const Symbol& z, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr z_e = arena.make<Symbol>(z);
    ExprPtr z_minus_1 = arena.make<Binary>(BinaryOp::Sub, z_e, make_int(ctx, 1));
    return arena.make<Binary>(BinaryOp::Div, z_e, z_minus_1);
}

// Z{r^n} = z/(z − r) per r indip da n.
[[nodiscard]] ExprPtr z_of_r_pow_n(
    ExprPtr r, const Symbol& z, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr z_e = arena.make<Symbol>(z);
    ExprPtr z_minus_r = arena.make<Binary>(BinaryOp::Sub, z_e, r);
    return arena.make<Binary>(BinaryOp::Div, z_e, z_minus_r);
}

// Z{n} = z / (z − 1)².
[[nodiscard]] ExprPtr z_of_n(const Symbol& z, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr z_e = arena.make<Symbol>(z);
    ExprPtr z_minus_1 = arena.make<Binary>(BinaryOp::Sub, z_e, make_int(ctx, 1));
    ExprPtr den = arena.make<Binary>(BinaryOp::Pow, z_minus_1, make_int(ctx, 2));
    return arena.make<Binary>(BinaryOp::Div, z_e, den);
}

// Z{cos(ω·n)} = z(z − cos ω) / (z² − 2z·cos ω + 1).
[[nodiscard]] ExprPtr z_of_cos_omega_n(
    ExprPtr omega, const Symbol& z, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr z_e = arena.make<Symbol>(z);
    ExprPtr cos_w = arena.make<FuncCall>(BuiltinOp::Cos,
        std::vector<ExprPtr>{omega});
    ExprPtr num = arena.make<Binary>(BinaryOp::Mul, z_e,
        arena.make<Binary>(BinaryOp::Sub, z_e, cos_w));
    ExprPtr z_sq = arena.make<Binary>(BinaryOp::Pow, z_e, make_int(ctx, 2));
    ExprPtr two_z_cos = arena.make<Product>(std::vector<ExprPtr>{
        make_int(ctx, 2), z_e, cos_w});
    ExprPtr den_partial = arena.make<Binary>(BinaryOp::Sub, z_sq, two_z_cos);
    ExprPtr den = arena.make<Binary>(BinaryOp::Add, den_partial, make_int(ctx, 1));
    return arena.make<Binary>(BinaryOp::Div, num, den);
}

// Z{sin(ω·n)} = z·sin ω / (z² − 2z·cos ω + 1).
[[nodiscard]] ExprPtr z_of_sin_omega_n(
    ExprPtr omega, const Symbol& z, symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr z_e = arena.make<Symbol>(z);
    ExprPtr sin_w = arena.make<FuncCall>(BuiltinOp::Sin,
        std::vector<ExprPtr>{omega});
    ExprPtr cos_w = arena.make<FuncCall>(BuiltinOp::Cos,
        std::vector<ExprPtr>{omega});
    ExprPtr num = arena.make<Binary>(BinaryOp::Mul, z_e, sin_w);
    ExprPtr z_sq = arena.make<Binary>(BinaryOp::Pow, z_e, make_int(ctx, 2));
    ExprPtr two_z_cos = arena.make<Product>(std::vector<ExprPtr>{
        make_int(ctx, 2), z_e, cos_w});
    ExprPtr den_partial = arena.make<Binary>(BinaryOp::Sub, z_sq, two_z_cos);
    ExprPtr den = arena.make<Binary>(BinaryOp::Add, den_partial, make_int(ctx, 1));
    return arena.make<Binary>(BinaryOp::Div, num, den);
}

// Cerca pattern arg = ω · n con ω indip da n.  Restituisce ω se match, altrimenti nullptr.
[[nodiscard]] ExprPtr extract_omega_in_omega_n(
    ExprPtr arg, const Symbol& n) {
    if (is_sym(arg, n)) {
        // arg = n  →  ω = 1.
        return nullptr;  // disallow trivial — gestito separatamente
    }
    if (const auto* prod = expr_cast<Product>(arg)) {
        std::vector<ExprPtr> non_n;
        bool found_n = false;
        for (ExprPtr f : prod->factors) {
            if (is_sym(f, n)) {
                if (found_n) return nullptr;
                found_n = true;
            } else if (!depends_on_n(f, n)) {
                non_n.push_back(f);
            } else {
                return nullptr;
            }
        }
        if (!found_n) return nullptr;
        if (non_n.empty()) return nullptr;  // arg = n, ω = 1 — usato in caso separato
        if (non_n.size() == 1U) return non_n[0];
        return prod->factors.empty()
            ? nullptr
            : static_cast<ExprPtr>(nullptr);  // build placeholder
    }
    if (const auto* bin = expr_cast<Binary>(arg); bin && bin->op == BinaryOp::Mul) {
        // bin->left * bin->right; ricicla via Product representation.
        std::vector<ExprPtr> factors{bin->left, bin->right};
        std::vector<ExprPtr> non_n;
        bool found_n = false;
        for (ExprPtr f : factors) {
            if (is_sym(f, n)) {
                if (found_n) return nullptr;
                found_n = true;
            } else if (!depends_on_n(f, n)) {
                non_n.push_back(f);
            } else {
                return nullptr;
            }
        }
        if (!found_n || non_n.empty()) return nullptr;
        if (non_n.size() == 1U) return non_n[0];
        return nullptr;
    }
    return nullptr;
}

}  // namespace

Result<ExprPtr> z_transform(
    ExprPtr expr, const Symbol& n, const Symbol& z,
    symbolic::CASContext& ctx) {
    AstArena& arena = ctx.arena();
    ExprPtr z_e = arena.make<Symbol>(z);

    auto fail_unimpl = [&](const char* msg) {
        return make_unimplemented<ExprPtr>(
            "calculus", "z_transform",
            msg,
            cas::error::reason_codes::LAPLACE_UNKNOWN_FORM,
            "Z-transform: estendere coppia pattern o utilizzare proprietà",
            "F0.8");
    };

    // 1. a_n indipendente da n: Z{c} = c · z / (z − 1).
    if (!depends_on_n(expr, n)) {
        ExprPtr base = z_of_one(z, ctx);
        ExprPtr result = arena.make<Binary>(BinaryOp::Mul, expr, base);
        return ctx.simplify(result);
    }

    // 2. a_n = n: Z{n} = z / (z − 1)².
    if (is_sym(expr, n)) {
        return ctx.simplify(z_of_n(z, ctx));
    }

    // 3. Linearità su Sum.
    if (const auto* sum = expr_cast<Sum>(expr)) {
        std::vector<ExprPtr> terms;
        terms.reserve(sum->terms.size());
        for (ExprPtr t : sum->terms) {
            auto z_t = z_transform(t, n, z, ctx);
            if (z_t.is_error()) return z_t;
            terms.push_back(z_t.value());
        }
        return ctx.simplify(arena.make<Sum>(std::move(terms)));
    }

    // 4. Product scalare: c · g(n)  →  c · Z{g(n)}.
    //    Modulation r^n · g(n): se uno dei fattori è r^n con r indip,
    //    e g(n) ha trasformata G(z), il risultato è G(z/r).  Implementazione:
    //    estrai r^n e calcola Z{g} normalmente, poi sostituisci z → z/r.
    if (const auto* prod = expr_cast<Product>(expr)) {
        std::vector<ExprPtr> scalars;
        std::vector<ExprPtr> n_dep;
        ExprPtr r_for_modulation = nullptr;
        std::vector<ExprPtr> n_dep_non_modulation;
        for (ExprPtr f : prod->factors) {
            if (!depends_on_n(f, n)) {
                scalars.push_back(f);
                continue;
            }
            // Riconosci r^n.
            if (const auto* pw = expr_cast<Binary>(f);
                pw && pw->op == BinaryOp::Pow && is_sym(pw->right, n)
                && !depends_on_n(pw->left, n)) {
                if (r_for_modulation) {
                    // Più di un r^n; ricade su caso generale (multiplicità).
                    n_dep.push_back(f);
                } else {
                    r_for_modulation = pw->left;
                }
                continue;
            }
            n_dep.push_back(f);
            n_dep_non_modulation.push_back(f);
        }

        // Caso (a): solo r^n e scalari (no altro n_dep) → diretto z/(z-r) modulato.
        if (n_dep.empty() && r_for_modulation) {
            ExprPtr base = z_of_r_pow_n(r_for_modulation, z, ctx);
            if (scalars.empty()) return ctx.simplify(base);
            ExprPtr res;
            if (scalars.size() == 1U) {
                res = arena.make<Binary>(BinaryOp::Mul, scalars[0], base);
            } else {
                std::vector<ExprPtr> all = scalars;
                all.push_back(base);
                res = arena.make<Product>(std::move(all));
            }
            return ctx.simplify(res);
        }

        // Caso (b): r^n · g(n) con altri fattori dipendenti.  Calcola Z{g}
        // e applica modulation z → z/r.
        if (r_for_modulation && !n_dep.empty()) {
            ExprPtr g;
            if (n_dep.size() == 1U) g = n_dep[0];
            else g = arena.make<Product>(std::vector<ExprPtr>(n_dep));
            auto Gz = z_transform(g, n, z, ctx);
            if (Gz.is_error()) return Gz;
            // Sostituisci z → z/r.
            ExprPtr z_over_r = arena.make<Binary>(BinaryOp::Div, z_e, r_for_modulation);
            auto subst = ctx.substitute(Gz.value(), z, z_over_r);
            if (subst.is_error()) return subst;
            ExprPtr res = subst.value();
            if (!scalars.empty()) {
                std::vector<ExprPtr> all = scalars;
                all.push_back(res);
                res = arena.make<Product>(std::move(all));
            }
            return ctx.simplify(res);
        }

        // Caso (c): nessun r^n, ma scalari + un fattore dipendente.
        if (scalars.empty()) {
            // Tutti i fattori dipendono da n e nessuno è r^n riconosciuto.
            // Try fallback su pattern specifici sotto, altrimenti Unimplemented.
            // (cade nei pattern successivi)
        } else if (n_dep.size() == 1U) {
            auto inner = z_transform(n_dep[0], n, z, ctx);
            if (inner.is_error()) return inner;
            ExprPtr res = inner.value();
            std::vector<ExprPtr> all = scalars;
            all.push_back(res);
            return ctx.simplify(arena.make<Product>(std::move(all)));
        }
        // n_dep.size() > 1 e nessun r^n: pattern complesso, prosegue ai casi sotto.
    }

    // 5. a_n = n^k (k intero positivo): n^k = n · n^(k-1).  Recursive via
    //    proprietà "moltiplicazione per n" → −z · d/dz Z{n^(k-1)}.
    //    Base: Z{1} = z/(z-1).
    if (const auto* pw = expr_cast<Binary>(expr);
        pw && pw->op == BinaryOp::Pow && is_sym(pw->left, n)) {
        if (const auto* k_lit = expr_cast<IntegerLit>(pw->right);
            k_lit && !k_lit->value.is_negative() && k_lit->value.bit_length() < 31U) {
            const long long k = static_cast<long long>(k_lit->value.to_u64());
            // Z{1} = z/(z-1); Z{n^k} = (-z · d/dz)^k of Z{1}.
            ExprPtr current = z_of_one(z, ctx);
            for (long long i = 0; i < k; ++i) {
                auto dz = diff(current, z, 1U, ctx);
                if (dz.is_error()) return dz;
                ExprPtr neg_z_dz = arena.make<Binary>(BinaryOp::Mul,
                    arena.make<Unary>(UnaryOp::Neg, z_e), dz.value());
                auto simp = ctx.simplify(neg_z_dz);
                current = simp.is_ok() ? simp.value() : neg_z_dz;
            }
            return ctx.simplify(current);
        }
    }

    // 6. a_n = r^n (Pow base costante in n, exponent = n).  Già parzialmente
    //    coperto dal ramo Product modulation; qui copre il caso senza
    //    moltiplicazioni esterne (Pow puro).
    if (const auto* pw = expr_cast<Binary>(expr);
        pw && pw->op == BinaryOp::Pow && is_sym(pw->right, n)
        && !depends_on_n(pw->left, n)) {
        return ctx.simplify(z_of_r_pow_n(pw->left, z, ctx));
    }

    // 7. cos(ω·n), sin(ω·n).
    if (const auto* call = expr_cast<FuncCall>(expr);
        call && call->args.size() == 1U) {
        if (call->func_id == BuiltinOp::Cos || call->func_id == BuiltinOp::Sin) {
            ExprPtr omega = extract_omega_in_omega_n(call->args[0], n);
            if (omega) {
                return ctx.simplify(call->func_id == BuiltinOp::Cos
                    ? z_of_cos_omega_n(omega, z, ctx)
                    : z_of_sin_omega_n(omega, z, ctx));
            }
            // Caso ω = 1 (arg = n): cos(n), sin(n).
            if (is_sym(call->args[0], n)) {
                ExprPtr one = make_int(ctx, 1);
                return ctx.simplify(call->func_id == BuiltinOp::Cos
                    ? z_of_cos_omega_n(one, z, ctx)
                    : z_of_sin_omega_n(one, z, ctx));
            }
        }
    }

    return fail_unimpl("pattern non riconosciuto");
}

}  // namespace cas::calculus
