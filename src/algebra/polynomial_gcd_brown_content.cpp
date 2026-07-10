// polynomial_gcd_brown_content.cpp — main-var polynomial-content pre-extraction
// for Brown's modular multivariate GCD (Geddes-Czapor-Labahn §7.4 setup).
//
// Split out of polynomial_gcd_brown_modular.cpp (anti-monolith): decompose
// ppP = cont_main(ppP) · pp_main(ppP), similarly ppQ, where cont_main is the
// (n-1)-variable gcd of the main_var-layer coefficients. Then
//   gcd(ppP, ppQ) = gcd(cont_main_P, cont_main_Q) · gcd(pp_main_P, pp_main_Q).
// The modular machinery handles only the primitive part w.r.t. main_var;
// without this pre-extraction, content factors common to every main_var layer
// (e.g. P = x·(yz+1)(w+...), gcd_content_z = x) get stripped by
// remove_spurious_main_var_factor leading to under-divided results.
//
// Reference: Geddes-Czapor-Labahn "Algorithms for Computer Algebra" §7.4.

#include "cas/algebra.hpp"
#include "cas/symbolic.hpp"
#include "algebra_internal.hpp"
#include "polynomial_gcd_fp_internal.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <utility>
#include <vector>

namespace cas::algebra {

using fp_helpers::BMonomial;
using fp_helpers::BSparsePoly;

MainVarContentGcd extract_main_var_content_gcd(
    BSparsePoly& ppP, BSparsePoly& ppQ,
    const std::vector<Symbol>& vars, std::size_t n_vars, std::size_t main_var,
    symbolic::CASContext& ctx) {

    MainVarContentGcd result;
    result.content_gcd =
        MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}});
    if (n_vars < 2U) return result;

    auto layers_main = [&](const BSparsePoly& sp) {
        std::map<std::size_t, BSparsePoly> out;
        for (const auto& [m, c] : sp) {
            std::size_t k = (main_var < m.size()) ? m[main_var] : 0U;
            BMonomial m2 = m;
            if (main_var < m2.size()) m2[main_var] = 0U;
            out[k][m2] = c;
        }
        return out;
    };
    auto content_main_z = [&](const BSparsePoly& sp) -> MultivariatePolynomial {
        auto layers = layers_main(sp);
        if (layers.size() <= 1U) {
            return MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}});
        }
        MultivariatePolynomial acc;
        bool first = true;
        for (auto& [k, layer] : layers) {
            if (layer.empty()) continue;
            MultivariatePolynomial cur = sub_sparse_to_mv(layer, vars, main_var);
            if (first) { acc = std::move(cur); first = false; continue; }
            auto rec = gcd_brown_modular(acc, cur, ctx);
            if (rec.is_error()) {
                return MultivariatePolynomial({{.coefficient = BigInt(1), .factors = {}}});
            }
            acc = std::move(rec.value());
            if (acc.terms().size() == 1U && acc.terms()[0].factors.empty()
                && acc.terms()[0].coefficient == BigInt(1)) {
                return acc;
            }
        }
        return acc;
    };

    const MultivariatePolynomial cont_main_P_mv = content_main_z(ppP);
    const MultivariatePolynomial cont_main_Q_mv = content_main_z(ppQ);
    auto is_one = [](const MultivariatePolynomial& m){
        return m.terms().size() == 1U && m.terms()[0].factors.empty()
            && m.terms()[0].coefficient == BigInt(1);
    };
    if (is_one(cont_main_P_mv) && is_one(cont_main_Q_mv)) return result;

    auto cg = gcd_brown_modular(cont_main_P_mv, cont_main_Q_mv, ctx);
    if (cg.is_error() || is_one(cg.value())) return result;

    MultivariatePolynomial content_gcd = std::move(cg.value());
    BSparsePoly cont_sp = mv_to_sub_sparse(content_gcd, vars, main_var);
    if (cont_sp.empty()) return result;

    auto strip = [&](const BSparsePoly& X) -> std::optional<BSparsePoly> {
        auto layers = layers_main(X);
        std::map<std::size_t, BSparsePoly> out_layers;
        for (const auto& [k, layer] : layers) {
            if (layer.empty()) { out_layers[k] = {}; continue; }
            BSparsePoly quo;
            if (!exact_divide_sparse_z(layer, cont_sp, n_vars, quo))
                return std::nullopt;
            out_layers[k] = std::move(quo);
        }
        BSparsePoly res;
        for (const auto& [k, layer] : out_layers) {
            for (const auto& [m, c] : layer) {
                BMonomial nm = m;
                if (nm.size() <= main_var) nm.resize(main_var + 1U, 0U);
                nm[main_var] = static_cast<unsigned int>(k);
                res[nm] = c;
            }
        }
        return res;
    };
    auto sp_pp_P = strip(ppP);
    auto sp_pp_Q = strip(ppQ);
    if (sp_pp_P.has_value() && sp_pp_Q.has_value()) {
        ppP = std::move(*sp_pp_P);
        ppQ = std::move(*sp_pp_Q);
        result.content_gcd = std::move(content_gcd);
        result.stripped = true;
    }
    return result;
}

}  // namespace cas::algebra
