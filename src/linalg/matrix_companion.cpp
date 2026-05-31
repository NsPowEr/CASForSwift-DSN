// CAS-F4.2d — Companion matrix di un polinomio univariato.
//
// Per p(x) monico di grado n:
//   p(x) = x^n + c_{n-1} x^{n-1} + ... + c_1 x + c_0
// la companion matrix in forma di Frobenius è:
//   C =  [ 0  0  0 ... 0  -c_0     ]
//        [ 1  0  0 ... 0  -c_1     ]
//        [ 0  1  0 ... 0  -c_2     ]
//        [ . . . . . . . . . . . . ]
//        [ 0  0  0 ... 1  -c_{n-1} ]
// con la proprietà det(λI − C) = p(λ).  Usata per ridurre lo studio degli
// autovalori di polinomi qualsiasi a un problema di eigenvalues matriciali.
//
// Pre-condizione: `polynomial` univariato in `var` con grado ≥ 1.  Se non
// monico (lc ≠ 1), si normalizza dividendo tutti i coefficienti per lc:
// la companion del polinomio normalizzato ha le stesse radici di p(x).

#include "cas/linalg/Matrix.hpp"
#include "cas/linalg/matrix_expr_helpers.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/error.hpp"
#include "cas/error_helpers.hpp"
#include "cas/symbolic.hpp"

#include <cstddef>
#include <vector>

namespace cas::linalg {

Result<MatrixExpr> companion_matrix(ExprPtr polynomial, const Symbol& var,
                                      symbolic::CASContext& ctx) {
    if (!polynomial) {
        return fail<MatrixExpr>(CASError{
            CASErrorKind::InvalidArgument,
            "companion_matrix: null polynomial expression", std::nullopt});
    }

    auto coeffs_res = algebra::univariate_coefficients(polynomial, var, ctx);
    if (coeffs_res.is_error()) {
        return fail<MatrixExpr>(coeffs_res.error());
    }
    std::vector<ExprPtr>& coeffs = coeffs_res.value();
    // Rimuovi zeri leading per determinare n effettivo.
    while (coeffs.size() > 1 && is_zero_expr(coeffs.back())) coeffs.pop_back();

    if (coeffs.size() < 2U) {
        return fail<MatrixExpr>(CASError{
            CASErrorKind::InvalidArgument,
            "companion_matrix: polynomial degree must be >= 1", std::nullopt});
    }

    const std::size_t n = coeffs.size() - 1U;
    ExprPtr lc = coeffs[n];

    // Normalizza al monico se lc ≠ 1: dividi tutti i coefficienti per lc.
    if (!is_one_expr(lc)) {
        // Verifica lc ≠ 0 simbolicamente: il loop pop_back rimuove solo
        // zeri strutturali; un simbolo con assumption zero passerebbe.
        // Senza questo check, div_expr(c, lc) potrebbe produrre c/0 silente.
        if (!is_known_nonzero(lc, ctx)) {
            return make_unimplemented<MatrixExpr>(
                "linalg", "companion_matrix",
                "leading coefficient not known to be non-zero",
                error::reason_codes::LINALG_ZERO_PIVOT,
                "Add assumption that leading coefficient is non-zero, or normalize polynomial manually",
                "F4.2d");
        }
        for (std::size_t i = 0; i < n; ++i) {
            auto q_res = div_expr(ctx, coeffs[i], lc);
            if (q_res.is_error()) return fail<MatrixExpr>(q_res.error());

            auto t = algebra::together(q_res.value(), ctx);
            ExprPtr norm = t.is_ok() ? t.value() : q_res.value();
            coeffs[i] = norm;
        }
        coeffs[n] = integer(ctx, 1);
    }

    // Costruisci C in forma Frobenius: sub-diagonale a 1, ultima colonna = -c_i.
    MatrixExpr C(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) C(i, j) = integer(ctx, 0);
    }
    for (std::size_t i = 1; i < n; ++i) C(i, i - 1) = integer(ctx, 1);
    for (std::size_t i = 0; i < n; ++i) {
        if (is_zero_expr(coeffs[i])) {
            C(i, n - 1) = integer(ctx, 0);
        } else {
            auto neg_res = negate_expr(ctx, coeffs[i]);
            if (neg_res.is_error()) return fail<MatrixExpr>(neg_res.error());
            C(i, n - 1) = neg_res.value();
        }
    }

    return ok(std::move(C));
}

}  // namespace cas::linalg
