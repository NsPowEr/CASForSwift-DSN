#pragma once
// polynomial_groebner_f5.hpp
//
// Signature-based Groebner computation (F5C criterion subset).
//
// HONEST SCOPE:
//   Implements F5C — the F5 signature criterion + Rewritten criterion applied
//   inside a Buchberger/Sugar framework. Full matrix-F5 (parallel reduction
//   with global signature bookkeeping) is deferred as OPEN entry F3.3-F5-FULL.
//
// What is guaranteed:
//   - f5c_groebner() returns a correct reduced Groebner basis.
//   - F5 criterion provably eliminates ALL reductions-to-zero caused by
//     syzygies between the annotated polynomials (Faugère 2002 Thm 1).
//   - Zero-reduction count with F5C ≤ zero-reduction count without (verified
//     by direct comparison probe in tests).
//
// All arithmetic is exact (Rational/BigInt). No double/float.

#include "polynomial_groebner_f4.hpp"
#include "cas/result.hpp"

#include <vector>

namespace cas::algebra {

struct F5Result {
    std::vector<PolyF4> basis;
    std::size_t zero_reductions_f5{0};       // zero reductions performed by F5C
    std::size_t zero_reductions_baseline{0}; // corresponding count for comparison
};

struct BuchbergerCountResult {
    std::vector<PolyF4> basis;
    std::size_t zero_reductions{0};
};

// F5C Groebner computation:
// Applies F5 signature criterion + Rewritten criterion to eliminate
// provably-redundant S-pairs (syzygies) during Buchberger computation.
// Returns the reduced GRevLex Groebner basis + zero-reduction counts.
[[nodiscard]] F5Result f5c_groebner(
    std::vector<PolyF4> F,
    MonomialOrder order = MonomialOrder::GRevLex);

// Plain Buchberger with zero-reduction counting, for baseline comparison.
[[nodiscard]] BuchbergerCountResult buchberger_with_zero_count(
    std::vector<PolyF4> F,
    MonomialOrder order = MonomialOrder::GRevLex);

} // namespace cas::algebra
