// F8.0-6.2 / Task 20 (BC-1..BC-3) — Branch-cut propagation helpers.
//
// Implements the unwinding-number correction terms required by
// Branch_Cut_Propagation.md when ctx.strict_branch_cuts() prevents the
// legacy positivity-assuming reductions.
//
// All corrections are built as raw AST nodes — no simplification is
// performed here. Downstream rules in the simplifier are free to fold
// (-1)^K(...) when K(...) is provably integer, or 2πi·K(z) when K(z) is
// provably zero (i.e. -π < Im(z) ≤ π).
//
// Reference: see header (Kahan 1987; Corless-Davenport-Jeffrey 2000).

#include "simplify_branch_cut.hpp"

#include <vector>

namespace cas::symbolic::branch_cut {

namespace {

[[nodiscard]] ExprPtr make_int(AstArena& arena, long long v) {
    return arena.make<IntegerLit>(BigInt(v));
}

[[nodiscard]] ExprPtr make_pi(AstArena& arena) {
    return arena.make<Constant>(MathConstant::Pi);
}

[[nodiscard]] ExprPtr make_imag_unit(AstArena& arena) {
    return arena.make<Constant>(MathConstant::I);
}

// 2·π·i — used in every correction.
[[nodiscard]] ExprPtr make_two_pi_i(AstArena& arena) {
    return arena.make<Product>(std::vector<ExprPtr>{
        make_int(arena, 2),
        make_pi(arena),
        make_imag_unit(arena),
    });
}

// ln(z) as a Ln FuncCall — left unsimplified so callers downstream can
// canonicalise according to whatever context they live in.
[[nodiscard]] ExprPtr make_ln(ExprPtr z, AstArena& arena) {
    return arena.make<FuncCall>(BuiltinOp::Ln, std::vector<ExprPtr>{z});
}

// K(z) as an UnwindingNumber FuncCall.
[[nodiscard]] ExprPtr make_K(ExprPtr z, AstArena& arena) {
    return arena.make<FuncCall>(BuiltinOp::UnwindingNumber,
        std::vector<ExprPtr>{z});
}

}  // namespace

ExprPtr make_pow_of_pow_correction(
    ExprPtr z,
    ExprPtr a,
    ExprPtr b,
    AstArena& arena) {

    // K(a · ln(z))
    ExprPtr a_ln_z = arena.make<Binary>(BinaryOp::Mul, a, make_ln(z, arena));
    ExprPtr K_inner = make_K(a_ln_z, arena);

    // exponent = 2πi · b · K(a · ln(z))
    ExprPtr exponent = arena.make<Product>(std::vector<ExprPtr>{
        make_two_pi_i(arena),
        b,
        K_inner,
    });

    // e^(exponent)
    return arena.make<FuncCall>(BuiltinOp::Exp,
        std::vector<ExprPtr>{exponent});
}

ExprPtr make_log_product_correction(
    ExprPtr z1,
    ExprPtr z2,
    AstArena& arena) {

    // ln(z1) + ln(z2)
    ExprPtr sum_logs = arena.make<Binary>(BinaryOp::Add,
        make_ln(z1, arena),
        make_ln(z2, arena));

    // K(ln(z1) + ln(z2))
    ExprPtr K = make_K(sum_logs, arena);

    // -2πi · K(...)
    ExprPtr neg_two_pi_i = arena.make<Unary>(UnaryOp::Neg,
        make_two_pi_i(arena));
    return arena.make<Binary>(BinaryOp::Mul, neg_two_pi_i, K);
}

ExprPtr make_log_quotient_correction(
    ExprPtr z1,
    ExprPtr z2,
    AstArena& arena) {

    // ln(z1) - ln(z2)
    ExprPtr diff_logs = arena.make<Binary>(BinaryOp::Sub,
        make_ln(z1, arena),
        make_ln(z2, arena));

    // K(ln(z1) - ln(z2))
    ExprPtr K = make_K(diff_logs, arena);

    // -2πi · K(...)
    ExprPtr neg_two_pi_i = arena.make<Unary>(UnaryOp::Neg,
        make_two_pi_i(arena));
    return arena.make<Binary>(BinaryOp::Mul, neg_two_pi_i, K);
}

ExprPtr make_sqrt_of_square_correction(
    ExprPtr z,
    AstArena& arena) {

    // 2·ln(z)
    ExprPtr two_ln_z = arena.make<Binary>(BinaryOp::Mul,
        make_int(arena, 2),
        make_ln(z, arena));

    // K(2·ln(z))
    ExprPtr K = make_K(two_ln_z, arena);

    // (-1)^K(2·ln(z))
    ExprPtr neg_one = make_int(arena, -1);
    return arena.make<Binary>(BinaryOp::Pow, neg_one, K);
}

}  // namespace cas::symbolic::branch_cut
