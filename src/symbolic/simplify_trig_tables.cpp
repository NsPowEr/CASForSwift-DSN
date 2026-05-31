#include "simplify_trig_tables_impl.hpp"

// F1.4 — Fermat-prime base angle constructors and exact sin/cos tables.
//
// Provides the exact radical forms for the "primitive" constructible angles
// at p = 3 (already given by 1/2 and √3/2) and p = 5 (the pentagon).
// These are the building blocks from which cos_ref_value_half_angle and
// try_angle_combination derive all other constructible angles.
//
// Mathematical basis (Gauss, Disquisitiones Arithmeticae §VII):
//   For the regular pentagon (q=5), the exact values follow from:
//     cos(π/5) = (1+√5)/4,   sin(π/5) = √(10−2√5)/4
//     cos(2π/5) = (√5−1)/4,  sin(2π/5) = √(10+2√5)/4
//   These are derivable via the quadratic period equation for the 5th cyclotomic
//   field (the "Gaussian period" ω + ω⁴ = (−1+√5)/2 where ω = e^{2πi/5}).

namespace cas::symbolic::detail {

ExprPtr cos_pi_over_5(AstArena& arena) {
    ExprPtr sqrt5 = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_integer(arena, BigInt(5))});
    ExprPtr sum = arena.make<Sum>(std::vector<ExprPtr>{make_integer(arena, BigInt(1)), sqrt5});
    return arena.make<Binary>(BinaryOp::Div, sum, make_integer(arena, BigInt(4)));
}

ExprPtr cos_2pi_over_5(AstArena& arena) {
    ExprPtr sqrt5 = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_integer(arena, BigInt(5))});
    ExprPtr neg_one = arena.make<Unary>(UnaryOp::Neg, make_integer(arena, BigInt(1)));
    ExprPtr sum = arena.make<Sum>(std::vector<ExprPtr>{sqrt5, neg_one});
    return arena.make<Binary>(BinaryOp::Div, sum, make_integer(arena, BigInt(4)));
}

ExprPtr sin_pi_over_5(AstArena& arena) {
    ExprPtr sqrt5 = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_integer(arena, BigInt(5))});
    ExprPtr two_sqrt5 = arena.make<Product>(std::vector<ExprPtr>{
        make_integer(arena, BigInt(2)), sqrt5});
    ExprPtr neg_two_sqrt5 = arena.make<Unary>(UnaryOp::Neg, two_sqrt5);
    ExprPtr inner_sum = arena.make<Sum>(std::vector<ExprPtr>{
        make_integer(arena, BigInt(10)), neg_two_sqrt5});
    ExprPtr outer = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{inner_sum});
    return arena.make<Binary>(BinaryOp::Div, outer, make_integer(arena, BigInt(4)));
}

ExprPtr sin_2pi_over_5(AstArena& arena) {
    ExprPtr sqrt5 = arena.make<FuncCall>(BuiltinOp::Sqrt,
        std::vector<ExprPtr>{make_integer(arena, BigInt(5))});
    ExprPtr two_sqrt5 = arena.make<Product>(std::vector<ExprPtr>{
        make_integer(arena, BigInt(2)), sqrt5});
    ExprPtr inner_sum = arena.make<Sum>(std::vector<ExprPtr>{
        make_integer(arena, BigInt(10)), two_sqrt5});
    ExprPtr outer = arena.make<FuncCall>(BuiltinOp::Sqrt, std::vector<ExprPtr>{inner_sum});
    return arena.make<Binary>(BinaryOp::Div, outer, make_integer(arena, BigInt(4)));
}

ExprPtr sin_ref_value_table(Rational ref, AstArena& arena) {
    const Rational zero(BigInt(0));
    const Rational one_tenth(BigInt(1), BigInt(10));
    const Rational one_sixth(BigInt(1), BigInt(6));
    const Rational one_fifth(BigInt(1), BigInt(5));
    const Rational one_quarter(BigInt(1), BigInt(4));
    const Rational three_tenths(BigInt(3), BigInt(10));
    const Rational one_third(BigInt(1), BigInt(3));
    const Rational two_fifths(BigInt(2), BigInt(5));
    const Rational one_half(BigInt(1), BigInt(2));

    if (ref == zero)       return make_integer(arena, BigInt(0));
    if (ref == one_tenth)  return cos_2pi_over_5(arena);
    if (ref == one_sixth)  return make_rational(arena, Rational(BigInt(1), BigInt(2)));
    if (ref == one_fifth)  return sin_pi_over_5(arena);
    if (ref == one_quarter)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{make_integer(arena, BigInt(2))}),
            make_integer(arena, BigInt(2)));
    if (ref == three_tenths) return cos_pi_over_5(arena);
    if (ref == one_third)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{make_integer(arena, BigInt(3))}),
            make_integer(arena, BigInt(2)));
    if (ref == two_fifths) return sin_2pi_over_5(arena);
    if (ref == one_half)   return make_integer(arena, BigInt(1));
    return nullptr;
}

ExprPtr cos_ref_value_table(Rational ref, AstArena& arena) {
    const Rational zero(BigInt(0));
    const Rational one_tenth(BigInt(1), BigInt(10));
    const Rational one_sixth(BigInt(1), BigInt(6));
    const Rational one_fifth(BigInt(1), BigInt(5));
    const Rational one_quarter(BigInt(1), BigInt(4));
    const Rational three_tenths(BigInt(3), BigInt(10));
    const Rational one_third(BigInt(1), BigInt(3));
    const Rational two_fifths(BigInt(2), BigInt(5));
    const Rational one_half(BigInt(1), BigInt(2));

    if (ref == zero)       return make_integer(arena, BigInt(1));
    if (ref == one_tenth)  return sin_2pi_over_5(arena);
    if (ref == one_sixth)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{make_integer(arena, BigInt(3))}),
            make_integer(arena, BigInt(2)));
    if (ref == one_fifth)  return cos_pi_over_5(arena);
    if (ref == one_quarter)
        return arena.make<Binary>(BinaryOp::Div,
            arena.make<FuncCall>(BuiltinOp::Sqrt,
                std::vector<ExprPtr>{make_integer(arena, BigInt(2))}),
            make_integer(arena, BigInt(2)));
    if (ref == three_tenths) return sin_pi_over_5(arena);
    if (ref == one_third)  return make_rational(arena, Rational(BigInt(1), BigInt(2)));
    if (ref == two_fifths) return cos_2pi_over_5(arena);
    if (ref == one_half)   return make_integer(arena, BigInt(0));
    return nullptr;
}

} // namespace cas::symbolic::detail
