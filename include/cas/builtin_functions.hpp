#pragma once

#include <cstdint>
#include <string_view>

namespace cas {

enum class BuiltinOp : std::uint16_t {
    Unknown = 0,
    Sin,
    Cos,
    Tan,
    Exp,
    Ln,
    Log,
    Log10,
    Sqrt,
    Asin,
    Acos,
    Atan,
    Cot,
    Sec,
    Csc,
    Sinh,
    Cosh,
    Tanh,
    Coth,
    Abs,
    Sign,
    Gamma,
    Zeta,
    Erf,
    Floor,
    Ceil,
    Round,
    Min,
    Max,
    Det,
    Rank,
    Trace,
    Inv,
    Transpose,
    N,
    BesselJ,
    BesselY,
    BesselZero,
    BesselI,
    BesselK,
    Binomial,
    SumFunc,
    RootSum,
    Re,
    Im,
    Conj,
    Arg,
    LegendreP,
    Beta,
    Pochhammer,
    Digamma,
    Polygamma,
    ChebyshevT,
    ChebyshevU,
    HermiteH,
    HermiteHe,
    Piecewise,
    JacobiP,
    LaguerreL,
    LambertW,
    DiracDelta,
    Heaviside,
    Hypergeometric0F1,
    Hypergeometric1F1,
    Hypergeometric2F1,
    EllipticK,
    EllipticE,
    EllipticPi,
    EllipticF,
    Factorial,
    Erfc,
    // F8.0-6.1: Unwinding number  K(z) = ⌈(Im(z) − π) / (2π)⌉.
    // Encodes branch-cut information so that  ln(exp(z)) = z + 2πi·K(z).
    UnwindingNumber,
    // A7 Brick 1 — Meijer G-function G^{m,n}_{p,q}(z | a_1..a_p; b_1..b_q).
    // FuncCall::args layout: [m, n, p, q, a_1..a_p, b_1..b_q, z] (flat, size
    // = 4 + p + q + 1). m,n,p,q are IntegerLit; construction/validation goes
    // through the make_meijerg factory (A7 Brick 2), not here — this enum
    // entry only gives the node a name for parse/print round-trip.
    // Spec: MISSING_FEATURES_SPECS/Meijer_G_Slater.md §7.
    MeijerG,
};

[[nodiscard]] constexpr BuiltinOp get_builtin_op(std::string_view name) noexcept {
    if (name == "sin") return BuiltinOp::Sin;
    if (name == "cos") return BuiltinOp::Cos;
    if (name == "tan") return BuiltinOp::Tan;
    if (name == "exp") return BuiltinOp::Exp;
    if (name == "ln") return BuiltinOp::Ln;
    if (name == "log") return BuiltinOp::Log;
    if (name == "log10") return BuiltinOp::Log10;
    if (name == "sqrt") return BuiltinOp::Sqrt;
    if (name == "asin" || name == "arcsin") return BuiltinOp::Asin;
    if (name == "acos" || name == "arccos") return BuiltinOp::Acos;
    if (name == "atan" || name == "arctan") return BuiltinOp::Atan;
    if (name == "atan2") return BuiltinOp::Atan;
    if (name == "cot") return BuiltinOp::Cot;
    if (name == "sec") return BuiltinOp::Sec;
    if (name == "csc") return BuiltinOp::Csc;
    if (name == "sinh") return BuiltinOp::Sinh;
    if (name == "cosh") return BuiltinOp::Cosh;
    if (name == "tanh") return BuiltinOp::Tanh;
    if (name == "coth") return BuiltinOp::Coth;
    if (name == "abs") return BuiltinOp::Abs;
    if (name == "sign") return BuiltinOp::Sign;
    if (name == "gamma") return BuiltinOp::Gamma;
    if (name == "zeta") return BuiltinOp::Zeta;
    if (name == "erf") return BuiltinOp::Erf;
    if (name == "erfc") return BuiltinOp::Erfc;
    if (name == "factorial") return BuiltinOp::Factorial;
    if (name == "floor") return BuiltinOp::Floor;
    if (name == "ceil") return BuiltinOp::Ceil;
    if (name == "round") return BuiltinOp::Round;
    if (name == "min") return BuiltinOp::Min;
    if (name == "max") return BuiltinOp::Max;
    if (name == "det") return BuiltinOp::Det;
    if (name == "rank") return BuiltinOp::Rank;
    if (name == "trace") return BuiltinOp::Trace;
    if (name == "inv") return BuiltinOp::Inv;
    if (name == "transpose") return BuiltinOp::Transpose;
    if (name == "N") return BuiltinOp::N;
    if (name == "BesselJ" || name == "bessel_j") return BuiltinOp::BesselJ;
    if (name == "BesselY" || name == "bessel_y") return BuiltinOp::BesselY;
    if (name == "BesselZero" || name == "bessel_zero") return BuiltinOp::BesselZero;
    if (name == "BesselI" || name == "bessel_i") return BuiltinOp::BesselI;
    if (name == "BesselK" || name == "bessel_k") return BuiltinOp::BesselK;
    if (name == "binomial") return BuiltinOp::Binomial;
    if (name == "sum") return BuiltinOp::SumFunc;
    if (name == "RootSum" || name == "rootsum") return BuiltinOp::RootSum;
    if (name == "re" || name == "Re") return BuiltinOp::Re;
    if (name == "im" || name == "Im") return BuiltinOp::Im;
    if (name == "conj" || name == "conjugate") return BuiltinOp::Conj;
    if (name == "arg") return BuiltinOp::Arg;
    if (name == "LegendreP" || name == "legendreP") return BuiltinOp::LegendreP;
    if (name == "Beta" || name == "beta") return BuiltinOp::Beta;
    if (name == "Pochhammer" || name == "pochhammer") return BuiltinOp::Pochhammer;
    if (name == "digamma" || name == "Digamma" || name == "psi") return BuiltinOp::Digamma;
    if (name == "polygamma" || name == "Polygamma") return BuiltinOp::Polygamma;
    if (name == "ChebyshevT" || name == "chebyshevT") return BuiltinOp::ChebyshevT;
    if (name == "ChebyshevU" || name == "chebyshevU") return BuiltinOp::ChebyshevU;
    if (name == "HermiteH" || name == "hermiteH") return BuiltinOp::HermiteH;
    if (name == "HermiteHe" || name == "hermiteHe") return BuiltinOp::HermiteHe;
    if (name == "piecewise" || name == "Piecewise") return BuiltinOp::Piecewise;
    if (name == "JacobiP" || name == "jacobiP") return BuiltinOp::JacobiP;
    if (name == "LaguerreL" || name == "laguerreL") return BuiltinOp::LaguerreL;
    if (name == "LambertW" || name == "lambertW" || name == "W") return BuiltinOp::LambertW;
    if (name == "DiracDelta" || name == "delta") return BuiltinOp::DiracDelta;
    if (name == "Heaviside" || name == "heaviside" || name == "theta" || name == "Theta") return BuiltinOp::Heaviside;
    if (name == "Hypergeometric0F1" || name == "hyp0F1") return BuiltinOp::Hypergeometric0F1;
    if (name == "Hypergeometric1F1" || name == "hyp1F1" || name == "Kummer1F1" || name == "KummerM") return BuiltinOp::Hypergeometric1F1;
    if (name == "Hypergeometric2F1" || name == "hyp2F1") return BuiltinOp::Hypergeometric2F1;
    if (name == "EllipticK" || name == "ellipticK") return BuiltinOp::EllipticK;
    if (name == "EllipticE" || name == "ellipticE") return BuiltinOp::EllipticE;
    if (name == "EllipticPi" || name == "ellipticPi") return BuiltinOp::EllipticPi;
    if (name == "EllipticF" || name == "ellipticF") return BuiltinOp::EllipticF;
    if (name == "UnwindingNumber" || name == "K") return BuiltinOp::UnwindingNumber;
    if (name == "MeijerG") return BuiltinOp::MeijerG;
    return BuiltinOp::Unknown;
}

[[nodiscard]] constexpr std::string_view builtin_op_name(BuiltinOp op) noexcept {
    switch (op) {
        case BuiltinOp::Sin: return "sin";
        case BuiltinOp::Cos: return "cos";
        case BuiltinOp::Tan: return "tan";
        case BuiltinOp::Exp: return "exp";
        case BuiltinOp::Ln: return "ln";
        case BuiltinOp::Log: return "log";
        case BuiltinOp::Log10: return "log10";
        case BuiltinOp::Sqrt: return "sqrt";
        case BuiltinOp::Asin: return "arcsin";
        case BuiltinOp::Acos: return "arccos";
        case BuiltinOp::Atan: return "arctan";
        case BuiltinOp::Cot: return "cot";
        case BuiltinOp::Sec: return "sec";
        case BuiltinOp::Csc: return "csc";
        case BuiltinOp::Sinh: return "sinh";
        case BuiltinOp::Cosh: return "cosh";
        case BuiltinOp::Tanh: return "tanh";
        case BuiltinOp::Coth: return "coth";
        case BuiltinOp::Abs: return "abs";
        case BuiltinOp::Sign: return "sign";
        case BuiltinOp::Gamma: return "gamma";
        case BuiltinOp::Zeta: return "zeta";
        case BuiltinOp::Erf: return "erf";
        case BuiltinOp::Floor: return "floor";
        case BuiltinOp::Ceil: return "ceil";
        case BuiltinOp::Round: return "round";
        case BuiltinOp::Min: return "min";
        case BuiltinOp::Max: return "max";
        case BuiltinOp::Det: return "det";
        case BuiltinOp::Rank: return "rank";
        case BuiltinOp::Trace: return "trace";
        case BuiltinOp::Inv: return "inv";
        case BuiltinOp::Transpose: return "transpose";
        case BuiltinOp::N: return "N";
        case BuiltinOp::BesselJ: return "BesselJ";
        case BuiltinOp::BesselY: return "BesselY";
        case BuiltinOp::BesselZero: return "BesselZero";
        case BuiltinOp::BesselI: return "BesselI";
        case BuiltinOp::BesselK: return "BesselK";
        case BuiltinOp::Binomial: return "binomial";
        case BuiltinOp::SumFunc: return "sum";
        case BuiltinOp::RootSum: return "RootSum";
        case BuiltinOp::Re: return "re";
        case BuiltinOp::Im: return "im";
        case BuiltinOp::Conj: return "conj";
        case BuiltinOp::Arg: return "arg";
        case BuiltinOp::LegendreP: return "LegendreP";
        case BuiltinOp::Beta: return "Beta";
        case BuiltinOp::Pochhammer: return "Pochhammer";
        case BuiltinOp::Digamma: return "digamma";
        case BuiltinOp::Polygamma: return "polygamma";
        case BuiltinOp::ChebyshevT: return "ChebyshevT";
        case BuiltinOp::ChebyshevU: return "ChebyshevU";
        case BuiltinOp::HermiteH: return "HermiteH";
        case BuiltinOp::HermiteHe: return "HermiteHe";
        case BuiltinOp::Piecewise: return "piecewise";
        case BuiltinOp::JacobiP: return "JacobiP";
        case BuiltinOp::LaguerreL: return "LaguerreL";
        case BuiltinOp::LambertW: return "LambertW";
        case BuiltinOp::DiracDelta: return "DiracDelta";
        case BuiltinOp::Heaviside: return "Heaviside";
        case BuiltinOp::Hypergeometric0F1: return "Hypergeometric0F1";
        case BuiltinOp::Hypergeometric1F1: return "Hypergeometric1F1";
        case BuiltinOp::Hypergeometric2F1: return "Hypergeometric2F1";
        case BuiltinOp::EllipticK: return "EllipticK";
        case BuiltinOp::EllipticE: return "EllipticE";
        case BuiltinOp::EllipticPi: return "EllipticPi";
        case BuiltinOp::EllipticF: return "EllipticF";
        case BuiltinOp::Factorial: return "factorial";
        case BuiltinOp::Erfc: return "erfc";
        case BuiltinOp::UnwindingNumber: return "UnwindingNumber";
        case BuiltinOp::MeijerG: return "MeijerG";
        case BuiltinOp::Unknown: return "unknown";
    }
    return "unknown";
}

} // namespace cas
