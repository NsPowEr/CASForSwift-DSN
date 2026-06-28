#include "differentiate_internal.hpp"

namespace cas::calculus {

struct Differentiator::Visitor {
    Differentiator& self;
    const Symbol& var;

    [[nodiscard]] Result<ExprPtr> operator()(const IntegerLit&) const { return ok(make_integer(self.arena_, 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const RationalLit&) const { return ok(make_integer(self.arena_, 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const Constant&) const { return ok(make_integer(self.arena_, 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const DecimalLit&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Decimal literals are not supported in symbolic differentiation")); }
    // A ComplexLit is a numeric constant (e.g. the imaginary unit −i emitted by
    // the rational integrator); it carries no variable, so d/dx ≡ 0 like any
    // other literal. Bailing here blocked verifying complex-form antiderivatives
    // such as ∫x²/(x²−1) = x − i·arctan(−i·x).
    [[nodiscard]] Result<ExprPtr> operator()(const ComplexLit&) const { return ok(make_integer(self.arena_, 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const Symbol& symbol) const { return ok(make_integer(self.arena_, symbol.name == var.name ? 1 : 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const Unary& unary) const { return self.differentiate_unary(unary, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Binary& binary) const { return self.differentiate_binary(binary, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const FuncCall& call) const { return self.differentiate_function(call, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Sum& sum) const { return self.differentiate_sum(sum, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Product& product) const { return self.differentiate_product(product, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Derivative& derivative) const { return self.differentiate(derivative.expression, derivative.variable, derivative.order + 1U); }
    [[nodiscard]] Result<ExprPtr> operator()(const Integral& integral) const { return self.differentiate_integral(integral, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Limit&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of limit not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const RootOf&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of RootOf not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const Matrix&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of matrix not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const SeriesExp&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of SeriesExp not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const Quantity&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of Quantity not implemented")); }
};

// Fundamental Theorem of Calculus + Leibniz integral rule.
Result<ExprPtr> Differentiator::differentiate_integral(const Integral& integral, const Symbol& var) {
    const bool has_lo = integral.lower.has_value();
    const bool has_hi = integral.upper.has_value();

    // Indefinite integral  ∫ f d(t).
    if (!has_lo && !has_hi) {
        // d/dx ∫ f dx = f  (FTC, integration variable matches).
        if (integral.variable.name == var.name) {
            return ok(integral.integrand);
        }
        // Otherwise differentiate under the integral sign:
        //   d/dx ∫ f(t,x) dt = ∫ (∂f/∂x) dt.
        auto df = differentiate(integral.integrand, var, 1U);
        if (df.is_error()) return df;
        if (is_exact_zero(df.value())) return ok(make_integer(arena_, 0));
        return ok(arena_.make<Integral>(df.value(), integral.variable,
                                        std::nullopt, std::nullopt));
    }

    // Definite integral  ∫_a^b f(t,x) dt — Leibniz rule:
    //   f(b,x)·b'(x) − f(a,x)·a'(x) + ∫_a^b (∂f/∂x) dt.
    if (has_lo && has_hi) {
        ExprPtr a = *integral.lower;
        ExprPtr b = *integral.upper;

        auto f_at_b = context_.substitute(integral.integrand, integral.variable, b);
        if (f_at_b.is_error()) return f_at_b;
        auto f_at_a = context_.substitute(integral.integrand, integral.variable, a);
        if (f_at_a.is_error()) return f_at_a;
        auto db = differentiate(b, var, 1U);
        if (db.is_error()) return db;
        auto da = differentiate(a, var, 1U);
        if (da.is_error()) return da;

        std::vector<ExprPtr> terms;
        terms.push_back(make_product(arena_, {f_at_b.value(), db.value()}));
        terms.push_back(make_unary(arena_, UnaryOp::Neg,
            make_product(arena_, {f_at_a.value(), da.value()})));

        auto df = differentiate(integral.integrand, var, 1U);
        if (df.is_error()) return df;
        if (!is_exact_zero(df.value())) {
            terms.push_back(arena_.make<Integral>(df.value(), integral.variable,
                                                  integral.lower, integral.upper));
        }
        return ok(make_sum(arena_, std::move(terms)));
    }

    return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented,
        "Differentiation of a one-sided improper integral is not implemented"));
}

Result<ExprPtr> Differentiator::differentiate_once(ExprPtr expr, const Symbol& var) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot differentiate a null expression"));
    }
    return visit_expr(expr, Visitor{*this, var});
}

}  // namespace cas::calculus
