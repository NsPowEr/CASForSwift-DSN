#include "differentiate_internal.hpp"

namespace cas::calculus {

struct Differentiator::Visitor {
    Differentiator& self;
    const Symbol& var;

    [[nodiscard]] Result<ExprPtr> operator()(const IntegerLit&) const { return ok(make_integer(self.arena_, 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const RationalLit&) const { return ok(make_integer(self.arena_, 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const Constant&) const { return ok(make_integer(self.arena_, 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const DecimalLit&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Decimal literals are not supported in symbolic differentiation")); }
    [[nodiscard]] Result<ExprPtr> operator()(const ComplexLit&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Complex literals differentiation not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const Symbol& symbol) const { return ok(make_integer(self.arena_, symbol.name == var.name ? 1 : 0)); }
    [[nodiscard]] Result<ExprPtr> operator()(const Unary& unary) const { return self.differentiate_unary(unary, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Binary& binary) const { return self.differentiate_binary(binary, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const FuncCall& call) const { return self.differentiate_function(call, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Sum& sum) const { return self.differentiate_sum(sum, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Product& product) const { return self.differentiate_product(product, var); }
    [[nodiscard]] Result<ExprPtr> operator()(const Derivative& derivative) const { return self.differentiate(derivative.expression, derivative.variable, derivative.order + 1U); }
    [[nodiscard]] Result<ExprPtr> operator()(const Integral&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of integral not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const Limit&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of limit not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const RootOf&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of RootOf not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const Matrix&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of matrix not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const SeriesExp&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of SeriesExp not implemented")); }
    [[nodiscard]] Result<ExprPtr> operator()(const Quantity&) const { return fail<ExprPtr>(make_error(CASErrorKind::Unimplemented, "Differentiation of Quantity not implemented")); }
};

Result<ExprPtr> Differentiator::differentiate_once(ExprPtr expr, const Symbol& var) {
    if (!expr) {
        return fail<ExprPtr>(make_error(CASErrorKind::InvalidArgument, "Cannot differentiate a null expression"));
    }
    return visit_expr(expr, Visitor{*this, var});
}

}  // namespace cas::calculus
