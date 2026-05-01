#include "cas/ast.hpp"
#include "cas/rational.hpp"
#include "cas/symbolic.hpp"
#include "cas/algebra.hpp"
#include "algebra_internal.hpp"
#include "polynomial_internal.hpp"
#include <algorithm>
#include <vector>
#include <numeric>

namespace cas::algebra {

// --- 2. Polynomial Conversions & Parsing ---

[[nodiscard]] Result<PolyExpr> parse_polynomial(ExprPtr expr_in, const Symbol& var, symbolic::CASContext& ctx) {
    auto expanded_res = expand(expr_in, ctx);
    if (expanded_res.is_error()) return fail<PolyExpr>(expanded_res.error());
    ExprPtr expr = expanded_res.value();
    
    auto simp_res = ctx.simplify(expr);
    if (simp_res.is_ok()) expr = simp_res.value();

    std::unordered_map<std::size_t, ExprPtr> coeffs;
    std::size_t max_deg = 0;

    std::function<Result<std::pair<std::size_t, ExprPtr>>(ExprPtr, bool)> extract_deg_coeff;
    extract_deg_coeff = [&](ExprPtr term, bool negate) -> Result<std::pair<std::size_t, ExprPtr>> {
        if (const auto* u = expr_cast<Unary>(term); u && u->op == UnaryOp::Neg) {
            return extract_deg_coeff(u->operand, !negate);
        }
        if (const auto* b_pow = expr_cast<Binary>(term); b_pow && b_pow->op == BinaryOp::Pow) {
            if (const auto* s = expr_cast<Symbol>(b_pow->left); s && s->name == var.name) {
                auto expo = poly_parse_nonnegative_integer_exponent(b_pow->right);
                if (expo.is_error()) return fail<std::pair<std::size_t, ExprPtr>>(expo.error());
                ExprPtr c = poly_make_integer(ctx.arena(), negate ? -1 : 1);
                return ok(std::make_pair(expo.value(), c));
            }
        }
        if (const auto* b_mul = expr_cast<Binary>(term); b_mul && b_mul->op == BinaryOp::Mul) {
            if (const auto* bp = expr_cast<Binary>(b_mul->right); bp && bp->op == BinaryOp::Pow) {
                if (const auto* sym = expr_cast<Symbol>(bp->left); sym && sym->name == var.name) {
                    auto expo = poly_parse_nonnegative_integer_exponent(bp->right);
                    if (expo.is_error()) return fail<std::pair<std::size_t, ExprPtr>>(expo.error());
                    ExprPtr c = negate ? ctx.arena().make<Unary>(UnaryOp::Neg, b_mul->left) : b_mul->left;
                    auto cs = ctx.simplify(c);
                    if (cs.is_error()) return fail<std::pair<std::size_t, ExprPtr>>(cs.error());
                    return ok(std::make_pair(expo.value(), cs.value()));
                }
            } else if (const auto* sym = expr_cast<Symbol>(b_mul->right); sym && sym->name == var.name) {
                ExprPtr c = negate ? ctx.arena().make<Unary>(UnaryOp::Neg, b_mul->left) : b_mul->left;
                auto cs = ctx.simplify(c);
                if (cs.is_error()) return fail<std::pair<std::size_t, ExprPtr>>(cs.error());
                return ok(std::make_pair(1UL, cs.value()));
            }
        }
        if (const auto* prod = expr_cast<Product>(term)) {
            std::size_t deg = 0;
            std::vector<ExprPtr> coeff_factors;
            for (ExprPtr f : prod->factors) {
                bool is_var_pow = false;
                if (const auto* s = expr_cast<Symbol>(f); s && s->name == var.name) {
                    deg = 1; is_var_pow = true;
                } else if (const auto* bp = expr_cast<Binary>(f); bp && bp->op == BinaryOp::Pow) {
                    if (const auto* s = expr_cast<Symbol>(bp->left); s && s->name == var.name) {
                        auto expo = poly_parse_nonnegative_integer_exponent(bp->right);
                        if (!expo.is_error()) { deg = expo.value(); is_var_pow = true; }
                    }
                }
                if (!is_var_pow) coeff_factors.push_back(f);
            }
            if (!coeff_factors.empty() && deg > 0) {
                ExprPtr c = coeff_factors.size() == 1U
                    ? coeff_factors.front()
                    : static_cast<ExprPtr>(ctx.arena().make<Product>(std::move(coeff_factors)));
                if (negate) c = ctx.arena().make<Unary>(UnaryOp::Neg, c);
                auto cs = ctx.simplify(c);
                if (cs.is_error()) return fail<std::pair<std::size_t, ExprPtr>>(cs.error());
                return ok(std::make_pair(deg, cs.value()));
            }
        }
        if (const auto* sym = expr_cast<Symbol>(term); sym && sym->name == var.name) {
            ExprPtr c = poly_make_integer(ctx.arena(), negate ? -1 : 1);
            return ok(std::make_pair(1UL, c));
        }
        ExprPtr c = negate ? ctx.arena().make<Unary>(UnaryOp::Neg, term) : term;
        auto cs = ctx.simplify(c);
        if (cs.is_error()) return fail<std::pair<std::size_t, ExprPtr>>(cs.error());
        return ok(std::make_pair(0UL, cs.value()));
    };

    auto process_term = [&](ExprPtr term) -> Result<void> {
        auto dc = extract_deg_coeff(term, false);
        if (dc.is_error()) return fail<void>(dc.error());
        auto [degree, coeff] = dc.value();
        auto it = coeffs.find(degree);
        if (it == coeffs.end()) {
            coeffs[degree] = coeff;
        } else {
            auto sum = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Add, it->second, coeff));
            if (sum.is_error()) return fail<void>(sum.error());
            coeffs[degree] = sum.value();
        }
        if (degree > max_deg) max_deg = degree;
        return ok();
    };

    if (const auto* s = expr_cast<Sum>(expr)) {
        for (auto t : s->terms) {
            auto res = process_term(t);
            if (res.is_error()) return fail<PolyExpr>(res.error());
        }
    } else {
        auto res = process_term(expr);
        if (res.is_error()) return fail<PolyExpr>(res.error());
    }

    PolyExpr poly;
    poly.resize(max_deg + 1, poly_make_integer(ctx.arena(), 0));
    for (auto const& [d, c] : coeffs) {
        poly[d] = c;
    }
    normalize_poly(poly);
    return ok(poly);
}

[[nodiscard]] Result<ExprPtr> polynomial_to_expr(const PolyExpr& poly, const Symbol& var, symbolic::CASContext& ctx) {
    if (poly.is_zero()) return ok(poly_make_integer(ctx.arena(), 0));
    std::vector<ExprPtr> terms;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        if (poly_is_zero_expr(poly[i])) continue;
        ExprPtr term;
        if (i == 0) {
            term = poly[i];
        } else {
            ExprPtr var_part = (i == 1) ? ctx.arena().make<Symbol>(var.name) 
                                        : ctx.arena().make<Binary>(BinaryOp::Pow, ctx.arena().make<Symbol>(var.name), poly_make_integer(ctx.arena(), i));
            if (poly_is_one_expr(poly[i])) term = var_part;
            else if (poly_is_minus_one_expr(poly[i])) term = ctx.arena().make<Unary>(UnaryOp::Neg, var_part);
            else term = ctx.arena().make<Binary>(BinaryOp::Mul, poly[i], var_part);
        }
        terms.push_back(term);
    }
    if (terms.empty()) return ok(poly_make_integer(ctx.arena(), 0));
    if (terms.size() == 1) return ok(terms[0]);
    std::reverse(terms.begin(), terms.end());
    return ok(ctx.arena().make<Sum>(std::move(terms)));
}

// --- 3. Arithmetic (Expr based) ---

void normalize_poly(PolyExpr& poly) {
    poly.normalize([](ExprPtr e) { return !e || poly_is_zero_expr(e); });
}

[[nodiscard]] Result<PolyExpr> poly_add(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx) {
    PolyExpr res;
    std::size_t n = std::max(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < n; ++i) {
        ExprPtr l = (i < lhs.size()) ? lhs[i] : poly_make_integer(ctx.arena(), 0);
        ExprPtr r = (i < rhs.size()) ? rhs[i] : poly_make_integer(ctx.arena(), 0);
        auto sum = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Add, l, r));
        if (sum.is_error()) return fail<PolyExpr>(sum.error());
        res.push_back(sum.value());
    }
    normalize_poly(res);
    return ok(res);
}

[[nodiscard]] Result<PolyExpr> poly_subtract(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx) {
    PolyExpr res;
    std::size_t n = std::max(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < n; ++i) {
        ExprPtr l = (i < lhs.size()) ? lhs[i] : poly_make_integer(ctx.arena(), 0);
        ExprPtr r = (i < rhs.size()) ? rhs[i] : poly_make_integer(ctx.arena(), 0);
        auto sub = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Sub, l, r));
        if (sub.is_error()) return fail<PolyExpr>(sub.error());
        res.push_back(sub.value());
    }
    normalize_poly(res);
    return ok(res);
}

[[nodiscard]] Result<PolyExpr> poly_negate(const PolyExpr& poly, symbolic::CASContext& ctx) {
    PolyExpr res;
    for (auto c : poly.coefficients()) {
        auto neg = ctx.simplify(ctx.arena().make<Unary>(UnaryOp::Neg, c));
        if (neg.is_error()) return fail<PolyExpr>(neg.error());
        res.push_back(neg.value());
    }
    return ok(res);
}

[[nodiscard]] Result<PolyExpr> poly_multiply(const PolyExpr& lhs, const PolyExpr& rhs, symbolic::CASContext& ctx) {
    if (lhs.empty() || rhs.empty()) return ok(PolyExpr{});
    PolyExpr res;
    res.resize(lhs.size() + rhs.size() - 1, poly_make_integer(ctx.arena(), 0));
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        if (poly_is_zero_expr(lhs[i])) continue;
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            if (poly_is_zero_expr(rhs[j])) continue;
            auto prod = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Mul, lhs[i], rhs[j]));
            if (prod.is_error()) return fail<PolyExpr>(prod.error());
            auto sum = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Add, res[i + j], prod.value()));
            if (sum.is_error()) return fail<PolyExpr>(sum.error());
            res[i + j] = sum.value();
        }
    }
    normalize_poly(res);
    return ok(res);
}

[[nodiscard]] Result<PolyExpr> poly_divide_by_scalar(const PolyExpr& poly, ExprPtr scalar, symbolic::CASContext& ctx) {
    PolyExpr res;
    for (auto c : poly.coefficients()) {
        auto div = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Div, c, scalar));
        if (div.is_error()) return fail<PolyExpr>(div.error());
        res.push_back(div.value());
    }
    normalize_poly(res);
    return ok(res);
}

[[nodiscard]] Result<PolyDivisionResult> divide_poly_with_remainder(
    const PolyExpr& dividend,
    const PolyExpr& divisor,
    symbolic::CASContext& ctx) {
    
    if (divisor.is_zero()) return fail<PolyDivisionResult>(make_error(CASErrorKind::Undefined, "Divisione polinomiale per zero"));
    PolyExpr q, r = dividend;
    if (r.size() < divisor.size()) return ok(PolyDivisionResult{q, r});
    q.resize(dividend.size() - divisor.size() + 1, poly_make_integer(ctx.arena(), 0));
    auto coeff_or_zero = [&](ExprPtr c) -> ExprPtr { return c ? c : poly_make_integer(ctx.arena(), 0); };

    for (int i = static_cast<int>(dividend.size() - divisor.size()); i >= 0; --i) {
        if (r.size() < divisor.size() + static_cast<std::size_t>(i)) continue;
        auto term_q = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Div, coeff_or_zero(r.leading_coeff()), coeff_or_zero(divisor.leading_coeff())));
        if (term_q.is_error()) return fail<PolyDivisionResult>(term_q.error());
        q[i] = term_q.value();
        for (std::size_t j = 0; j < divisor.size(); ++j) {
            auto prod = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Mul, q[i], coeff_or_zero(divisor[j])));
            if (prod.is_error()) return fail<PolyDivisionResult>(prod.error());
            auto sub = ctx.simplify(ctx.arena().make<Binary>(BinaryOp::Sub, coeff_or_zero(r[i + j]), prod.value()));
            if (sub.is_error()) return fail<PolyDivisionResult>(sub.error());
            r[i + j] = sub.value();
        }
        normalize_poly(r);
    }
    return ok(PolyDivisionResult{q, r});
}

// --- 4. Rational Poly Logic ---

void normalize_rational_coefficients(RatPoly& coefficients) {
    coefficients.normalize([](const Rational& r) { return r.numerator().is_zero(); });
}

[[nodiscard]] RatPoly add_rational_poly(const RatPoly& a, const RatPoly& b) {
    RatPoly res;
    std::size_t n = std::max(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) {
        Rational ca = (i < a.size()) ? a[i] : Rational(0);
        Rational cb = (i < b.size()) ? b[i] : Rational(0);
        res.push_back(ca + cb);
    }
    res.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return res;
}

[[nodiscard]] RatPoly sub_rational_poly(const RatPoly& a, const RatPoly& b) {
    RatPoly res;
    std::size_t n = std::max(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) {
        Rational ca = (i < a.size()) ? a[i] : Rational(0);
        Rational cb = (i < b.size()) ? b[i] : Rational(0);
        res.push_back(ca - cb);
    }
    res.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return res;
}

[[nodiscard]] RatPoly mul_rational_poly(const RatPoly& a, const RatPoly& b) {
    if (a.empty() || b.empty()) return RatPoly{};
    RatPoly res;
    res.resize(a.size() + b.size() - 1, Rational(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].numerator().is_zero()) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            if (b[j].numerator().is_zero()) continue;
            res[i + j] = res[i + j] + (a[i] * b[j]);
        }
    }
    res.normalize([](const Rational& r) { return r.numerator().is_zero(); });
    return res;
}

[[nodiscard]] std::pair<RatPoly, RatPoly> div_rem_rational_poly(const RatPoly& a, const RatPoly& b) {
    if (b.is_zero()) return {RatPoly{}, a};
    RatPoly q, r = a;
    if (r.size() < b.size()) return {q, r};
    q.resize(a.size() - b.size() + 1, Rational(0));
    Rational lb_inv = Rational(b.leading_coeff().denominator(), b.leading_coeff().numerator());
    for (int i = static_cast<int>(a.size() - b.size()); i >= 0; --i) {
        if (r.size() < b.size() + i) continue;
        Rational term_q = r.leading_coeff() * lb_inv;
        q[i] = term_q;
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] = r[i + j] - (term_q * b[j]);
        }
        r.normalize([](const Rational& rat) { return rat.numerator().is_zero(); });
    }
    return {q, r};
}

[[nodiscard]] std::tuple<RatPoly, RatPoly, RatPoly> extended_gcd_rational_poly(const RatPoly& a, const RatPoly& b) {
    if (a.is_zero()) return {b, RatPoly{}, RatPoly{{Rational(1)}}};
    if (b.is_zero()) return {a, RatPoly{{Rational(1)}}, RatPoly{}};
    auto [q, r] = div_rem_rational_poly(a, b);
    if (r.is_zero()) return {b, RatPoly{}, RatPoly{{Rational(1)}}};
    auto [g, s1, t1] = extended_gcd_rational_poly(b, r);
    return {g, t1, sub_rational_poly(s1, mul_rational_poly(q, t1))};
}
// --- 5. Advanced & Bridge Functions ---

[[nodiscard]] std::size_t poly_degree(const PolyExpr& poly) { return poly.degree(); }
[[nodiscard]] ExprPtr leading_coefficient(const PolyExpr& poly) { return poly.empty() ? nullptr : poly.leading_coeff(); }

[[nodiscard]] bool is_zero_poly(const PolyExpr& poly) {
    for (auto c : poly.coefficients()) if (!poly_is_zero_expr(c)) return false;
    return true;
}

[[nodiscard]] PolyExpr poly_make_monomial(ExprPtr coefficient, std::size_t degree) {
    PolyExpr res;
    res.resize(degree + 1, nullptr);
    res[degree] = coefficient;
    return res;
}

[[nodiscard]] PolyExpr integer_coefficients_to_poly(const IntPoly& poly, AstArena& arena) {
    std::vector<ExprPtr> expr_coeffs;
    for(auto c : poly.coefficients()) expr_coeffs.push_back(arena.make<IntegerLit>(c));
    return PolyExpr(expr_coeffs);
}

[[nodiscard]] Result<PolyXGCDResult> poly_extended_gcd(const PolyExpr& a, const PolyExpr& b, symbolic::CASContext& ctx) {
    if (a.is_zero()) {
        std::vector<ExprPtr> s_coeffs = {poly_make_integer(ctx.arena(), 0)};
        std::vector<ExprPtr> t_coeffs = {poly_make_integer(ctx.arena(), 1)};
        return ok(PolyXGCDResult{b, PolyExpr{s_coeffs}, PolyExpr{t_coeffs}});
    }
    if (b.is_zero()) {
        std::vector<ExprPtr> s_coeffs = {poly_make_integer(ctx.arena(), 1)};
        std::vector<ExprPtr> t_coeffs = {poly_make_integer(ctx.arena(), 0)};
        return ok(PolyXGCDResult{a, PolyExpr{s_coeffs}, PolyExpr{t_coeffs}});
    }

    auto div_res = divide_poly_with_remainder(a, b, ctx);
    if (div_res.is_error()) return fail<PolyXGCDResult>(div_res.error());
    
    auto q = div_res.value().quotient;
    auto r = div_res.value().remainder;

    if (r.is_zero()) {
        std::vector<ExprPtr> s_coeffs = {poly_make_integer(ctx.arena(), 0)};
        std::vector<ExprPtr> t_coeffs = {poly_make_integer(ctx.arena(), 1)};
        return ok(PolyXGCDResult{b, PolyExpr{s_coeffs}, PolyExpr{t_coeffs}});
    }

    auto ext_res = poly_extended_gcd(b, r, ctx);
    if (ext_res.is_error()) return fail<PolyXGCDResult>(ext_res.error());
    
    auto g = ext_res.value().gcd;
    auto s1 = ext_res.value().s;
    auto t1 = ext_res.value().t;

    auto q_t1_res = poly_multiply(q, t1, ctx);
    if (q_t1_res.is_error()) return fail<PolyXGCDResult>(q_t1_res.error());
    
    auto s_res = poly_subtract(s1, q_t1_res.value(), ctx);
    if (s_res.is_error()) return fail<PolyXGCDResult>(s_res.error());

    return ok(PolyXGCDResult{g, t1, s_res.value()});
}

[[nodiscard]] Result<PolyExpr> normalize_poly_monic(const PolyExpr& poly, symbolic::CASContext& ctx) {
    if (poly.is_zero()) return ok(poly);
    return poly_divide_by_scalar(poly, poly.leading_coeff(), ctx);
}

[[nodiscard]] Result<std::size_t> poly_parse_nonnegative_integer_exponent(ExprPtr expr) {
    if (const auto* i = expr_cast<IntegerLit>(expr)) {
        if (!i->value.is_negative()) return ok(static_cast<std::size_t>(i->value.to_u64()));
    }
    return fail<std::size_t>(make_error(CASErrorKind::InvalidArgument, "Esponente non valido"));
}

[[nodiscard]] std::optional<RationalRootCandidate> find_rational_root_candidate(const IntPoly& coefficients) {
    if (coefficients.size() <= 1U) return std::nullopt;
    const std::vector<BigInt> p_divs = positive_divisors_or_one(coefficients[0]);
    const std::vector<BigInt> q_divs = positive_divisors_or_one(coefficients.leading_coeff());
    for (const BigInt& q : q_divs) {
        for (const BigInt& p : p_divs) {
            for (int sign : {1, -1}) {
                BigInt num = sign == 1 ? p : -p;
                Rational candidate(num, q);
                if (evaluate_integer_polynomial_at(coefficients, candidate).numerator().is_zero()) {
                    return RationalRootCandidate{num, q};
                }
            }
        }
    }
    return std::nullopt;
}

} // namespace cas::algebra
