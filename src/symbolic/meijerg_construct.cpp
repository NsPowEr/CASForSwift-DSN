// A7 Brick 2 — make_meijerg factory + view_meijerg accessor.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Meijer_G_Slater.md §7.

#include "cas/meijerg.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>

namespace cas::symbolic {

namespace {

[[nodiscard]] ExprPtr make_size_literal(AstArena& arena, std::size_t n) {
    return arena.make<IntegerLit>(BigInt(static_cast<std::int64_t>(n)));
}

// §2.2: "IF is_positive_integer(a[j] - b[k])" -- decided via the engine
// itself (ctx.simplify), so this is exact for every case the engine can
// reduce to a literal, and correctly inert (no false rejection) when the
// difference stays symbolic.
[[nodiscard]] Result<bool> is_decidably_positive_integer_diff(
    CASContext& ctx, ExprPtr a_j, ExprPtr b_k) {
    ExprPtr diff = ctx.arena().make<Sum>(std::vector<ExprPtr>{
        a_j, ctx.arena().make<Unary>(UnaryOp::Neg, b_k)});
    auto simplified = ctx.simplify(diff);
    if (simplified.is_error()) return fail<bool>(simplified.error());
    const auto* lit = expr_cast<IntegerLit>(simplified.value());
    return ok(lit != nullptr && lit->value.is_positive());
}

// §7.2: within one sub-group, G is invariant under permutation (each
// sub-group is a plain Gamma product) -- sort by the engine's own
// structural order for canonical form / pointer-equality sharing.
void canonicalize_group(std::vector<ExprPtr>& group) {
    std::sort(group.begin(), group.end(),
        [](ExprPtr lhs, ExprPtr rhs) { return canonical_compare(lhs, rhs) < 0; });
}

}  // namespace

Result<ExprPtr> make_meijerg(
    CASContext& ctx, std::size_t m, std::size_t n,
    std::vector<ExprPtr> a, std::vector<ExprPtr> b, ExprPtr z) {
    const std::size_t p = a.size();
    const std::size_t q = b.size();

    if (m > q || n > p) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "make_meijerg: requires 0<=m<=q and 0<=n<=p (got m="
                + std::to_string(m) + ", n=" + std::to_string(n)
                + ", p=" + std::to_string(p) + ", q=" + std::to_string(q) + ")"});
    }
    if (z == nullptr) {
        return fail<ExprPtr>(CASError{
            .kind = CASErrorKind::InvalidArgument,
            .message = "make_meijerg: argument z must not be null"});
    }
    if (p + q > ctx.meijerg_max_param_count()) {
        return fail<ExprPtr>(make_unimplemented_error(
            UnimplementedInfo{
                .module = "symbolic",
                .function = "make_meijerg",
                .input_shape = "p+q = " + std::to_string(p + q)
                    + " parameters",
                .reason = error::reason_codes::GENERIC,
                .suggestion = "Increase meijerg_max_param_count in CASContext",
                .ticket = "A7"},
            "Meijer G parameter count exceeds meijerg_max_param_count"));
    }

    // §2.2 pole-overlap: j in 1..n (a's first sub-group), k in 1..m (b's
    // first sub-group) -- both are numerator Gamma factors of the integrand.
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t k = 0; k < m; ++k) {
            auto overlap = is_decidably_positive_integer_diff(ctx, a[j], b[k]);
            if (overlap.is_error()) return fail<ExprPtr>(overlap.error());
            if (overlap.value()) {
                return fail<ExprPtr>(CASError{
                    .kind = CASErrorKind::InvalidArgument,
                    .message = "make_meijerg: pole overlap, a[" + std::to_string(j)
                        + "] - b[" + std::to_string(k)
                        + "] is a positive integer"});
            }
        }
    }

    // §7.2 canonicalization: the four sub-groups independently, order
    // between them left untouched.
    std::vector<ExprPtr> a_n(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(n));
    std::vector<ExprPtr> a_rest(a.begin() + static_cast<std::ptrdiff_t>(n), a.end());
    std::vector<ExprPtr> b_m(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(m));
    std::vector<ExprPtr> b_rest(b.begin() + static_cast<std::ptrdiff_t>(m), b.end());
    canonicalize_group(a_n);
    canonicalize_group(a_rest);
    canonicalize_group(b_m);
    canonicalize_group(b_rest);

    AstArena& arena = ctx.arena();
    std::vector<ExprPtr> args;
    args.reserve(4U + p + q + 1U);
    args.push_back(make_size_literal(arena, m));
    args.push_back(make_size_literal(arena, n));
    args.push_back(make_size_literal(arena, p));
    args.push_back(make_size_literal(arena, q));
    args.insert(args.end(), a_n.begin(), a_n.end());
    args.insert(args.end(), a_rest.begin(), a_rest.end());
    args.insert(args.end(), b_m.begin(), b_m.end());
    args.insert(args.end(), b_rest.begin(), b_rest.end());
    args.push_back(z);

    return ok(arena.make<FuncCall>(BuiltinOp::MeijerG, std::move(args)));
}

Result<MeijerGView> view_meijerg(const FuncCall& call) {
    if (call.func_id != BuiltinOp::MeijerG) {
        return fail<MeijerGView>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "view_meijerg: not a MeijerG FuncCall"});
    }
    if (call.args.size() < 5U) {
        return fail<MeijerGView>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "view_meijerg: malformed MeijerG args (too few)"});
    }
    auto read_size = [&](std::size_t idx, const char* field) -> Result<std::size_t> {
        const auto* lit = expr_cast<IntegerLit>(call.args[idx]);
        if (lit == nullptr || lit->value.is_negative()) {
            return fail<std::size_t>(CASError{
                .kind = CASErrorKind::InternalError,
                .message = std::string("view_meijerg: ") + field
                    + " is not a non-negative IntegerLit"});
        }
        return ok(static_cast<std::size_t>(lit->value.to_u64()));
    };
    auto m = read_size(0, "m");
    if (m.is_error()) return fail<MeijerGView>(m.error());
    auto n = read_size(1, "n");
    if (n.is_error()) return fail<MeijerGView>(n.error());
    auto p = read_size(2, "p");
    if (p.is_error()) return fail<MeijerGView>(p.error());
    auto q = read_size(3, "q");
    if (q.is_error()) return fail<MeijerGView>(q.error());

    if (call.args.size() != 4U + p.value() + q.value() + 1U) {
        return fail<MeijerGView>(CASError{
            .kind = CASErrorKind::InternalError,
            .message = "view_meijerg: args.size() does not match 4+p+q+1"});
    }

    MeijerGView view;
    view.m = m.value();
    view.n = n.value();
    view.p = p.value();
    view.q = q.value();
    view.a.assign(call.args.begin() + 4,
                  call.args.begin() + static_cast<std::ptrdiff_t>(4U + p.value()));
    view.b.assign(call.args.begin() + static_cast<std::ptrdiff_t>(4U + p.value()),
                  call.args.begin() + static_cast<std::ptrdiff_t>(4U + p.value() + q.value()));
    view.z = call.args.back();
    return ok(view);
}

}  // namespace cas::symbolic
