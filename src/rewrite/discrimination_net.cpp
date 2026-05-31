#include "discrimination_net.hpp"
#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "../symbolic/symbolic_internal.hpp"

// F1.5 — Discrimination Net implementation.
//
// Algorithm reference:
//   Forgy (1982) "RETE: A Fast Algorithm for the Many Pattern / Many Object
//   Pattern Match Problem", Artificial Intelligence 19:17-37.
//   Hoffmann & O'Donnell (1982) "Pattern Matching in Trees",
//   JACM 29(1):68-95 — typed-wildcard tree discrimination.
//
// Two-level trie:
//   Root → ExprKind bucket → optional FuncCall sub-bucket.
//   Universal (wildcard root) → checked for every expression.
//
// Pattern type filters (Level 1 annotation):
//   Symbol name "x_"         → WildcardType::Any
//   Symbol name "x_Integer"  → WildcardType::Integer  (IntegerLit only)
//   Symbol name "x_Positive" → WildcardType::Positive (positive literal only)
//   Symbol name "x_Symbol"   → WildcardType::Symbol   (Symbol node only)

namespace cas::rewrite {

// ── WildcardType helpers ──────────────────────────────────────────────────────

WildcardType wildcard_type_from_name(const std::string& name) noexcept {
    // Typed wildcard naming convention (two forms accepted):
    //   Form A (trailing underscore, preferred): "x_Integer_", "n_Positive_"
    //     → ends in '_'; type annotation is the segment between last two underscores.
    //   Form B (no trailing underscore): "x_Integer", "n_Positive"
    //     → does not end in '_'; type annotation is suffix after last underscore.
    //
    // In both cases the minimum for WildcardType::Any is "x_" (ends in '_',
    // no type annotation between the final two underscores).
    //
    // Extract suffix:
    //   If name ends in '_': find second-to-last '_', suffix = chars between them.
    //   Otherwise: suffix = chars after last '_'.
    if (name.size() < 2U) return WildcardType::Any;

    std::string_view suffix;
    if (name.back() == '_') {
        // Form A: "x_Integer_"
        const auto pos = name.rfind('_', name.size() - 2U);
        if (pos == std::string::npos || pos == 0U) return WildcardType::Any;
        suffix = std::string_view(name.data() + pos + 1U, name.size() - pos - 2U);
    } else {
        // Form B: "x_Integer"
        const auto pos = name.rfind('_');
        if (pos == std::string::npos) return WildcardType::Any;
        suffix = std::string_view(name.data() + pos + 1U, name.size() - pos - 1U);
    }

    if (suffix == "Integer")  return WildcardType::Integer;
    if (suffix == "Positive") return WildcardType::Positive;
    if (suffix == "Symbol")   return WildcardType::Symbol;
    return WildcardType::Any;
}

bool expr_satisfies_wildcard_type(ExprPtr expr, WildcardType wt) noexcept {
    if (!expr) return false;
    switch (wt) {
        case WildcardType::Any:
            return true;
        case WildcardType::Integer:
            return expr_kind(expr) == ExprKind::IntegerLit;
        case WildcardType::Positive: {
            if (const auto* il = expr_cast<IntegerLit>(expr))
                return !il->value.is_zero() && !il->value.is_negative();
            if (const auto* rl = expr_cast<RationalLit>(expr))
                return !rl->numerator.is_zero() && !rl->numerator.is_negative();
            return false;
        }
        case WildcardType::Symbol:
            return expr_kind(expr) == ExprKind::Symbol;
    }
    return false;  // unreachable
}

// ── DiscriminationNet ─────────────────────────────────────────────────────────

// Helper: is this pattern a plain wildcard (root is a Symbol ending in '_')?
bool DiscriminationNet::is_universal_wildcard(ExprPtr pattern) noexcept {
    if (!pattern) return false;
    if (expr_kind(pattern) != ExprKind::Symbol) return false;
    const auto& name = expr_ref<Symbol>(pattern).name;
    return name.size() >= 2U && name.back() == '_';
}

void DiscriminationNet::insert(const symbolic::RewriteRule& rule) {
    if (!rule.pattern) return;

    if (is_universal_wildcard(rule.pattern)) {
        universal_bucket_.push_back(&rule);
        ++total_rules_;
        return;
    }

    const ExprKind k = expr_kind(rule.pattern);
    const auto idx = static_cast<std::size_t>(k);

    if (k == ExprKind::FuncCall) {
        // Index by function name for FuncCall.
        const auto& fc = expr_ref<FuncCall>(rule.pattern);
        funcall_buckets_[fc.name].push_back(&rule);
    } else {
        if (idx < kNumKinds) {
            kind_buckets_[idx].push_back(&rule);
        }
    }
    ++total_rules_;
}

const std::vector<const symbolic::RewriteRule*>&
DiscriminationNet::lookup(ExprPtr expr) const noexcept {
    if (!expr) return empty_;

    const ExprKind k = expr_kind(expr);
    const auto idx = static_cast<std::size_t>(k);

    // Determine whether we need to merge universal + kind-specific candidates
    // or can return a single bucket directly.
    const bool has_universal = !universal_bucket_.empty();

    if (k == ExprKind::FuncCall) {
        const auto& fc = expr_ref<FuncCall>(expr);
        const auto it = funcall_buckets_.find(fc.name);
        const bool has_specific = (it != funcall_buckets_.end() && !it->second.empty());

        if (!has_universal && !has_specific) return empty_;
        if (!has_universal && has_specific)  return it->second;
        if (has_universal && !has_specific)  return universal_bucket_;

        // Merge into scratch.
        scratch_.clear();
        scratch_.insert(scratch_.end(), it->second.begin(), it->second.end());
        scratch_.insert(scratch_.end(), universal_bucket_.begin(), universal_bucket_.end());
        return scratch_;
    }

    // Non-FuncCall: kind bucket.
    if (idx >= kNumKinds) return empty_;
    const auto& bucket = kind_buckets_[idx];

    if (!has_universal && bucket.empty()) return empty_;
    if (!has_universal) return bucket;
    if (bucket.empty()) return universal_bucket_;

    // Merge.
    scratch_.clear();
    scratch_.insert(scratch_.end(), bucket.begin(), bucket.end());
    scratch_.insert(scratch_.end(), universal_bucket_.begin(), universal_bucket_.end());
    return scratch_;
}

void DiscriminationNet::clear() noexcept {
    for (auto& b : kind_buckets_) b.clear();
    funcall_buckets_.clear();
    universal_bucket_.clear();
    scratch_.clear();
    total_rules_ = 0U;
}

} // namespace cas::rewrite
