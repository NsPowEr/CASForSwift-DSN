// A31 fase 1 — SideConditionSet implementation + CASContext::emit_side_condition.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Domain_Conditions_Propagation.md

#include "cas/side_conditions.hpp"
#include "cas/error.hpp"
#include "cas/symbolic.hpp"

#include <algorithm>

namespace cas::symbolic {

namespace {

[[nodiscard]] bool is_stronger_and_covers(const DomainCondition& stronger,
                                          const DomainCondition& c) noexcept {
    // Positive(e) subsumes NonZero(e) and NonNegative(e) — spec §3.4.
    if (stronger.kind != DomainConditionKind::Positive) return false;
    if (c.kind != DomainConditionKind::NonZero && c.kind != DomainConditionKind::NonNegative) {
        return false;
    }
    return ExprEqual{}(stronger.subject, c.subject);
}

}  // namespace

bool SideConditionSet::contains(const DomainCondition& c) const noexcept {
    return std::find(items_.begin(), items_.end(), c) != items_.end();
}

bool SideConditionSet::add(DomainCondition c) {
    if (contains(c)) return true;  // exact duplicate: no-op, not a failure.

    for (const auto& existing : items_) {
        if (is_stronger_and_covers(existing, c)) return true;  // already subsumed.
    }

    if (c.kind == DomainConditionKind::Positive) {
        items_.erase(
            std::remove_if(items_.begin(), items_.end(),
                [&](const DomainCondition& existing) {
                    return is_stronger_and_covers(c, existing);
                }),
            items_.end());
    }

    if (items_.size() >= max_size_) return false;
    items_.push_back(c);
    return true;
}

bool SideConditionSet::merge(const SideConditionSet& other) {
    for (const auto& c : other.items_) {
        if (!add(c)) return false;
    }
    return true;
}

SideConditionSet SideConditionSet::since(const SideConditionSet& mark) const {
    SideConditionSet delta;
    delta.set_max_size(max_size_);
    for (const auto& c : items_) {
        if (!mark.contains(c)) delta.add(c);
    }
    return delta;
}

[[nodiscard]] bool is_condition_already_proven(
    const Assumptions& assumptions, const DomainCondition& c) {
    switch (c.kind) {
        case DomainConditionKind::NonZero:      return assumptions.is_nonzero(c.subject);
        case DomainConditionKind::Positive:     return assumptions.is_positive(c.subject);
        case DomainConditionKind::NonNegative:  return assumptions.is_nonnegative(c.subject);
        case DomainConditionKind::Real:         return assumptions.is_real(c.subject);
        case DomainConditionKind::IntegerVal:   return assumptions.is_integer(c.subject);
        // No Assumptions predicate corresponds to "principal branch was used"
        // -- it is not a provable fact about the subject's value, but a
        // statement about which analytic continuation the rewrite picked.
        case DomainConditionKind::PrincipalBranch: return false;
    }
    return false;
}

Result<void> CASContext::emit_side_condition(DomainConditionKind kind, ExprPtr subject) {
    DomainCondition c{kind, subject};
    if (is_condition_already_proven(assumptions_, c)) {
        return ok();  // fact of the input, not an assumption taken.
    }
    side_conditions_.set_max_size(max_side_conditions());
    if (side_conditions_.add(c)) {
        return ok();
    }
    return Result<void>(make_unimplemented_error(
        UnimplementedInfo{
            .module = "symbolic",
            .function = "CASContext::emit_side_condition",
            .input_shape = "more than " + std::to_string(max_side_conditions())
                + " distinct domain conditions in one simplify() call",
            .reason = error::reason_codes::SIDE_CONDITION_BUDGET_EXCEEDED,
            .suggestion = "Increase max_side_conditions in CASContext",
            .ticket = "A31"},
        "Side-condition budget exceeded"));
}

const SideConditionSet& CASContext::last_side_conditions() const noexcept {
    return side_conditions_;
}

}  // namespace cas::symbolic
