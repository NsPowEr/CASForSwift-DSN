#pragma once

// A31 fase 1 — side-conditions accumulated by simplify() rewrite steps that
// are only valid under a domain assumption (e.g. x/x -> 1 requires x != 0).
// simplify() itself is UNCHANGED (still generic-point in its ExprPtr output);
// this is a parallel channel exposed via CASContext::last_side_conditions()
// / emit_side_condition(), consumed opt-in.
// Spec: .APROJECT_REFERENCES/MISSING_FEATURES_SPECS/Domain_Conditions_Propagation.md

#include "cas/ast.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cas::symbolic {

class Assumptions;

// Closed, extendable vocabulary (CLAUDE.md hardcode category 3: exhaustive
// enum, never a free-form string subject/predicate pair).
enum class DomainConditionKind : std::uint8_t {
    NonZero,
    Positive,
    NonNegative,
    Real,
    IntegerVal,
    PrincipalBranch,
};

// P(subject): the simplification that emitted this is valid only where the
// predicate `kind` holds of `subject` (spec §3.1). An assumption taken, not
// a fact proved — facts already provable via Assumptions are never emitted
// (see is_condition_already_proven below).
struct DomainCondition {
    DomainConditionKind kind;
    ExprPtr subject;

    [[nodiscard]] bool operator==(const DomainCondition& other) const noexcept {
        return kind == other.kind && ExprEqual{}(subject, other.subject);
    }
};

// Ordered (insertion order), deduplicated collection of DomainCondition, with
// the subsumption rule of the spec (§3.4): Positive(e) subsumes NonZero(e)
// and NonNegative(e) on the same subject; exact duplicates collapse to one.
class SideConditionSet {
public:
    // Adds `c` unless a stronger existing condition on the same subject
    // already covers it (subsumption); removes any existing weaker condition
    // that `c` itself subsumes. Returns false iff adding would exceed
    // max_size() — the caller (CASContext::emit_side_condition) turns that
    // into a structured Unimplemented; SideConditionSet itself never drops a
    // condition silently, it just refuses and reports failure.
    bool add(DomainCondition c);

    // add()s every item of `other`; stops (returns false) at the first
    // budget failure, having applied everything before it.
    bool merge(const SideConditionSet& other);

    // Items present in `*this` but not in `mark` (set difference by ==).
    // Used to capture exactly the conditions emitted during one specific
    // simplify() call, for precise per-key cache attribution (spec §4.2).
    [[nodiscard]] SideConditionSet since(const SideConditionSet& mark) const;

    [[nodiscard]] bool contains(const DomainCondition& c) const noexcept;
    [[nodiscard]] const std::vector<DomainCondition>& items() const noexcept { return items_; }
    [[nodiscard]] bool empty() const noexcept { return items_.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return items_.size(); }
    void clear() noexcept { items_.clear(); }

    void set_max_size(std::size_t max_size) noexcept { max_size_ = max_size; }
    [[nodiscard]] std::size_t max_size() const noexcept { return max_size_; }

private:
    std::vector<DomainCondition> items_;
    std::size_t max_size_{256U};
};

// True iff `c` is already a PROVEN fact under `assumptions` (not merely
// assumed) — in which case simplify() need not register it as a side
// condition; it is a fact of the input, not a restriction the rewrite
// introduced (spec §3.3). Dispatches on DomainConditionKind exhaustively;
// PrincipalBranch has no Assumptions-provable counterpart today and is
// therefore never "already proven" — always registered when emitted.
[[nodiscard]] bool is_condition_already_proven(
    const Assumptions& assumptions, const DomainCondition& c);

}  // namespace cas::symbolic
