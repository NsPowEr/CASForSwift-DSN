#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "term_order_internal.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <vector>

namespace cas::symbolic {

TermOrderRelation relation_from_compare(int cmp) noexcept {
    if (cmp < 0) {
        return TermOrderRelation::Less;
    }
    if (cmp > 0) {
        return TermOrderRelation::Greater;
    }
    return TermOrderRelation::Equivalent;
}

TermOrderRelation compare_knuth_bendix_weight_order(ExprPtr lhs, ExprPtr rhs) {
    if (!lhs || !rhs) {
        return lhs == rhs ? TermOrderRelation::Equivalent : TermOrderRelation::Incomparable;
    }
    if (lhs == rhs) {
        return TermOrderRelation::Equivalent;
    }

    const std::size_t lhs_weight = expr_weight(lhs);
    const std::size_t rhs_weight = expr_weight(rhs);
    if (lhs_weight < rhs_weight) {
        return TermOrderRelation::Less;
    }
    if (lhs_weight > rhs_weight) {
        return TermOrderRelation::Greater;
    }

    const int head_cmp = compare_head_precedence(lhs, rhs);
    if (head_cmp != 0) {
        return relation_from_compare(head_cmp);
    }

    const std::vector<ExprPtr> lhs_children = term_order_children(lhs);
    const std::vector<ExprPtr> rhs_children = term_order_children(rhs);
    const std::size_t shared_size = std::min(lhs_children.size(), rhs_children.size());
    for (std::size_t index = 0; index < shared_size; ++index) {
        const TermOrderRelation child_cmp = compare_knuth_bendix_weight_order(lhs_children[index], rhs_children[index]);
        if (child_cmp != TermOrderRelation::Equivalent) {
            return child_cmp;
        }
    }
    if (lhs_children.size() < rhs_children.size()) {
        return TermOrderRelation::Less;
    }
    if (lhs_children.size() > rhs_children.size()) {
        return TermOrderRelation::Greater;
    }

    return relation_from_compare(canonical_compare(lhs, rhs));
}

}  // namespace cas::symbolic
