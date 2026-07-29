#include "cas/symbolic.hpp"
#include "cas/rational.hpp"
#include "term_order_internal.hpp"
#include "symbolic_internal.hpp"
#include <algorithm>
#include <vector>

namespace cas::symbolic {

bool term_order_gt(ExprPtr lhs, ExprPtr rhs) {
    if (!lhs || !rhs || lhs == rhs) {
        return false;
    }

    const std::vector<ExprPtr> lhs_children = term_order_children(lhs);
    for (ExprPtr child : lhs_children) {
        if (term_order_ge(child, rhs)) {
            return true;
        }
    }

    const std::vector<ExprPtr> rhs_children = term_order_children(rhs);
    const auto dominates_rhs_children = [&]() {
        return std::all_of(
            rhs_children.begin(),
            rhs_children.end(),
            [&](ExprPtr child) {
                return term_order_gt(lhs, child);
            });
    };

    const int head_cmp = compare_head_precedence(lhs, rhs);
    if (head_cmp > 0 && dominates_rhs_children()) {
        return true;
    }

    if (head_cmp == 0 && lhs_children.size() == rhs_children.size() && dominates_rhs_children()) {
        for (std::size_t index = 0; index < lhs_children.size(); ++index) {
            if (structural_equal(lhs_children[index], rhs_children[index])) {
                continue;
            }
            return term_order_gt(lhs_children[index], rhs_children[index]);
        }
    }

    return false;
}

bool term_order_ge(ExprPtr lhs, ExprPtr rhs) {
    return lhs == rhs || term_order_gt(lhs, rhs);
}

}  // namespace cas::symbolic
