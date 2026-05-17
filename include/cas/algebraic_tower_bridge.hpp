#pragma once

#include "cas/algebraic_number.hpp"
#include "cas/algebraic_tower.hpp"
#include "cas/ast.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <optional>
#include <vector>

namespace cas {
namespace algebra {

struct TowerGenerators {
    ExprPtr alpha_1;
    AlgebraicNumber::CoeffVec min_poly_1;

    ExprPtr alpha_2;
    std::vector<AlgebraicNumber> min_poly_2;
};

[[nodiscard]] Result<std::optional<TowerGenerators>> detect_two_level_tower(
    ExprPtr expr,
    symbolic::CASContext& ctx);

[[nodiscard]] Result<std::optional<AlgebraicTowerTwoLevel>> try_express_in_tower_two_level(
    ExprPtr expr,
    const TowerGenerators& gens,
    symbolic::CASContext& ctx);

[[nodiscard]] ExprPtr tower_to_expr(
    const AlgebraicTowerTwoLevel& value,
    const TowerGenerators& gens,
    AstArena& arena);

}  // namespace algebra
}  // namespace cas
