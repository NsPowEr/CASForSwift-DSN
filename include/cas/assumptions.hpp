#pragma once

// Assumptions — fatti PROVATI sui simboli (dominio, segno, relazioni d'ordine),
// distinti dalle side-conditions di A31, che sono assunzioni PRESE da un
// rewrite (cas/side_conditions.hpp). La differenza e' operativa: cio' che
// `Assumptions` prova non va registrato come condizione
// (`is_condition_already_proven`, spec A31 §3.3).
//
// Estratto da cas/symbolic.hpp (anti-monolito, A53): blocco coeso e con vita
// propria — implementazione in src/symbolic/assumptions_core.cpp e
// assumptions_predicates.cpp.

#include "cas/ast.hpp"
#include "cas/result.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cas::symbolic {

struct RangeAssumption {
    ExprPtr lower;
    ExprPtr upper;
};

enum class RelType : std::uint8_t {
    Less,          // <
    LessEqual      // <=
};

struct Relation {
    ExprPtr target;
    RelType type;
};

enum class Domain : std::uint8_t {
    Complex,
    Real,
    Integer,
    Rational,
    Natural,     // >= 0
    Positive,    // > 0
    Negative,    // < 0
    NonZero,
};

class Assumptions {
public:
    void assume_real(const Symbol& symbol);
    void assume_positive(const Symbol& symbol);
    void assume_integer(const Symbol& symbol);
    void assume_nonzero(const Symbol& symbol);
    void assume_in_range(const Symbol& symbol, ExprPtr lower, ExprPtr upper);
    void assume_domain(const Symbol& symbol, Domain domain);

    // Advanced Assumptions (F9-A4)
    void assume_greater(ExprPtr lhs, ExprPtr rhs);
    void assume_greater_equal(ExprPtr lhs, ExprPtr rhs);
    void assume(ExprPtr condition);

    [[nodiscard]] bool is_real(const Symbol& symbol) const;
    [[nodiscard]] bool is_real(ExprPtr expr) const;
    [[nodiscard]] bool is_positive(const Symbol& symbol) const;
    [[nodiscard]] bool is_positive(ExprPtr expr) const;
    [[nodiscard]] bool is_nonnegative(ExprPtr expr) const;
    [[nodiscard]] bool is_negative(ExprPtr expr) const;
    [[nodiscard]] bool is_greater(ExprPtr lhs, ExprPtr rhs) const;
    [[nodiscard]] bool is_greater_equal(ExprPtr lhs, ExprPtr rhs) const;
    [[nodiscard]] bool is_nonzero(const Symbol& symbol) const;
    [[nodiscard]] bool is_nonzero(ExprPtr expr) const;
    [[nodiscard]] bool could_be_zero(const Symbol& symbol) const;
    [[nodiscard]] bool could_be_zero(ExprPtr expr) const;
    [[nodiscard]] bool is_integer(const Symbol& symbol) const;
    [[nodiscard]] bool is_integer(ExprPtr expr) const;
    [[nodiscard]] Domain get_domain(const Symbol& symbol) const;
    [[nodiscard]] std::optional<RangeAssumption> get_range(const Symbol& symbol) const;
    [[nodiscard]] Result<void> check_consistency() const;

    void update_roots(AstArena& target, std::unordered_map<ExprPtr, ExprPtr>& cache);

    // F7.0-A4.1: monotonically increasing revision counter, incremented by
    // every mutator (assume_*, assume()). Used by CASContext::simplify() to
    // detect assumption changes since the last cache fill, and invalidate
    // simplify_cache_ / integrate_cache_ on mismatch.
    //
    // Mathematical correctness: simplify(abs(x)) returns x when x is known
    // positive, but -x when x is known negative — without cache invalidation
    // on assumption change, a stale entry would corrupt the session.
    [[nodiscard]] std::uint64_t revision() const noexcept { return revision_; }

private:
    // F7.0-A4.1: bump in every mutator above (assume_*, assume).
    std::uint64_t revision_{0};

    [[nodiscard]] bool prove_relation(ExprPtr start, ExprPtr end, bool strict, std::unordered_set<const ExprNode*>& visited) const;
    [[nodiscard]] bool prove_positive_linear(ExprPtr expr) const;
    [[nodiscard]] bool prove_positive_product(const Product& prod) const;

    std::unordered_set<std::string> real_symbols_;
    std::unordered_set<std::string> positive_symbols_;
    std::unordered_set<std::string> negative_symbols_;
    std::unordered_set<std::string> integer_symbols_;
    std::unordered_set<std::string> nonzero_symbols_;
    std::unordered_map<std::string, RangeAssumption> range_symbols_;
    std::unordered_map<std::string, Domain> symbol_domains_;

    // Relation graph for deduction chains
    std::unordered_map<ExprPtr, std::vector<Relation>, ExprHash, ExprEqual> relations_;
};

}  // namespace cas::symbolic
