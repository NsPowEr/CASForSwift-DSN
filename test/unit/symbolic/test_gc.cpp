#include "cas/ast.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"
#include <gtest/gtest.h>
#include <vector>

namespace cas::symbolic {
namespace {

Result<ExprPtr> parse_expr(const std::string& input, AstArena& arena) {
    auto tokens = Lexer(input).tokenize();
    if (tokens.is_error()) return fail<ExprPtr>(tokens.error());
    Parser parser(tokens.value(), arena);
    return parser.parse();
}

TEST(CASContextGCTest, CollectsGarbageAndPreservesRoots) {
    CASContext ctx;
    
    // 1. Definiamo una variabile (radice interna)
    auto x_expr = parse_expr("x^2 + y", ctx.arena()).value();
    ctx.define(Symbol("f"), x_expr);
    
    // 2. Creiamo un root esterno
    auto external_root = parse_expr("z + 1", ctx.arena()).value();
    auto old_external_ptr = external_root.get();
    
    // 3. Eseguiamo GC
    std::vector<ExprPtr*> roots = { &external_root };
    ctx.collect_garbage(roots);
    
    // 4. Verifichiamo che i puntatori siano cambiati (arena nuova)
    EXPECT_NE(external_root.get(), old_external_ptr);
    
    // 5. Verifichiamo che il contenuto sia preservato
    auto f_lookup = ctx.lookup(Symbol("f"));
    ASSERT_TRUE(f_lookup.has_value());
    EXPECT_NE(f_lookup->get(), x_expr.get()); // Deve essere una copia
    
    // 6. Verifichiamo l'interning nella nuova arena
    auto one1 = ctx.arena().make<IntegerLit>(BigInt(1));
    auto one2 = ctx.arena().make<IntegerLit>(BigInt(1));
    EXPECT_EQ(one1.get(), one2.get());
}

TEST(CASContextGCTest, MemoryPressureStressTest) {
    CASContext ctx;
    
    // Eseguiamo 10.000 semplificazioni che creano molti nodi temporanei
    for (int i = 0; i < 10000; ++i) {
        // Espressione che genera molti nodi intermedi
        auto expr_res = parse_expr("x + 1 - x + 2 - 2 + y - y", ctx.arena());
        ASSERT_TRUE(expr_res.is_ok());
        auto res = ctx.simplify(expr_res.value());
        ASSERT_TRUE(res.is_ok());
        
        // Ogni 100 iterazioni puliamo
        if (i % 100 == 0) {
            ctx.collect_garbage();
        }
    }
    
    // Se la GC non funzionasse, l'arena avrebbe > 10.000 nodi.
    // Con la GC, deve contenere solo eventuali costanti internate e poco altro.
    // x + 1 - x + 2 - 2 + y - y si semplifica a 0.
    EXPECT_LT(ctx.arena().size(), 100U); 
}

} // namespace
} // namespace cas::symbolic
