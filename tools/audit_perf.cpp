
#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>
#include "cas/symbolic.hpp"
#include "cas/ast.hpp"

using namespace cas;
using namespace cas::symbolic;

void benchmark_arena_gc() {
    std::cout << "--- Arena Compaction Benchmark ---" << std::endl;
    CASContext ctx;
    auto& arena = ctx.arena();

    // 1. Generate 1,000,000 transient nodes
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000000; ++i) {
        // Create transient nodes (not stored anywhere)
        (void)arena.make<Binary>(BinaryOp::Add, 
            arena.make<IntegerLit>(BigInt(i)), 
            arena.make<IntegerLit>(BigInt(i + 1)));
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time to generate 1,000,000 transient nodes: " << duration.count() << "ms" << std::endl;
    std::cout << "Nodes in arena: " << arena.size() << std::endl;

    // 2. Define some live roots
    ExprPtr live1 = arena.make<Symbol>("x");
    ExprPtr live2 = arena.make<Binary>(BinaryOp::Add, live1, arena.make<IntegerLit>(BigInt(42)));
    ctx.define(Symbol("x_val"), live2);

    // 3. Perform Garbage Collection
    start = std::chrono::high_resolution_clock::now();
    ctx.collect_garbage({});
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "GC (collect_garbage) duration: " << duration.count() << "ms" << std::endl;
    std::cout << "Nodes in arena after GC: " << ctx.arena().size() << std::endl;

    if (duration.count() > 100) {
        std::cout << "[CRITICAL] GC latency exceeded 100ms!" << std::endl;
    } else {
        std::cout << "[OK] GC latency within limits." << std::endl;
    }
}

void benchmark_hash_consing() {
    std::cout << "\n--- Hash-Consing Scalability Benchmark ---" << std::endl;
    CASContext ctx;
    auto& arena = ctx.arena();

    const int N = 100000;
    std::cout << "Inserting " << N << " unique small expressions..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        (void)arena.make<Binary>(BinaryOp::Mul, 
            arena.make<Symbol>("v" + std::to_string(i)), 
            arena.make<IntegerLit>(BigInt(i)));
        
        if (i > 0 && i % 20000 == 0) {
             auto now = std::chrono::high_resolution_clock::now();
             auto lap = std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
             std::cout << "  " << i << " nodes: " << lap.count() << "ms" << std::endl;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto total = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Total time: " << total.count() << "ms" << std::endl;
    std::cout << "Average time per node: " << (double)total.count() / N << "ms" << std::endl;

    // Measure time for collision-prone or complex expressions
    std::cout << "Testing interning of deeply nested expressions..." << std::endl;
    ExprPtr deep = arena.make<Symbol>("x");
    for (int i = 0; i < 500; ++i) {
        deep = arena.make<Binary>(BinaryOp::Add, deep, arena.make<IntegerLit>(BigInt(1)));
    }

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000; ++i) {
        (void)arena.make<Binary>(BinaryOp::Add, deep, arena.make<IntegerLit>(BigInt(1)));
    }
    end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time to re-intern 1,000 deep expressions: " << duration.count() << "ms" << std::endl;
}

int main() {
    try {
        benchmark_arena_gc();
        benchmark_hash_consing();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
