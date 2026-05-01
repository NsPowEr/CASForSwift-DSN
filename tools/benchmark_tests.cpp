#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/ast_debug.hpp"
#include "cas/calculus.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/symbolic.hpp"

#include <chrono>
#include <iostream>
#include <string>

using namespace cas;
using namespace cas::symbolic;

static Result<ExprPtr> parse(const std::string& s, AstArena& arena) {
    auto tok = Lexer(s).tokenize();
    if (tok.is_error()) return fail<ExprPtr>(tok.error());
    return Parser(tok.value(), arena).parse();
}

static std::string fmt(ExprPtr e) {
    if (!e) return "(null)";
    auto s = cas::formatter::TextFormatter{}.format(e);
    if (s.empty()) s = cas::debug_print(e);
    return s;
}

template<typename Fn>
static void run_test(int n, const char* label, Fn fn) {
    auto t0 = std::chrono::steady_clock::now();
    fn();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - t0).count();
    std::cout << "\n=== TEST " << n << ": " << label << " [" << ms << "ms] ===" << std::endl;
}

int main() {
    // ── T1: integrate(1/(x^8+1), x) ──
    {
        auto t0 = std::chrono::steady_clock::now();
        CASContext ctx;
        auto e = parse("1/(x^8+1)", ctx.arena());
        auto r = e.is_ok() ? calculus::integrate(e.value(), Symbol("x"), ctx) : fail<ExprPtr>(e.error());
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-t0).count();
        std::cout << "\n=== TEST 1: integrate(1/(x^8+1), x) [" << ms << "ms] ===\n";
        if (r.is_error()) std::cout << "FAIL: " << r.error().message << "\n";
        else std::cout << "OK: " << fmt(r.value()) << "\n";
    }

    // ── T2: factor_over_integers(x^4+1, x) ──
    {
        auto t0 = std::chrono::steady_clock::now();
        CASContext ctx;
        auto e = parse("x^4+1", ctx.arena());
        auto ms_fn = [&]() { return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-t0).count(); };
        std::cout << "\n=== TEST 2: factor_over_integers(x^4+1, x) [";
        if (e.is_error()) { std::cout << "0ms] ===\nFAIL parse\n"; }
        else {
            auto r = algebra::factor_over_integers(e.value(), Symbol("x"), ctx);
            std::cout << ms_fn() << "ms] ===\n";
            if (r.is_error()) {
                std::cout << "FAIL: " << r.error().message << "\n";
            } else {
                const auto& f = r.value();
                std::cout << "content=" << fmt(f.content) << " factors=" << f.factors.size() << "\n";
                for (const auto& fac : f.factors)
                    std::cout << "  mult=" << fac.multiplicity << " poly=" << fmt(fac.factor) << "\n";
                if (f.factors.empty()) std::cout << "  (irreducible over Z — expected for x^4+1)\n";
            }
        }
    }

    // ── T3: limit((cos(x)-exp(-x^2/2))/x^4, x->0) ──
    {
        auto t0 = std::chrono::steady_clock::now();
        CASContext ctx;
        auto e = parse("(cos(x) - exp(-x^2/2)) / x^4", ctx.arena());
        auto pt = parse("0", ctx.arena());
        std::cout << "\n=== TEST 3: limit((cos(x)-exp(-x^2/2))/x^4, x, 0) [";
        if (e.is_error() || pt.is_error()) { std::cout << "0ms] ===\nFAIL parse\n"; }
        else {
            auto r = calculus::limit(e.value(), Symbol("x"), pt.value(), LimitDirection::Both, ctx);
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-t0).count();
            std::cout << ms << "ms] ===\n";
            if (r.is_error()) std::cout << "FAIL: " << r.error().message << "\n";
            else std::cout << "OK: " << fmt(r.value()) << "\n";
        }
    }

    // ── T4: sqrt(-1)*sqrt(-1) - sqrt(1*1) ──
    {
        auto t0 = std::chrono::steady_clock::now();
        CASContext ctx;
        // Step by step
        auto e1 = parse("sqrt(-1)*sqrt(-1)", ctx.arena());
        auto e2 = parse("sqrt((-1)*(-1))", ctx.arena());
        auto e_full = parse("sqrt(-1)*sqrt(-1) - sqrt((-1)*(-1))", ctx.arena());
        std::cout << "\n=== TEST 4: branch cuts sqrt(-1)*sqrt(-1) - sqrt((-1)*(-1)) ===\n";
        if (e1.is_ok()) {
            auto s1 = ctx.simplify(e1.value());
            std::cout << "  sqrt(-1)*sqrt(-1) = " << (s1.is_ok() ? fmt(s1.value()) : "ERR:"+s1.error().message) << "\n";
        }
        if (e2.is_ok()) {
            auto s2 = ctx.simplify(e2.value());
            std::cout << "  sqrt(1) = " << (s2.is_ok() ? fmt(s2.value()) : "ERR:"+s2.error().message) << "\n";
        }
        if (e_full.is_ok()) {
            auto sf = ctx.simplify(e_full.value());
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-t0).count();
            std::cout << "  full [" << ms << "ms]: " << (sf.is_ok() ? fmt(sf.value()) : "ERR:"+sf.error().message) << "\n";
            std::cout << "  debug: " << (sf.is_ok() ? debug_print(sf.value()) : "err") << "\n";
        }
    }

    // ── T5: det(Vandermonde 4x4) ──
    {
        auto t0 = std::chrono::steady_clock::now();
        CASContext ctx;
        auto e = parse("det([[1,a,a^2,a^3],[1,b,b^2,b^3],[1,c,c^2,c^3],[1,d,d^2,d^3]])", ctx.arena());
        std::cout << "\n=== TEST 5: det(Vandermonde 4x4) ===\n";
        if (e.is_error()) { std::cout << "FAIL parse: " << e.error().message << "\n"; }
        else {
            auto s = ctx.simplify(e.value());
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-t0).count();
            std::cout << "  simplify [" << ms << "ms]: ";
            if (s.is_error()) std::cout << "ERR: " << s.error().message << "\n";
            else {
                std::cout << fmt(s.value()) << "\n";
                std::cout << "  debug: " << debug_print(s.value()) << "\n";
            }
        }
    }

    return 0;
}
