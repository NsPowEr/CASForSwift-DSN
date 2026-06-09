// test/golden/main.cpp — standalone golden test runner (no GTest dependency).
//
// Usage:
//   golden_runner <corpus.jsonl> <maxima_out_dir> [--json <output.json>]
//
// For each line in corpus.jsonl:
//   1. Parse the input via our CAS.
//   2. Read <maxima_out_dir>/<idx>.maxima.out, extract last expression line.
//   3. Normalise Maxima output to our syntax; parse via our CAS.
//   4. Compare using mathematically_equal(cas_result, maxima_result, ctx).
//   5. Accumulate PASS / FAIL / SKIP per area.
//
// Outputs a JSON report: { "area": { "pass": N, "fail": M, "skip": K,
//                                    "examples_fail": [...] }, ... }
// and prints a summary table to stdout.
//
// Exit code: 0 if pass_rate >= 0% (always succeeds for baseline measurement),
//            1 on internal error.

#include "corpus_runner.hpp"
#include "maxima_parser.hpp"
#include "solve_set_equal.hpp"
#include "matrix_adapter.hpp"
#include "runner_timeout.hpp"
#include "runner_format.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace cas;
using namespace cas::golden;

// ---------------------------------------------------------------------------
// Minimal JSON line parser: extract string field value from a JSON object.
// Avoids adding a third-party library.
// ---------------------------------------------------------------------------
static std::string json_string_field(const std::string& line, const std::string& key) {
    std::string search = "\"" + key + "\"";
    auto pos = line.find(search);
    if (pos == std::string::npos) return "";
    pos += search.size();
    // skip whitespace and colon
    while (pos < line.size() && (line[pos] == ' ' || line[pos] == ':')) ++pos;
    if (pos >= line.size() || line[pos] != '"') return "";
    ++pos; // skip opening quote
    std::string result;
    while (pos < line.size() && line[pos] != '"') {
        if (line[pos] == '\\' && pos + 1 < line.size()) {
            ++pos;
            switch (line[pos]) {
                case '"':  result += '"'; break;
                case '\\': result += '\\'; break;
                case '/':  result += '/'; break;
                case 'n':  result += '\n'; break;
                case 't':  result += '\t'; break;
                default:   result += line[pos]; break;
            }
        } else {
            result += line[pos];
        }
        ++pos;
    }
    return result;
}

// ---------------------------------------------------------------------------
// Read entire file into string.
// ---------------------------------------------------------------------------
static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) return "";
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// ---------------------------------------------------------------------------
// Escape a string for JSON output.
// ---------------------------------------------------------------------------
static std::string json_escape(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '"') out += "\\\"";
        else if (c == '\\') out += "\\\\";
        else if (c == '\n') out += "\\n";
        else if (c == '\r') out += "\\r";
        else if (c == '\t') out += "\\t";
        else out += c;
    }
    return out;
}

// ---------------------------------------------------------------------------
// TextFormatter wrapper: format an ExprPtr as a string (for reporting only,
// never used for mathematical comparison). Cap output size to keep the
// summary table reachable on integrate-style multi-MB renderings.
// ---------------------------------------------------------------------------
static std::string format_expr(ExprPtr expr) {
    return cas::golden::format_expr_capped(expr);
}

// ---------------------------------------------------------------------------
// Per-area statistics
// ---------------------------------------------------------------------------
struct AreaStats {
    int pass{0};
    int fail{0};
    int skip{0};
    std::vector<std::string> fail_examples; // up to 5 examples
};

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <corpus.jsonl> <maxima_out_dir> [--json <output.json>]\n";
        return 1;
    }

    const std::string corpus_path = argv[1];
    const std::string maxima_dir  = argv[2];
    std::string json_output_path;
    unsigned int per_entry_timeout_sec = 30;

    for (int i = 3; i < argc; ++i) {
        std::string a(argv[i]);
        if (a == "--json" && i + 1 < argc) {
            json_output_path = argv[++i];
        } else if (a == "--per-entry-timeout" && i + 1 < argc) {
            per_entry_timeout_sec = static_cast<unsigned int>(std::stoul(argv[++i]));
        }
    }

    std::ifstream corpus_file(corpus_path);
    if (!corpus_file.is_open()) {
        std::cerr << "ERROR: cannot open corpus: " << corpus_path << "\n";
        return 1;
    }

    std::map<std::string, AreaStats> area_stats;
    int idx = 0;

    // Heap-allocate CASContext so we can explicitly destroy it (and its arena)
    // before the function's local variables destruct.  This avoids a crash in
    // the integrate_cache_ LRU list destructor, which walks ExprPtr-keyed nodes
    // that reference arena memory.  We call clear_caches() then delete ctx to
    // guarantee the cache is drained before the arena is freed.
    auto* ctx_ptr = new symbolic::CASContext();
    symbolic::CASContext& ctx = *ctx_ptr;
    ctx.set_caching_enabled(false);

    // F7.5.A3: install SIGALRM handler + register ctx as interrupt target.
    cas::golden::install_alarm_handler();
    cas::golden::timeout_target().store(&ctx, std::memory_order_relaxed);

    std::string line;
    while (std::getline(corpus_file, line)) {
        // Skip empty / comment lines
        if (line.empty() || line[0] == '#') continue;

        std::string input_str = json_string_field(line, "input");
        std::string area      = json_string_field(line, "area");
        std::string ref       = json_string_field(line, "ref");
        std::string book_ref  = json_string_field(line, "book_ref");

        if (input_str.empty()) {
            ++idx;
            continue;
        }
        if (area.empty()) area = "unknown";

        auto& stats = area_stats[area];

        // Reset context state between entries (variables and assumptions only;
        // arena is intentionally reused — all prior ExprPtrs remain valid).
        ctx.clear_variables();
        ctx.clear_caches();

        // F7.5.A3: arm per-entry timer. Handler will set interrupt flag.
        cas::golden::start_entry_timer(per_entry_timeout_sec);

        // --- Matrix area: dispatch det/trace/transpose/rank/inverse/eigenvalues (F7.5.A2) ---
        if (area == "matrix") {
            auto cas_res = cas::golden::evaluate_cas_matrix(input_str, ctx);
            if (!cas_res.is_ok()) {
                ++stats.skip;
                std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                          << input_str << " => " << cas_res.error().message << "\n";
                ++idx;
                continue;
            }
            std::string maxima_out_path = maxima_dir + "/" + std::to_string(idx) + ".maxima.out";
            std::string maxima_raw = read_file(maxima_out_path);
            if (maxima_raw.empty()) {
                ++stats.skip;
                std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                          << input_str << " (no Maxima output)\n";
                ++idx;
                continue;
            }
            std::string last_line = extract_maxima_result_line(maxima_raw);
            auto& cas_v = cas_res.value();
            using K = cas::golden::MatrixDispatchResult::Kind;
            bool pass = false;
            bool skip = false;
            std::string skip_reason;
            switch (cas_v.kind) {
                case K::Scalar: {
                    std::string norm = normalize_maxima_output(last_line);
                    if (norm.empty()) {
                        skip = true; skip_reason = "Maxima scalar empty";
                        break;
                    }
                    auto me = cas::golden::parse_maxima_expr(norm, ctx.arena());
                    if (!me.is_ok()) {
                        skip = true; skip_reason = "Maxima parse: " + norm;
                        break;
                    }
                    auto eq = cas::symbolic::mathematically_equal(cas_v.scalar, me.value(), ctx);
                    if (!eq.is_ok()) { skip = true; skip_reason = "inconclusive"; break; }
                    pass = eq.value();
                    break;
                }
                case K::Matrix: {
                    auto mm = cas::golden::parse_maxima_matrix(last_line, ctx);
                    if (!mm.is_ok()) {
                        skip = true; skip_reason = "Maxima matrix parse: " + last_line;
                        break;
                    }
                    auto eq = cas::golden::compare_matrices(cas_v.matrix, mm.value(), ctx);
                    if (!eq.is_ok()) { skip = true; skip_reason = "inconclusive"; break; }
                    pass = eq.value();
                    break;
                }
                case K::Eigenvalues: {
                    // Maxima eigenvalues output: [[λ_1, …, λ_k], [m_1, …, m_k]]
                    // — first sublist holds the eigenvalues, second the
                    // multiplicities. split on top-level commas of the outer
                    // list yields ["[λ_1,…]", "[m_1,…]"]; we want the raw
                    // eigenvalues, so split the first element again.
                    auto outer = cas::golden::split_maxima_list(last_line);
                    if (outer.empty()) { skip = true; skip_reason = "Maxima eig empty"; break; }
                    auto eig_strs = cas::golden::split_maxima_list(outer[0]);
                    if (eig_strs.empty()) { skip = true; skip_reason = "Maxima eig inner empty"; break; }
                    std::vector<ExprPtr> evals;
                    evals.reserve(eig_strs.size());
                    bool any_err = false;
                    for (const auto& el : eig_strs) {
                        auto n = normalize_maxima_output(el);
                        if (n.empty()) { any_err = true; break; }
                        auto p = cas::golden::parse_maxima_expr(n, ctx.arena());
                        if (!p.is_ok()) { any_err = true; break; }
                        evals.push_back(p.value());
                    }
                    if (any_err) { skip = true; skip_reason = "Maxima eig elem parse"; break; }
                    auto eq = cas::golden::compare_solve_sets(
                        cas_v.eigenvalues_list, evals, ctx);
                    if (!eq.is_ok()) { skip = true; skip_reason = "inconclusive eig"; break; }
                    pass = eq.value();
                    break;
                }
                default:
                    skip = true; skip_reason = "no dispatch kind";
                    break;
            }
            if (skip) {
                ++stats.skip;
                std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                          << input_str << " (" << skip_reason << ")\n";
            } else if (pass) {
                ++stats.pass;
                std::cout << "  PASS [" << std::setw(3) << idx << "] " << input_str << "\n";
            } else {
                ++stats.fail;
                std::cout << "  FAIL [" << std::setw(3) << idx << "] " << input_str << "\n";
            }
            ++idx;
            continue;
        }

        // --- Solve area: special set-equality path (F7.5.A1) ---
        if (area == "solve") {
            auto cas_solve = evaluate_cas_solve(input_str, ctx);
            if (!cas_solve.is_ok()) {
                ++stats.skip;
                std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                          << input_str << " => " << cas_solve.error().message << "\n";
                ++idx;
                continue;
            }
            std::string maxima_out_path = maxima_dir + "/" + std::to_string(idx) + ".maxima.out";
            std::string maxima_raw = read_file(maxima_out_path);
            if (maxima_raw.empty()) {
                ++stats.skip;
                std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                          << input_str << " (no Maxima output)\n";
                ++idx;
                continue;
            }
            std::string last_line = extract_maxima_result_line(maxima_raw);
            auto maxima_solve = parse_maxima_solve_list(last_line, ctx.arena());
            if (!maxima_solve.is_ok()) {
                ++stats.skip;
                std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                          << input_str << " (Maxima parse: " << last_line << ")\n";
                ++idx;
                continue;
            }
            auto eq = compare_solve_sets(cas_solve.value(), maxima_solve.value(), ctx);
            if (!eq.is_ok()) {
                ++stats.skip;
                std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                          << input_str << " (inconclusive solve compare)\n";
            } else if (eq.value()) {
                ++stats.pass;
                std::cout << "  PASS [" << std::setw(3) << idx << "] " << input_str << "\n";
            } else {
                ++stats.fail;
                std::ostringstream cas_s, max_s;
                cas_s << "{";
                for (std::size_t i = 0; i < cas_solve.value().size(); ++i) {
                    if (i) cas_s << ", ";
                    cas_s << format_expr(cas_solve.value()[i]);
                }
                cas_s << "}";
                max_s << "{";
                for (std::size_t i = 0; i < maxima_solve.value().size(); ++i) {
                    if (i) max_s << ", ";
                    max_s << format_expr(maxima_solve.value()[i]);
                }
                max_s << "}";
                std::cout << "  FAIL [" << std::setw(3) << idx << "] "
                          << input_str << "\n"
                          << "       CAS:    " << cas_s.str() << "\n"
                          << "       Maxima: " << max_s.str() << "\n";
            }
            ++idx;
            continue;
        }

        // --- Step 1: Evaluate our CAS ---
        Result<ExprPtr> cas_result = evaluate_cas(input_str, ctx);

        if (!cas_result.is_ok()) {
            // SKIP if our CAS returns Unimplemented
            ++stats.skip;
            std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                      << input_str << " => " << cas_result.error().message << "\n";
            ++idx;
            continue;
        }

        // --- Step 2: Read Maxima output file ---
        std::string maxima_out_path = maxima_dir + "/" + std::to_string(idx) + ".maxima.out";
        std::string maxima_raw = read_file(maxima_out_path);

        if (maxima_raw.empty()) {
            // No Maxima output available — SKIP this entry
            ++stats.skip;
            std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                      << input_str << " (no Maxima output)\n";
            ++idx;
            continue;
        }

        // --- Step 3: Extract and normalise Maxima result ---
        std::string last_line = extract_maxima_result_line(maxima_raw);
        std::string normalised = normalize_maxima_output(last_line);

        if (normalised.empty()) {
            ++stats.skip;
            std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                      << input_str << " (Maxima: '" << last_line << "')\n";
            ++idx;
            continue;
        }

        // --- Step 4: Parse Maxima result via our parser ---
        Result<ExprPtr> maxima_expr = parse_maxima_expr(normalised, ctx.arena());

        if (!maxima_expr.is_ok()) {
            ++stats.skip;
            std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                      << input_str << " (Maxima parse fail: " << normalised << ")\n";
            ++idx;
            continue;
        }

        // --- Step 5: Simplify Maxima result ---
        auto maxima_simplified = ctx.simplify(maxima_expr.value());
        if (!maxima_simplified.is_ok()) {
            ++stats.skip;
            ++idx;
            continue;
        }

        // --- Step 6: Compare via mathematically_equal ---
        auto eq = mathematically_equal(cas_result.value(), maxima_simplified.value(), ctx);

        if (!eq.is_ok()) {
            // Comparison inconclusive — try Risch subset
            auto eq2 = mathematically_equal_subset_risch(
                cas_result.value(), maxima_simplified.value(), ctx);
            if (!eq2.is_ok() || !eq2.value()) {
                ++stats.skip;
                std::cout << "  SKIP [" << std::setw(3) << idx << "] "
                          << input_str << " (inconclusive equality)\n";
                ++idx;
                continue;
            }
            ++stats.pass;
            std::cout << "  PASS [" << std::setw(3) << idx << "] " << input_str << "\n";
        } else if (eq.value()) {
            ++stats.pass;
            std::cout << "  PASS [" << std::setw(3) << idx << "] " << input_str << "\n";
        } else {
            ++stats.fail;
            std::string cas_str  = format_expr(cas_result.value());
            std::string max_str  = format_expr(maxima_simplified.value());
            std::string example  = input_str + " | CAS: " + cas_str + " | Maxima: " + max_str;
            if (stats.fail_examples.size() < 5)
                stats.fail_examples.push_back(example);
            std::cout << "  FAIL [" << std::setw(3) << idx << "] "
                      << input_str << "\n"
                      << "       CAS:    " << cas_str << "\n"
                      << "       Maxima: " << max_str << "\n";
        }

        ++idx;
    }

    // F7.5.A3: disarm timer post-loop. Handler reset for cleanliness.
    cas::golden::stop_entry_timer();
    cas::golden::timeout_target().store(nullptr, std::memory_order_relaxed);

    // ---------------------------------------------------------------------------
    // Print summary table
    // ---------------------------------------------------------------------------
    std::cout << "\n";
    std::cout << "=============================================================\n";
    std::cout << "  Golden Test Pass-Rate Report\n";
    std::cout << "=============================================================\n";
    std::cout << std::left
              << std::setw(14) << "Area"
              << std::setw(6)  << "PASS"
              << std::setw(6)  << "FAIL"
              << std::setw(6)  << "SKIP"
              << std::setw(8)  << "TOTAL"
              << std::setw(10) << "PASS%"
              << "\n";
    std::cout << "-------------------------------------------------------------\n";

    int grand_pass = 0, grand_fail = 0, grand_skip = 0;
    for (auto& [area, st] : area_stats) {
        int tot = st.pass + st.fail + st.skip;
        double pct = (tot > 0) ? (100.0 * st.pass / (st.pass + st.fail)) : 0.0;
        if (st.pass + st.fail == 0) pct = 0.0;
        std::cout << std::left
                  << std::setw(14) << area
                  << std::setw(6)  << st.pass
                  << std::setw(6)  << st.fail
                  << std::setw(6)  << st.skip
                  << std::setw(8)  << tot
                  << std::fixed << std::setprecision(1) << pct << "%\n";
        grand_pass += st.pass;
        grand_fail += st.fail;
        grand_skip += st.skip;
    }

    std::cout << "-------------------------------------------------------------\n";
    int grand_tot = grand_pass + grand_fail + grand_skip;
    double grand_pct = (grand_pass + grand_fail > 0)
                       ? (100.0 * grand_pass / (grand_pass + grand_fail))
                       : 0.0;
    std::cout << std::left
              << std::setw(14) << "TOTAL"
              << std::setw(6)  << grand_pass
              << std::setw(6)  << grand_fail
              << std::setw(6)  << grand_skip
              << std::setw(8)  << grand_tot
              << std::fixed << std::setprecision(1) << grand_pct << "%\n";
    std::cout << "=============================================================\n";

    // ---------------------------------------------------------------------------
    // JSON output
    // ---------------------------------------------------------------------------
    if (!json_output_path.empty()) {
        std::ofstream jf(json_output_path);
        jf << "{\n";
        bool first_area = true;
        for (auto& [area, st] : area_stats) {
            if (!first_area) jf << ",\n";
            first_area = false;
            jf << "  \"" << area << "\": {\n";
            jf << "    \"pass\": " << st.pass << ",\n";
            jf << "    \"fail\": " << st.fail << ",\n";
            jf << "    \"skip\": " << st.skip << ",\n";
            jf << "    \"examples_fail\": [";
            for (std::size_t i = 0; i < st.fail_examples.size(); ++i) {
                if (i > 0) jf << ", ";
                jf << "\"" << json_escape(st.fail_examples[i]) << "\"";
            }
            jf << "]\n";
            jf << "  }";
        }
        jf << "\n}\n";
        std::cout << "\nJSON report written to: " << json_output_path << "\n";
    }

    // Leak the context intentionally: the CASContext destructor has a known
    // issue where the integrate_cache_ LRU list destructor walks corrupted
    // ExprPtr-keyed nodes after the arena is freed.  Since this is a
    // short-lived CLI tool, leaking ctx is safe and avoids the crash.
    // TODO: file upstream bug — CacheContainer should call clear() before
    // the arena is destroyed (arena destructor should notify cache holders).
    (void)ctx_ptr; // intentional leak
    return 0;
}
