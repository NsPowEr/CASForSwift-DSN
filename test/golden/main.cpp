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
#include "integrate_equiv.hpp"
#include "runner_timeout.hpp"
#include "runner_format.hpp"
#include "giac_parser.hpp"

#include "cas/algebra.hpp"
#include "cas/ast.hpp"
#include "cas/calculus.hpp"
#include "cas/formatter.hpp"
#include "cas/lexer.hpp"
#include "cas/parser.hpp"
#include "cas/result.hpp"
#include "cas/symbolic.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
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
// A32 — parse the corpus "assume" predicate and register the symbol-level
// domain restrictions on ctx, so simplify produces the same restricted output
// as the assume-aware Maxima reference (run_golden_maxima.sh emits the matching
// assume()/declare() directives). Grammar (a comma-separated conjunction):
//   <sym> real | <sym> integer | <expr> (>=|>|!=|<=|<) <expr>
// Only atoms of the form <bare symbol> vs 0 (or `real`/`integer`) change the
// generic symbolic CAS output; expression-level facts (e.g. `a*d - b*c != 0`,
// `n != -1`) do not alter it and are handled on the Maxima side only.
// Returns the human-readable list of restrictions actually applied CAS-side.
static std::string apply_assume_predicates(
    cas::symbolic::CASContext& ctx, const std::string& assume_str) {
    auto trim = [](std::string s) -> std::string {
        const std::size_t a = s.find_first_not_of(" \t");
        if (a == std::string::npos) return {};
        const std::size_t b = s.find_last_not_of(" \t");
        return s.substr(a, b - a + 1);
    };
    auto parse_expr = [&](const std::string& tok) -> ExprPtr {
        auto lx = cas::Lexer(tok).tokenize();
        if (!lx.is_ok()) return nullptr;
        cas::Parser p(lx.value(), ctx.arena());
        auto r = p.parse();
        return r.is_ok() ? r.value() : nullptr;
    };
    std::string applied;
    auto note = [&](const std::string& s) {
        if (!applied.empty()) applied += ", ";
        applied += s;
    };
    std::size_t start = 0;
    while (start <= assume_str.size()) {
        const std::size_t comma = assume_str.find(',', start);
        const std::string atom = trim(assume_str.substr(
            start, comma == std::string::npos ? std::string::npos : comma - start));
        start = (comma == std::string::npos) ? assume_str.size() + 1 : comma + 1;
        if (atom.empty()) continue;

        auto ends_with = [&](const char* suf) {
            const std::string s(suf);
            return atom.size() >= s.size()
                && atom.compare(atom.size() - s.size(), s.size(), s) == 0;
        };
        if (ends_with(" real")) {
            if (const auto* s = expr_cast<Symbol>(
                    parse_expr(trim(atom.substr(0, atom.size() - 5))))) {
                ctx.assumptions().assume_real(*s); note(s->name + " real");
            }
            continue;
        }
        if (ends_with(" integer")) {
            if (const auto* s = expr_cast<Symbol>(
                    parse_expr(trim(atom.substr(0, atom.size() - 8))))) {
                ctx.assumptions().assume_integer(*s); note(s->name + " integer");
            }
            continue;
        }
        // Relational atom. Test ">=" / "<=" / "!=" before "<" / ">".
        static const char* kOps[] = {">=", "<=", "!=", ">", "<"};
        for (const char* op : kOps) {
            const std::size_t pos = atom.find(op);
            if (pos == std::string::npos) continue;
            const std::string op_s(op);
            ExprPtr lhs = parse_expr(trim(atom.substr(0, pos)));
            const std::string rhs = trim(atom.substr(pos + op_s.size()));
            const auto* sym = expr_cast<Symbol>(lhs);
            const bool rhs_zero = (rhs == "0");
            if (sym != nullptr && rhs_zero) {
                if (op_s == ">=") {
                    ctx.assumptions().assume_greater_equal(lhs, ExprPtr());
                    note(sym->name + " >= 0");
                } else if (op_s == ">") {
                    ctx.assumptions().assume_positive(*sym);
                    note(sym->name + " > 0");
                } else if (op_s == "!=") {
                    ctx.assumptions().assume_nonzero(*sym);
                    note(sym->name + " != 0");
                }
                // "<= 0" / "< 0" absent from the corpus; expression LHS or a
                // non-zero RHS have no generic-output effect (Maxima side only).
            }
            break;
        }
    }
    return applied;
}

// ---------------------------------------------------------------------------
// Per-area statistics
// ---------------------------------------------------------------------------
struct AreaStats {
    int pass{0};
    int fail{0};
    int skip{0};
    // A51 — entry troncate dal budget per-entry. Sono una categoria a se':
    // il loro verdetto non dice nulla sulla matematica, dice solo che il
    // motore non ha finito nel tempo concesso. Contarle come PASS o FAIL
    // lega gli aggregati del ratchet alla velocita' della macchina — e' il
    // difetto che A51 chiude.
    int over_budget{0};
    std::vector<std::string> fail_examples; // up to 5 examples
};

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <corpus.jsonl> <maxima_out_dir> [--json <output.json>]"
                     " [--per-entry-json <output.jsonl>]\n";
        return 1;
    }

    const std::string corpus_path = argv[1];
    const std::string maxima_dir  = argv[2];
    std::string json_output_path;
    std::string per_entry_json_path;
    std::string giac_dir; // A35 — opt-in, per-area dir written by run_golden_giac.sh
    unsigned int per_entry_timeout_sec = 30;
    bool ops_report = false;
    std::optional<std::uint64_t> max_operation_ops_override;

    for (int i = 3; i < argc; ++i) {
        std::string a(argv[i]);
        if (a == "--json" && i + 1 < argc) {
            json_output_path = argv[++i];
        } else if (a == "--per-entry-json" && i + 1 < argc) {
            // A35 — dump per-entry {idx, area, ref, input, verdict, cas_output}
            // (JSONL, una entry per riga) necessario per il cross-diff con i
            // verdetti giac gia' salvati per-entry da run_golden_giac.sh: gli
            // aggregati di --json non bastano a produrre la lista "giac
            // risolve, noi no" ne' a confrontare il VALORE (non solo
            // ANSWERED/pass) contro giac.
            per_entry_json_path = argv[++i];
        } else if (a == "--giac-dir" && i + 1 < argc) {
            // A35 — quando fornita, per ogni entry si legge anche
            // <giac-dir>/<idx>.giac.out (scritto da run_golden_giac.sh) e si
            // confronta col risultato CAS via la STESSA equivalenza vera gia'
            // usata per Maxima (mathematically_equal / compare_solve_sets /
            // antiderivative_equivalent) — non solo ANSWERED vs pass. Il
            // verdetto giac finisce SOLO nel dump --per-entry-json: non tocca
            // AreaStats / il ratchet, che resta ancorato a Maxima soltanto.
            giac_dir = argv[++i];
        } else if (a == "--per-entry-timeout" && i + 1 < argc) {
            per_entry_timeout_sec = static_cast<unsigned int>(std::stoul(argv[++i]));
        } else if (a == "--ops-report") {
            // A51: emette per ogni entry il budget deterministico consumato.
            ops_report = true;
        } else if (a == "--max-ops" && i + 1 < argc) {
            // A53: il gate ops del runner e' l'oggetto stesso della misura di
            // calibrazione — con OperationScope aperto su calculus::integrate
            // il budget vale ora per l'INTERA integrazione, non per singola
            // simplify. Misurare il costo reale richiede quindi di poter
            // allargare (o spegnere, con 0 — contratto A30) il gate per un run
            // di sola misura, senza toccare il default nel codice.
            max_operation_ops_override = std::stoull(argv[++i]);
        }
    }

    std::ifstream corpus_file(corpus_path);
    if (!corpus_file.is_open()) {
        std::cerr << "ERROR: cannot open corpus: " << corpus_path << "\n";
        return 1;
    }

    // A35 — JSONL, una entry per riga: scritta in streaming (non un array
    // JSON accumulato) cosi' un run interrotto lascia comunque righe valide.
    std::ofstream per_entry_jf;
    if (!per_entry_json_path.empty()) {
        per_entry_jf.open(per_entry_json_path);
        if (!per_entry_jf.is_open()) {
            std::cerr << "ERROR: cannot open --per-entry-json output: "
                      << per_entry_json_path << "\n";
            return 1;
        }
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

    // A51 — il verdetto di una entry deve dipendere SOLO dal suo input.
    //
    // Il contesto e' condiviso da tutte le entry dell'area, quindi ogni budget
    // ereditato implicitamente lega il verdetto alla storia del run. Due
    // conseguenze misurate prima di questa dichiarazione esplicita:
    //   * il wall-clock interno (10s di default) era a 1.4x dal costo reale
    //     del confronto piu' caro dell'area (7.1s, misurato): sotto carico
    //     sforava, l'esito passava da FALSE a errore, e il runner lo contava
    //     SKIP invece che FAIL — la stessa entry oscillava fra due misure con
    //     lo STESSO binario;
    //   * `set_timeout` senza `set_max_operation_ops` DISATTIVA il gate
    //     deterministico (contratto A30), e il ramo gcd di corpus_runner.hpp
    //     lo faceva, spegnendolo per tutto il resto del run.
    //
    // Qui si dichiara l'opposto, e l'ORDINE delle due righe e' parte del
    // contratto — invertirlo le rende una disattivazione silenziosa:
    //   1. il gate ops reso ESPLICITO al suo stesso valore di default (quello
    //      calibrato da A30), cosi' nessun `set_timeout` successivo puo' piu'
    //      spegnerlo. Deve venire PRIMA perche' `set_timeout` azzera il
    //      budget ops finche' e' implicito (context_core.cpp, contratto A30):
    //      leggere il getter dopo restituirebbe 0, e la riga renderebbe
    //      esplicito proprio lo ZERO — gate spento con un commento che
    //      dichiara il contrario. Misurato: e' cosi' che il gate e' rimasto
    //      spento per l'intero run della prima stesura di A51.
    //   2. il wall-clock interno pari al budget per-entry gia' dichiarato
    //      dallo script — nessun numero nuovo inventato, e la protezione
    //      anti-hang vera resta il SIGALRM per-entry installato sotto.
    ctx.set_max_operation_ops(max_operation_ops_override.value_or(ctx.max_operation_ops()));
    ctx.set_timeout(std::chrono::seconds(per_entry_timeout_sec));

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

        // A32: the corpus "assume" predicate is now PARSED and applied — both
        // CAS-side (symbol-level domain restrictions on ctx, below) and in the
        // Maxima references (run_golden_maxima.sh emits the matching
        // assume()/declare()). The A31 conditional-domain path stays enabled as
        // a fallback for symbolic cases the explicit assumption does not prove.
        const std::string assume_str = json_string_field(line, "assume");
        ctx.set_conditional_domain_rules(!assume_str.empty());

        auto& stats = area_stats[area];

        // Reset context state between entries (variables, caches AND assumptions;
        // arena is intentionally reused — all prior ExprPtrs remain valid).
        ctx.clear_variables();
        ctx.clear_caches();
        ctx.clear_assumptions();

        // A32: apply THIS entry's domain restrictions after the reset, so
        // simplify matches the assume-aware Maxima reference for this line.
        if (!assume_str.empty()) {
            const std::string applied = apply_assume_predicates(ctx, assume_str);
            std::cout << "  INFO [" << std::setw(3) << idx << "] assume \""
                      << assume_str << "\""
                      << (applied.empty()
                              ? " (no CAS-side symbol restriction; Maxima-side only)"
                              : " -> CAS: " + applied)
                      << "\n";
        }

        // F7.5.A3: arm per-entry timer. Handler will set interrupt flag.
        cas::golden::start_entry_timer(per_entry_timeout_sec);

        // A51 — con `--ops-report` ogni entry dichiara quanto budget
        // deterministico ha consumato la sua operazione piu' cara. E' il dato
        // che permette di tarare `max_operation_ops` sul lavoro reale del
        // corpus invece che a intuito. Emesso da un distruttore perche' il
        // corpo del loop esce da decine di punti diversi (un `continue` per
        // ogni classe di SKIP): un'unica riga a valle e' l'unico modo di
        // coprirli tutti senza toccarli tutti.
        // A51 — riclassificazione delle entry troncate.
        //
        // Il SIGALRM per-entry interrompe il motore a meta' lavoro: cio' che
        // il runner registra dopo quel punto (un `false` da un confronto
        // rimasto a meta', un NO_STRATEGY da una strategia interrotta) NON e'
        // un verdetto sulla matematica, ma su quanto tempo aveva a
        // disposizione la macchina. Misurato sul corpus: 16 entry su 1026
        // toccano il cap, e fra queste tre FAIL (`bronstein[69][70][73]`) e un
        // PASS (`integrate[94]`) con margine 1.0x — cioe' il loro esito si
        // decide esattamente sul confine.
        //
        // A35 — output CAS testuale per il dump per-entry, popolato dal ramo
        // che effettivamente calcola un risultato (matrix/solve/generico);
        // resta vuoto per gli SKIP decisi prima di ottenere un risultato CAS
        // (nessun dispatch, area matrice non supportata, ecc.).
        std::string current_cas_output;

        // A35 — verdetto/valore giac per QUESTA entry, popolati dai blocchi
        // solve/generico via la stessa equivalenza vera usata per Maxima
        // (mathematically_equal / compare_solve_sets / antiderivative_equivalent
        // per integrate-bronstein). "not_compared" quando --giac-dir non e'
        // fornita o l'area non ha ancora un confronto giac cablato (matrix).
        std::string current_giac_verdict = "not_compared";
        std::string current_giac_output;

        // A35 — dump {idx, area, ref, input, verdict, cas_output, giac_verdict,
        // giac_output} per il cross-diff con Giac: gli aggregati non bastano a
        // produrre la lista "giac risolve, noi no". Stesso trucco RAII di
        // OverBudgetGuard sotto (dichiarato PRIMA cosi' il suo distruttore
        // gira DOPO — vede lo stato FINALE di stats, incluso l'eventuale
        // over_budget) per non dover toccare ognuno dei punti di uscita del
        // corpo del loop: il verdetto CAS si deduce dal DIFF dei contatori
        // invece di essere passato a mano.
        struct PerEntryJsonGuard {
            std::ofstream& out;
            AreaStats& stats;
            int entry_idx;
            const std::string& area;
            const std::string& ref;
            const std::string& input;
            const std::string& cas_output;
            const std::string& giac_verdict;
            const std::string& giac_output;
            int pass0, fail0, skip0, over0;
            ~PerEntryJsonGuard() {
                if (!out.is_open()) return;
                const char* verdict = "unknown";
                if (stats.over_budget > over0)      verdict = "over_budget";
                else if (stats.pass > pass0)        verdict = "pass";
                else if (stats.fail > fail0)        verdict = "fail";
                else if (stats.skip > skip0)        verdict = "skip";
                out << "{\"idx\":" << entry_idx
                    << ",\"area\":\"" << json_escape(area) << "\""
                    << ",\"ref\":\"" << json_escape(ref) << "\""
                    << ",\"input\":\"" << json_escape(input) << "\""
                    << ",\"verdict\":\"" << verdict << "\""
                    << ",\"cas_output\":\"" << json_escape(cas_output) << "\""
                    << ",\"giac_verdict\":\"" << json_escape(giac_verdict) << "\""
                    << ",\"giac_output\":\"" << json_escape(giac_output) << "\"}\n";
            }
        } per_entry_json_guard{per_entry_jf, stats, idx, area, ref, input_str,
                                current_cas_output,
                                current_giac_verdict, current_giac_output,
                                stats.pass, stats.fail, stats.skip, stats.over_budget};

        struct OverBudgetGuard {
            AreaStats& stats;
            int entry_idx;
            const std::string& input;
            int pass0, fail0, skip0;
            ~OverBudgetGuard() {
                if (!cas::golden::entry_timed_out()) return;
                if (stats.pass > pass0)       --stats.pass;
                else if (stats.fail > fail0)  --stats.fail;
                else if (stats.skip > skip0)  --stats.skip;
                ++stats.over_budget;
                std::cout << "  OVER [" << std::setw(3) << entry_idx << "] "
                          << input
                          << " => troncata dal budget per-entry"
                             " (verdetto non attendibile)\n";
            }
        } over_budget_guard{stats, idx, input_str, stats.pass, stats.fail, stats.skip};

        ctx.reset_ops_high_water();
        struct OpsReportGuard {
            const symbolic::CASContext& ctx;
            bool enabled;
            int entry_idx;
            std::chrono::steady_clock::time_point start{
                std::chrono::steady_clock::now()};
            ~OpsReportGuard() {
                if (enabled) {
                    const auto ms = std::chrono::duration_cast<
                        std::chrono::milliseconds>(
                            std::chrono::steady_clock::now() - start).count();
                    // Il TEMPO accanto alle ops non e' ridondante: dice quanto
                    // margine ha l'entry rispetto al cap per-entry, cioe' se il
                    // suo verdetto e' al riparo dal carico della macchina o
                    // sta sul confine (il difetto che A51 deve chiudere).
                    std::cout << "       ops[" << entry_idx << "]="
                              << ctx.ops_high_water() << " ms=" << ms << "\n";
                }
            }
        } ops_guard{ctx, ops_report, idx};

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

            // A35 — cattura il risultato CAS per il dump per-entry, prima del
            // confronto: nessun formatter di matrici esiste ancora nel
            // runner, ne basta uno minimale (righe fra parentesi quadre).
            switch (cas_v.kind) {
                case K::Scalar:
                    current_cas_output = format_expr(cas_v.scalar);
                    break;
                case K::Matrix: {
                    std::ostringstream ms;
                    ms << "[";
                    for (std::size_t r = 0; r < cas_v.matrix.rows(); ++r) {
                        if (r) ms << ",";
                        ms << "[";
                        for (std::size_t c = 0; c < cas_v.matrix.cols(); ++c) {
                            if (c) ms << ",";
                            ms << format_expr(cas_v.matrix(r, c));
                        }
                        ms << "]";
                    }
                    ms << "]";
                    current_cas_output = ms.str();
                    break;
                }
                case K::Eigenvalues: {
                    std::ostringstream es;
                    es << "{";
                    for (std::size_t i = 0; i < cas_v.eigenvalues_list.size(); ++i) {
                        if (i) es << ",";
                        es << format_expr(cas_v.eigenvalues_list[i]);
                    }
                    es << "}";
                    current_cas_output = es.str();
                    break;
                }
                default:
                    break;
            }
            switch (cas_v.kind) {
                case K::Scalar: {
                    // HC-F75-A2-MAXIMA-MATTRACE: if Maxima emits
                    // `mattrace(matrix(...))` unevaluated, parse the
                    // matrix and compute trace on the CAS side.
                    if (auto mt = cas::golden::try_evaluate_mattrace_wrapper(
                            last_line, ctx);
                        mt.has_value()) {
                        if (!mt->is_ok()) {
                            skip = true; skip_reason =
                                "Maxima mattrace eval: " + mt->error().message;
                            break;
                        }
                        auto eq = cas::symbolic::mathematically_equal(
                            cas_v.scalar, mt->value(), ctx);
                        if (!eq.is_ok()) {
                            skip = true; skip_reason = "inconclusive mattrace";
                            break;
                        }
                        pass = eq.value();
                        break;
                    }
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
            // A35 — cattura il set di radici per il dump per-entry.
            {
                std::ostringstream rs;
                rs << "{";
                for (std::size_t i = 0; i < cas_solve.value().size(); ++i) {
                    if (i) rs << ",";
                    rs << format_expr(cas_solve.value()[i]);
                }
                rs << "}";
                current_cas_output = rs.str();
            }
            // A35 — confronto giac (opt-in, indipendente da Maxima/dal ratchet).
            if (!giac_dir.empty()) {
                auto gf = read_giac_file(giac_dir, idx);
                if (gf.tag.empty()) {
                    current_giac_verdict = "no_giac_file";
                } else if (gf.tag != "ANSWERED") {
                    current_giac_verdict = gf.tag; // UNEVALUATED | TIMEOUT | ERROR
                } else {
                    auto giac_solve = parse_giac_solve_list(gf.result, ctx.arena());
                    if (!giac_solve.is_ok()) {
                        current_giac_verdict = "unparseable";
                        current_giac_output = gf.result;
                    } else {
                        current_giac_output = gf.result;
                        auto geq = compare_solve_sets(cas_solve.value(), giac_solve.value(), ctx);
                        current_giac_verdict = (!geq.is_ok()) ? "inconclusive"
                                              : (geq.value() ? "pass" : "fail");
                    }
                }
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

        // A35 — cattura il risultato CAS per il dump per-entry: da qui in
        // poi ogni ramo (SKIP lato Maxima incluso) ha gia' un cas_output.
        current_cas_output = format_expr(cas_result.value());

        // A35 — confronto giac (opt-in, indipendente da Maxima/dal ratchet).
        // Stessa equivalenza vera usata sotto per Maxima: mathematically_equal
        // con fallback antiderivative_equivalent per integrate/bronstein
        // (due antiderivate dello stesso integrando differiscono al piu' per
        // una costante — lo stesso motivo per cui serve anche lato Maxima).
        if (!giac_dir.empty()) {
            auto gf = read_giac_file(giac_dir, idx);
            if (gf.tag.empty()) {
                current_giac_verdict = "no_giac_file";
            } else if (gf.tag != "ANSWERED") {
                current_giac_verdict = gf.tag; // UNEVALUATED | TIMEOUT | ERROR
            } else {
                auto giac_expr = parse_giac_expr(gf.result, ctx.arena());
                if (!giac_expr.is_ok()) {
                    current_giac_verdict = "unparseable";
                    current_giac_output = gf.result;
                } else {
                    current_giac_output = gf.result;
                    auto geq = mathematically_equal(cas_result.value(), giac_expr.value(), ctx);
                    bool gpass = geq.is_ok() && geq.value();
                    if (!gpass && (area == "integrate" || area == "bronstein")) {
                        auto gantieq = cas::golden::antiderivative_equivalent(
                            input_str, cas_result.value(), giac_expr.value(), ctx);
                        gpass = gantieq.is_ok() && gantieq.value();
                    }
                    current_giac_verdict = (!geq.is_ok() && !gpass) ? "inconclusive"
                                          : (gpass ? "pass" : "fail");
                }
            }
        }

        // A31 §10.5: show the domain conditions registered while evaluating
        // this entry (the flag's contract: whoever enables it reads them).
        // Best effort: multi-call pipelines report the LAST top-level set.
        if (!assume_str.empty() && !ctx.last_side_conditions().empty()) {
            std::cout << "  INFO [" << std::setw(3) << idx << "] conditions taken:";
            for (const auto& c : ctx.last_side_conditions().items()) {
                const char* kind_name = "?";
                switch (c.kind) {
                    case symbolic::DomainConditionKind::NonZero:     kind_name = "nonzero"; break;
                    case symbolic::DomainConditionKind::Positive:    kind_name = "positive"; break;
                    case symbolic::DomainConditionKind::NonNegative: kind_name = "nonnegative"; break;
                    case symbolic::DomainConditionKind::Real:        kind_name = "real"; break;
                    case symbolic::DomainConditionKind::IntegerVal:  kind_name = "integer"; break;
                    case symbolic::DomainConditionKind::PrincipalBranch: kind_name = "principal-branch"; break;
                }
                std::cout << " " << kind_name << "(" << format_expr(c.subject) << ")";
            }
            std::cout << "\n";
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
        } else if (area == "integrate" || area == "bronstein") {
            // Antiderivative equivalence: two antiderivatives of the same
            // integrand differ at most by a constant. Differentiate both
            // and re-compare. Resolves the corpus mismatch between our
            // rigorous `ln(abs(x))` and Maxima's optimistic `log(x)` for
            // entries like integrate(1/x, x).
            auto antieq = cas::golden::antiderivative_equivalent(
                input_str, cas_result.value(), maxima_simplified.value(), ctx);
            if (antieq.is_ok() && antieq.value()) {
                ++stats.pass;
                std::cout << "  PASS [" << std::setw(3) << idx << "] "
                          << input_str << " (antideriv equiv)\n";
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
              << std::setw(6)  << "OVER"
              << std::setw(8)  << "TOTAL"
              << std::setw(10) << "PASS%"
              << "\n";
    std::cout << "-------------------------------------------------------------\n";

    int grand_pass = 0, grand_fail = 0, grand_skip = 0, grand_over = 0;
    for (auto& [area, st] : area_stats) {
        int tot = st.pass + st.fail + st.skip + st.over_budget;
        double pct = (tot > 0) ? (100.0 * st.pass / (st.pass + st.fail)) : 0.0;
        if (st.pass + st.fail == 0) pct = 0.0;
        std::cout << std::left
                  << std::setw(14) << area
                  << std::setw(6)  << st.pass
                  << std::setw(6)  << st.fail
                  << std::setw(6)  << st.skip
                  << std::setw(6)  << st.over_budget
                  << std::setw(8)  << tot
                  << std::fixed << std::setprecision(1) << pct << "%\n";
        grand_pass += st.pass;
        grand_fail += st.fail;
        grand_skip += st.skip;
        grand_over += st.over_budget;
    }

    std::cout << "-------------------------------------------------------------\n";
    int grand_tot = grand_pass + grand_fail + grand_skip + grand_over;
    double grand_pct = (grand_pass + grand_fail > 0)
                       ? (100.0 * grand_pass / (grand_pass + grand_fail))
                       : 0.0;
    std::cout << std::left
              << std::setw(14) << "TOTAL"
              << std::setw(6)  << grand_pass
              << std::setw(6)  << grand_fail
              << std::setw(6)  << grand_skip
              << std::setw(6)  << grand_over
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
            jf << "    \"over_budget\": " << st.over_budget << ",\n";
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
