#pragma once
// bipartite_matching.hpp — Hopcroft-Karp maximum bipartite matching.
//
// Reference: Hopcroft & Karp (1973), "An n^{5/2} algorithm for maximum
// matchings in bipartite graphs". Standard layered-BFS + DFS-augment
// formulation, complexity O(E · sqrt(V)).
//
// Designed for the golden-runner solve adapter: |V_left| ≤ 20,
// |V_right| ≤ 20 (corpus has ≤ 6 roots per entry). Edge construction
// cost dominates (each edge requires a `mathematically_equal` call) —
// the matching algorithm itself is negligible.
//
// Generic over edge predicate via std::function so callers don't have
// to materialise the adjacency list eagerly. Edges are memoised
// internally to avoid redundant predicate evaluations during DFS retries.

#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cas::golden {

// Find a maximum matching in a bipartite graph with `n_left` vertices
// on the left and `n_right` vertices on the right.
//
// `has_edge(i, j)` must return true iff left-vertex `i` matches
// right-vertex `j`. The function is queried lazily and each (i, j)
// pair is queried at most once.
//
// Returns the size of the maximum matching.
inline std::size_t hopcroft_karp_max_matching(
    std::size_t n_left,
    std::size_t n_right,
    const std::function<bool(std::size_t, std::size_t)>& has_edge) {

    constexpr std::size_t kNil = std::numeric_limits<std::size_t>::max();
    constexpr std::size_t kInf = std::numeric_limits<std::size_t>::max();

    // pair_l[i] = j matched to left i, or kNil.
    // pair_r[j] = i matched to right j, or kNil.
    std::vector<std::size_t> pair_l(n_left, kNil);
    std::vector<std::size_t> pair_r(n_right, kNil);
    std::vector<std::size_t> dist(n_left + 1, kInf);  // dist[n_left] = NIL distance

    // Edge memoisation: pair (i, j) → bool. Saves redundant calls
    // during DFS path attempts.
    std::unordered_map<std::size_t, bool> edge_cache;
    auto edge = [&](std::size_t i, std::size_t j) -> bool {
        std::size_t key = i * n_right + j;
        auto it = edge_cache.find(key);
        if (it != edge_cache.end()) return it->second;
        bool v = has_edge(i, j);
        edge_cache.emplace(key, v);
        return v;
    };

    auto bfs = [&]() -> bool {
        std::queue<std::size_t> q;
        for (std::size_t i = 0; i < n_left; ++i) {
            if (pair_l[i] == kNil) {
                dist[i] = 0;
                q.push(i);
            } else {
                dist[i] = kInf;
            }
        }
        dist[n_left] = kInf;  // sentinel layer

        while (!q.empty()) {
            std::size_t u = q.front();
            q.pop();
            if (dist[u] >= dist[n_left]) continue;
            for (std::size_t v = 0; v < n_right; ++v) {
                if (!edge(u, v)) continue;
                std::size_t next = (pair_r[v] == kNil) ? n_left : pair_r[v];
                if (dist[next] == kInf) {
                    dist[next] = dist[u] + 1;
                    if (next != n_left) q.push(next);
                }
            }
        }
        return dist[n_left] != kInf;
    };

    std::function<bool(std::size_t)> dfs = [&](std::size_t u) -> bool {
        if (u == n_left) return true;
        for (std::size_t v = 0; v < n_right; ++v) {
            if (!edge(u, v)) continue;
            std::size_t next = (pair_r[v] == kNil) ? n_left : pair_r[v];
            if (dist[next] == dist[u] + 1) {
                if (dfs(next)) {
                    pair_r[v] = u;
                    pair_l[u] = v;
                    return true;
                }
            }
        }
        dist[u] = kInf;
        return false;
    };

    std::size_t matching = 0;
    while (bfs()) {
        for (std::size_t u = 0; u < n_left; ++u) {
            if (pair_l[u] == kNil && dfs(u)) {
                ++matching;
            }
        }
    }
    return matching;
}

}  // namespace cas::golden
