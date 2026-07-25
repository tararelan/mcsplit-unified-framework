// #define _GNU_SOURCE

#include "graph.h"
#include "common.h"
#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <limits.h>

// Forward declarations
static int calc_bound(const std::vector<Bidomain>& domains);
int select_bidomain(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size);
int find_min_value(const std::vector<int>& arr, int start_idx, int len);
int partition(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
void remove_vtx_from_left_domain(std::vector<int>& left, Bidomain& bd, int v);
static void remove_bidomain(std::vector<Bidomain>& domains, int idx);
std::vector<Bidomain> filter_domains(const std::vector<Bidomain>& d, std::vector<int>& left, std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool multiway);
void solve(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, bool multiway, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs(const Graph& g, const Graph& h, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout);

// McSplit's upper bound.
// For each bidomain, at most min(|G-side|, |H-side|) further pairs can be matched
// since each vertex can only be used once. Summing over all bidomains gives the
// maximum number of additional pairs that could possibly be added to the current
// partial solution.
static int calc_bound(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains)
        bound += std::min(bd.left_len, bd.right_len);
    return bound;
}

// Finds the smallest vertex index in left[start_idx .. start_idx+len-1].
// Used as a tiebreaker in bidomain selection - when two bidomains have the
// same max(left_len, right_len), we prefer the one whose left set contains
// the smallest-indexed vertex. This makes the search deterministic and
// consistent with the reference implementation.
int find_min_value(const std::vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i = 0; i < len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

// Selects which bidomain to branch on next.
// McSplit's heuristic: pick the bidomain with the smallest max(left_len, right_len).
// This is a fail-first strategy - branch on the most constrained choice first,
// which tends to find contradictions earlier and prune more of the search tree.
// Ties broken by smallest vertex index in the left set.
int select_bidomain(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

// Partitions the slice all_vv[start .. start+len-1] in-place so that vertices
// adjacent to the vertex represented by adjrow come first, followed by
// non-adjacent vertices. Returns the count of adjacent vertices (i.e. the
// length of the adjacent prefix).
// This is the core operation of filter_domains - it splits one side of a
// bidomain into adjacent and non-adjacent parts in O(len) time without
// allocating any new memory.
int partition(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow) {
    int i = 0;
    for (int j = 0; j < len; j++) {
        if (adjrow[all_vv[start + j]]) {
            std::swap(all_vv[start + i], all_vv[start + j]);
            i++;
        }
    }
    return i;
}

// Removes vertex v from the left side of its bidomain by swapping it with
// the last element and decrementing left_len. This is O(n) in the bidomain
// size but avoids shifting elements. v must be in the bidomain.
void remove_vtx_from_left_domain(std::vector<int>& left, Bidomain& bd, int v) {
    int i = 0;
    while (left[bd.l + i] != v) i++;
    std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
    bd.left_len--;
}

// Removes a bidomain from the list by replacing it with the last element
// and popping the back. O(1) - order within the list does not matter.
static void remove_bidomain(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

// Computes the new set of bidomains after committing to the match (v, w).
// For each existing bidomain (L, R), splits it into:
//   - vertices in L adjacent to v, paired with vertices in R adjacent to w
//   - vertices in L not adjacent to v, paired with vertices in R not adjacent to w
// Only non-empty splits are kept. This enforces the induced subgraph constraint:
// any future match must preserve adjacency/non-adjacency relative to (v,w).
//
// In multiway mode (labelled or directed graphs), the adjacent part is further
// split by edge label value - vertices with different edge labels to v (or w)
// go into separate bidomains, since they can only be matched with vertices
// carrying the same label.
std::vector<Bidomain> filter_domains(const std::vector<Bidomain>& d, std::vector<int>& left,
        std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool multiway) {
    std::vector<Bidomain> new_d;
    new_d.reserve(d.size() * 2);
    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        // Split left side: adjacent to v goes to front, non-adjacent to back
        int left_len = partition(left, l, old_bd.left_len, g.adjmat[v]);
        // Split right side: adjacent to w goes to front, non-adjacent to back
        int right_len = partition(right, r, old_bd.right_len, h.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        // Non-adjacent part: only keep if both sides are non-empty
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
            // Labelled/directed: further split adjacent part by edge label value.
            // Sort both sides by their edge label to v/w respectively, then
            // merge-scan to pair up vertices with matching labels.
            auto& adjrow_v = g.adjmat[v];
            auto& adjrow_w = h.adjmat[w];
            std::sort(left.begin() + l, left.begin() + l + left_len,
                    [&](int a, int b) { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(right.begin() + r, right.begin() + r + right_len,
                    [&](int a, int b) { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            int li = l, ri = r;
            while (li < l_top && ri < r_top) {
                unsigned int left_label = adjrow_v[left[li]];
                unsigned int right_label = adjrow_w[right[ri]];
                if (left_label < right_label) { li++; }
                else if (left_label > right_label) { ri++; }
                else {
                    // Both sides have the same edge label - create one bidomain
                    // for this label value, spanning all consecutive vertices with it
                    int lmin = li, rmin = ri;
                    do { li++; } while (li < l_top && adjrow_v[left[li]] == left_label);
                    do { ri++; } while (ri < r_top && adjrow_w[right[ri]] == left_label);
                    new_d.push_back({lmin, rmin, li - lmin, ri - rmin, true});
                }
            }
        } else if (left_len && right_len) {
            // Unlabelled/undirected: adjacent part stays as one bidomain
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }
    return new_d;
}

// Core branch-and-bound recursive search.
// At each node:
//   1. Check timeout
//   2. Update incumbent if current solution is larger
//   3. Compute upper bound - prune if it cannot beat incumbent
//   4. Select the most constrained bidomain (smallest max(|L|,|R|))
//   5. Select the smallest-indexed vertex v from its left side
//   6. Try matching v with each w in the right side in ascending order:
//      - compute new domains after committing to (v,w)
//      - recurse
//      - backtrack
//   7. Also try not matching v at all (the "exclude v" branch)
void solve(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent,
        std::vector<VtxPair>& current, std::vector<Bidomain>& domains,
        std::vector<int>& left, std::vector<int>& right, unsigned int goal,
        bool multiway, Stats& stats,
        std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout) {

    if (abort_due_to_timeout) return;

    stats.nodes++;

    // Update incumbent if current partial solution is the best seen so far
    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best = stats.nodes;
        stats.time_to_best = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    // Compute upper bound: current size + max possible additional matches
    // Prune if bound cannot beat incumbent or reach goal
    int bound = current.size() + calc_bound(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    // Select bidomain to branch on: smallest max(left_len, right_len)
    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1) return; // No bidomains left - search complete at this node

    Bidomain& bd = domains[bd_idx];

    // Always branch on the smallest-indexed vertex in the left set
    // This makes the search deterministic and matches the reference implementation
    int v = find_min_value(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain(left, domains[bd_idx], v);

    int w = -1;
    bd.right_len--;

    // Branch 1: try matching v with each w in the right set, in ascending order
    for (int i = 0; i <= bd.right_len; i++) {
        // Find next smallest w greater than the previous w tried
        int idx = -1;
        int smallest = INT_MAX;
        for (int j = 0; j <= bd.right_len; j++) {
            if (right[bd.r + j] > w && right[bd.r + j] < smallest) {
                smallest = right[bd.r + j];
                idx = j;
            }
        }
        if (idx == -1) break;
        w = right[bd.r + idx];
        // Move w to end of right slice so filter_domains doesn't see it
        std::swap(right[bd.r + idx], right[bd.r + bd.right_len]);

        // Compute new bidomains after committing to match (v, w)
        auto new_domains = filter_domains(domains, left, right, g, h, v, w, multiway);
        current.push_back(VtxPair(v, w));
        solve(g, h, incumbent, current, new_domains, left, right, goal,
                multiway, stats, start_time, abort_due_to_timeout);
        current.pop_back(); // Backtrack
    }

    bd.right_len++;
    // If v was the last vertex in its bidomain's left side, remove the bidomain
    if (bd.left_len == 0)
        remove_bidomain(domains, bd_idx);

    // Branch 2: try not matching v at all - recurse without v in any solution
    solve(g, h, incumbent, current, domains, left, right, goal,
            multiway, stats, start_time, abort_due_to_timeout);
}

// Entry point for McSplit search.
// Sorts vertices by degree, builds initial bidomains, runs solve(), and
// converts the solution back to original vertex indices.
std::vector<VtxPair> mcs(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout) {

    // Calculate degree of each vertex for sorting.
    // Only counts outgoing edges (adjmat[v][w] & 1) - consistent with
    // the reference implementation's degree-based vertex ordering.
    auto calc_degrees = [](const Graph& g) {
        std::vector<int> degree(g.n, 0);
        for (int v = 0; v < g.n; v++)
            for (int w = 0; w < g.n; w++)
                if (g.adjmat[v][w] & 1) degree[v]++;
        return degree;
    };

    std::vector<int> g_deg = calc_degrees(g);
    std::vector<int> h_deg = calc_degrees(h);

    // Determine graph density to decide sort direction.
    // Dense graphs: sort ascending (low-degree vertices first).
    // Sparse graphs: sort descending (high-degree vertices first).
    // The sort is applied to G using H's density and vice versa - this
    // matches the reference implementation's sorting logic.
    int g_edges = 0, h_edges = 0;
    for (int d : g_deg) g_edges += d;
    for (int d : h_deg) h_edges += d;
    bool h_dense = h_edges > h.n * (h.n - 1);
    bool g_dense = g_edges > g.n * (g.n - 1);

    std::vector<int> vv0(g.n), vv1(h.n);
    std::iota(vv0.begin(), vv0.end(), 0);
    std::iota(vv1.begin(), vv1.end(), 0);
    std::stable_sort(vv0.begin(), vv0.end(), [&](int a, int b) {
        return h_dense ? g_deg[a] < g_deg[b] : g_deg[a] > g_deg[b];
    });
    std::stable_sort(vv1.begin(), vv1.end(), [&](int a, int b) {
        return g_dense ? h_deg[a] < h_deg[b] : h_deg[a] > h_deg[b];
    });

    // Create sorted copies of both graphs - search runs on these,
    // solution is converted back to original indices at the end
    Graph g_sorted = induced_subgraph(const_cast<Graph&>(g), vv0);
    Graph h_sorted = induced_subgraph(const_cast<Graph&>(h), vv1);

    // Build initial bidomains - one per label value shared by both graphs.
    // Each bidomain pairs all G-vertices with a given label against all
    // H-vertices with the same label. At the start, no vertices are matched
    // so label compatibility is the only constraint.
    std::vector<int> left, right;
    std::vector<Bidomain> domains;

    std::set<unsigned int> left_labels, right_labels;
    for (unsigned int label : g_sorted.label) left_labels.insert(label);
    for (unsigned int label : h_sorted.label) right_labels.insert(label);
    std::set<unsigned int> labels; // Labels appearing in both graphs
    std::set_intersection(left_labels.begin(), left_labels.end(),
                          right_labels.begin(), right_labels.end(),
                          std::inserter(labels, labels.begin()));

    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();
        for (int i = 0; i < g_sorted.n; i++)
            if (g_sorted.label[i] == label)
                left.push_back(i);
        for (int i = 0; i < h_sorted.n; i++)
            if (h_sorted.label[i] == label)
                right.push_back(i);
        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        if (left_len && right_len)
            domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    std::vector<VtxPair> incumbent, current;
    auto start_time = std::chrono::steady_clock::now();
    solve(g_sorted, h_sorted, incumbent, current, domains, left, right, 1,
            multiway, stats, start_time, abort_due_to_timeout);

    // Convert solution indices from sorted graph back to original graph indices
    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}