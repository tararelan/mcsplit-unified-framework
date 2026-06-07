// #define _GNU_SOURCE

#include "graph.h"
#include "common.h"
#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <chrono>
#include <atomic>
#include <limits.h>

using gtype = double;
const int short_memory_threshold_ll = 1e5;
const long long int long_memory_threshold_ll = 1e9;

// Forward declarations
static int calc_bound_ll(const std::vector<Bidomain>& domains);
static void remove_bidomain_ll(std::vector<Bidomain>& domains, int idx);
int selectV_index_ll(const std::vector<int>& arr, const std::vector<gtype>& grade, int start_idx, int len);
int selectW_index_ll(const std::vector<int>& arr, const std::vector<gtype>& grade, int start_idx, int len, const std::vector<int>& wselected);
int select_bidomain_ll(const std::vector<Bidomain>& domains, const std::vector<int>& left, const std::vector<gtype>& grade, int current_matching_size);
int partition_ll(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
int remove_matched_vertex_ll(std::vector<int>& arr, int start, int len, const std::vector<int>& matched);
void remove_vtx_from_array_ll(std::vector<int>& arr, int start_idx, int& len, int remove_idx);
std::vector<Bidomain> rewardfeed_ll(const std::vector<Bidomain>& d, std::vector<VtxPair>& current, std::vector<int>& g_matched, std::vector<int>& h_matched, std::vector<int>& left, std::vector<int>& right, std::vector<gtype>& lgrade, std::vector<gtype>& rgrade, const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats);
void solve_ll(const Graph& g, const Graph& h, std::vector<gtype>& lgrade, std::vector<gtype>& rgrade, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<int>& g_matched, std::vector<int>& h_matched, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, bool multiway, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_ll(const Graph& g, const Graph& h, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout);

static int calc_bound_ll(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

static void remove_bidomain_ll(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

// Selects vertex with highest score, ties broken on smallest index.
int selectV_index_ll(const std::vector<int>& arr, const std::vector<gtype>& grade,
        int start_idx, int len) {
    int idx = -1;
    gtype max_g = -1;
    int best_vtx = INT_MAX;
    for (int i = 0; i < len; i++) {
        int vtx = arr[start_idx + i];
        if (grade[vtx] > max_g || (grade[vtx] == max_g && vtx < best_vtx)) {
            idx = i;
            best_vtx = vtx;
            max_g = grade[vtx];
        }
    }
    return idx;
}

// Selects w with highest score, skipping already-tried vertices.
int selectW_index_ll(const std::vector<int>& arr, const std::vector<gtype>& grade,
        int start_idx, int len, const std::vector<int>& wselected) {
    int idx = -1;
    gtype max_g = -1;
    int best_vtx = INT_MAX;
    for (int i = 0; i < len; i++) {
        int vtx = arr[start_idx + i];
        if (wselected[vtx] == 0) {
            if (grade[vtx] > max_g || (grade[vtx] == max_g && vtx < best_vtx)) {
                idx = i;
                best_vtx = vtx;
                max_g = grade[vtx];
            }
        }
    }
    return idx;
}

// Selects bidomain with smallest max(left_len, right_len),
// ties broken by highest-scored vertex in left set.
int select_bidomain_ll(const std::vector<Bidomain>& domains, const std::vector<int>& left,
        const std::vector<gtype>& grade, int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        int tie_breaker = left[bd.l + selectV_index_ll(left, grade, bd.l, bd.left_len)];
        if (len < min_size || (len == min_size && tie_breaker < min_tie_breaker)) {
            min_size = len;
            min_tie_breaker = tie_breaker;
            best = i;
        }
    }
    return best;
}

// Partitions arr[start..start+len-1] so adjacent vertices come first.
int partition_ll(std::vector<int>& all_vv, int start, int len,
        const std::vector<unsigned int>& adjrow) {
    int i = 0;
    for (int j = 0; j < len; j++) {
        if (adjrow[all_vv[start + j]]) {
            std::swap(all_vv[start + i], all_vv[start + j]);
            i++;
        }
    }
    return i;
}

// Removes already-matched vertices from a domain side in-place.
// Used by LUM to clean up non-adjacent bidomains after leaf matching.
int remove_matched_vertex_ll(std::vector<int>& arr, int start, int len,
        const std::vector<int>& matched) {
    int p = 0;
    for (int i = 0; i < len; i++) {
        if (!matched[arr[start + i]]) {
            std::swap(arr[start + i], arr[start + p]);
            p++;
        }
    }
    return p;
}

// Removes a vertex from an array by swapping with last element. O(1).
void remove_vtx_from_array_ll(std::vector<int>& arr, int start_idx, int& len, int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

// Combined LUM + RL reward function.
// After matching (v, w):
//   1. Immediately match all compatible leaf pairs (LUM) — leaf vertices have
//      degree 1 so they always end up in the same bidomain; matching them in
//      bulk avoids branching on trivial choices.
//   2. Compute the RL reward = reduction in the upper bound caused by the split.
//   3. Update lgrade[v] and rgrade[w] by the reward. Decay all scores by half
//      if either score exceeds short_memory_threshold.
// This is the key difference from DAL: LL uses only the RL reward (no V/Q
// scores, no domain count term, no policy alternation).
std::vector<Bidomain> rewardfeed_ll(const std::vector<Bidomain>& d,
        std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype>& lgrade, std::vector<gtype>& rgrade,
        const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats) {
    current.push_back(VtxPair(v, w));
    g_matched[v] = 1;
    h_matched[w] = 1;

    // LUM: match leaf vertices of v and w in bulk.
    // g.leaves[v] is a sorted list of (edge_label, vertex_label) → [leaf vertices].
    // Merge-scan to find matching label groups and greedily match unmatched pairs.
    int leaves_match_size = 0;
    for (unsigned int i = 0, j = 0; i < g.leaves[v].size() && j < h.leaves[w].size(); ) {
        if (g.leaves[v][i].first < h.leaves[w][j].first) { i++; }
        else if (g.leaves[v][i].first > h.leaves[w][j].first) { j++; }
        else {
            const std::vector<int>& leaf_g = g.leaves[v][i].second;
            const std::vector<int>& leaf_h = h.leaves[w][j].second;
            for (unsigned int p = 0, q = 0; p < leaf_g.size() && q < leaf_h.size(); ) {
                if (g_matched[leaf_g[p]]) { p++; }
                else if (h_matched[leaf_h[q]]) { q++; }
                else {
                    int v_leaf = leaf_g[p], w_leaf = leaf_h[q];
                    p++; q++;
                    current.push_back(VtxPair(v_leaf, w_leaf));
                    g_matched[v_leaf] = 1;
                    h_matched[w_leaf] = 1;
                    leaves_match_size++;
                }
            }
            i++; j++;
        }
    }

    std::vector<Bidomain> new_d;
    new_d.reserve(d.size());
    int temp = 0, total = 0;
    int unmatched_left_len, unmatched_right_len;

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;

        // If LUM matched leaves, remove them from non-adjacent bidomains
        if (leaves_match_size > 0 && !old_bd.is_adjacent) {
            unmatched_left_len = remove_matched_vertex_ll(left, l, old_bd.left_len, g_matched);
            unmatched_right_len = remove_matched_vertex_ll(right, r, old_bd.right_len, h_matched);
        } else {
            unmatched_left_len = old_bd.left_len;
            unmatched_right_len = old_bd.right_len;
        }

        int left_len = partition_ll(left, l, unmatched_left_len, g.adjmat[v]);
        int right_len = partition_ll(right, r, unmatched_right_len, h.adjmat[w]);
        int left_len_noedge = unmatched_left_len - left_len;
        int right_len_noedge = unmatched_right_len - right_len;

        // RL reward: how much the bound dropped due to this split
        temp = std::min(old_bd.left_len, old_bd.right_len)
             - std::min(left_len, right_len)
             - std::min(left_len_noedge, right_len_noedge);
        total += temp;

        if (left_len_noedge && right_len_noedge) {
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        }
        if (multiway && left_len && right_len) {
            auto& adjrow_v = g.adjmat[v];
            auto& adjrow_w = h.adjmat[w];
            std::sort(left.begin() + l, left.begin() + l + left_len,
                    [&](int a, int b) { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(right.begin() + r, right.begin() + r + right_len,
                    [&](int a, int b) { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len, r_top = r + right_len;
            int li = l, ri = r;
            while (li < l_top && ri < r_top) {
                unsigned int ll = adjrow_v[left[li]];
                unsigned int rl = adjrow_w[right[ri]];
                if (ll < rl) { li++; }
                else if (ll > rl) { ri++; }
                else {
                    int lmin = li, rmin = ri;
                    do { li++; } while (li < l_top && adjrow_v[left[li]] == ll);
                    do { ri++; } while (ri < r_top && adjrow_w[right[ri]] == ll);
                    new_d.push_back({lmin, rmin, li - lmin, ri - rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }

    // Update RL scores and decay if threshold exceeded
    if (total > 0) {
        stats.conflicts++;
        lgrade[v] += total;
        rgrade[w] += total;
        if (lgrade[v] > short_memory_threshold_ll) {
            for (int i = 0; i < g.n; i++) { lgrade[i] /= 2; }
        }
        if (rgrade[w] > long_memory_threshold_ll) {
            for (int i = 0; i < h.n; i++) { rgrade[i] /= 2; }
        }
    }
    return new_d;
}

// Core BnB recursive search for McSplit+LL.
// Identical structure to McSplit+DAL with M=1 (RL-only mode), but without
// policy alternation. Uses lgrade/rgrade for both bidomain selection and
// v/w selection throughout the entire search.
void solve_ll(const Graph& g, const Graph& h,
        std::vector<gtype>& lgrade, std::vector<gtype>& rgrade,
        std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right,
        unsigned int goal, bool multiway, Stats& stats,
        std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout) {

    if (abort_due_to_timeout) { return; }

    stats.nodes++;

    // Update incumbent if current solution is the best seen so far
    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best = stats.nodes;
        stats.time_to_best = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    // Prune if upper bound cannot beat incumbent
    int bound = current.size() + calc_bound_ll(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    // Select bidomain using lgrade scores
    int bd_idx = select_bidomain_ll(domains, left, lgrade, current.size());
    if (bd_idx == -1) { return; }

    Bidomain& bd = domains[bd_idx];

    // Select v using lgrade scores
    int tmp_idx = selectV_index_ll(left, lgrade, bd.l, bd.left_len);
    int v = left[bd.l + tmp_idx];
    remove_vtx_from_array_ll(left, bd.l, bd.left_len, tmp_idx);

    std::vector<int> wselected(h.n, 0);
    bd.right_len--;

    // Try matching v with each w, selected by rgrade score
    for (int i = 0; i <= bd.right_len; i++) {
        int w_idx = selectW_index_ll(right, rgrade, bd.r, bd.right_len + 1, wselected);
        if (w_idx == -1) { break; }
        int w = right[bd.r + w_idx];
        wselected[w] = 1;
        std::swap(right[bd.r + w_idx], right[bd.r + bd.right_len]);

        // Record current size before rewardfeed adds LUM matches
        unsigned int cur_len = current.size();

        auto new_domains = rewardfeed_ll(domains, current, g_matched, h_matched,
                left, right, lgrade, rgrade, g, h, v, w, multiway, stats);

        solve_ll(g, h, lgrade, rgrade, incumbent, current,
                g_matched, h_matched, new_domains, left, right,
                goal, multiway, stats, start_time, abort_due_to_timeout);

        // Undo all matches added by rewardfeed (v,w and any LUM leaves)
        while (current.size() > cur_len) {
            VtxPair pr = current.back();
            current.pop_back();
            g_matched[pr.v] = 0;
            h_matched[pr.w] = 0;
        }
    }

    bd.right_len++;
    if (bd.left_len == 0) { remove_bidomain_ll(domains, bd_idx); }

    // Branch where v is not matched
    solve_ll(g, h, lgrade, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right,
            goal, multiway, stats, start_time, abort_due_to_timeout);
}

// Entry point for McSplit+LL search.
// Sorts vertices by degree, packs leaf lists for LUM, builds per-label
// initial bidomains, initialises RL score vectors, runs solve_ll().
std::vector<VtxPair> mcs_ll(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout) {

    auto calc_degrees = [](const Graph& g) {
        std::vector<int> degree(g.n, 0);
        for (int v = 0; v < g.n; v++) {
            for (int w = 0; w < g.n; w++) {
                if (g.adjmat[v][w] & 1) { degree[v]++; }
            }
        }
        return degree;
    };

    std::vector<int> g_deg = calc_degrees(g);
    std::vector<int> h_deg = calc_degrees(h);

    int g_edges = 0, h_edges = 0;
    for (int d : g_deg) { g_edges += d; }
    for (int d : h_deg) { h_edges += d; }
    bool h_dense = h_edges > h.n * (h.n - 1) / 2;
    bool g_dense = g_edges > g.n * (g.n - 1) / 2;

    // Sort vertices: ascending degree if sparse, descending if dense.
    // Note: LL reference uses > n*(n-1) threshold (same as RRSplit/SymSplit)
    // not the McSplit /2 formulation.
    std::vector<int> vv0(g.n), vv1(h.n);
    std::iota(vv0.begin(), vv0.end(), 0);
    std::iota(vv1.begin(), vv1.end(), 0);
    std::stable_sort(vv0.begin(), vv0.end(), [&](int a, int b) {
        return h_dense ? g_deg[a] < g_deg[b] : g_deg[a] > g_deg[b];
    });
    std::stable_sort(vv1.begin(), vv1.end(), [&](int a, int b) {
        return g_dense ? h_deg[a] < h_deg[b] : h_deg[a] > h_deg[b];
    });

    Graph g_sorted = induced_subgraph(const_cast<Graph&>(g), vv0);
    Graph h_sorted = induced_subgraph(const_cast<Graph&>(h), vv1);

    // Precompute leaf lists for LUM
    pack_leaves(g_sorted);
    pack_leaves(h_sorted);

    // Build per-label bidomains — LL uses label splitting unlike RRSplit/SymSplit
    std::vector<int> left, right;
    std::vector<Bidomain> domains;

    std::set<unsigned int> left_labels, right_labels;
    for (unsigned int label : g_sorted.label) { left_labels.insert(label); }
    for (unsigned int label : h_sorted.label) { right_labels.insert(label); }
    std::set<unsigned int> labels;
    std::set_intersection(left_labels.begin(), left_labels.end(),
                          right_labels.begin(), right_labels.end(),
                          std::inserter(labels, labels.begin()));

    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();
        for (int i = 0; i < g_sorted.n; i++) {
            if (g_sorted.label[i] == label) { left.push_back(i); }
        }
        for (int i = 0; i < h_sorted.n; i++) {
            if (h_sorted.label[i] == label) { right.push_back(i); }
        }
        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        if (left_len && right_len) {
            domains.push_back({start_l, start_r, left_len, right_len, false});
        }
    }

    stats.root_upper_bound = calc_bound_ll(domains);

    // Initialise RL score vectors — all zero at start of search
    std::vector<gtype> lgrade(g_sorted.n, 0);
    std::vector<gtype> rgrade(h_sorted.n, 0);
    std::vector<int> g_matched(g_sorted.n, 0);
    std::vector<int> h_matched(h_sorted.n, 0);

    std::vector<VtxPair> incumbent, current;
    auto start_time = std::chrono::steady_clock::now();
    solve_ll(g_sorted, h_sorted, lgrade, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right, 1,
            multiway, stats, start_time, abort_due_to_timeout);

    // Convert solution indices back to original vertex indices
    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}