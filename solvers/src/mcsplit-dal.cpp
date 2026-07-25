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
#include <iterator>

using gtype = double;
const int short_memory_threshold = 1e5;       // decay threshold for short-term RL scores
const long long int long_memory_threshold = 1e9; // decay threshold for long-term DAL scores

// Forward declarations
static int calc_bound(const std::vector<Bidomain>& domains);
int select_bidomain_dal(const std::vector<Bidomain>& domains, const std::vector<int>& left, const std::vector<gtype>& grade, int current_matching_size);
int selectV_index(const std::vector<int>& arr, const std::vector<gtype>& grade, int start_idx, int len);
int selectW_index(const std::vector<int>& arr, const std::vector<gtype>& grade, int start_idx, int len, const std::vector<int>& wselected);
int partition_dal(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
int remove_matched_vertex(std::vector<int>& arr, int start, int len, const std::vector<int>& matched);
void remove_vtx_from_array(std::vector<int>& arr, int start_idx, int& len, int remove_idx);
static void remove_bidomain(std::vector<Bidomain>& domains, int idx);
std::vector<Bidomain> rewardfeed_RL(const std::vector<Bidomain>& d, std::vector<VtxPair>& current, std::vector<int>& left, std::vector<int>& right, std::vector<gtype>& lgrade, std::vector<gtype>& rgrade, std::vector<int>& g_matched, std::vector<int>& h_matched, const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats);
std::vector<Bidomain> rewardfeed_DAL(const std::vector<Bidomain>& d, std::vector<VtxPair>& current, std::vector<int>& g_matched, std::vector<int>& h_matched, std::vector<int>& left, std::vector<int>& right, std::vector<gtype>& V, std::vector<gtype>& Qv, const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats);
void solve_dal(const Graph& g, const Graph& h, std::vector<gtype>& V, std::vector<gtype>& lgrade, std::vector<gtype>& rgrade, std::vector<std::vector<gtype>>& Q, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<int>& g_matched, std::vector<int>& h_matched, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, bool multiway, int& M, int& num, int Maxnum, int policy_mode, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_dal(const Graph& g, const Graph& h, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout, int policy_mode = 0);  // 0=hybrid, 1=RL_only, 2=DAL_only

// McSplit's upper bound: sum of min(left_len, right_len) over all bidomains.
// From each bidomain at most min(|L|, |R|) further pairs can be matched.
static int calc_bound(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

// Partitions arr[start..start+len-1] so adjacent vertices come first.
// Returns count of adjacent vertices.
// Renamed to partition_dal to avoid conflict with std::partition from <algorithm>.
int partition_dal(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow) {
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
// Returns the new length after removal.
// Used by LUM to clean up non-adjacent bidomains after leaf matching.
int remove_matched_vertex(std::vector<int>& arr, int start, int len, const std::vector<int>& matched) {
    int p = 0;
    for (int i = 0; i < len; i++) {
        if (!matched[arr[start + i]]) {
            std::swap(arr[start + i], arr[start + p]);
            p++;
        }
    }
    return p;
}

// Removes a vertex from an array by swapping with last element and decrementing len.
// O(1) removal - used to remove v from the left side of its bidomain before branching.
void remove_vtx_from_array(std::vector<int>& arr, int start_idx, int& len, int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

// Removes a bidomain by replacing it with the last element and popping.
// O(1) - order within the list does not matter.
static void remove_bidomain(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

// Selects the vertex with the highest accumulated score in grade[],
// breaking ties on smallest vertex index.
// Used for both v selection (G-side) and bidomain tiebreaking.
int selectV_index(const std::vector<int>& arr, const std::vector<gtype>& grade, int start_idx, int len) {
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

// Selects the w vertex (H-side) with the highest score, skipping already-selected vertices.
// wselected tracks which w values have already been tried in the current branching loop.
int selectW_index(const std::vector<int>& arr, const std::vector<gtype>& grade, int start_idx, int len, const std::vector<int>& wselected) {
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

// Selects bidomain with smallest max(left_len, right_len).
// Ties broken by the highest-scored vertex in the left set (not smallest index as in McSplit).
// This is the key difference from McSplit's select_bidomain - DAL/RL use scores for tiebreaking.
int select_bidomain_dal(const std::vector<Bidomain>& domains, const std::vector<int>& left, const std::vector<gtype>& grade, int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        int tie_breaker = left[bd.l + selectV_index(left, grade, bd.l, bd.left_len)];
        if (len < min_size || (len == min_size && tie_breaker < min_tie_breaker)) {
            min_size = len;
            min_tie_breaker = tie_breaker;
            best = i;
        }
    }
    return best;
}

// RL reward function: splits bidomains after matching (v, w) and updates lgrade/rgrade.
// The RL reward is the reduction in the upper bound caused by the match:
//   reward = sum over bidomains of (old_min - new_adj_min - new_nadj_min)
// If reward > 0, add it to lgrade[v] and rgrade[w].
// Decay all scores by half if any score exceeds short_memory_threshold,
// preventing historical scores from dominating future decisions.
std::vector<Bidomain> rewardfeed_RL(const std::vector<Bidomain>& d, std::vector<VtxPair>& current,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype>& lgrade, std::vector<gtype>& rgrade,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats) {
    std::vector<Bidomain> new_d;
    new_d.reserve(d.size() * 2);
    current.push_back(VtxPair(v, w));
    g_matched[v] = 1;
    h_matched[w] = 1;

    int temp = 0, total = 0;
    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        int left_len = partition_dal(left, l, old_bd.left_len, g.adjmat[v]);
        int right_len = partition_dal(right, r, old_bd.right_len, h.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;

        // RL reward contribution from this bidomain:
        // how much the bound dropped due to the split
        temp = std::min(old_bd.left_len, old_bd.right_len)
             - std::min(left_len, right_len)
             - std::min(left_len_noedge, right_len_noedge);
        total += temp;

        if (left_len_noedge && right_len_noedge) {
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        }
        if (multiway && left_len && right_len) {
            // Labelled/directed: further split adjacent part by edge label
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
        lgrade[v] += total;
        if (lgrade[v] > short_memory_threshold) {
            for (int i = 0; i < g.n; i++) { lgrade[i] /= 2; }
        }
        rgrade[w] += total;
        if (rgrade[w] > short_memory_threshold) {
            for (int i = 0; i < h.n; i++) { rgrade[i] /= 2; }
        }
    }
    return new_d;
}

// DAL reward function: splits bidomains after matching (v, w) and updates V/Q scores.
// DAL reward = bound reduction + number of new domains after the split.
// The domain count term rewards matches that fragment the problem into more independent
// subproblems, addressing the limitation of RL which only rewards bound reduction.
//
// Also performs LUM (Leaf Vertex Union Match): when (v, w) is matched, immediately
// match all compatible leaf pairs in bulk without branching. A leaf vertex has exactly
// one neighbour. Leaf pairs always end up in the same bidomain regardless of other
// decisions, so matching them all at once collapses many branching decisions into one.
std::vector<Bidomain> rewardfeed_DAL(const std::vector<Bidomain>& d, std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype>& V, std::vector<gtype>& Qv,
        const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats) {
    current.push_back(VtxPair(v, w));
    g_matched[v] = 1;
    h_matched[w] = 1;

    // LUM: match leaf vertices of v and w in bulk.
    // g.leaves[v] and h.leaves[w] are precomputed sorted lists of leaf neighbours
    // grouped by (edge_label, vertex_label). We merge-scan the two lists and
    // greedily match unmatched leaf pairs with matching labels.
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
    new_d.reserve(d.size() * 2);
    int temp = 0, total = 0;

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        int unmatched_left_len, unmatched_right_len;

        // If LUM matched any leaves, remove them from non-adjacent bidomains
        // since they are no longer available for future matches
        if (leaves_match_size > 0 && !old_bd.is_adjacent) {
            unmatched_left_len = remove_matched_vertex(left, l, old_bd.left_len, g_matched);
            unmatched_right_len = remove_matched_vertex(right, r, old_bd.right_len, h_matched);
        } else {
            unmatched_left_len = old_bd.left_len;
            unmatched_right_len = old_bd.right_len;
        }

        int left_len = partition_dal(left, l, unmatched_left_len, g.adjmat[v]);
        int right_len = partition_dal(right, r, unmatched_right_len, h.adjmat[w]);
        int left_len_noedge = unmatched_left_len - left_len;
        int right_len_noedge = unmatched_right_len - right_len;

        // RL component of DAL reward: bound reduction from this bidomain
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

    // DAL reward = bound reduction + domain count (number of new bidomains).
    // The domain count rewards fragmentation of the residual problem.
    int domgrade = new_d.size();
    total = total + domgrade;
    if (total > 0) {
        V[v] += total;
        Qv[w] += total;
        // Decay V scores (short-term) if threshold exceeded
        if (V[v] > short_memory_threshold) {
            for (int i = 0; i < g.n; i++) { V[i] /= 2; }
        }
        // Decay Q scores (long-term, per-pair) if threshold exceeded
        if (Qv[w] > long_memory_threshold) {
            for (int i = 0; i < h.n; i++) { Qv[i] /= 2; }
        }
    }
    return new_d;
}

// Core BnB recursive search for McSplit+DAL.
// Alternates between two branching policies on a fixed schedule:
//   M=1: RL policy - select v and w by lgrade/rgrade (bound-reduction scores)
//   M=2: DAL policy - select v and w by V/Q[v] (bound-reduction + domain-count scores)
// Policy switches every Maxnum = 2*min(|G|,|H|) branching decisions.
// This hybrid approach avoids the Matthew effect where a single policy
// concentrates on a small subset of high-score vertices.
void solve_dal(const Graph& g, const Graph& h,
        std::vector<gtype>& V, std::vector<gtype>& lgrade, std::vector<gtype>& rgrade,
        std::vector<std::vector<gtype>>& Q,
        std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right,
        unsigned int goal, bool multiway, int& M, int& num, int Maxnum,
        int policy_mode, Stats& stats,
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
		num = 0;
    }

    // Prune if upper bound cannot beat incumbent
    int bound = current.size() + calc_bound(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    // Increment branching counter and switch policy if threshold reached
    num++;
    if (num > Maxnum) {
		if (policy_mode == 0) { M = M % 2 + 1; }
		else { M = policy_mode; }
		num = 0;
	}

    // Select bidomain using current policy's scores
    int bd_idx = -1;
    if (M == 1) { bd_idx = select_bidomain_dal(domains, left, lgrade, current.size()); }
    else        { bd_idx = select_bidomain_dal(domains, left, V, current.size()); }
    if (bd_idx == -1) { return; }

    Bidomain& bd = domains[bd_idx];

    // Select v using current policy's scores
    int tmp_idx = -1;
    if (M == 1) { tmp_idx = selectV_index(left, lgrade, bd.l, bd.left_len); }
    else        { tmp_idx = selectV_index(left, V, bd.l, bd.left_len); }

    int v = left[bd.l + tmp_idx];
    remove_vtx_from_array(left, bd.l, bd.left_len, tmp_idx);

    std::vector<int> wselected(h.n, 0);
    bd.right_len--;

    // Try matching v with each w in the right side of the bidomain
    for (int i = 0; i <= bd.right_len; i++) {
        // Update policy counter for each w tried (not just each v)
        if (i != 0) {
            num++;
            if (num > Maxnum) {
				if (policy_mode == 0) { M = M % 2 + 1; }
				else { M = policy_mode; }
				num = 0;
			}
        }

        // Select w using current policy's scores, skipping already-tried vertices
        int w_idx = -1;
        if (M == 1) { w_idx = selectW_index(right, rgrade, bd.r, bd.right_len + 1, wselected); }
        else        { w_idx = selectW_index(right, Q[v], bd.r, bd.right_len + 1, wselected); }
        if (w_idx == -1) { break; }

        int w = right[bd.r + w_idx];
        wselected[w] = 1;
        std::swap(right[bd.r + w_idx], right[bd.r + bd.right_len]);

        // Record current solution length before reward functions add LUM matches
        unsigned int cur_len = current.size();
        std::vector<Bidomain> new_domains;

        // Apply reward function for current policy
        // rewardfeed_RL: adds (v,w) to current, splits domains, updates RL scores
        // rewardfeed_DAL: adds (v,w) + leaf matches to current, splits domains, updates DAL scores
        if (M == 1) {
            new_domains = rewardfeed_RL(domains, current, left, right, lgrade, rgrade, g_matched, h_matched, g, h, v, w, multiway, stats);
        } else {
            new_domains = rewardfeed_DAL(domains, current, g_matched, h_matched, left, right, V, Q[v], g, h, v, w, multiway, stats);
        }

        solve_dal(g, h, V, lgrade, rgrade, Q, incumbent, current,
			g_matched, h_matched, new_domains, left, right, goal,
			multiway, M, num, Maxnum, policy_mode, stats, start_time, abort_due_to_timeout);

        // Undo all matches added by the reward function (v,w and any LUM leaves)
        while (current.size() > cur_len) {
            VtxPair pr = current.back();
            current.pop_back();
            g_matched[pr.v] = 0;
            h_matched[pr.w] = 0;
        }
    }

    bd.right_len++;
    if (bd.left_len == 0) { remove_bidomain(domains, bd_idx); }

    // Branch where v is not matched at all
    solve_dal(g, h, V, lgrade, rgrade, Q, incumbent, current,
        g_matched, h_matched, domains, left, right,
        goal, multiway, M, num, Maxnum, policy_mode, stats, start_time, abort_due_to_timeout);
}

// Entry point for McSplit+DAL search.
// Sorts vertices by degree, builds initial bidomains, initialises score vectors,
// packs leaf lists for LUM, runs solve_dal(), and converts solution back to
// original vertex indices.
std::vector<VtxPair> mcs_dal(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout,
        int policy_mode) {

    // Calculate degrees for vertex sorting
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

    // Determine graph density to decide sort direction
    int g_edges = 0, h_edges = 0;
    for (int d : g_deg) { g_edges += d; }
    for (int d : h_deg) { h_edges += d; }
    bool h_dense = h_edges > h.n * (h.n - 1);
    bool g_dense = g_edges > g.n * (g.n - 1);

    // Sort vertices: ascending if sparse, descending if dense
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

    // Precompute leaf vertex lists for LUM optimisation
    pack_leaves(g_sorted);
    pack_leaves(h_sorted);

    // Build initial label class domains
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

    // Initialise score vectors
    // lgrade/rgrade: per-vertex RL scores (short-term memory)
    // V: per-vertex DAL scores (short-term memory)
    // Q: per-pair DAL scores (long-term memory)
    // Maxnum: how many branching decisions before policy switches
    int Maxnum = 2 * std::min(g.n, h.n);
    int M = 1, num = 0;
    std::vector<gtype> lgrade(g_sorted.n, 0);
    std::vector<gtype> V(g_sorted.n, 0);
    std::vector<gtype> rgrade(h_sorted.n, 0);
    std::vector<std::vector<gtype>> Q(g_sorted.n, std::vector<gtype>(h_sorted.n, 0));
    std::vector<int> g_matched(g_sorted.n, 0);
    std::vector<int> h_matched(h_sorted.n, 0);

    std::vector<VtxPair> incumbent, current;
    auto start_time = std::chrono::steady_clock::now();
    solve_dal(g, h, V, lgrade, rgrade, Q, incumbent, current,
        g_matched, h_matched, domains, left, right,
        1, multiway, M, num, Maxnum, policy_mode, stats, start_time, abort_due_to_timeout);

    // Convert solution indices back to original (unsorted) vertex indices
    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}