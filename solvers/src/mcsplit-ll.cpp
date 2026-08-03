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

// Forward declarations for McSplit+LL (LSM scoring + LUM leaf matching,
// both independently toggleable via use_lsm/use_lum for the ablation
// study - see Section 4.3).
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

// NOTE: Qv and rgrade are alternatives, not both used simultaneously -
// Qv (per-pair Q[v] row) is read when use_lsm=true, rgrade (per-vertex)
// is read when use_lsm=false. Both are passed so the caller doesn't need
// two separate rewardfeed variants for the ablation's LSM-on/off modes.
std::vector<Bidomain> rewardfeed_ll(const std::vector<Bidomain>& d,
        std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype>& lgrade,
        std::vector<gtype>& Qv,        // per-pair row Q[v] - used when use_lsm=true
        std::vector<gtype>& rgrade,    // per-vertex rgrade - used when use_lsm=false
        const Graph& g, const Graph& h,
        int v, int w, bool multiway,
        bool use_lum, bool use_lsm,
        Stats& stats);

void solve_ll(const Graph& g, const Graph& h,
        std::vector<gtype>& lgrade,
        std::vector<std::vector<gtype>>& Q,
        std::vector<gtype>& rgrade,
        std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right,
        unsigned int goal, bool multiway,
        bool use_lum, bool use_lsm,
        Stats& stats,
        std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout);

std::vector<VtxPair> mcs_ll(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout,
        bool use_lum = true, bool use_lsm = true);


// McSplit's upper bound: sum of min(left_len, right_len) over all bidomains
static int calc_bound_ll(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains)
        bound += std::min(bd.left_len, bd.right_len);
    return bound;
}

// Removes a bidomain by replacing it with the last element and popping;
// order within the list does not matter
static void remove_bidomain_ll(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

// Selects the vertex with the highest accumulated score in grade[],
// breaking ties on smallest vertex index. max_g starts at -1 so the
// first candidate is always accepted even if all scores are still 0.
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

// Same as selectV_index_ll but skips vertices already tried in the
// current branching loop (wselected). Returns -1 if all candidates
// in range have already been tried - callers must check for this.
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

// Selects bidomain with smallest max(left_len, right_len). Ties broken
// by comparing each domain's highest-scored left vertex (via
// selectV_index_ll), then taking the domain whose best vertex has the
// smallest index.
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

// Boolean pre-partition (adjacent vs non-adjacent); direction handled
// correctly by the multiway split downstream, same as elsewhere
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
// Returns the new length after removal. Used by LUM to clean up
// non-adjacent bidomains after leaf matching.
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

// O(1) removal by swap-with-last-and-shrink; side-agnostic (works for
// either the G-side or H-side array)
void remove_vtx_from_array_ll(std::vector<int>& arr, int start_idx, int& len, int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

// use_lum: if true, bulk-match leaf pairs after matching (v,w)
// use_lsm: if true, use per-pair Q[v][w] for w-scoring; if false, use per-vertex rgrade[w]
//
// NOTE: correctly marks each LUM-matched leaf pair as matched immediately
// (g_matched[v_leaf]=1, h_matched[w_leaf]=1 below), same as
// rewardfeed_RL/rewardfeed_DAL. Does not have the reference DAL binary's
// missing-update bug documented in Section 4.1.4.
std::vector<Bidomain> rewardfeed_ll(
        const std::vector<Bidomain>& d,
        std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype>& lgrade,
        std::vector<gtype>& Qv,
        std::vector<gtype>& rgrade,
        const Graph& g, const Graph& h,
        int v, int w, bool multiway,
        bool use_lum, bool use_lsm,
        Stats& stats) {

    current.push_back(VtxPair(v, w));
    g_matched[v] = 1;
    h_matched[w] = 1;

    // LUM: bulk-match compatible leaf pairs
    int leaves_match_size = 0;
    if (use_lum) {
        for (unsigned int i = 0, j = 0;
             i < g.leaves[v].size() && j < h.leaves[w].size(); ) {
            if (g.leaves[v][i].first < h.leaves[w][j].first) { i++; }
            else if (g.leaves[v][i].first > h.leaves[w][j].first) { j++; }
            else {
                const std::vector<int>& leaf_g = g.leaves[v][i].second;
                const std::vector<int>& leaf_h = h.leaves[w][j].second;
                for (unsigned int p = 0, q = 0;
                     p < leaf_g.size() && q < leaf_h.size(); ) {
                    if (g_matched[leaf_g[p]]) { p++; }
                    else if (h_matched[leaf_h[q]]) { q++; }
                    else {
                        int v_leaf = leaf_g[p], w_leaf = leaf_h[q];
                        p++; q++;
                        current.push_back(VtxPair(v_leaf, w_leaf));
                        g_matched[v_leaf] = 1; // matched immediately - correct, unlike reference
                        h_matched[w_leaf] = 1;
                        leaves_match_size++;
                    }
                }
                i++; j++;
            }
        }
    }

    // Domain split
    std::vector<Bidomain> new_d;
    new_d.reserve(d.size() * 2);
    int temp = 0, total = 0;
    int unmatched_left_len, unmatched_right_len;

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;

        if (leaves_match_size > 0 && !old_bd.is_adjacent) {
            unmatched_left_len  = remove_matched_vertex_ll(left,  l, old_bd.left_len,  g_matched);
            unmatched_right_len = remove_matched_vertex_ll(right, r, old_bd.right_len, h_matched);
        } else {
            unmatched_left_len  = old_bd.left_len;
            unmatched_right_len = old_bd.right_len;
        }

        int left_len  = partition_ll(left,  l, unmatched_left_len,  g.adjmat[v]);
        int right_len = partition_ll(right, r, unmatched_right_len, h.adjmat[w]);
        int left_len_noedge  = unmatched_left_len  - left_len;
        int right_len_noedge = unmatched_right_len - right_len;

        temp = std::min(old_bd.left_len, old_bd.right_len)
             - std::min(left_len, right_len)
             - std::min(left_len_noedge, right_len_noedge);
        total += temp;

        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l + left_len, r + right_len,
                             left_len_noedge, right_len_noedge, old_bd.is_adjacent});

        if (multiway && left_len && right_len) {
            // direction-aware label split (see rewardfeed_RL in mcsplit-dal.cpp)
            auto& adjrow_v = g.adjmat[v];
            auto& adjrow_w = h.adjmat[w];
            std::sort(left.begin()  + l, left.begin()  + l + left_len,
                    [&](int a, int b) { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(right.begin() + r, right.begin() + r + right_len,
                    [&](int a, int b) { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len, r_top = r + right_len;
            int li = l, ri = r;
            while (li < l_top && ri < r_top) {
                unsigned int ll = adjrow_v[left[li]];
                unsigned int rl = adjrow_w[right[ri]];
                if      (ll < rl) { li++; }
                else if (ll > rl) { ri++; }
                else {
                    int lmin = li, rmin = ri;
                    do { li++; } while (li < l_top && adjrow_v[left[li]]  == ll);
                    do { ri++; } while (ri < r_top && adjrow_w[right[ri]] == ll);
                    new_d.push_back({lmin, rmin, li - lmin, ri - rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }

    // LSM score update
    // use_lsm=true:  update lgrade[v] (S0, short-term) and Qv[w] (St, long-term per-pair)
    // use_lsm=false: update lgrade[v] and rgrade[w] (both per-vertex, McSplit+RL behaviour)
    if (total > 0) {
        lgrade[v] += total;
        if (use_lsm) {
            Qv[w] += total;
            if (lgrade[v] > short_memory_threshold_ll)
                for (int i = 0; i < (int)lgrade.size(); i++) lgrade[i] /= 2;
            if (Qv[w] > long_memory_threshold_ll)
                for (int i = 0; i < (int)Qv.size(); i++) Qv[i] /= 2;
        } else {
            rgrade[w] += total;
            if (lgrade[v] > short_memory_threshold_ll)
                for (int i = 0; i < (int)lgrade.size(); i++) lgrade[i] /= 2;
            if (rgrade[w] > short_memory_threshold_ll)
                for (int i = 0; i < (int)rgrade.size(); i++) rgrade[i] /= 2;
        }
    }
    return new_d;
}

// BnB search

void solve_ll(const Graph& g, const Graph& h,
        std::vector<gtype>& lgrade,
        std::vector<std::vector<gtype>>& Q,
        std::vector<gtype>& rgrade,
        std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right,
        unsigned int goal, bool multiway,
        bool use_lum, bool use_lsm,
        Stats& stats,
        std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout) {

    if (abort_due_to_timeout) return;
    stats.nodes++;

    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best  = stats.nodes;
        stats.time_to_best   = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    int bound = current.size() + calc_bound_ll(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    // v selection always uses lgrade (S0 / per-vertex short-term score)
    int bd_idx = select_bidomain_ll(domains, left, lgrade, current.size());
    if (bd_idx == -1) return;

    Bidomain& bd = domains[bd_idx];
    int tmp_idx = selectV_index_ll(left, lgrade, bd.l, bd.left_len);
    int v = left[bd.l + tmp_idx];
    remove_vtx_from_array_ll(left, bd.l, bd.left_len, tmp_idx);

    std::vector<int> wselected(h.n, 0);
    bd.right_len--;

    // visits each right-side candidate once (same guaranteed-non-negative
    // idx pattern as solve_dal/solve_dsb, but here uses an explicit
    // wselected mask rather than the increasing-order trick)
    for (int i = 0; i <= bd.right_len; i++) {
        // w selection: per-pair Q[v] if use_lsm, per-vertex rgrade otherwise
        const std::vector<gtype>& w_scores = use_lsm ? Q[v] : rgrade;
        int w_idx = selectW_index_ll(right, w_scores, bd.r, bd.right_len + 1, wselected);
        if (w_idx == -1) break;
        int w = right[bd.r + w_idx];
        wselected[w] = 1;
        std::swap(right[bd.r + w_idx], right[bd.r + bd.right_len]);

        unsigned int cur_len = current.size();

        auto new_domains = rewardfeed_ll(domains, current, g_matched, h_matched,
                left, right, lgrade, Q[v], rgrade, g, h, v, w,
                multiway, use_lum, use_lsm, stats);

        solve_ll(g, h, lgrade, Q, rgrade, incumbent, current,
                g_matched, h_matched, new_domains, left, right,
                goal, multiway, use_lum, use_lsm,
                stats, start_time, abort_due_to_timeout);

        while (current.size() > cur_len) {
            VtxPair pr = current.back();
            current.pop_back();
            g_matched[pr.v] = 0;
            h_matched[pr.w] = 0;
        }
    }

    bd.right_len++;
    if (bd.left_len == 0) remove_bidomain_ll(domains, bd_idx);

    // branch where v is not matched at all
    solve_ll(g, h, lgrade, Q, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right,
            goal, multiway, use_lum, use_lsm,
            stats, start_time, abort_due_to_timeout);
}

// Entry point
// use_lum=true,  use_lsm=true  → full McSplit+LL
// use_lum=false, use_lsm=true  → LSM only (no leaf matching)
// use_lum=true,  use_lsm=false → LUM only (per-vertex w-scoring, McSplit+RL style)
// use_lum=false, use_lsm=false → plain McSplit+RL (sanity check - should match rl results)

std::vector<VtxPair> mcs_ll(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout,
        bool use_lum, bool use_lsm) {

    // total degree (out + in), matching the reference's calculate_degrees
    auto calc_degrees = [](const Graph& g) {
		std::vector<int> degree(g.n, 0);
		for (int v = 0; v < g.n; v++)
			for (int w = 0; w < g.n; w++) {
				unsigned int mask = 0xFFFFu;
				if (g.adjmat[v][w] & mask)  degree[v]++;
				if (g.adjmat[v][w] & ~mask) degree[v]++;
			}
		return degree;
	};

    std::vector<int> g_deg = calc_degrees(g);
    std::vector<int> h_deg = calc_degrees(h);

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

    Graph g_sorted = induced_subgraph(const_cast<Graph&>(g), vv0);
    Graph h_sorted = induced_subgraph(const_cast<Graph&>(h), vv1);

    // leaf lists only needed when LUM is enabled - skipped otherwise to
    // avoid unnecessary O(n^2) work in the LSM-only ablation variant
    if (use_lum) {
        pack_leaves(g_sorted);
        pack_leaves(h_sorted);
    }

    std::vector<int> left, right;
    std::vector<Bidomain> domains;

    std::set<unsigned int> left_labels, right_labels;
    for (unsigned int label : g_sorted.label) left_labels.insert(label);
    for (unsigned int label : h_sorted.label) right_labels.insert(label);
    std::set<unsigned int> labels;
    std::set_intersection(left_labels.begin(), left_labels.end(),
                          right_labels.begin(), right_labels.end(),
                          std::inserter(labels, labels.begin()));

    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();
        for (int i = 0; i < g_sorted.n; i++)
            if (g_sorted.label[i] == label) left.push_back(i);
        for (int i = 0; i < h_sorted.n; i++)
            if (h_sorted.label[i] == label) right.push_back(i);
        int left_len  = left.size()  - start_l;
        int right_len = right.size() - start_r;
        if (left_len && right_len)
            domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    // S0: short-term per-vertex score for v (and w when use_lsm=false)
    std::vector<gtype> lgrade(g_sorted.n, 0);
    // St: long-term per-pair score table for w (used when use_lsm=true)
    std::vector<std::vector<gtype>> Q(g_sorted.n, std::vector<gtype>(h_sorted.n, 0));
    // per-vertex w score (used when use_lsm=false, i.e. McSplit+RL style)
    std::vector<gtype> rgrade(h_sorted.n, 0);

    std::vector<int> g_matched(g_sorted.n, 0);
    std::vector<int> h_matched(h_sorted.n, 0);
    std::vector<VtxPair> incumbent, current;

    auto start_time = std::chrono::steady_clock::now();
    // correctly uses g_sorted/h_sorted, matching domains/left/right
    solve_ll(g_sorted, h_sorted, lgrade, Q, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right, 1,
            multiway, use_lum, use_lsm,
            stats, start_time, abort_due_to_timeout);

    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }
    return incumbent;
}