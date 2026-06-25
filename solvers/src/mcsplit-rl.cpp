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

using gtype_rl = double;

// Forward declarations
static int calc_bound_rl(const std::vector<Bidomain>& domains);
static void remove_bidomain_rl(std::vector<Bidomain>& domains, int idx);
static int selectV_index_rl(const std::vector<int>& arr, const std::vector<gtype_rl>& grade, int start_idx, int len);
static int selectW_index_rl(const std::vector<int>& arr, const std::vector<gtype_rl>& grade, int start_idx, int len, const std::vector<int>& wselected);
static int select_bidomain_rl(const std::vector<Bidomain>& domains, const std::vector<int>& left, const std::vector<gtype_rl>& grade, int current_matching_size);
static int partition_rl(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
static int remove_matched_vertex_rl(std::vector<int>& arr, int start, int len, const std::vector<int>& matched);
static void remove_vtx_from_array_rl(std::vector<int>& arr, int start_idx, int& len, int remove_idx);
static std::vector<Bidomain> rewardfeed_rl(const std::vector<Bidomain>& d, std::vector<VtxPair>& current, std::vector<int>& g_matched, std::vector<int>& h_matched, std::vector<int>& left, std::vector<int>& right, std::vector<gtype_rl>& lgrade, std::vector<gtype_rl>& rgrade, const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats);
static std::vector<Bidomain> rewardfeed_rllum(const std::vector<Bidomain>& d, std::vector<VtxPair>& current, std::vector<int>& g_matched, std::vector<int>& h_matched, std::vector<int>& left, std::vector<int>& right, std::vector<gtype_rl>& lgrade, std::vector<gtype_rl>& rgrade, const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats);
static void solve_rl(const Graph& g, const Graph& h, std::vector<gtype_rl>& lgrade, std::vector<gtype_rl>& rgrade, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<int>& g_matched, std::vector<int>& h_matched, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, bool multiway, bool use_lum, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_rl(const Graph& g, const Graph& h, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_rllum(const Graph& g, const Graph& h, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout);

static int calc_bound_rl(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

static void remove_bidomain_rl(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

// Selects vertex with highest RL score, ties broken on smallest index.
static int selectV_index_rl(const std::vector<int>& arr, const std::vector<gtype_rl>& grade,
        int start_idx, int len) {
    int idx = -1;
    gtype_rl max_g = -1;
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

// Selects w with highest RL score, skipping already-tried vertices.
static int selectW_index_rl(const std::vector<int>& arr, const std::vector<gtype_rl>& grade,
        int start_idx, int len, const std::vector<int>& wselected) {
    int idx = -1;
    gtype_rl max_g = -1;
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
static int select_bidomain_rl(const std::vector<Bidomain>& domains,
        const std::vector<int>& left, const std::vector<gtype_rl>& grade,
        int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        int tie_breaker = left[bd.l + selectV_index_rl(left, grade, bd.l, bd.left_len)];
        if (len < min_size || (len == min_size && tie_breaker < min_tie_breaker)) {
            min_size = len;
            min_tie_breaker = tie_breaker;
            best = i;
        }
    }
    return best;
}

static int partition_rl(std::vector<int>& all_vv, int start, int len,
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
// Used by LUM variant to clean up non-adjacent bidomains after leaf matching.
static int remove_matched_vertex_rl(std::vector<int>& arr, int start, int len,
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

static void remove_vtx_from_array_rl(std::vector<int>& arr, int start_idx, int& len,
        int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

// RL reward function without LUM.
// Matches (v, w), computes the RL reward (bound reduction caused by the split),
// and updates lgrade[v] and rgrade[w]. Both scores use the same threshold (1e9)
// — this is the key difference from McSplit+LL which uses asymmetric thresholds
// (1e5 for lgrade, 1e9 for rgrade).
static std::vector<Bidomain> rewardfeed_rl(const std::vector<Bidomain>& d,
        std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype_rl>& lgrade, std::vector<gtype_rl>& rgrade,
        const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats) {
    current.push_back(VtxPair(v, w));
    g_matched[v] = 1;
    h_matched[w] = 1;

    std::vector<Bidomain> new_d;
    new_d.reserve(d.size());
    int temp = 0, total = 0;

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;

        int left_len = partition_rl(left, l, old_bd.left_len, g.adjmat[v]);
        int right_len = partition_rl(right, r, old_bd.right_len, h.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;

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

    // Update RL scores — both use 1e9 threshold (symmetric, unlike LL)
    if (total > 0) {
        stats.conflicts++;
        lgrade[v] += total;
        if (lgrade[v] > 1e9) {
            for (int i = 0; i < g.n; i++) { lgrade[i] /= 1e9; }
        }
        rgrade[w] += total;
        if (rgrade[w] > 1e9) {
            for (int i = 0; i < h.n; i++) { rgrade[i] /= 1e9; }
        }
    }
    return new_d;
}

// RL reward function with LUM.
// Extends rewardfeed_rl with leaf union matching: after matching (v, w),
// immediately match all compatible leaf pairs. Leaf vertices have degree 1
// so they always end up in the same bidomain — matching them in bulk avoids
// branching on trivial choices.
// Score update is identical to rewardfeed_rl (symmetric 1e9 threshold).
static std::vector<Bidomain> rewardfeed_rllum(const std::vector<Bidomain>& d,
        std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype_rl>& lgrade, std::vector<gtype_rl>& rgrade,
        const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats) {
    current.push_back(VtxPair(v, w));
    g_matched[v] = 1;
    h_matched[w] = 1;

    // LUM: merge-scan leaf lists to find matching label groups,
    // greedily match unmatched leaf pairs.
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
            unmatched_left_len = remove_matched_vertex_rl(left, l, old_bd.left_len, g_matched);
            unmatched_right_len = remove_matched_vertex_rl(right, r, old_bd.right_len, h_matched);
        } else {
            unmatched_left_len = old_bd.left_len;
            unmatched_right_len = old_bd.right_len;
        }

        int left_len = partition_rl(left, l, unmatched_left_len, g.adjmat[v]);
        int right_len = partition_rl(right, r, unmatched_right_len, h.adjmat[w]);
        int left_len_noedge = unmatched_left_len - left_len;
        int right_len_noedge = unmatched_right_len - right_len;

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

    // Score update identical to RL — symmetric 1e9 threshold on both sides
    if (total > 0) {
        stats.conflicts++;
        lgrade[v] += total;
        if (lgrade[v] > 1e9) {
            for (int i = 0; i < g.n; i++) { lgrade[i] /= 1e9; }
        }
        rgrade[w] += total;
        if (rgrade[w] > 1e9) {
            for (int i = 0; i < h.n; i++) { rgrade[i] /= 1e9; }
        }
    }
    return new_d;
}

// Core BnB search shared by RL and RL+LUM.
// use_lum flag selects which reward function is called.
// Structure is identical to McSplit+DAL with M fixed to 1 (RL-only),
// but without the policy alternation machinery (no M, num, Maxnum).
static void solve_rl(const Graph& g, const Graph& h,
        std::vector<gtype_rl>& lgrade, std::vector<gtype_rl>& rgrade,
        std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right,
        unsigned int goal, bool multiway, bool use_lum, Stats& stats,
        std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout) {

    if (abort_due_to_timeout) { return; }

    stats.nodes++;

    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best = stats.nodes;
        stats.time_to_best = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    int bound = current.size() + calc_bound_rl(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    // Select bidomain and v using lgrade scores
    int bd_idx = select_bidomain_rl(domains, left, lgrade, current.size());
    if (bd_idx == -1) { return; }

    Bidomain& bd = domains[bd_idx];

    int tmp_idx = selectV_index_rl(left, lgrade, bd.l, bd.left_len);
    int v = left[bd.l + tmp_idx];
    remove_vtx_from_array_rl(left, bd.l, bd.left_len, tmp_idx);

    std::vector<int> wselected(h.n, 0);
    bd.right_len--;

    // Try matching v with each w, selected by rgrade score
    for (int i = 0; i <= bd.right_len; i++) {
        int w_idx = selectW_index_rl(right, rgrade, bd.r, bd.right_len + 1, wselected);
        if (w_idx == -1) { break; }
        int w = right[bd.r + w_idx];
        wselected[w] = 1;
        std::swap(right[bd.r + w_idx], right[bd.r + bd.right_len]);

        // Record size before rewardfeed (LUM may add extra pairs)
        unsigned int cur_len = current.size();

        std::vector<Bidomain> new_domains;
        if (use_lum) {
            new_domains = rewardfeed_rllum(domains, current, g_matched, h_matched,
                    left, right, lgrade, rgrade, g, h, v, w, multiway, stats);
        } else {
            new_domains = rewardfeed_rl(domains, current, g_matched, h_matched,
                    left, right, lgrade, rgrade, g, h, v, w, multiway, stats);
        }

        solve_rl(g, h, lgrade, rgrade, incumbent, current,
                g_matched, h_matched, new_domains, left, right,
                goal, multiway, use_lum, stats, start_time, abort_due_to_timeout);

        // Undo all matches added by rewardfeed (v,w and any LUM leaves)
        while (current.size() > cur_len) {
            VtxPair pr = current.back();
            current.pop_back();
            g_matched[pr.v] = 0;
            h_matched[pr.w] = 0;
        }
    }

    bd.right_len++;
    if (bd.left_len == 0) { remove_bidomain_rl(domains, bd_idx); }

    solve_rl(g, h, lgrade, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right,
            goal, multiway, use_lum, stats, start_time, abort_due_to_timeout);
}

// Shared setup for both RL and RL+LUM.
static std::vector<VtxPair> mcs_rl_common(const Graph& g, const Graph& h,
        bool multiway, bool use_lum, Stats& stats,
        std::atomic<bool>& abort_due_to_timeout) {

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

    if (use_lum) {
        pack_leaves(g_sorted);
        pack_leaves(h_sorted);
    }

    // Per-label bidomains — same as McSplit, DAL, LL
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

    stats.root_upper_bound = calc_bound_rl(domains);

    std::vector<gtype_rl> lgrade(g_sorted.n, 0);
    std::vector<gtype_rl> rgrade(h_sorted.n, 0);
    std::vector<int> g_matched(g_sorted.n, 0);
    std::vector<int> h_matched(h_sorted.n, 0);

    std::vector<VtxPair> incumbent, current;
    auto start_time = std::chrono::steady_clock::now();
    solve_rl(g_sorted, h_sorted, lgrade, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right, 1,
            multiway, use_lum, stats, start_time, abort_due_to_timeout);

    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}

std::vector<VtxPair> mcs_rl(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout) {
    return mcs_rl_common(g, h, multiway, false, stats, abort_due_to_timeout);
}

std::vector<VtxPair> mcs_rllum(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout) {
    return mcs_rl_common(g, h, multiway, true, stats, abort_due_to_timeout);
}