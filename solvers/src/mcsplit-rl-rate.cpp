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

// McSplit+RL Rate (Section 2.9, Liu et al. AAAI 2020's ablation of their own
// McSplit+RL). Identical to plain McSplit+RL (mcsplit-rl.cpp) except for the
// reward itself: base RL accumulates the raw bound-drop R(v,w); RL Rate
// instead accumulates a rate-normalised reward
//     R_rate(v, w) = R(v, w) / overEstimate(D)
// where overEstimate(D) = sum_{bidomain in D} min(|G'|,|H'|) is McSplit's
// upper bound (Eq. 1) evaluated on the label classes D as they stood before
// this match - i.e. calc_bound_rlr(d), the same quantity solve_rlr's own
// pruning check already computes, just captured before the split rather
// than after. Both S_G(v)/lgrade and S_H(w)/rgrade are updated with this
// normalised value; vertex/bidomain selection, the split logic, and the
// 1e9 rescale-on-overflow mechanism are otherwise unchanged from base RL.
//
// There is no public reference implementation of this ablation (see
// mcsplit-rl.cpp's header); this follows Section 2.9's formula directly.
// overEstimate(D) is always >= 1 whenever there is a bidomain left to
// branch on (both sides non-empty), so no zero-division guard is needed on
// the path that reaches this function.

using gtype_rl_rate = double;

// Forward declarations
static int calc_bound_rlr(const std::vector<Bidomain>& domains);
static void remove_bidomain_rlr(std::vector<Bidomain>& domains, int idx);
static int selectV_index_rlr(const std::vector<int>& arr, const std::vector<gtype_rl_rate>& grade, int start_idx, int len);
static int selectW_index_rlr(const std::vector<int>& arr, const std::vector<gtype_rl_rate>& grade, int start_idx, int len, const std::vector<int>& wselected);
static int select_bidomain_rlr(const std::vector<Bidomain>& domains, const std::vector<int>& left, const std::vector<gtype_rl_rate>& grade, int current_matching_size);
static int partition_rlr(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
static void remove_vtx_from_array_rlr(std::vector<int>& arr, int start_idx, int& len, int remove_idx);

static int calc_bound_rlr(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains) bound += std::min(bd.left_len, bd.right_len);
    return bound;
}

static void remove_bidomain_rlr(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

static int selectV_index_rlr(const std::vector<int>& arr, const std::vector<gtype_rl_rate>& grade,
        int start_idx, int len) {
    int idx = -1;
    gtype_rl_rate max_g = -1;
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

static int selectW_index_rlr(const std::vector<int>& arr, const std::vector<gtype_rl_rate>& grade,
        int start_idx, int len, const std::vector<int>& wselected) {
    int idx = -1;
    gtype_rl_rate max_g = -1;
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

static int select_bidomain_rlr(const std::vector<Bidomain>& domains,
        const std::vector<int>& left, const std::vector<gtype_rl_rate>& grade,
        int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        int tie_breaker = left[bd.l + selectV_index_rlr(left, grade, bd.l, bd.left_len)];
        if (len < min_size || (len == min_size && tie_breaker < min_tie_breaker)) {
            min_size = len;
            min_tie_breaker = tie_breaker;
            best = i;
        }
    }
    return best;
}

static int partition_rlr(std::vector<int>& all_vv, int start, int len,
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

static void remove_vtx_from_array_rlr(std::vector<int>& arr, int start_idx, int& len,
        int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

/**
 * RL Rate's reward step. `total` is R(v,w) = overEstimate(D)-overEstimate(D'),
 * computed exactly as in base RL's rewardfeed_rl. The only change: before
 * accumulating into lgrade/rgrade, it is divided by overEstimate(D) - the
 * bound on the *pre-match* domain list `d`, i.e. calc_bound_rlr(d) - per
 * Section 2.9's R_rate(v,w) = R(v,w) / overEstimate(D).
 */
static std::vector<Bidomain> rewardfeed_rlr(const std::vector<Bidomain>& d,
        std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype_rl_rate>& lgrade, std::vector<gtype_rl_rate>& rgrade,
        const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats) {
    current.push_back(VtxPair(v, w));
    g_matched[v] = 1;
    h_matched[w] = 1;

    // overEstimate(D): bound on the domain list as it stood before this
    // match, captured before filtering/splitting mutates left/right.
    int over_estimate_d = calc_bound_rlr(d);

    std::vector<Bidomain> new_d;
    new_d.reserve(d.size() * 2);
    int temp = 0, total = 0;

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;

        int left_len = partition_rlr(left, l, old_bd.left_len, g.adjmat[v]);
        int right_len = partition_rlr(right, r, old_bd.right_len, h.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;

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

    if (total > 0) {
        // over_estimate_d is >= 1 here: d had at least the bidomain we just
        // branched from, with both sides non-empty (solve_rlr only reaches
        // this function when select_bidomain_rlr found something to branch
        // on), so min(left_len,right_len) for that bidomain alone is >= 1.
        gtype_rl_rate rate_reward = static_cast<gtype_rl_rate>(total) / over_estimate_d;

        lgrade[v] += rate_reward;
        if (lgrade[v] > 1e9) {
            for (int i = 0; i < g.n; i++) { lgrade[i] /= 1e9; }
        }
        rgrade[w] += rate_reward;
        if (rgrade[w] > 1e9) {
            for (int i = 0; i < h.n; i++) { rgrade[i] /= 1e9; }
        }
    }
    return new_d;
}

static void solve_rlr(const Graph& g, const Graph& h,
        std::vector<gtype_rl_rate>& lgrade, std::vector<gtype_rl_rate>& rgrade,
        std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right,
        unsigned int goal, bool multiway, Stats& stats,
        std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout) {

    if (abort_due_to_timeout) { return; }

    stats.nodes++;
    stats.bidomain_counts.push_back(domains.size());

    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best = stats.nodes;
        stats.time_to_best = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    int bound = current.size() + calc_bound_rlr(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    int bd_idx = select_bidomain_rlr(domains, left, lgrade, current.size());
    if (bd_idx == -1) { return; }

    Bidomain& bd = domains[bd_idx];

    int tmp_idx = selectV_index_rlr(left, lgrade, bd.l, bd.left_len);
    int v = left[bd.l + tmp_idx];
    remove_vtx_from_array_rlr(left, bd.l, bd.left_len, tmp_idx);

    std::vector<int> wselected(h.n, 0);
    bd.right_len--;

    for (int i = 0; i <= bd.right_len; i++) {
        int w_idx = selectW_index_rlr(right, rgrade, bd.r, bd.right_len + 1, wselected);
        if (w_idx == -1) { break; }
        int w = right[bd.r + w_idx];
        wselected[w] = 1;
        std::swap(right[bd.r + w_idx], right[bd.r + bd.right_len]);

        unsigned int cur_len = current.size();

        std::vector<Bidomain> new_domains = rewardfeed_rlr(domains, current, g_matched, h_matched,
                left, right, lgrade, rgrade, g, h, v, w, multiway, stats);

        solve_rlr(g, h, lgrade, rgrade, incumbent, current,
                g_matched, h_matched, new_domains, left, right,
                goal, multiway, stats, start_time, abort_due_to_timeout);

        while (current.size() > cur_len) {
            VtxPair pr = current.back();
            current.pop_back();
            g_matched[pr.v] = 0;
            h_matched[pr.w] = 0;
        }
    }

    bd.right_len++;
    if (bd.left_len == 0) { remove_bidomain_rlr(domains, bd_idx); }

    solve_rlr(g, h, lgrade, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right,
            goal, multiway, stats, start_time, abort_due_to_timeout);
}

std::vector<VtxPair> mcs_rl_rate(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout) {

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
        for (int i = 0; i < g_sorted.n; i++) { if (g_sorted.label[i] == label) left.push_back(i); }
        for (int i = 0; i < h_sorted.n; i++) { if (h_sorted.label[i] == label) right.push_back(i); }
        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        if (left_len && right_len) domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    std::vector<gtype_rl_rate> lgrade(g_sorted.n, 0);
    std::vector<gtype_rl_rate> rgrade(h_sorted.n, 0);
    std::vector<int> g_matched(g_sorted.n, 0);
    std::vector<int> h_matched(h_sorted.n, 0);

    std::vector<VtxPair> incumbent, current;
    auto start_time = std::chrono::steady_clock::now();
    solve_rlr(g_sorted, h_sorted, lgrade, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right, 1,
            multiway, stats, start_time, abort_due_to_timeout);

    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }
    return incumbent;
}
