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

// McSplit+RL-down (Section 2.10): McSplit+RL's learned branching heuristic
// combined with McSplit-down's top-down search strategy (mcsplit-down.cpp).
// The two mechanisms are, per the survey, complementary and independent:
// McSplit-down changes the search *direction* (ask "does a mapping of size
// k exist?", decrementing k, instead of building up from empty and
// maximising); McSplit+RL changes *which vertices get branched on* within
// whatever search is run. This file is therefore base RL's branching/reward
// machinery (mcsplit-rl.cpp) unchanged, driven by McSplit-down's outer
// goal-size loop (mcsplit-down.cpp) instead of a single top-level solve
// call.
//
// One design decision without a stated answer in the survey or any public
// reference (neither variant's combination has a released implementation -
// see mcsplit-rl.cpp's and mcsplit-down.cpp's headers): whether the RL
// score tables (lgrade/rgrade) persist across goal-size attempts, or reset
// for each. This resets them every attempt: each goal value poses a
// differently-scoped question ("does a mapping of size `goal` exist?"),
// so a vertex that scored well while (unsuccessfully) chasing a larger
// goal is not assumed to be the right vertex to prioritise once the
// search has moved on to a smaller one, and keeping RL-down's branching
// heuristic and down's search strategy as cleanly separable components
// (matching this project's whole framework philosophy of isolating
// algorithms' distinguishing pieces) argues against carrying learned
// state across what are, by design, independent attempts - k-down's own
// static structural filters have no equivalent persistence question at
// all, having nothing that accumulates across k. g_matched/h_matched and
// the domains/left/right arrays also reset every attempt, same as
// McSplit-down, since those describe the current partial mapping rather
// than accumulated experience.

using gtype_rl_down = double;

// Forward declarations
static int calc_bound_rld(const std::vector<Bidomain>& domains);
static void remove_bidomain_rld(std::vector<Bidomain>& domains, int idx);
static int selectV_index_rld(const std::vector<int>& arr, const std::vector<gtype_rl_down>& grade, int start_idx, int len);
static int selectW_index_rld(const std::vector<int>& arr, const std::vector<gtype_rl_down>& grade, int start_idx, int len, const std::vector<int>& wselected);
static int select_bidomain_rld(const std::vector<Bidomain>& domains, const std::vector<int>& left, const std::vector<gtype_rl_down>& grade, int current_matching_size);
static int partition_rld(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
static void remove_vtx_from_array_rld(std::vector<int>& arr, int start_idx, int& len, int remove_idx);
static std::vector<Bidomain> build_domains_rld(const Graph& g, const Graph& h, std::vector<int>& left, std::vector<int>& right);

static int calc_bound_rld(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains) bound += std::min(bd.left_len, bd.right_len);
    return bound;
}

static void remove_bidomain_rld(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

static int selectV_index_rld(const std::vector<int>& arr, const std::vector<gtype_rl_down>& grade,
        int start_idx, int len) {
    int idx = -1;
    gtype_rl_down max_g = -1;
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

static int selectW_index_rld(const std::vector<int>& arr, const std::vector<gtype_rl_down>& grade,
        int start_idx, int len, const std::vector<int>& wselected) {
    int idx = -1;
    gtype_rl_down max_g = -1;
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

static int select_bidomain_rld(const std::vector<Bidomain>& domains,
        const std::vector<int>& left, const std::vector<gtype_rl_down>& grade,
        int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        int tie_breaker = left[bd.l + selectV_index_rld(left, grade, bd.l, bd.left_len)];
        if (len < min_size || (len == min_size && tie_breaker < min_tie_breaker)) {
            min_size = len;
            min_tie_breaker = tie_breaker;
            best = i;
        }
    }
    return best;
}

static int partition_rld(std::vector<int>& all_vv, int start, int len,
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

static void remove_vtx_from_array_rld(std::vector<int>& arr, int start_idx, int& len,
        int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

/**
 * Identical reward computation to mcsplit-rl.cpp's rewardfeed_rl (raw
 * bound-drop, accumulated into per-vertex lgrade/rgrade with the same 1e9
 * rescale). RL-down does not change the reward mechanism at all - only the
 * pruning condition and driving loop differ, exactly as McSplit-down only
 * changes those relative to plain McSplit.
 */
static std::vector<Bidomain> rewardfeed_rld(const std::vector<Bidomain>& d,
        std::vector<VtxPair>& current,
        std::vector<int>& g_matched, std::vector<int>& h_matched,
        std::vector<int>& left, std::vector<int>& right,
        std::vector<gtype_rl_down>& lgrade, std::vector<gtype_rl_down>& rgrade,
        const Graph& g, const Graph& h, int v, int w, bool multiway, Stats& stats) {
    current.push_back(VtxPair(v, w));
    g_matched[v] = 1;
    h_matched[w] = 1;

    std::vector<Bidomain> new_d;
    new_d.reserve(d.size() * 2);
    int temp = 0, total = 0;

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;

        int left_len = partition_rld(left, l, old_bd.left_len, g.adjmat[v]);
        int right_len = partition_rld(right, r, old_bd.right_len, h.adjmat[w]);
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

/**
 * Core BnB search: base RL's branching/reward logic, plus McSplit-down's
 * goal-reached cut (see mcsplit-down.cpp's solve_down for the identical
 * check and its justification).
 */
static void solve_rld(const Graph& g, const Graph& h,
        std::vector<gtype_rl_down>& lgrade, std::vector<gtype_rl_down>& rgrade,
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

    int bound = current.size() + calc_bound_rld(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    // McSplit-down's goal-reached cut (see mcsplit-down.cpp).
    if (incumbent.size() == goal) {
        stats.cut_branches++;
        return;
    }

    int bd_idx = select_bidomain_rld(domains, left, lgrade, current.size());
    if (bd_idx == -1) { return; }

    Bidomain& bd = domains[bd_idx];

    int tmp_idx = selectV_index_rld(left, lgrade, bd.l, bd.left_len);
    int v = left[bd.l + tmp_idx];
    remove_vtx_from_array_rld(left, bd.l, bd.left_len, tmp_idx);

    std::vector<int> wselected(h.n, 0);
    bd.right_len--;

    for (int i = 0; i <= bd.right_len; i++) {
        int w_idx = selectW_index_rld(right, rgrade, bd.r, bd.right_len + 1, wselected);
        if (w_idx == -1) { break; }
        int w = right[bd.r + w_idx];
        wselected[w] = 1;
        std::swap(right[bd.r + w_idx], right[bd.r + bd.right_len]);

        unsigned int cur_len = current.size();

        std::vector<Bidomain> new_domains = rewardfeed_rld(domains, current, g_matched, h_matched,
                left, right, lgrade, rgrade, g, h, v, w, multiway, stats);

        solve_rld(g, h, lgrade, rgrade, incumbent, current,
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
    if (bd.left_len == 0) { remove_bidomain_rld(domains, bd_idx); }

    solve_rld(g, h, lgrade, rgrade, incumbent, current,
            g_matched, h_matched, domains, left, right,
            goal, multiway, stats, start_time, abort_due_to_timeout);
}

/**
 * Fresh initial bidomains for one goal-size attempt (see mcsplit-down.cpp's
 * build_domains_down - identical, duplicated per this project's
 * one-file-per-algorithm convention).
 */
static std::vector<Bidomain> build_domains_rld(const Graph& g, const Graph& h,
        std::vector<int>& left, std::vector<int>& right) {
    left.clear();
    right.clear();
    std::vector<Bidomain> domains;

    std::set<unsigned int> left_labels, right_labels;
    for (unsigned int label : g.label) left_labels.insert(label);
    for (unsigned int label : h.label) right_labels.insert(label);
    std::set<unsigned int> labels;
    std::set_intersection(left_labels.begin(), left_labels.end(),
                          right_labels.begin(), right_labels.end(),
                          std::inserter(labels, labels.begin()));

    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();
        for (int i = 0; i < g.n; i++) if (g.label[i] == label) left.push_back(i);
        for (int i = 0; i < h.n; i++) if (h.label[i] == label) right.push_back(i);
        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        if (left_len && right_len) domains.push_back({start_l, start_r, left_len, right_len, false});
    }
    return domains;
}

/**
 * Entry point for McSplit+RL-down (Section 2.10). Sorting/orientation setup
 * is identical to mcs_rl_common in mcsplit-rl.cpp; the difference is the
 * top-down driving loop (see mcs_down in mcsplit-down.cpp), decrementing
 * the goal size from |V(g)| until a mapping of that size is found.
 * lgrade/rgrade persist across goal attempts - see this file's header for
 * why.
 *
 * @param g Graph G (must satisfy |V(g)| <= |V(h)|; enforced by main.cpp).
 */
std::vector<VtxPair> mcs_rl_down(const Graph& g, const Graph& h, bool multiway,
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

    std::vector<VtxPair> incumbent;
    auto start_time = std::chrono::steady_clock::now();

    for (int k = 0; k < g_sorted.n; k++) {
        if (abort_due_to_timeout) break;

        unsigned int goal = g_sorted.n - k;
        std::vector<int> left, right;
        std::vector<Bidomain> domains = build_domains_rld(g_sorted, h_sorted, left, right);
        std::vector<VtxPair> current;
        std::vector<int> g_matched(g_sorted.n, 0);
        std::vector<int> h_matched(h_sorted.n, 0);
        // Fresh RL score tables every goal attempt - each attempt is a
        // differently-scoped question ("does a mapping of size `goal`
        // exist?"), so scores learned while chasing one goal are not
        // assumed to transfer to the next (see this file's header).
        std::vector<gtype_rl_down> lgrade(g_sorted.n, 0);
        std::vector<gtype_rl_down> rgrade(h_sorted.n, 0);

        solve_rld(g_sorted, h_sorted, lgrade, rgrade, incumbent, current,
                g_matched, h_matched, domains, left, right,
                goal, multiway, stats, start_time, abort_due_to_timeout);

        if (incumbent.size() == goal) break;
    }

    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }
    return incumbent;
}
