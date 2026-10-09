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

// McSplit-down ("McSplit↓" / "big-first" search, Section 2.6): rather than
// building a mapping from empty and maximising, this variant asks "does a
// mapping of size k exist?" starting at k = |V(G)| (assuming G is the
// smaller graph, per main.cpp's convention) and decrementing k until one is
// found. Everything else - bidomain representation, branching heuristic,
// bound, filter_domains' multiway split - is identical to plain McSplit
// (mcsplit.cpp); only the pruning condition and the top-level driving loop
// differ, matching reference/mcsplit/mcsp.c's --big-first (-b) flag, which
// is what this file re-implements as a standalone algorithm rather than a
// runtime option.

// Forward declarations
static int calc_bound_down(const std::vector<Bidomain>& domains);
static int select_bidomain_down(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size);
static int find_min_value_down(const std::vector<int>& arr, int start_idx, int len);
static int partition_down(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
static void remove_vtx_from_left_domain_down(std::vector<int>& left, Bidomain& bd, int v);
static void remove_bidomain_down(std::vector<Bidomain>& domains, int idx);
static std::vector<Bidomain> filter_domains_down(const std::vector<Bidomain>& d, std::vector<int>& left, std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool multiway);
static void solve_down(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, bool multiway, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout);
static std::vector<Bidomain> build_domains_down(const Graph& g, const Graph& h, std::vector<int>& left, std::vector<int>& right);
std::vector<VtxPair> mcs_down(const Graph& g, const Graph& h, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout);

/**
 * McSplit's upper bound: sum of min(left_len, right_len) over all bidomains.
 * Identical to mcsplit.cpp's calc_bound - duplicated here per this project's
 * convention of each algorithm file being self-contained (see mcsplit-rl.cpp
 * etc. for the same pattern).
 */
static int calc_bound_down(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains)
        bound += std::min(bd.left_len, bd.right_len);
    return bound;
}

/**
 * Finds the smallest vertex index in arr[start_idx .. start_idx+len-1].
 * Used as select_bidomain_down's tiebreaker, identical to mcsplit.cpp.
 */
static int find_min_value_down(const std::vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i = 0; i < len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

/**
 * Selects the bidomain with smallest max(left_len, right_len), tie-broken
 * by smallest vertex index in the left set. Identical to mcsplit.cpp.
 */
static int select_bidomain_down(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value_down(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value_down(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

/**
 * Partitions all_vv[start .. start+len-1] in-place: adjacent vertices first,
 * non-adjacent after. Returns the length of the adjacent prefix. Identical
 * to mcsplit.cpp's partition.
 */
static int partition_down(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow) {
    int i = 0;
    for (int j = 0; j < len; j++) {
        if (adjrow[all_vv[start + j]]) {
            std::swap(all_vv[start + i], all_vv[start + j]);
            i++;
        }
    }
    return i;
}

/**
 * Removes vertex v from the left side of its bidomain (swap-with-last).
 * Identical to mcsplit.cpp.
 */
static void remove_vtx_from_left_domain_down(std::vector<int>& left, Bidomain& bd, int v) {
    int i = 0;
    while (left[bd.l + i] != v) i++;
    std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
    bd.left_len--;
}

/**
 * Removes a bidomain (swap-with-last, order doesn't matter). Identical to
 * mcsplit.cpp.
 */
static void remove_bidomain_down(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

/**
 * Splits bidomains after matching (v, w), with the same multiway
 * (direction/label-aware) split as mcsplit.cpp's filter_domains. Identical
 * logic - McSplit-down changes only the goal/pruning strategy, not how
 * domains are split.
 */
static std::vector<Bidomain> filter_domains_down(const std::vector<Bidomain>& d, std::vector<int>& left,
        std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool multiway) {
    std::vector<Bidomain> new_d;
    new_d.reserve(d.size() * 2);
    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        int left_len = partition_down(left, l, old_bd.left_len, g.adjmat[v]);
        int right_len = partition_down(right, r, old_bd.right_len, h.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
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
                    int lmin = li, rmin = ri;
                    do { li++; } while (li < l_top && adjrow_v[left[li]] == left_label);
                    do { ri++; } while (ri < r_top && adjrow_w[right[ri]] == left_label);
                    new_d.push_back({lmin, rmin, li - lmin, ri - rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }
    return new_d;
}

/**
 * Core branch-and-bound recursive search - identical to mcsplit.cpp's solve,
 * with one addition: once the incumbent already matches this iteration's
 * goal size, there is nothing left to search for (McSplit-down only asks
 * "does a mapping of size `goal` exist?", not "what is the largest
 * mapping?"), so the branch is cut immediately. This mirrors
 * reference/mcsplit/mcsp.c's `if (arguments.big_first && d.incumbent->len
 * == d.matching_size_goal) return;`.
 *
 * @param goal Exact solution size mcs_down is currently testing for.
 */
static void solve_down(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent,
        std::vector<VtxPair>& current, std::vector<Bidomain>& domains,
        std::vector<int>& left, std::vector<int>& right, unsigned int goal,
        bool multiway, Stats& stats,
        std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout) {

    if (abort_due_to_timeout) return;

    stats.nodes++;
    stats.bidomain_counts.push_back(domains.size());

    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best = stats.nodes;
        stats.time_to_best = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    int bound = current.size() + calc_bound_down(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    // McSplit-down-specific cut: this iteration only asks whether a mapping
    // of size `goal` exists. Once the incumbent has reached it, any further
    // search below this node can only reconfirm existence, not improve on
    // it, so there is nothing left worth exploring.
    if (incumbent.size() == goal) {
        stats.cut_branches++;
        return;
    }

    int bd_idx = select_bidomain_down(domains, left, current.size());
    if (bd_idx == -1) return;

    Bidomain& bd = domains[bd_idx];

    int v = find_min_value_down(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain_down(left, domains[bd_idx], v);

    int w = -1;
    bd.right_len--;

    for (int i = 0; i <= bd.right_len; i++) {
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
        std::swap(right[bd.r + idx], right[bd.r + bd.right_len]);

        auto new_domains = filter_domains_down(domains, left, right, g, h, v, w, multiway);
        current.push_back(VtxPair(v, w));
        solve_down(g, h, incumbent, current, new_domains, left, right, goal,
                multiway, stats, start_time, abort_due_to_timeout);
        current.pop_back();
    }

    bd.right_len++;
    if (bd.left_len == 0)
        remove_bidomain_down(domains, bd_idx);

    solve_down(g, h, incumbent, current, domains, left, right, goal,
            multiway, stats, start_time, abort_due_to_timeout);
}

/**
 * Builds the initial bidomain list for a fresh search at some goal size.
 * Split out from mcs_down (rather than inlined once, as in mcsplit.cpp's
 * mcs) because McSplit-down calls this once per goal value tried - left,
 * right and domains are all consumed/reordered during solve_down, so each
 * goal-size attempt needs its own fresh copy.
 *
 * @param g Graph G (sorted/relabelled copy).
 * @param h Graph H (sorted/relabelled copy).
 * @param left  Output: G-side vertex array, populated fresh (overwritten).
 * @param right Output: H-side vertex array, populated fresh (overwritten).
 * @return Initial list of bidomains, one per label shared by both graphs.
 */
static std::vector<Bidomain> build_domains_down(const Graph& g, const Graph& h,
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
        for (int i = 0; i < g.n; i++)
            if (g.label[i] == label)
                left.push_back(i);
        for (int i = 0; i < h.n; i++)
            if (h.label[i] == label)
                right.push_back(i);
        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        if (left_len && right_len)
            domains.push_back({start_l, start_r, left_len, right_len, false});
    }
    return domains;
}

/**
 * Entry point for McSplit-down search (Section 2.6).
 *
 * Sorts vertices by degree exactly as mcsplit.cpp's mcs() does, then - since
 * main.cpp's convention guarantees |V(g)| <= |V(h)| - tries goal sizes
 * k = |V(g)|, |V(g)|-1, |V(g)|-2, ... in turn, rebuilding fresh bidomains
 * each time (see build_domains_down), stopping as soon as a mapping of the
 * current goal size is found. This matches reference/mcsplit/mcsp.c's
 * --big-first driving loop in mcs().
 *
 * @param g Graph G (must satisfy |V(g)| <= |V(h)|; enforced by main.cpp).
 * @param h Graph H.
 * @param multiway Whether to use label/direction-aware domain splitting.
 * @param stats Search statistics, updated in-place across all goal attempts.
 * @param abort_due_to_timeout Flag checked to abort the search early.
 * @return The largest common subgraph solution found, in original graph indices.
 */
std::vector<VtxPair> mcs_down(const Graph& g, const Graph& h, bool multiway,
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

    std::vector<VtxPair> incumbent;
    auto start_time = std::chrono::steady_clock::now();

    // Try goal sizes |V(g)|, |V(g)|-1, ... down to 1, stopping as soon as a
    // mapping of that size is found. g_sorted.n == g.n (sorting only
    // reorders vertices, it doesn't drop any), so this still relies on
    // main.cpp's g-is-the-smaller-graph convention.
    for (int k = 0; k < g_sorted.n; k++) {
        if (abort_due_to_timeout) break;

        unsigned int goal = g_sorted.n - k;
        std::vector<int> left, right;
        std::vector<Bidomain> domains = build_domains_down(g_sorted, h_sorted, left, right);
        std::vector<VtxPair> current;

        solve_down(g_sorted, h_sorted, incumbent, current, domains, left, right,
                goal, multiway, stats, start_time, abort_due_to_timeout);

        if (incumbent.size() == goal) break;
    }

    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}
