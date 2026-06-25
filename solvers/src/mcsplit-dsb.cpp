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

// DSB constants
static const int DSB_THRESHOLD = 16;
static const int DSB_INDEX_PARTITION_SIZE = 2;
static const int DSB_BEST_PARTITION_UB = DSB_THRESHOLD;
static const int DSB_BEST_PARTITION_DIFF = 10;
static const int DSB_ROWS = DSB_THRESHOLD - DSB_INDEX_PARTITION_SIZE + 2;
static const int DSB_BOXES = 2;

// Forward declarations
static void dsb_init_array(long double* arr, int size);
static unsigned int calc_bound_dsb(const Graph& g, const Graph& h, std::vector<Bidomain>& domains, const std::vector<int>& left, const std::vector<int>& right, const std::vector<VtxPair>& incumbent, const std::vector<VtxPair>& current, long double* array, bool& bound_enabled, bool& count_enabled, int& called_count, int& success_count, int max_count, int max_success, int* s_degrees);
static int find_min_value_dsb(const std::vector<int>& arr, int start_idx, int len);
static int select_bidomain_dsb(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size);
static int partition_dsb(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
static std::vector<Bidomain> filter_domains_dsb(const std::vector<Bidomain>& d, std::vector<int>& left, std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool multiway);
static int index_of_next_smallest_dsb(const std::vector<int>& arr, int start_idx, int len, int w);
static void remove_vtx_from_left_domain_dsb(std::vector<int>& left, Bidomain& bd, int v);
static void remove_bidomain_dsb(std::vector<Bidomain>& domains, int idx);
void solve_dsb(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, bool multiway, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout, long double* array, bool& bound_enabled, bool& count_enabled, int& called_count, int& success_count, int max_count, int max_success, int* s_degrees);
std::vector<VtxPair> mcs_dsb(const Graph& g, const Graph& h, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout);

static void dsb_init_array(long double* arr, int size) {
    for (int i = 0; i < size; i++) { arr[i] = 1.0; }
}

// DSB's tighter upper bound.
// For small bidomains (size <= THRESHOLD, diff <= BEST_PARTITION_DIFF),
// analyses internal edge structure to find a tighter bound than min(left, right).
// Uses a learned probability table to decide when to spend time on the tighter computation.
// For large or unbalanced bidomains, falls back to McSplit's min(left, right) bound.
static unsigned int calc_bound_dsb(const Graph& g, const Graph& h,
        std::vector<Bidomain>& domains,
        const std::vector<int>& left, const std::vector<int>& right,
        const std::vector<VtxPair>& incumbent, const std::vector<VtxPair>& current,
        long double* array, bool& bound_enabled, bool& count_enabled,
        int& called_count, int& success_count, int max_count, int max_success,
        int* s_degrees) {

    unsigned int incumbent_size = incumbent.size();
    unsigned int current_size = current.size();
    unsigned int mcsplit_bound = current_size;
    long double possible_reductions = 0;

    if (!bound_enabled) {
        for (const Bidomain& bd : domains) {
            mcsplit_bound += std::min(bd.left_len, bd.right_len);
        }
        return mcsplit_bound;
    }

    for (const Bidomain& bd : domains) {
        const auto [min_len, max_len] = std::minmax(bd.left_len, bd.right_len);
        const int size_diff = std::abs(bd.left_len - bd.right_len);
        if (min_len > 2 && max_len <= DSB_BEST_PARTITION_UB && size_diff <= DSB_BEST_PARTITION_DIFF) {
            if (bd.is_valid) {
                possible_reductions = min_len - bd.size;
            } else {
                possible_reductions += array[((max_len - DSB_INDEX_PARTITION_SIZE) * DSB_INDEX_PARTITION_SIZE) + (size_diff > 2 ? 1 : 0)];
            }
        }
        mcsplit_bound += min_len;
    }

    if (mcsplit_bound <= incumbent_size ||
        (possible_reductions < 0.5) ||
        ((mcsplit_bound - possible_reductions) > incumbent_size)) {
        return mcsplit_bound;
    }

    if (count_enabled && called_count == max_count) {
        count_enabled = false;
        if (success_count < max_success) {
            bound_enabled = false;
        }
        return mcsplit_bound;
    }

    unsigned int bound = current_size;
    for (Bidomain& bd : domains) {
        const auto [min_len, max_len] = std::minmax(bd.left_len, bd.right_len);
        const int size_diff = std::abs(bd.left_len - bd.right_len);

        if (bd.is_valid) {
            bound += bd.size;
            continue;
        }

        // Special case: bidomain of size 2x2 — check if the two G-vertices
        // are adjacent and if the two H-vertices are adjacent. If edge structure
        // differs, at most 1 can be matched.
        if (bd.left_len == 2 && bd.right_len == 2) {
            if (g.adjmat[left[bd.l]][left[bd.l + 1]] != h.adjmat[right[bd.r]][right[bd.r + 1]]) {
                bound += 1;
                bd.size = 1;
                int index = (max_len - DSB_INDEX_PARTITION_SIZE) * DSB_INDEX_PARTITION_SIZE;
                array[index] = (array[index] + 1) / 2;
            } else {
                bd.size = 2;
                int index = (max_len - DSB_INDEX_PARTITION_SIZE) * DSB_INDEX_PARTITION_SIZE;
                array[index] = array[index] / 2;
                bound += 2;
            }
            bd.is_valid = true;
            continue;
        }

        // Graph partition bound for medium bidomains.
        // Compares internal edge counts of both sides. If one side has more
        // edges than the other, vertices with the highest degree must be
        // excluded to make a valid matching possible, reducing the bound.
        if (min_len > 2 && max_len <= DSB_BEST_PARTITION_UB && size_diff <= DSB_BEST_PARTITION_DIFF) {
            int s_len, b_len;
            int val = min_len;
            int index = (max_len - DSB_INDEX_PARTITION_SIZE) * DSB_INDEX_PARTITION_SIZE + (size_diff > 2 ? 1 : 0);
            int s_connections = 0, b_connections = 0;

            std::fill(s_degrees, s_degrees + min_len, 0);

            if (bd.left_len <= bd.right_len) {
                s_len = bd.left_len;
                b_len = bd.right_len;
                for (int i = 0; i < s_len; i++) {
                    const int v1 = left[bd.l + i];
                    for (int j = i + 1; j < s_len; j++) {
                        const int v2 = left[bd.l + j];
                        if (g.adjmat[v1][v2]) {
                            s_degrees[i]++;
                            s_degrees[j]++;
                            s_connections++;
                        }
                    }
                }
                for (int i = 0; i < b_len; i++) {
                    const int u1 = right[bd.r + i];
                    for (int j = i + 1; j < b_len; j++) {
                        const int u2 = right[bd.r + j];
                        if (h.adjmat[u1][u2]) { b_connections++; }
                    }
                }
            } else {
                s_len = bd.right_len;
                b_len = bd.left_len;
                for (int i = 0; i < s_len; i++) {
                    const int v1 = right[bd.r + i];
                    for (int j = i + 1; j < s_len; j++) {
                        const int v2 = right[bd.r + j];
                        if (h.adjmat[v1][v2]) {
                            s_degrees[i]++;
                            s_degrees[j]++;
                            s_connections++;
                        }
                    }
                }
                for (int i = 0; i < b_len; i++) {
                    const int u1 = left[bd.l + i];
                    for (int j = i + 1; j < b_len; j++) {
                        const int u2 = left[bd.l + j];
                        if (g.adjmat[u1][u2]) { b_connections++; }
                    }
                }
            }

            int s_non_connections = s_len * (s_len - 1) / 2 - s_connections;
            int b_non_connections = b_len * (b_len - 1) / 2 - b_connections;

            int diff = s_connections - b_connections;
            if (diff > 0) {
                std::sort(s_degrees, s_degrees + s_len);
                diff -= s_degrees[--val];
                while (diff > 0) { diff -= s_degrees[--val]; }
            } else {
                diff = s_non_connections - b_non_connections;
                if (diff > 0) {
                    std::sort(s_degrees, s_degrees + s_len, std::greater<int>());
                    diff -= s_len - s_degrees[--val];
                    while (diff > 0) { diff -= s_len - s_degrees[--val]; }
                }
            }

            array[index] = (array[index] + min_len - val) / 2;
            bound += val;
            bd.size = val;
        } else {
            bound += min_len;
            bd.size = min_len;
        }
        bd.is_valid = true;
    }

    if (count_enabled) {
        if (bound <= incumbent_size) { success_count++; }
        called_count++;
    }
    return bound;
}

static int find_min_value_dsb(const std::vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i = 0; i < len; i++) {
        if (arr[start_idx + i] < min_v) { min_v = arr[start_idx + i]; }
    }
    return min_v;
}

static int select_bidomain_dsb(const std::vector<Bidomain>& domains,
        const std::vector<int>& left, int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value_dsb(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value_dsb(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

static int partition_dsb(std::vector<int>& all_vv, int start, int len,
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

// Splits bidomains after matching (v, w).
// Key difference from McSplit: when a bidomain splits cleanly (all vertices
// go to one side), the cached bound (size, is_valid) can be carried over to
// the child bidomain since the internal structure is unchanged.
// When the split mixes vertices, the cache is invalidated.
static std::vector<Bidomain> filter_domains_dsb(const std::vector<Bidomain>& d,
        std::vector<int>& left, std::vector<int>& right,
        const Graph& g, const Graph& h, int v, int w, bool multiway) {
    std::vector<Bidomain> new_d;
    new_d.reserve(d.size());

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        int left_len = partition_dsb(left, l, old_bd.left_len, g.adjmat[v]);
        int right_len = partition_dsb(right, r, old_bd.right_len, h.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;

        if (left_len_noedge && right_len_noedge) {
            if (left_len != 0 || right_len != 0) {
                // Split mixed — cache invalid
                new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent, -1, false});
            } else {
                // All vertices went to non-adjacent side — carry cache over
                new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent, old_bd.size, old_bd.is_valid});
            }
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
                unsigned int left_label = adjrow_v[left[li]];
                unsigned int right_label = adjrow_w[right[ri]];
                if (left_label < right_label) { li++; }
                else if (left_label > right_label) { ri++; }
                else {
                    int lmin = li, rmin = ri;
                    do { li++; } while (li < l_top && adjrow_v[left[li]] == left_label);
                    do { ri++; } while (ri < r_top && adjrow_w[right[ri]] == left_label);
                    new_d.push_back({lmin, rmin, li - lmin, ri - rmin, true, -1, false});
                }
            }
        } else if (left_len && right_len) {
            if (left_len_noedge != 0 || right_len_noedge != 0) {
                // Split mixed — cache invalid
                new_d.push_back({l, r, left_len, right_len, true, -1, false});
            } else {
                // All vertices went to adjacent side — carry cache over
                new_d.push_back({l, r, left_len, right_len, true, old_bd.size, old_bd.is_valid});
            }
        }
    }
    return new_d;
}

static int index_of_next_smallest_dsb(const std::vector<int>& arr, int start_idx, int len, int w) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i = 0; i < len; i++) {
        if (arr[start_idx + i] > w && arr[start_idx + i] < smallest) {
            smallest = arr[start_idx + i];
            idx = i;
        }
    }
    return idx;
}

static void remove_vtx_from_left_domain_dsb(std::vector<int>& left, Bidomain& bd, int v) {
    int i = 0;
    while (left[bd.l + i] != v) { i++; }
    std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
    bd.left_len--;
}

static void remove_bidomain_dsb(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

// Core BnB search for McSplit+DSB.
// Identical structure to McSplit except calc_bound_dsb replaces calc_bound.
// The DSB bound is tighter on small balanced bidomains — it can prune branches
// that McSplit's bound cannot.
void solve_dsb(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent,
        std::vector<VtxPair>& current, std::vector<Bidomain>& domains,
        std::vector<int>& left, std::vector<int>& right, unsigned int goal,
        bool multiway, Stats& stats,
        std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout,
        long double* array, bool& bound_enabled, bool& count_enabled,
        int& called_count, int& success_count, int max_count, int max_success,
        int* s_degrees) {

    if (abort_due_to_timeout) { return; }

    stats.nodes++;

    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best = stats.nodes;
        stats.time_to_best = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    unsigned int bound = calc_bound_dsb(g, h, domains, left, right, incumbent, current,
            array, bound_enabled, count_enabled, called_count, success_count,
            max_count, max_success, s_degrees);

    if (bound <= (unsigned int)incumbent.size() || bound < goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    int bd_idx = select_bidomain_dsb(domains, left, current.size());
    if (bd_idx == -1) { return; }

    Bidomain& bd = domains[bd_idx];

    int v = find_min_value_dsb(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain_dsb(left, domains[bd_idx], v);

    int w = -1;
    bd.right_len--;

    for (int i = 0; i <= bd.right_len; i++) {
        int idx = index_of_next_smallest_dsb(right, bd.r, bd.right_len + 1, w);
        w = right[bd.r + idx];

        right[bd.r + idx] = right[bd.r + bd.right_len];
        right[bd.r + bd.right_len] = w;

        auto new_domains = filter_domains_dsb(domains, left, right, g, h, v, w, multiway);
        current.push_back(VtxPair(v, w));
        solve_dsb(g, h, incumbent, current, new_domains, left, right, goal,
                multiway, stats, start_time, abort_due_to_timeout,
                array, bound_enabled, count_enabled, called_count, success_count,
                max_count, max_success, s_degrees);
        current.pop_back();
    }

    bd.right_len++;
    if (bd.left_len == 0) { remove_bidomain_dsb(domains, bd_idx); }

    solve_dsb(g, h, incumbent, current, domains, left, right, goal,
            multiway, stats, start_time, abort_due_to_timeout,
            array, bound_enabled, count_enabled, called_count, success_count,
            max_count, max_success, s_degrees);
}

// Entry point for McSplit+DSB.
// Initialises the probability table, sorts vertices, builds per-label bidomains,
// runs solve_dsb().
std::vector<VtxPair> mcs_dsb(const Graph& g, const Graph& h, bool multiway,
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

    // Build per-label bidomains — DSB uses label splitting like McSplit
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
            domains.push_back({start_l, start_r, left_len, right_len, false, -1, false});
        }
    }

    // Compute root upper bound using McSplit's simple bound (before DSB tightening)
    int root_ub = 0;
    for (const Bidomain& bd : domains) {
        root_ub += std::min(bd.left_len, bd.right_len);
    }
    stats.root_upper_bound = root_ub;

    // Initialise DSB probability table and adaptive bound control
    long double array[DSB_ROWS * DSB_BOXES];
    dsb_init_array(array, DSB_ROWS * DSB_BOXES);
    bool bound_enabled = true;
    bool count_enabled = true;
    int called_count = 0;
    int success_count = 0;
    const int max_count = 100;
    const int max_success = 30;
    int s_degrees[DSB_BEST_PARTITION_UB];

    std::vector<VtxPair> incumbent, current;
    auto start_time = std::chrono::steady_clock::now();
    solve_dsb(g_sorted, h_sorted, incumbent, current, domains, left, right, 1,
            multiway, stats, start_time, abort_due_to_timeout,
            array, bound_enabled, count_enabled, called_count, success_count,
            max_count, max_success, s_degrees);

    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}