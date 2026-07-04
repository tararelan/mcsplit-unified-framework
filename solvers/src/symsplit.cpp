// #define _GNU_SOURCE

#include "graph.h"
#include "common.h"
#include <algorithm>
#include <numeric>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <atomic>
#include <limits.h>
#include <boost/functional/hash.hpp>

// Forward declarations
static int calc_bound_sym(const std::vector<Bidomain>& domains);
static void remove_bidomain_sym(std::vector<Bidomain>& domains, int idx);
static int select_bidomain_sym(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size);
static int find_min_value_sym(const std::vector<int>& arr, int start_idx, int len);
static int index_of_next_smallest_sym(const std::vector<int>& arr, int start_idx, int len, int w);
int partition_sym(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
int partition_right_sym(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow, std::vector<int>& index_right);
int partition_sparse_sym(std::vector<int>& all_vv, int start, int len, int degree, const unsigned int* adjlist, std::vector<int>& index_right);
std::vector<Bidomain> filter_domains_sym(const std::vector<Bidomain>& d, std::vector<int>& left, std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool& best_match, std::vector<int>& index_right);
int find_vertices_with_common_neighbors(const Graph& g, std::vector<int>& eqn_classes);
bool break_h_sym(const std::vector<int>& arr, int start_idx, int len, int w, const std::vector<int>& h_eqn_classes);
void solve_sym(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout, const std::vector<int>& g_eqn_classes, const std::vector<int>& h_eqn_classes, std::vector<int>& index_right);
std::vector<VtxPair> mcs_sym(const Graph& g, const Graph& h, Stats& stats, std::atomic<bool>& abort_due_to_timeout);

static int calc_bound_sym(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

static void remove_bidomain_sym(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

static int find_min_value_sym(const std::vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i = 0; i < len; i++) {
        if (arr[start_idx + i] < min_v) {
            min_v = arr[start_idx + i];
        }
    }
    return min_v;
}

static int select_bidomain_sym(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size) {
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i = 0; i < domains.size(); i++) {
        const Bidomain& bd = domains[i];
        int len = std::max(bd.left_len, bd.right_len);
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value_sym(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value_sym(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

static int index_of_next_smallest_sym(const std::vector<int>& arr, int start_idx, int len, int w) {
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

int partition_sym(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow) {
    int i = 0;
    for (int j = 0; j < len; j++) {
        if (adjrow[all_vv[start + j]]) {
            std::swap(all_vv[start + i], all_vv[start + j]);
            i++;
        }
    }
    return i;
}

int partition_right_sym(std::vector<int>& all_vv, int start, int len,
        const std::vector<unsigned int>& adjrow, std::vector<int>& index_right) {
    int i = 0;
    for (int j = 0; j < len; j++) {
        if (adjrow[all_vv[start + j]]) {
            std::swap(index_right[all_vv[start + i]], index_right[all_vv[start + j]]);
            std::swap(all_vv[start + i], all_vv[start + j]);
            i++;
        }
    }
    return i;
}

int partition_sparse_sym(std::vector<int>& all_vv, int start, int len,
        int degree, const unsigned int* adjlist, std::vector<int>& index_right) {
    int j = 0;
    for (int i = 0; i < degree; i++) {
        int pos = index_right[adjlist[i]];
        if (pos >= start && pos < start + len) {
            std::swap(index_right[all_vv[start + j]], index_right[all_vv[pos]]);
            std::swap(all_vv[start + j], all_vv[pos]);
            j++;
        }
    }
    return j;
}

std::vector<Bidomain> filter_domains_sym(const std::vector<Bidomain>& d, std::vector<int>& left,
        std::vector<int>& right, const Graph& g, const Graph& h, int v, int w,
        bool& best_match, std::vector<int>& index_right) {
    std::vector<Bidomain> new_d;
    new_d.reserve(d.size() * 2);
    unsigned int ccount = 0;

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;

        int left_len = partition_sym(left, l, old_bd.left_len, g.adjmat[v]);

        int right_len;
        if (old_bd.right_len > (int)h.degree[w]) {
            right_len = partition_sparse_sym(right, r, old_bd.right_len, h.degree[w], h.adjlist[w], index_right);
        } else {
            right_len = partition_right_sym(right, r, old_bd.right_len, h.adjmat[w], index_right);
        }

        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;

        if ((left_len == 0 && right_len == 0) ||
            (left_len_noedge == 0 && right_len_noedge == 0) ||
            old_bd.left_len == 0) {
            ccount++;
        }

        if (left_len_noedge && right_len_noedge) {
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        }
        if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }

    best_match = (ccount == d.size());
    return new_d;
}

// Replicates boost::hash_combine exactly.
static void hash_combine(std::size_t& seed, std::size_t v) {
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Computes symmetry equivalence classes using neighbourhood hashing.
// Replicates the reference SymSplit implementation exactly, replacing
// boost::hash_combine with an equivalent inline implementation.
// Each vertex v contributes two hashes to the map:
//   n_hash: hash of v's neighbours only
//   p_hash: hash of v's neighbours plus v itself
// Vertices sharing either hash are placed in the same group.
// eqn_classes[v] = -1 if v has no detected symmetry partner.
int find_vertices_with_common_neighbors(const Graph& g, std::vector<int>& eqn_classes) {
    std::unordered_map<std::size_t, std::vector<int>> n_groups;
    n_groups.reserve(g.n);
    eqn_classes.assign(g.n, -1);

    for (int v = 0; v < g.n; v++) {
        std::size_t n_hash = 0, p_hash = 0;
        for (int u = 0; u < g.n; u++) {
            if (v != u && g.adjmat[v][u]) {
                boost::hash_combine(n_hash, (unsigned int)u);
                boost::hash_combine(p_hash, (unsigned int)u);
            }
            if (v == u) {
                boost::hash_combine(p_hash, (unsigned int)v);
            }
        }
        n_groups[p_hash].push_back(v);
        n_groups[n_hash].push_back(v);
    }

    int n_syms = 0;
    int label = 0;
    for (const auto& batch : n_groups) {
        if (batch.second.size() > 1) {
            for (int v : batch.second) {
                eqn_classes[v] = label;
                n_syms++;
            }
        }
        label++;
    }
    return n_syms;
}

// Value symmetry check for H.
// Returns true if a smaller equivalent H-vertex exists in the right domain,
// meaning matching v to w would produce an isomorphic solution.
bool break_h_sym(const std::vector<int>& arr, int start_idx, int len, int w,
        const std::vector<int>& h_eqn_classes) {
    int w_class = h_eqn_classes[w];
    for (int i = 0; i < len; i++) {
        if (h_eqn_classes[arr[start_idx + i]] == w_class && arr[start_idx + i] < w) {
            return true;
        }
    }
    return false;
}

// Core BnB recursive search for SymSplit.
// Extends RRSplit with dual symmetry breaking:
// 1. Variable symmetry (G-side): skip w <= largest w matched to equivalent G-vertex,
//    remove equivalent G-vertices from left domain in exclude branch.
// 2. Value symmetry (H-side): skip w if smaller equivalent H-vertex exists in domain.
// 3. Maximality reduction: prune if filter_domains_sym sets best_match.
// nodes++ at top before bound check, consistent with SymSplit reference.
void solve_sym(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent,
        std::vector<VtxPair>& current, std::vector<Bidomain>& domains,
        std::vector<int>& left, std::vector<int>& right, unsigned int goal,
        Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout,
        const std::vector<int>& g_eqn_classes, const std::vector<int>& h_eqn_classes,
        std::vector<int>& index_right) {

    if (abort_due_to_timeout) { return; }

    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best = stats.nodes;
        stats.time_to_best = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    stats.nodes++;

    int bound = current.size() + calc_bound_sym(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    int bd_idx = select_bidomain_sym(domains, left, current.size());
    if (bd_idx == -1) { return; }

    Bidomain& bd = domains[bd_idx];

    int v = find_min_value_sym(left, bd.l, bd.left_len);
    {
        int i = 0;
        while (left[bd.l + i] != v) { i++; }
        std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
        bd.left_len--;
    }

    int w_lower = -1;
    if (g_eqn_classes[v] != -1) {
        for (const VtxPair& a : current) {
            if (g_eqn_classes[a.v] == g_eqn_classes[v] && w_lower < a.w) {
                w_lower = a.w;
            }
        }
    }

    bd.right_len--;
    int w = w_lower;

    bool best_match = false;
    bool skip_exclude = false;  // ← NEW

    for (int i = bd.right_len; i >= 0; --i) {
        int idx = index_of_next_smallest_sym(right, bd.r, bd.right_len + 1, w);
        if (idx == -1) { break; }
        w = right[bd.r + idx];

        if (h_eqn_classes[w] != -1 && break_h_sym(right, bd.r, bd.right_len + 1, w, h_eqn_classes)) {
            stats.sym_pruned++;
            continue;
        }

        std::swap(index_right[w], index_right[right[bd.r + bd.right_len]]);
        right[bd.r + idx] = right[bd.r + bd.right_len];
        right[bd.r + bd.right_len] = w;

        auto new_domains = filter_domains_sym(domains, left, right, g, h, v, w,
                best_match, index_right);
        current.push_back(VtxPair(v, w));
        solve_sym(g, h, incumbent, current, new_domains, left, right, goal,
                stats, start_time, abort_due_to_timeout,
                g_eqn_classes, h_eqn_classes, index_right);
        current.pop_back();

        if (best_match || bound <= (int)incumbent.size()) {
            stats.sym_pruned++;
            skip_exclude = true;  // ← NEW: skip exclude branch, matching reference
            break;
        }
    }

    bd.right_len++;

    if (g_eqn_classes[v] != -1) {
        for (int i = 0; i < bd.left_len; i++) {
            if (g_eqn_classes[left[bd.l + i]] == g_eqn_classes[v]) {
                std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
                bd.left_len--;
                i--;
                stats.sym_pruned++;
            }
        }
    }

    if (bd.left_len == 0) { remove_bidomain_sym(domains, bd_idx); }

    if (!skip_exclude) {  // ← NEW: only explore exclude branch if not pruned
        solve_sym(g, h, incumbent, current, domains, left, right, goal,
                stats, start_time, abort_due_to_timeout,
                g_eqn_classes, h_eqn_classes, index_right);
    }
}

// Entry point for SymSplit search.
// Computes symmetry classes on unsorted graphs (matching reference),
// sorts vertices by degree, builds single initial bidomain, runs solve_sym().
std::vector<VtxPair> mcs_sym(const Graph& g, const Graph& h,
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

    auto sum_vec = [](const std::vector<int>& v) {
        int s = 0; for (int x : v) { s += x; } return s;
    };

    bool g1_dense = sum_vec(h_deg) < h.n * (h.n - 1);
    bool g0_dense = sum_vec(g_deg) < g.n * (g.n - 1);

    std::vector<int> vv0(g.n), vv1(h.n);
    std::iota(vv0.begin(), vv0.end(), 0);
    std::iota(vv1.begin(), vv1.end(), 0);
    std::stable_sort(vv0.begin(), vv0.end(), [&](int a, int b) {
        return !g1_dense ? (g_deg[a] < g_deg[b]) : (g_deg[a] > g_deg[b]);
    });
    std::stable_sort(vv1.begin(), vv1.end(), [&](int a, int b) {
        return !g0_dense ? (h_deg[a] < h_deg[b]) : (h_deg[a] > h_deg[b]);
    });

    // Compute symmetry classes on unsorted graphs — matches reference which
    // calls find_vertices_with_common_neighbors before induced_subgraph
    std::vector<int> g_eqn_classes, h_eqn_classes;
    find_vertices_with_common_neighbors(g, g_eqn_classes);
    find_vertices_with_common_neighbors(h, h_eqn_classes);

    // Remap equivalence classes through the sort permutation so they align
    // with the sorted vertex indices used during search
    std::vector<int> g_eqn_sorted(g.n), h_eqn_sorted(h.n);
    for (int i = 0; i < g.n; i++) { g_eqn_sorted[i] = g_eqn_classes[vv0[i]]; }
    for (int i = 0; i < h.n; i++) { h_eqn_sorted[i] = h_eqn_classes[vv1[i]]; }

    Graph g_sorted = induced_subgraph(const_cast<Graph&>(g), vv0);
    Graph h_sorted = induced_subgraph(const_cast<Graph&>(h), vv1);

    set_adjlist(g_sorted);
    set_adjlist(h_sorted);

    std::vector<int> index_right(h_sorted.n);
    for (int i = 0; i < h_sorted.n; i++) { index_right[i] = i; }

    std::vector<int> left, right;
    for (int i = 0; i < g_sorted.n; i++) { left.push_back(i); }
    for (int i = 0; i < h_sorted.n; i++) { right.push_back(i); }
    std::vector<Bidomain> domains;
    domains.push_back({0, 0, g_sorted.n, h_sorted.n, false});

    stats.root_upper_bound = calc_bound_sym(domains);

    std::vector<VtxPair> incumbent, current;
    auto start_time = std::chrono::steady_clock::now();
    solve_sym(g_sorted, h_sorted, incumbent, current, domains, left, right, 1,
            stats, start_time, abort_due_to_timeout,
            g_eqn_sorted, h_eqn_sorted, index_right);

    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}