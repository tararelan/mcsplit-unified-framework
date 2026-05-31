#define _GNU_SOURCE

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
int calc_bound(const std::vector<Bidomain>& domains);
int select_bidomain(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size);
int find_min_value(const std::vector<int>& arr, int start_idx, int len);
int partition(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
void remove_vtx_from_left_domain(std::vector<int>& left, Bidomain& bd, int v);
void remove_bidomain(std::vector<Bidomain>& domains, int idx);
std::vector<Bidomain> filter_domains(const std::vector<Bidomain>& d, std::vector<int>& left, std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool multiway);
void solve(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs(const Graph& g, const Graph& h, bool multiway, Stats& stats, std::atomic<bool>& abort_due_to_timeout);

// McSplit's upper bound, summing min(left_len, right_len) over all label classes.
int calc_bound(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains)
        bound += std::min(bd.left_len, bd.right_len);
    return bound;
}

// Finds smallest vertex index in a section of the left buffer, used as tiebreaker in bidomain selection.
int find_min_value(const std::vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i = 0; i < len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

// Selects label class with smallest max(left_len, right_len), breaking ties on smallest vertex index.
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

// Partitions arr[start..start+len] into adjacent (front) and non-adjacent (back). Returns adjacent count.
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

// Removes vertex v from left side of its bidomain.
void remove_vtx_from_left_domain(std::vector<int>& left, Bidomain& bd, int v) {
    int i = 0;
    while (left[bd.l + i] != v) i++;
    std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
    bd.left_len--;
}

// Removes bidomain from list by replacing with last element.
void remove_bidomain(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

// Splits each bidomain into adjacent and non-adjacent parts after matching (v,w).
// In multiway mode (labelled/directed), further splits adjacent part by edge label.
std::vector<Bidomain> filter_domains(const std::vector<Bidomain>& d, std::vector<int>& left,
        std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool multiway) {
    std::vector<Bidomain> new_d;
    new_d.reserve(d.size());
    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        int left_len = partition(left, l, old_bd.left_len, g.adjmat[v]);
        int right_len = partition(right, r, old_bd.right_len, h.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
            // Sort by edge label and create separate bidomains per distinct label value
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

// Core BnB recursive search.
void solve(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, bool multiway, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout) {

    if (abort_due_to_timeout) return;

    stats.nodes++;

    if (current.size() > incumbent.size()) {
        incumbent = current;
        stats.incumbent_size = incumbent.size();
        stats.nodes_to_best = stats.nodes;
        stats.time_to_best = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
    }

    int bound = current.size() + calc_bound(domains);
    if (bound <= (int)incumbent.size() || bound < (int)goal) {
        stats.cut_branches++;
        stats.bound_pruned++;
        return;
    }

    int bd_idx = select_bidomain(domains, left, current.size());
    if (bd_idx == -1) return;

    Bidomain& bd = domains[bd_idx];

    int v = find_min_value(left, bd.l, bd.left_len);
    remove_vtx_from_left_domain(left, domains[bd_idx], v);

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

        auto new_domains = filter_domains(domains, left, right, g, h, v, w, multiway);
        current.push_back(VtxPair(v, w));
        solve(g, h, incumbent, current, new_domains, left, right, goal,
                multiway, stats, start_time, abort_due_to_timeout);
        current.pop_back();
    }

    bd.right_len++;
    if (bd.left_len == 0)
        remove_bidomain(domains, bd_idx);

    solve(g, h, incumbent, current, domains, left, right, goal,
            multiway, stats, start_time, abort_due_to_timeout);
}

std::vector<VtxPair> mcs(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout) {

    // Calculate degrees for vertex sorting
    auto calc_degrees = [](const Graph& g) {
        std::vector<int> degree(g.n, 0);
        for (int v = 0; v < g.n; v++)
            for (int w = 0; w < g.n; w++)
                if (g.adjmat[v][w] & 1) degree[v]++;
        return degree;
    };

    std::vector<int> g_deg = calc_degrees(g);
    std::vector<int> h_deg = calc_degrees(h);

    int g_edges = 0, h_edges = 0;
    for (int d : g_deg) g_edges += d;
    for (int d : h_deg) h_edges += d;
    bool h_dense = h_edges > h.n * (h.n - 1) / 2;
    bool g_dense = g_edges > g.n * (g.n - 1) / 2;

    // Sort vertices: ascending degree if graph is sparse, descending if dense
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

    // Build initial label class domains
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

    stats.root_upper_bound = calc_bound(domains);

    std::vector<VtxPair> incumbent, current;
	auto start_time = std::chrono::steady_clock::now();
    solve(g_sorted, h_sorted, incumbent, current, domains, left, right, 1, multiway, stats, start_time, abort_due_to_timeout);

    // Convert back to original vertex indices
    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}