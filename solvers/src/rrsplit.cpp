// #define _GNU_SOURCE

#include "graph.h"
#include "common.h"
#include <algorithm>
#include <numeric>
#include <vector>
#include <chrono>
#include <atomic>
#include <limits.h>
#include <iostream>

// NOTE: RRSplit does not implement the multiway (direction/label-aware)
// split used by the other six algorithms - see Section 4.1.2.2. A fix was
// attempted (patching filter_domains_rr with the same split used
// elsewhere) and confirmed sound via brute-force testing, but it
// surfaced a deeper issue: EqClass (structural equivalence, see graph.cpp)
// is itself computed without regard to edge direction, so the equivalence
// reduction can discard vertices that are not actually redundant on
// directed graphs, producing valid-but-suboptimal solutions. Fixing this
// properly would require making EqClass direction-aware, which was out
// of scope for this project's timeline. RRSplit is validated and used
// for undirected datasets (MIVIA, LV) only; it is excluded from BI.

// Forward declarations
static int calc_bound(const std::vector<Bidomain>& domains);
static void remove_bidomain(std::vector<Bidomain>& domains, int idx);
static int select_bidomain(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size);
static int find_min_value(const std::vector<int>& arr, int start_idx, int len);
static int index_of_next_smallest(const std::vector<int>& arr, int start_idx, int len, int w);
int partition_rr(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow);
int partition_right_rr(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow, std::vector<int>& index_right);
int partition_sparse_rr(std::vector<int>& all_vv, int start, int len, int degree, const unsigned int* adjlist, std::vector<int>& index_right);
std::vector<Bidomain> filter_domains_rr(const std::vector<Bidomain>& d, std::vector<int>& left, std::vector<int>& right, const Graph& g, const Graph& h, int v, int w, bool& best_match, std::vector<int>& index_right);
void solve_rr(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent, std::vector<VtxPair>& current, std::vector<Bidomain>& domains, std::vector<int>& left, std::vector<int>& right, unsigned int goal, Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool>& abort_due_to_timeout, const ui* EqClass, std::vector<int>& index_right, int red_mode);
std::vector<VtxPair> mcs_rr(const Graph& g, const Graph& h, Stats& stats, std::atomic<bool>& abort_due_to_timeout, int red_mode = 7);

/**
 * McSplit's upper bound: sum of min(left_len, right_len) over all bidomains
 *
 * @param domains Current list of bidomains.
 * @return Upper bound on the number of additional pairs that could be matched.
 */
static int calc_bound(const std::vector<Bidomain>& domains) {
    int bound = 0;
    for (const Bidomain& bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

/**
 * Removes a bidomain by replacing it with the last element and popping;
 * order within the list does not matter
 *
 * @param domains List of bidomains to remove from, modified in-place.
 * @param idx Index of the bidomain to remove.
 */
static void remove_bidomain(std::vector<Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

/**
 * Finds the smallest vertex index in arr[start_idx..start_idx+len-1];
 * used as select_bidomain's tie-break
 *
 * @param arr Array to search (typically the `left` vertex array).
 * @param start_idx Index of the first element of the slice to search.
 * @param len Number of elements in the slice.
 * @return The smallest value found in arr[start_idx .. start_idx+len-1].
 */
static int find_min_value(const std::vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i = 0; i < len; i++) {
        if (arr[start_idx + i] < min_v) {
            min_v = arr[start_idx + i];
        }
    }
    return min_v;
}

/**
 * Selects bidomain with smallest max(left_len, right_len), ties broken
 * by smallest vertex index in the left set - same rule as plain McSplit
 *
 * @param domains Current list of bidomains to choose from.
 * @param left G-side vertex array (used for the tie-break lookup).
 * @param current_matching_size Size of the partial solution so far (unused directly, kept for interface consistency).
 * @return Index into domains of the selected bidomain, or -1 if domains is empty.
 */
static int select_bidomain(const std::vector<Bidomain>& domains, const std::vector<int>& left, int current_matching_size) {
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

/**
 * Returns index of smallest value in arr[start_idx..start_idx+len-1] that is > w.
 * Used by solve_rr to visit right-side candidates in increasing order one at
 * a time without pre-sorting, same pattern as index_of_next_smallest_dsb in
 * mcsplit-dsb.cpp. Returns -1 if no such value exists (w is the maximum).
 *
 * @param arr Array to search (typically the `right` vertex array).
 * @param start_idx Index of the first element of the slice to search.
 * @param len Number of elements in the slice.
 * @param w Lower bound; the returned value must be strictly greater than w.
 * @return Index (relative to start_idx) of the smallest qualifying value, or -1 if none exists.
 */
static int index_of_next_smallest(const std::vector<int>& arr, int start_idx, int len, int w) {
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

/**
 * Boolean pre-partition (adjacent vs non-adjacent), G-side only - no
 * direction awareness (see Section 4.1.2.2, and the note at the top of
 * this file re: RRSplit's EqClass limitation)
 *
 * @param all_vv Vertex array to partition in-place (the `left` array).
 * @param start Index of the first element of the slice to partition.
 * @param len Number of elements in the slice.
 * @param adjrow Adjacency row of the vertex being matched, indexed by vertex id.
 * @return Number of adjacent vertices, i.e. the length of the adjacent prefix.
 */
int partition_rr(std::vector<int>& all_vv, int start, int len, const std::vector<unsigned int>& adjrow) {
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
 * H-side partition, additionally maintaining index_right (a reverse
 * lookup: index_right[vertex] = its current position in the right array)
 * so that other RRSplit machinery (e.g. index_of_next_smallest, EqClass
 * lookups) can find a vertex's position in O(1) without a linear scan
 *
 * @param all_vv Vertex array to partition in-place (the `right` array).
 * @param start Index of the first element of the slice to partition.
 * @param len Number of elements in the slice.
 * @param adjrow Adjacency row of the vertex being matched, indexed by vertex id.
 * @param index_right Reverse lookup (vertex -> current position in `right`), updated in-place.
 * @return Number of adjacent vertices, i.e. the length of the adjacent prefix.
 */
int partition_right_rr(std::vector<int>& all_vv, int start, int len,
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

/**
 * Sparse variant of partition_right_rr: when w's degree is small relative
 * to the bidomain size, it is cheaper to iterate w's adjacency list
 * directly (via adjlist) and look up each neighbour's current position
 * (via index_right) than to scan the whole bidomain checking adjrow
 *
 * @param all_vv Vertex array to partition in-place (the `right` array).
 * @param start Index of the first element of the slice to partition.
 * @param len Number of elements in the slice.
 * @param degree Degree of w, i.e. the length of adjlist.
 * @param adjlist w's adjacency list (neighbour vertex ids).
 * @param index_right Reverse lookup (vertex -> current position in `right`), updated in-place.
 * @return Number of adjacent vertices, i.e. the length of the adjacent prefix.
 */
int partition_sparse_rr(std::vector<int>& all_vv, int start, int len,
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

/**
 * Splits bidomains after matching (v, w). No multiway/direction-aware
 * split - see the note at the top of this file re: RRSplit's
 * direction-blindness (Section 4.1.2.2).
 *
 * best_match (via ccount) implements RRSplit's maximality reduction:
 * ccount counts, across all bidomains, how many were either fully
 * consumed (both sides split entirely to adjacent or entirely to
 * non-adjacent) or were already empty before this match. If every
 * bidomain satisfies this, best_match is set true, signalling to the
 * caller (solve_rr) that (v,w) is guaranteed part of an optimal
 * solution here - no other branch at this node needs exploring.
 *
 * @param d Bidomains to split, as they stood before matching (v, w).
 * @param left G-side vertex array, reordered in-place during partitioning.
 * @param right H-side vertex array, reordered in-place during partitioning.
 * @param g Graph G.
 * @param h Graph H.
 * @param v G-side vertex just matched.
 * @param w H-side vertex just matched.
 * @param best_match Set to true if (v,w) is guaranteed part of an optimal solution here.
 * @param index_right Reverse lookup (vertex -> current position in `right`), updated in-place.
 * @return The new list of bidomains after the split.
 */
std::vector<Bidomain> filter_domains_rr(const std::vector<Bidomain>& d, std::vector<int>& left,
        std::vector<int>& right, const Graph& g, const Graph& h, int v, int w,
        bool& best_match, std::vector<int>& index_right) {
    std::vector<Bidomain> new_d;
    new_d.reserve(d.size() * 2);
    unsigned int ccount = 0;

    for (const Bidomain& old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;

        int left_len = partition_rr(left, l, old_bd.left_len, g.adjmat[v]);

        // choose the cheaper H-side partition strategy based on w's degree
        // relative to the bidomain size - see partition_sparse_rr/
        // partition_right_rr comments
        int right_len;
        if (old_bd.right_len > (int)h.degree[w]) {
            right_len = partition_sparse_rr(right, r, old_bd.right_len, h.degree[w], h.adjlist[w], index_right);
        } else {
            right_len = partition_right_rr(right, r, old_bd.right_len, h.adjmat[w], index_right);
        }

        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;

        // this bidomain contributes to a "clean" split (or was already
        // empty) - see docstring above for what this means for best_match
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

/**
 * red_mode: bitmask controlling which RRSplit reductions are enabled
 *   bit 0 (value 1) = vertex-equivalence reduction (w_lower skip + exclude-branch removal)
 *   bit 1 (value 2) = maximality reduction (best_match commit/prune)
 *   bit 2 (value 4) = tighter equivalence-class bound
 * Presets:
 *   7 (111) = full RRSplit - all reductions on
 *   6 (110) = no vertex-equivalence (max + bound only)
 *   5 (101) = no maximality (veq + bound only)
 *   3 (011) = no tighter bound (veq + max only)
 *
 * @param g Graph G.
 * @param h Graph H.
 * @param incumbent Best solution found so far, updated in-place.
 * @param current Partial solution being built along this search path.
 * @param domains Bidomains available at this node.
 * @param left G-side vertex array.
 * @param right H-side vertex array.
 * @param goal Minimum solution size required to keep searching.
 * @param stats Search statistics, updated in-place.
 * @param start_time Time the search began, used for timing incumbent updates.
 * @param abort_due_to_timeout Flag checked to abort the search early.
 * @param EqClass Structural equivalence class of each G-side vertex.
 * @param index_right Reverse lookup (vertex -> current position in `right`), updated in-place.
 * @param red_mode Bitmask controlling which RRSplit reductions are enabled.
 */
void solve_rr(const Graph& g, const Graph& h, std::vector<VtxPair>& incumbent,
        std::vector<VtxPair>& current, std::vector<Bidomain>& domains,
        std::vector<int>& left, std::vector<int>& right, unsigned int goal,
        Stats& stats, std::chrono::time_point<std::chrono::steady_clock> start_time,
        std::atomic<bool>& abort_due_to_timeout, const ui* EqClass,
        std::vector<int>& index_right, int red_mode) {

    if (abort_due_to_timeout) { return; }

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
    if (bd_idx == -1) { return; }

    Bidomain& bd = domains[bd_idx];

    int v = find_min_value(left, bd.l, bd.left_len);

    {
        int i = 0;
        while (left[bd.l + i] != v) { i++; }
        std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
        bd.left_len--;
    }

    // Vertex-equivalence reduction: find largest w already matched to
    // an equivalent vertex, skip all w <= that value
    int w_lower = -1;
    if (red_mode & 1) {   // VEQ enabled
        for (const VtxPair& a : current) {
            if (EqClass[a.v] == EqClass[v] && w_lower < a.w) {
                w_lower = a.w;
            }
        }
    }

    // reserve v's/w's slot: bd.left_len already shrunk above when v was
    // pulled out; bd.right_len is shrunk here so the later main loop
    // (right_len down to 0, index_of_next_smallest) sees exactly the
    // remaining unmatched candidates
    bd.right_len--;

    // Tighter equivalence-class bound check (from reference):
    if ((red_mode & 4) && w_lower > 0) {   // tighter BOUND enabled
        int count_left = 0, count_right = 0;
        for (int i = bd.left_len; i >= 0; --i) {
            if (EqClass[left[bd.l + i]] == EqClass[v]) { count_left++; }
        }
        for (int i = bd.right_len; i >= 0; --i) {
            if (right[bd.r + i] > w_lower) { count_right++; }
        }
        if (bd.left_len <= bd.right_len && count_left > count_right) {
            if (bound + count_right - count_left <= (int)incumbent.size()) {
                // full undo of both v-removal and right_len decrement
                // before pruning this branch
                bd.right_len++;
                left[bd.l + bd.left_len] = v;
                bd.left_len++;
                stats.sym_pruned++;
                return;
            }
        }
        if (bd.left_len > bd.right_len && (bd.right_len - count_right) > (bd.left_len - count_left)) {
            if (bound + (bd.left_len - count_left) - (bd.right_len - count_right) <= (int)incumbent.size()) {
                bd.right_len++;
                left[bd.l + bd.left_len] = v;
                bd.left_len++;
                stats.sym_pruned++;
                return;
            }
        }
    }

    int w = w_lower;

    stats.nodes++;

#ifdef DEBUG_RR
    if (stats.nodes <= 5) {
        std::cerr << "NODE " << stats.nodes
                  << " current=" << current.size()
                  << " bound=" << bound
                  << " v=" << v
                  << " EqClass[v]=" << EqClass[v]
                  << " w_lower=" << w_lower
                  << " bd.left_len=" << bd.left_len
                  << " bd.right_len=" << (bd.right_len + 1)
                  << std::endl;
    }
#endif

    bool best_match = false;
    for (int i = bd.right_len; i >= 0; --i) {
        int idx = index_of_next_smallest(right, bd.r, bd.right_len + 1, w);
        if (idx == -1) { break; }
        w = right[bd.r + idx];

        std::swap(index_right[w], index_right[right[bd.r + bd.right_len]]);
        right[bd.r + idx] = right[bd.r + bd.right_len];
        right[bd.r + bd.right_len] = w;

        auto new_domains = filter_domains_rr(domains, left, right, g, h, v, w,
                best_match, index_right);
        current.push_back(VtxPair(v, w));
        solve_rr(g, h, incumbent, current, new_domains, left, right, goal,
                stats, start_time, abort_due_to_timeout, EqClass, index_right, red_mode);
        current.pop_back();

        // Maximality reduction: if best_match, current solution is optimal here
        bool maximality_prune = (red_mode & 2) && best_match;   // MAX enabled
        if (maximality_prune || bound <= (int)incumbent.size()) {
            stats.sym_pruned++;
            break;
        }
    }

    bd.right_len++;

    // Vertex-equivalence reduction (exclude branch):
    // Remove all G-vertices equivalent to v - their branches are isomorphic
    if (red_mode & 1) {   // VEQ enabled
        for (int i = 0; i < bd.left_len; i++) {
            if (EqClass[left[bd.l + i]] == EqClass[v]) {
                std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
                bd.left_len--;
                i--;
                stats.sym_pruned++;
            }
        }
    }

    if (bd.left_len == 0) { remove_bidomain(domains, bd_idx); }

    // branch where v is excluded entirely from the mapping
    solve_rr(g, h, incumbent, current, domains, left, right, goal,
            stats, start_time, abort_due_to_timeout, EqClass, index_right, red_mode);
}

/**
 * Entry point for RRSplit.
 * Note the density comparison direction here (sum < n*(n-1), with the
 * sort lambdas using !h_dense/!g_dense) is the opposite convention from
 * mcsplit-dal.cpp/mcsplit-dsb.cpp's density check (sum > n*(n-1)) - this
 * matches the reference's own convention for RRSplit specifically, not
 * an inconsistency; see comment below.
 *
 * @param g Graph G.
 * @param h Graph H.
 * @param stats Search statistics, updated in-place.
 * @param abort_due_to_timeout Flag checked to abort the search early.
 * @param red_mode Bitmask controlling which RRSplit reductions are enabled.
 * @return The largest common subgraph solution found, as vertex pairs in the original graphs' indices.
 */
std::vector<VtxPair> mcs_rr(const Graph& g, const Graph& h,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout,
        int red_mode) {

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

    // Reference uses sum < n*(n-1) for density (not /2)
    auto sum_vec = [](const std::vector<int>& v) {
        int s = 0; for (int x : v) { s += x; } return s;
    };

    bool h_dense = sum_vec(h_deg) < h.n * (h.n - 1);
    bool g_dense = sum_vec(g_deg) < g.n * (g.n - 1);

    std::vector<int> vv0(g.n), vv1(h.n);
    std::iota(vv0.begin(), vv0.end(), 0);
    std::iota(vv1.begin(), vv1.end(), 0);
    std::stable_sort(vv0.begin(), vv0.end(), [&](int a, int b) {
        return !h_dense ? (g_deg[a] < g_deg[b]) : (g_deg[a] > g_deg[b]);
    });
    std::stable_sort(vv1.begin(), vv1.end(), [&](int a, int b) {
        return !g_dense ? (h_deg[a] < h_deg[b]) : (h_deg[a] > h_deg[b]);
    });

    Graph g_sorted = induced_subgraph(const_cast<Graph&>(g), vv0);
    Graph h_sorted = induced_subgraph(const_cast<Graph&>(h), vv1);

    // adjacency lists + structural equivalence classes, both required by
    // RRSplit's reductions (see the direction-blindness note at the top
    // of this file - both set_adjlist and GetEqClass are boolean-only)
    set_adjlist(g_sorted);
    set_adjlist(h_sorted);

    ui* EqClass = nullptr;
    GetEqClass(g_sorted, EqClass);

    std::vector<int> index_right(h_sorted.n);
    for (int i = 0; i < h_sorted.n; i++) { index_right[i] = i; }

    // Single bidomain - reference does not split by label
    std::vector<int> left, right;
    for (int i = 0; i < g_sorted.n; i++) { left.push_back(i); }
    for (int i = 0; i < h_sorted.n; i++) { right.push_back(i); }
    std::vector<Bidomain> domains;
    domains.push_back({0, 0, g_sorted.n, h_sorted.n, false});

    std::vector<VtxPair> incumbent, current;
    auto start_time = std::chrono::steady_clock::now();
    // correctly uses g_sorted/h_sorted, matching domains/left/right/EqClass
    solve_rr(g_sorted, h_sorted, incumbent, current, domains, left, right, 1,
        stats, start_time, abort_due_to_timeout, EqClass, index_right, red_mode);

    delete[] EqClass;

    for (auto& p : incumbent) {
        p.v = vv0[p.v];
        p.w = vv1[p.w];
    }

    return incumbent;
}