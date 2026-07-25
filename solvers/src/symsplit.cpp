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
static int calc_bound_sym(const std::vector<Bidomain> &domains);
static void remove_bidomain_sym(std::vector<Bidomain> &domains, int idx);
static int select_bidomain_sym(const std::vector<Bidomain> &domains, const std::vector<int> &left, int current_matching_size);
static int find_min_value_sym(const std::vector<int> &arr, int start_idx, int len);
static int index_of_next_smallest_sym(const std::vector<int> &arr, int start_idx, int len, int w);
int partition_sym(std::vector<int> &all_vv, int start, int len, const std::vector<unsigned int> &adjrow);
int partition_right_sym(std::vector<int> &all_vv, int start, int len, const std::vector<unsigned int> &adjrow, std::vector<int> &index_right);
int partition_sparse_sym(std::vector<int> &all_vv, int start, int len, int degree, const unsigned int *adjlist, std::vector<int> &index_right);
std::vector<Bidomain> filter_domains_sym(const std::vector<Bidomain> &d, std::vector<int> &left, std::vector<int> &right, const Graph &g, const Graph &h, int v, int w, bool &best_match, std::vector<int> &index_right);
int find_vertices_with_common_neighbors(const Graph &g, std::vector<int> &eqn_classes);
bool break_h_sym(const std::vector<int> &arr, int start_idx, int len, int w, const std::vector<int> &h_eqn_classes);
void solve_sym(const Graph &g, const Graph &h, std::vector<VtxPair> &incumbent, std::vector<VtxPair> &current, std::vector<Bidomain> &domains, std::vector<int> &left, std::vector<int> &right, unsigned int goal, Stats &stats, std::chrono::time_point<std::chrono::steady_clock> start_time, std::atomic<bool> &abort_due_to_timeout, const std::vector<int> &g_eqn_classes, const std::vector<int> &h_eqn_classes, std::vector<int> &index_right, int sym_mode);
std::vector<VtxPair> mcs_sym(const Graph &g, const Graph &h, Stats &stats, std::atomic<bool> &abort_due_to_timeout, int sym_mode = 3);

/**
 * Computes the McSplit upper bound on the number of additional vertex
 * pairs that can be added to the current mapping.
 *
 * Each bidomain can contribute at most min(left_len, right_len) pairs,
 * since vertices must be matched one-to-one. Summing across all
 * bidomains gives the standard McSplit bound as defined in
 * McCreesh et al. (2017).
 *
 * @param domains  Current list of bidomains
 * @return         Upper bound on remaining matchable vertex pairs
 */
static int calc_bound_sym(const std::vector<Bidomain> &domains)
{
	int bound = 0;
	for (const Bidomain &bd : domains)
	{
		bound += std::min(bd.left_len, bd.right_len);
	}
	return bound;
}

/**
 * Removes the bidomain at index idx from domains in O(1) by
 * swapping it with the last element before popping.
 * Does not preserve ordering of the remaining bidomains.
 *
 * @param domains  Bidomain list to modify
 * @param idx      Index of the bidomain to remove
 */
static void remove_bidomain_sym(std::vector<Bidomain> &domains, int idx)
{
	domains[idx] = domains[domains.size() - 1];
	domains.pop_back();
}

/**
 * Returns the minimum value in a subarray of length len starting at start_idx.
 *
 * Used during vertex selection to find the smallest available vertex
 * index within a bidomain's slice of the permutation array.
 *
 * @param arr        Permutation array (G-side or H-side)
 * @param start_idx  Start offset into arr
 * @param len        Number of elements to search
 * @return           Minimum value in arr[start_idx .. start_idx + len - 1]
 */
static int find_min_value_sym(const std::vector<int> &arr, int start_idx, int len)
{
	int min_v = INT_MAX;
	for (int i = 0; i < len; i++)
	{
		if (arr[start_idx + i] < min_v)
		{
			min_v = arr[start_idx + i];
		}
	}
	return min_v;
}

/**
 * Selects the branching bidomain using the McSplit heuristic.
 *
 * Chooses the bidomain with the smallest max(left_len, right_len),
 * breaking ties by selecting the bidomain whose G-side contains the
 * smallest vertex index. This corresponds to the SelectLabelClass
 * heuristic from McCreesh et al. (2017), which minimises the number
 * of branches created at the current search node while preferring
 * a consistent vertex ordering to reduce symmetry in practice.
 *
 * @param domains               Current list of bidomains
 * @param left                  G-side permutation array
 * @param current_matching_size Current mapping size (unused, reserved)
 * @return                      Index into domains of the selected bidomain,
 *                              or -1 if domains is empty
 */
static int select_bidomain_sym(const std::vector<Bidomain> &domains, const std::vector<int> &left, int current_matching_size)
{
	int min_size = INT_MAX;
	int min_tie_breaker = INT_MAX;
	int best = -1;
	for (unsigned int i = 0; i < domains.size(); i++)
	{
		const Bidomain &bd = domains[i];
		int len = std::max(bd.left_len, bd.right_len);
		if (len < min_size)
		{
			min_size = len;
			min_tie_breaker = find_min_value_sym(left, bd.l, bd.left_len);
			best = i;
		}
		else if (len == min_size)
		{
			int tie_breaker = find_min_value_sym(left, bd.l, bd.left_len);
			if (tie_breaker < min_tie_breaker)
			{
				min_tie_breaker = tie_breaker;
				best = i;
			}
		}
	}
	return best;
}

/**
 * Finds the index of the smallest value in a subarray that is strictly
 * greater than w.
 *
 * Used during value symmetry breaking to enumerate H-side vertices in
 * increasing order, skipping those already considered. By always picking
 * the next smallest value above the previous choice, the search visits
 * H-vertices in a fixed ascending order, which is the ordering required
 * by SymSplit's value symmetry breaking rule.
 *
 * @param arr        Permutation array (H-side)
 * @param start_idx  Start offset into arr
 * @param len        Number of elements to search
 * @param w          Lower bound - returned index must have arr value > w
 * @return           Offset from start_idx of the next smallest element,
 *                   or -1 if no such element exists
 */
static int index_of_next_smallest_sym(const std::vector<int> &arr, int start_idx, int len, int w)
{
	int idx = -1;
	int smallest = INT_MAX;
	for (int i = 0; i < len; i++)
	{
		if (arr[start_idx + i] > w && arr[start_idx + i] < smallest)
		{
			smallest = arr[start_idx + i];
			idx = i;
		}
	}
	return idx;
}

/**
 * Partitions a subarray of the permutation array in place by adjacency.
 *
 * Reorders all_vv[start .. start+len-1] so that vertices adjacent to
 * the most recently matched vertex (adjrow[v] != 0) come first, followed
 * by non-adjacent vertices. This is the bidomain splitting step performed
 * after each match: the adjacent partition becomes the new adjacent
 * bidomain and the non-adjacent partition becomes the new non-adjacent
 * bidomain for the next search level.
 *
 * Uses a two-pointer sweep (stable relative to the adjacent group,
 * unstable for the non-adjacent group) and runs in O(len).
 *
 * @param all_vv   Shared permutation array (G-side or H-side)
 * @param start    Start offset of the subarray to partition
 * @param len      Length of the subarray
 * @param adjrow   Adjacency row of the matched vertex; adjrow[v] != 0 means v is adjacent
 * @return         Number of adjacent vertices (size of the left partition)
 */
int partition_sym(std::vector<int> &all_vv, int start, int len, const std::vector<unsigned int> &adjrow)
{
	int i = 0;
	for (int j = 0; j < len; j++)
	{
		if (adjrow[all_vv[start + j]])
		{
			std::swap(all_vv[start + i], all_vv[start + j]);
			i++;
		}
	}
	return i;
}

/**
 * Partitions the H-side subarray by adjacency, maintaining a position
 * index for O(1) vertex lookups.
 *
 * Equivalent to partition_sym but additionally keeps index_right
 * consistent with all_vv after each swap. index_right[v] always holds
 * the current position of vertex v in all_vv, which SymSplit's value
 * symmetry breaking rule uses to check in O(1) whether a vertex
 * equivalent to the current candidate has already been placed in the
 * adjacent partition (and therefore already considered).
 *
 * @param all_vv       H-side permutation array
 * @param start        Start offset of the subarray to partition
 * @param len          Length of the subarray
 * @param adjrow       Adjacency row of the matched H-vertex; adjrow[v] != 0 means v is adjacent
 * @param index_right  Position index kept in sync with all_vv:
 *                     index_right[v] = current index of v in all_vv
 * @return             Number of adjacent vertices (size of the left partition)
 */
int partition_right_sym(std::vector<int> &all_vv, int start, int len,
                        const std::vector<unsigned int> &adjrow, std::vector<int> &index_right)
{
	int i = 0;
	for (int j = 0; j < len; j++)
	{
		if (adjrow[all_vv[start + j]])
		{
			std::swap(index_right[all_vv[start + i]], index_right[all_vv[start + j]]);
			std::swap(all_vv[start + i], all_vv[start + j]);
			i++;
		}
	}
	return i;
}

/**
 * Partitions the H-side subarray by adjacency using the adjacency list
 * rather than a full adjacency row scan.
 *
 * Functionally equivalent to partition_right_sym but iterates over the
 * neighbour list of the matched vertex instead of scanning all len
 * candidates. For sparse graphs this reduces the partition cost from
 * O(len) to O(degree), which is significant when bidomains are large
 * but the matched vertex has few neighbours.
 *
 * Uses index_right to locate each neighbour in all_vv in O(1), swapping
 * it into the adjacent partition if it falls within the current subarray
 * bounds [start, start+len). index_right is kept consistent with all_vv
 * after every swap.
 *
 * @param all_vv       H-side permutation array
 * @param start        Start offset of the subarray to partition
 * @param len          Length of the subarray
 * @param degree       Number of neighbours of the matched H-vertex
 * @param adjlist      Neighbour list of the matched H-vertex
 * @param index_right  Position index kept in sync with all_vv:
 *                     index_right[v] = current index of v in all_vv
 * @return             Number of adjacent vertices moved to the left partition
 */
int partition_sparse_sym(std::vector<int> &all_vv, int start, int len,
                         int degree, const unsigned int *adjlist, std::vector<int> &index_right)
{
	int j = 0;
	for (int i = 0; i < degree; i++)
	{
		int pos = index_right[adjlist[i]];
		if (pos >= start && pos < start + len)
		{
			std::swap(index_right[all_vv[start + j]], index_right[all_vv[pos]]);
			std::swap(all_vv[start + j], all_vv[pos]);
			j++;
		}
	}
	return j;
}

/**
 * Splits all bidomains after committing to the match (v, w).
 *
 * For each existing bidomain, partitions both the G-side and H-side
 * vertex sets by adjacency to v and w respectively, producing up to
 * two new bidomains per old one:
 *   - Adjacent partition: vertices adjacent to both v and w (is_adjacent = true)
 *   - Non-adjacent partition: vertices non-adjacent to both v and w
 *
 * The H-side partition uses the sparse adjacency list when
 * right_len > h.degree[w] (i.e. fewer neighbours than candidates),
 * falling back to the dense adjmat scan otherwise. index_right is
 * kept consistent with the right permutation array throughout.
 *
 * Bidomains where either partition is empty on both sides are counted
 * toward ccount. If every old bidomain collapses this way (ccount ==
 * d.size()), best_match is set to true, indicating that the current
 * mapping cannot be extended further - used by the maximality reduction
 * in SymSplit to detect and record maximal solutions early.
 *
 * @param d            Current bidomain list
 * @param left         G-side permutation array (modified in place)
 * @param right        H-side permutation array (modified in place)
 * @param g            Pattern graph
 * @param h            Target graph
 * @param v            Matched vertex in G
 * @param w            Matched vertex in H
 * @param best_match   Output: set to true if no bidomain can be extended
 * @param index_right  H-side position index, kept in sync with right
 * @return             New bidomain list after splitting
 */
std::vector<Bidomain> filter_domains_sym(const std::vector<Bidomain> &d, std::vector<int> &left,
                                         std::vector<int> &right, const Graph &g, const Graph &h,
                                         int v, int w, bool &best_match, std::vector<int> &index_right)
{
	std::vector<Bidomain> new_d;
	new_d.reserve(d.size() * 2);
	unsigned int ccount = 0;

	for (const Bidomain &old_bd : d)
	{
		int l = old_bd.l;
		int r = old_bd.r;

		int left_len = partition_sym(left, l, old_bd.left_len, g.adjmat[v]);

		int right_len;
		if (old_bd.right_len > (int)h.degree[w])
		{
			right_len = partition_sparse_sym(right, r, old_bd.right_len, h.degree[w], h.adjlist[w], index_right);
		}
		else
		{
			right_len = partition_right_sym(right, r, old_bd.right_len, h.adjmat[w], index_right);
		}

		int left_len_noedge = old_bd.left_len - left_len;
		int right_len_noedge = old_bd.right_len - right_len;

		if ((left_len == 0 && right_len == 0) ||
			(left_len_noedge == 0 && right_len_noedge == 0) ||
			old_bd.left_len == 0)
		{
			ccount++;
		}

		if (left_len_noedge && right_len_noedge)
		{
			new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
		}
		if (left_len && right_len)
		{
			new_d.push_back({l, r, left_len, right_len, true});
		}
	}

	best_match = (ccount == d.size());
	return new_d;
}

/**
 * Combines a hash value into an existing seed using the same mixing
 * function as boost::hash_combine.
 *
 *
 * Used by GetEqClass to hash vertex neighbourhoods for O(n) symmetry
 * class detection instead of O(n^2) pairwise comparison.
 *
 * @param seed  Hash accumulator, modified in place
 * @param v     Value to mix into the seed
 */
static void hash_combine(std::size_t &seed, std::size_t v)
{
	seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
 * Computes modular symmetry equivalence classes using neighbourhood hashing.
 *
 *
 * For each vertex v, two hashes are computed:
 *   - n_hash: hash of v's neighbours only - detects negative symmetry
 *     (N^-(u) = N^-(v), vertices share the same open neighbourhood)
 *   - p_hash: hash of v's neighbours plus v itself - detects positive
 *     symmetry (N^+(u) = N^+(v), vertices share the same closed neighbourhood)
 *
 * Both hashes are inserted into the same map, so vertices are grouped
 * together if they are symmetric under either criterion. Groups of size
 * greater than one are assigned a shared class label; singleton groups
 * (no symmetry partner) are left with eqn_classes[v] = -1.
 *
 * Time complexity is O(n^2) due to the full adjacency row scan per vertex,
 * consistent with the O(n^2) symmetry detection bound stated in the
 * SymSplit paper.
 *
 * @param g            Input graph with adjmat initialised
 * @param eqn_classes  Output vector of length g.n; eqn_classes[v] is the
 *                     symmetry class label of v, or -1 if v has no partner
 * @return             Total number of vertices belonging to a non-trivial
 *                     symmetry class (i.e. vertices where eqn_classes[v] != -1)
 */
int find_vertices_with_common_neighbors(const Graph &g, std::vector<int> &eqn_classes)
{
	std::unordered_map<std::size_t, std::vector<int>> n_groups;
	n_groups.reserve(g.n);
	eqn_classes.assign(g.n, -1);

	for (int v = 0; v < g.n; v++)
	{
		std::size_t n_hash = 0, p_hash = 0;
		for (int u = 0; u < g.n; u++)
		{
			if (v != u && g.adjmat[v][u])
			{
				boost::hash_combine(n_hash, (unsigned int)u);
				boost::hash_combine(p_hash, (unsigned int)u);
			}
			if (v == u)
			{
				boost::hash_combine(p_hash, (unsigned int)v);
			}
		}
		n_groups[p_hash].push_back(v);
		n_groups[n_hash].push_back(v);
	}

	int n_syms = 0;
	int label = 0;
	for (const auto &batch : n_groups)
	{
		if (batch.second.size() > 1)
		{
			for (int v : batch.second)
			{
				eqn_classes[v] = label;
				n_syms++;
			}
		}
		label++;
	}
	return n_syms;
}

/**
 * Value symmetry breaking check for SymSplit.
 *
 * Returns true if the branch matching v to w should be pruned because
 * a smaller H-vertex equivalent to w is still present in the current
 * bidomain. By Rule 2 (Value Symmetry Breaking) from Kothalawala et al.
 * (2026), when multiple symmetric H-vertices are available as candidates,
 * only the smallest under the fixed ordering need be considered - all
 * others produce isomorphic solutions and can be skipped.
 *
 * @param arr            H-side permutation array
 * @param start_idx      Start offset of the current bidomain in arr
 * @param len            Number of H-side candidates in the bidomain
 * @param w              Candidate H-vertex being considered
 * @param h_eqn_classes  Symmetry class labels for H; h_eqn_classes[v] = -1
 *                       if v has no symmetry partner
 * @return               True if a smaller equivalent H-vertex exists in
 *                       the bidomain, meaning w should be skipped
 */
bool break_h_sym(const std::vector<int> &arr, int start_idx, int len, int w,
                 const std::vector<int> &h_eqn_classes)
{
	int w_class = h_eqn_classes[w];
	for (int i = 0; i < len; i++)
	{
		if (h_eqn_classes[arr[start_idx + i]] == w_class && arr[start_idx + i] < w)
		{
			return true;
		}
	}
	return false;
}

/**
 * Core branch-and-bound recursive search for SymSplit.
 *
 * Extends RRSplit with a dual symmetry breaking framework applied to
 * both input graphs simultaneously. Three pruning mechanisms are active,
 * controlled by sym_mode (bitmask: bit 0 = variable symmetry, bit 1 = value symmetry):
 *
 * 1. Variable symmetry breaking (G-side, bit 0):
 *    When branching on v, computes w_lower = the largest w' already matched
 *    to a G-vertex equivalent to v. Any H-vertex w <= w_lower would produce
 *    a cs-isomorphic solution already explored, so the include branch only
 *    considers w > w_lower. In the exclude branch, all G-vertices equivalent
 *    to v are removed from the left domain, as their exclude branches are
 *    isomorphic to v's.
 *
 * 2. Value symmetry breaking (H-side, bit 1):
 *    For each candidate w, skips it if a smaller H-vertex equivalent to w
 *    is still present in the bidomain (break_h_sym). Only the smallest
 *    representative of each H-symmetry class needs to be explored.
 *
 * 3. Maximality reduction:
 *    If filter_domains_sym sets best_match (no bidomain can be extended),
 *    the current mapping is maximal and the exclude branch is skipped.
 *    This is also triggered if the incumbent has grown to match the bound,
 *    making further search futile.
 *
 * Node counting follows the SymSplit reference: stats.nodes is incremented
 * at the top of each call before the bound check, so pruned nodes are counted.
 *
 * @param g                    Pattern graph
 * @param h                    Target graph
 * @param incumbent            Best solution found so far (updated in place)
 * @param current              Current partial mapping being extended
 * @param domains              Current bidomain list
 * @param left                 G-side permutation array
 * @param right                H-side permutation array
 * @param goal                 Target solution size (for top-down search variant)
 * @param stats                Search statistics (updated in place)
 * @param start_time           Search start time for elapsed time computation
 * @param abort_due_to_timeout Atomic flag set by timeout thread; checked at entry
 * @param g_eqn_classes        Symmetry class labels for G vertices (-1 if none)
 * @param h_eqn_classes        Symmetry class labels for H vertices (-1 if none)
 * @param index_right          H-side position index kept in sync with right
 * @param sym_mode             Bitmask controlling active symmetry rules
 */
void solve_sym(const Graph &g, const Graph &h, std::vector<VtxPair> &incumbent,
               std::vector<VtxPair> &current, std::vector<Bidomain> &domains,
               std::vector<int> &left, std::vector<int> &right, unsigned int goal,
               Stats &stats, std::chrono::time_point<std::chrono::steady_clock> start_time,
               std::atomic<bool> &abort_due_to_timeout,
               const std::vector<int> &g_eqn_classes, const std::vector<int> &h_eqn_classes,
               std::vector<int> &index_right, int sym_mode)
{
	if (abort_due_to_timeout)
	{
		return;
	}

	if (current.size() > incumbent.size())
	{
		incumbent = current;
		stats.incumbent_size = incumbent.size();
		stats.nodes_to_best = stats.nodes;
		stats.time_to_best = std::chrono::duration<double>(
								 std::chrono::steady_clock::now() - start_time)
								 .count();
	}

	stats.nodes++;

	int bound = current.size() + calc_bound_sym(domains);
	if (bound <= (int)incumbent.size() || bound < (int)goal)
	{
		stats.cut_branches++;
		stats.bound_pruned++;
		return;
	}

	int bd_idx = select_bidomain_sym(domains, left, current.size());
	if (bd_idx == -1)
	{
		return;
	}

	Bidomain &bd = domains[bd_idx];

	int v = find_min_value_sym(left, bd.l, bd.left_len);
	{
		int i = 0;
		while (left[bd.l + i] != v)
		{
			i++;
		}
		std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
		bd.left_len--;
	}

	int w_lower = -1;
	if ((sym_mode & 1) && g_eqn_classes[v] != -1)
	{
		for (const VtxPair &a : current)
		{
			if (g_eqn_classes[a.v] == g_eqn_classes[v] && w_lower < a.w)
			{
				w_lower = a.w;
			}
		}
	}

	bd.right_len--;
	int w = w_lower;

	bool best_match = false;
	bool skip_exclude = false; // ← NEW

	for (int i = bd.right_len; i >= 0; --i)
	{
		int idx = index_of_next_smallest_sym(right, bd.r, bd.right_len + 1, w);
		if (idx == -1)
		{
			break;
		}
		w = right[bd.r + idx];

		if ((sym_mode & 2) && h_eqn_classes[w] != -1 && break_h_sym(right, bd.r, bd.right_len + 1, w, h_eqn_classes))
		{
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
				  g_eqn_classes, h_eqn_classes, index_right, sym_mode);
		current.pop_back();

		if (best_match || bound <= (int)incumbent.size())
		{
			stats.sym_pruned++;
			skip_exclude = true; // ← NEW: skip exclude branch, matching reference
			break;
		}
	}

	bd.right_len++;

	if (g_eqn_classes[v] != -1)
	{
		if ((sym_mode & 1) && g_eqn_classes[v] != -1)
		{
			for (int i = 0; i < bd.left_len; i++)
			{
				if (g_eqn_classes[left[bd.l + i]] == g_eqn_classes[v])
				{
					std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
					bd.left_len--;
					i--;
					stats.sym_pruned++;
				}
			}
		}
	}

	if (bd.left_len == 0)
	{
		remove_bidomain_sym(domains, bd_idx);
	}

	if (!skip_exclude)
	{ // ← NEW: only explore exclude branch if not pruned
		solve_sym(g, h, incumbent, current, domains, left, right, goal,
				  stats, start_time, abort_due_to_timeout,
				  g_eqn_classes, h_eqn_classes, index_right, sym_mode);
	}
}

/**
 * Entry point for SymSplit: sets up the search state and runs solve_sym().
 *
 * 1. Degree sorting: vertices are sorted by degree before search to improve
 *    branching quality. Sort direction depends on graph density - ascending
 *    for sparse graphs (fewer branches early), descending for dense graphs.
 *    Density is approximated by comparing total edge count to n*(n-1).
 *
 * 2. Symmetry class computation: find_vertices_with_common_neighbors is called
 *    on the original unsorted graphs, matching the reference implementation
 *    which computes symmetry before reindexing. Class labels are then remapped
 *    through the sort permutation so they align with the sorted vertex indices
 *    used during search.
 *
 * 3. Graph reindexing: induced_subgraph produces sorted copies of g and h.
 *    Adjacency lists are built on the sorted copies for use by the sparse
 *    partition function in filter_domains_sym.
 *
 * 4. Initial state: a single bidomain covering all vertices of both graphs
 *    is created, index_right is initialised to the identity permutation,
 *    and the root upper bound is recorded before branching begins.
 *
 * 5. Result remapping: vertex indices in the returned incumbent are remapped
 *    back to the original unsorted numbering via the sort permutations.
 *
 * @param g                    Pattern graph
 * @param h                    Target graph
 * @param stats                Search statistics (updated in place)
 * @param abort_due_to_timeout Atomic flag set by timeout thread
 * @param sym_mode             Bitmask: bit 0 = variable symmetry, bit 1 = value symmetry
 * @return                     Maximum common induced subgraph as matched vertex pairs,
 *                             indexed into the original unsorted graphs
 */
std::vector<VtxPair> mcs_sym(const Graph &g, const Graph &h,
                              Stats &stats, std::atomic<bool> &abort_due_to_timeout, int sym_mode)
{

	auto calc_degrees = [](const Graph &g)
	{
		std::vector<int> degree(g.n, 0);
		for (int v = 0; v < g.n; v++)
		{
			for (int w = 0; w < g.n; w++)
			{
				if (g.adjmat[v][w] & 1)
				{
					degree[v]++;
				}
			}
		}
		return degree;
	};

	std::vector<int> g_deg = calc_degrees(g);
	std::vector<int> h_deg = calc_degrees(h);

	auto sum_vec = [](const std::vector<int> &v)
	{
		int s = 0;
		for (int x : v)
		{
			s += x;
		}
		return s;
	};

	bool h_dense = sum_vec(h_deg) < h.n * (h.n - 1);
	bool g_dense = sum_vec(g_deg) < g.n * (g.n - 1);

	std::vector<int> vv0(g.n), vv1(h.n);
	std::iota(vv0.begin(), vv0.end(), 0);
	std::iota(vv1.begin(), vv1.end(), 0);
	std::stable_sort(vv0.begin(), vv0.end(), [&](int a, int b)
					 { return !h_dense ? (g_deg[a] < g_deg[b]) : (g_deg[a] > g_deg[b]); });
	std::stable_sort(vv1.begin(), vv1.end(), [&](int a, int b)
					 { return !g_dense ? (h_deg[a] < h_deg[b]) : (h_deg[a] > h_deg[b]); });

	// Compute symmetry classes on unsorted graphs - matches reference which
	// calls find_vertices_with_common_neighbors before induced_subgraph
	std::vector<int> g_eqn_classes, h_eqn_classes;
	find_vertices_with_common_neighbors(g, g_eqn_classes);
	find_vertices_with_common_neighbors(h, h_eqn_classes);

	// Remap equivalence classes through the sort permutation so they align
	// with the sorted vertex indices used during search
	std::vector<int> g_eqn_sorted(g.n), h_eqn_sorted(h.n);
	for (int i = 0; i < g.n; i++)
	{
		g_eqn_sorted[i] = g_eqn_classes[vv0[i]];
	}
	for (int i = 0; i < h.n; i++)
	{
		h_eqn_sorted[i] = h_eqn_classes[vv1[i]];
	}

	Graph g_sorted = induced_subgraph(const_cast<Graph &>(g), vv0);
	Graph h_sorted = induced_subgraph(const_cast<Graph &>(h), vv1);

	set_adjlist(g_sorted);
	set_adjlist(h_sorted);

	std::vector<int> index_right(h_sorted.n);
	for (int i = 0; i < h_sorted.n; i++)
	{
		index_right[i] = i;
	}

	std::vector<int> left, right;
	for (int i = 0; i < g_sorted.n; i++)
	{
		left.push_back(i);
	}
	for (int i = 0; i < h_sorted.n; i++)
	{
		right.push_back(i);
	}
	std::vector<Bidomain> domains;
	domains.push_back({0, 0, g_sorted.n, h_sorted.n, false});

	std::vector<VtxPair> incumbent, current;
	auto start_time = std::chrono::steady_clock::now();
	solve_sym(g_sorted, h_sorted, incumbent, current, domains, left, right, 1,
			  stats, start_time, abort_due_to_timeout,
			  g_eqn_sorted, h_eqn_sorted, index_right, sym_mode);

	for (auto &p : incumbent)
	{
		p.v = vv0[p.v];
		p.w = vv1[p.w];
	}

	return incumbent;
}