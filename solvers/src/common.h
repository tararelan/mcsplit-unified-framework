#pragma once

#include <vector>
#include <chrono>

/**
 * A matched pair of vertices, one from each input graph.
 * Represents a single assignment in the current partial mapping M,
 * where v in G is mapped to w in H.
 */
struct VtxPair
{
	int v; // vertex in G
	int w; // vertex in H
	VtxPair(int v, int w) : v(v), w(w) {}
};

/**
 * A bidomain (label class) in the McSplit search state.
 *
 * Represents a pair of vertex subsets - one from G and one from H -
 * whose members share the same adjacency pattern relative to all
 * currently matched pairs. Any vertex in the left set may be matched
 * with any vertex in the right set without violating the isomorphism
 * constraint. The upper bound contribution of this bidomain is
 * min(left_len, right_len).
 *
 * Vertices are not stored directly; instead, l and r are offsets into
 * a shared permutation array (one for G, one for H), with left_len and
 * right_len giving the number of valid entries from that offset.
 *
 * is_adjacent distinguishes bidomains whose vertices are adjacent to
 * the most recently matched pair (label bit = 1) from those that are
 * not (label bit = 0), which is used by the MCCS connected-subgraph
 * variant to restrict branching to adjacent vertices only.
 *
 * size and is_valid are used by McSplit+DSB to cache the degree-sequence
 * bound computation for this bidomain, avoiding redundant recalculation.
 */
struct Bidomain
{
	int l;				   // offset into the G permutation array
	int r;				   // offset into the H permutation array
	int left_len;		   // number of G vertices in this bidomain
	int right_len;		   // number of H vertices in this bidomain
	bool is_adjacent;	   // true if vertices are adjacent to the last matched pair
	int size = -1;		   // cached DSB size (-1 if not yet computed)
	bool is_valid = false; // true if the cached size is current
	Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent) : l(l), r(r), left_len(left_len), right_len(right_len), is_adjacent(is_adjacent) {}
	Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent, int size, bool is_valid) : l(l), r(r), left_len(left_len), right_len(right_len), is_adjacent(is_adjacent), size(size), is_valid(is_valid) {}
};

/**
 * Search statistics collected during a single algorithm run.
 *
 * Populated incrementally during search and written to the output CSV
 * after termination. All counts refer to the entire search, including
 * both the finding phase (reaching the incumbent) and the proving phase
 * (exhausting remaining branches to certify optimality).
 *
 * nodes_to_best and time_to_best isolate the finding phase, which is
 * useful for comparing how quickly algorithms reach the optimal solution
 * independently of how long they take to prove it.
 *
 * bound_pruned and sym_pruned are subsets of cut_branches, broken out
 * to allow independent analysis of bounding versus symmetry reduction
 * contributions to pruning.
 */
struct Stats
{
	unsigned long long nodes = 0;		  // total search tree nodes explored
	unsigned long long cut_branches = 0;  // total branches pruned
	unsigned long long bound_pruned = 0;  // branches pruned by the upper bound
	unsigned long long sym_pruned = 0;	  // branches pruned by symmetry or equivalence reduction
	double time_elapsed = 0.0;			  // total wall-clock time in seconds
	double time_to_best = 0.0;			  // wall-clock time when incumbent was last updated
	unsigned long long nodes_to_best = 0; // node count when incumbent was last updated
	int incumbent_size = 0;				  // number of vertices in the best solution found
	int solution_edges = 0;				  // number of edges in the best solution found
	bool aborted = false;				  // true if the timeout was reached before completion
};