#pragma once

// VeriPB 3 proof logging for the McSplit-family solvers (McSplit, RL, LL
// without LUM, DAL with the RL policy). The logger writes an OPB model of
// the maximum common induced subgraph instance and a proof log that
// certifies the optimal solution size found by the search.
//
// Model (written from the graphs as main() sees them, after the swap that
// makes g the smaller graph, using their original vertex indices):
//   x<v>_<w>  vertex v of G is mapped to vertex w of H (only label-compatible pairs)
//   x<v>_b    vertex v of G is left unmapped ("bottom")
//   minimise  sum_v x<v>_b
//   each v takes exactly one value (written as >= 1 and <= 1 constraints)
//   each w is used at most once
//   for every pair of assignments v1->w1, v2->w2 whose adjacency (including
//   direction and edge label) disagrees: ~x<v1>_<w1> + ~x<v2>_<w2> >= 1
//
// Proof: the search trail is a sequence of decisions, either v->w (a match)
// or v->bottom (the "exclude v" branch). Every node that is exhausted or
// pruned logs the negation of its trail as a RUP clause. A bound prune is
// justified by one pol that adds the objective-improving constraint to the
// at-least-one constraints of the left side and the injectivity constraints
// of the right side of every bidomain with |L| > |R|; unit propagation then
// sees the conflict. New incumbents are logged with soli.
//
// McSplit+DSB's tighter bound counts edges (or non-edges) inside a bidomain:
// if the smaller side S has more internal edges than the larger side B, the
// unmatched vertices of S must cover the difference in degree. The logger
// recomputes this bound itself and, when McSplit's bound is not enough,
// justifies the prune with a pbc subproof under the trail (see
// derive_degree_bound). When the H side is the smaller one, this needs the
// proof-only variables y<w> <=> "w is not the image of any vertex", which
// enable_unused_vertex_vars() introduces.
//
// LUM (McSplit+LL, McSplit+DAL) pairs leaves of v with leaves of w after
// matching v->w. Each pair is justified by redundance under the trail
// (see lum_match): the alternatives for the G leaf are ruled out one at a
// time, each by a red whose witness moves the leaf onto the H leaf, until
// unit propagation forces the pair. These clauses are guarded by the trail
// and deleted when the decision v->w is popped.

#include "graph.h"
#include "common.h"
#include <fstream>
#include <string>
#include <vector>

class ProofLogger
{
public:
	/**
	 * Opens <prefix>.opb and <prefix>.pbp for writing.
	 *
	 * @param prefix  Path prefix for the model and proof files
	 */
	explicit ProofLogger(const std::string &prefix);

	/**
	 * Writes the OPB model and the proof header. Must be called once, with the
	 * same graphs that are passed to the solver (g is the smaller graph).
	 * Also logs the empty mapping as the first solution, matching the solver's
	 * initial empty incumbent.
	 */
	void write_model(const Graph &g, const Graph &h);

	/**
	 * Tells the logger how the solver's internal (degree-sorted) vertex
	 * indices relate to the model's vertex indices: sorted vertex i of G is
	 * model vertex g_order[i], and likewise for H.
	 */
	void set_vertex_orders(const std::vector<int> &g_order, const std::vector<int> &h_order);

	/**
	 * Introduces y<w> <=> "w is unused" for every vertex w of H. Required
	 * before searching with the DSB bound; must be called at the top level.
	 */
	void enable_unused_vertex_vars();

	/// Logs a new incumbent (in sorted indices) with soli.
	void new_incumbent(const std::vector<VtxPair> &current);

	/// Pushes the decision v->w (sorted indices) onto the proof trail.
	void push_match(int v, int w);

	/**
	 * Justifies a LUM pair: under the current trail, whose last decision
	 * is v->w, the leaf vl of v is matched to the leaf wl of w (sorted
	 * indices). Both must be free, with the same edge and vertex labels.
	 */
	void lum_match(int vl, int wl);

	/// Pushes the decision v->bottom (sorted index) onto the proof trail.
	void push_exclude(int v);

	/// Pops the most recent decision, and any LUM pairs justified under it, from the proof trail.
	void pop();

	/**
	 * Justifies a bound prune at the current node and logs the backtrack
	 * clause. left and right are the solver's vertex arrays that the
	 * bidomains index into (sorted indices). Uses McSplit's bound when that
	 * suffices, and the DSB degree bounds otherwise; throws if neither
	 * justifies pruning against the current incumbent.
	 */
	void prune_by_bound(const std::vector<Bidomain> &domains,
						const std::vector<int> &left, const std::vector<int> &right);

	/// Logs the backtrack clause for a node whose children have all been refuted.
	void backtrack();

	/// Records that the search was interrupted; no further reasoning is logged.
	void abort();

	/**
	 * Writes the conclusion and footer. If the search completed, this
	 * certifies that best_size is the maximum common induced subgraph size;
	 * otherwise it only certifies that a solution of size best_size exists.
	 */
	void finish(int best_size);

private:
	// How the logger justifies the bound on one bidomain (L, R): the number
	// of matched vertices of L is at most `bound`.
	struct BidomainPlan
	{
		enum Kind
		{
			MinSize, // min(|L|, |R|), McSplit's bound
			Pair,	 // 2x2 bidomain whose internal adjacencies differ: at most 1
			Degree	 // DSB edge or non-edge counting on one side
		} kind = MinSize;
		int bound = 0;
		bool s_is_left = true;	 // Degree: is S (the counted side) the G side?
		bool complement = false; // Degree: count non-edges instead of edges
		int r = 0;				 // Degree: at least r vertices of S are unmatched
	};

	std::ofstream opb;
	std::ofstream pbp;

	int g_n = 0, h_n = 0;
	std::vector<std::vector<bool>> compatible; // compatible[v][w]: model has x<v>_<w>
	std::vector<std::vector<unsigned int>> g_adj, h_adj; // model adjacency matrices
	std::vector<int> users;					   // users[w]: number of v compatible with w
	std::vector<long long> inj_id;			   // constraint ID of "w used at most once", -1 if absent
	std::vector<long long> unused_alo_id;	   // ID of y<w> + sum_v x<v>_<w> >= 1, -1 if absent
	std::vector<long long> unused_amo_id;	   // ID of y<w> + sum_v x<v>_<w> <= 1, -1 if absent
	std::vector<int> g_order, h_order;
	int incumbent_size = 0;
	int n_matched = 0; // vertices matched by the trail and LUM

	// the current mapping, in model indices
	std::vector<int> g_image;		 // g_image[v]: w if v->w, -1 otherwise
	std::vector<bool> h_used;		 // h_used[w]: some v->w
	std::vector<bool> g_excluded;	 // g_excluded[v]: v->bottom is on the trail

	long long next_id = 1;	   // ID the next constraint-producing rule will receive
	long long obj_id = -1;	   // ID of the most recent objective-improving constraint
	long long contra_id = -1;  // ID of the root contradiction, once derived
	bool aborted = false;

	std::vector<std::string> trail;					 // decision literals, outermost first
	std::vector<bool> trail_is_match;				 // whether each decision is v->w
	std::vector<int> trail_vertex;					 // the G vertex of each decision (model index)
	std::vector<std::vector<int>> trail_lum;		 // G leaves paired by LUM under each decision
	std::vector<std::vector<long long>> trail_lum_ids; // their forcing clauses, deleted on pop
	std::vector<std::vector<long long>> child_nogoods; // per depth: backtrack clauses awaiting deletion

	long long alo_id(int v) const { return v + 1; }
	long long amo_id(int v) const { return g_n + v + 1; }

	void write_backtrack_clause();
	void record_backtrack_clause(long long id);
	void write_trail_clause();

	BidomainPlan plan_bidomain(const std::vector<int> &L, const std::vector<int> &R) const;
	long long derive_pair_bound(const std::vector<int> &L, const std::vector<int> &R);
	long long derive_degree_bound(const std::vector<int> &L, const std::vector<int> &R, const BidomainPlan &plan);
};
