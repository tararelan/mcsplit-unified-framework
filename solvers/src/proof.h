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

	/// Logs a new incumbent (in sorted indices) with soli.
	void new_incumbent(const std::vector<VtxPair> &current);

	/// Pushes the decision v->w (sorted indices) onto the proof trail.
	void push_match(int v, int w);

	/// Pushes the decision v->bottom (sorted index) onto the proof trail.
	void push_exclude(int v);

	/// Pops the most recent decision from the proof trail.
	void pop();

	/**
	 * Justifies a McSplit bound prune at the current node and logs the
	 * backtrack clause. left and right are the solver's vertex arrays that
	 * the bidomains index into (sorted indices).
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
	std::ofstream opb;
	std::ofstream pbp;

	int g_n = 0, h_n = 0;
	std::vector<std::vector<bool>> compatible; // compatible[v][w]: model has x<v>_<w>
	std::vector<long long> inj_id;			   // constraint ID of "w used at most once", -1 if absent
	std::vector<int> g_order, h_order;

	long long next_id = 1;	   // ID the next constraint-producing rule will receive
	long long obj_id = -1;	   // ID of the most recent objective-improving constraint
	long long contra_id = -1;  // ID of the root contradiction, once derived
	bool aborted = false;

	std::vector<std::string> trail;					 // decision literals, outermost first
	std::vector<std::vector<long long>> child_nogoods; // per depth: backtrack clauses awaiting deletion

	long long alo_id(int v) const { return v + 1; }
	long long amo_id(int v) const { return g_n + v + 1; }

	void write_backtrack_clause();
};
