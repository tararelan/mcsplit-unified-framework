#pragma once

#include <vector>
#include <string>

using ui = unsigned int;

/**
 * Represents an input graph for MCS search.
 *
 * Supports two representations simultaneously:
 *   - Adjacency matrix (adjmat): primary representation used by all
 *     algorithms for O(1) edge existence and label queries
 *   - Adjacency list (adjlist + degree): built on demand by set_adjlist()
 *     and used by RRSplit and SymSplit for efficient neighbourhood iteration
 *     during symmetry class computation
 *
 * Vertex labels are stored in label[v]. For graphs with self-loops,
 * the most significant bit of label[v] is set instead of adding a
 * self-entry to adjmat.
 *
 * leaves[u] groups the leaf neighbours of u by (edge label, vertex label)
 * pair, precomputed by pack_leaves() for use by the LUM strategy in
 * McSplit+LL and McSplit+DAL. Each entry is:
 *   { (edge_label, vertex_label) -> [leaf vertex indices] }
 *
 * degree and adjlist are owning raw pointers allocated by set_adjlist()
 * and freed by the destructor. They are nullptr until set_adjlist() is called
 * (assuming the constructor initialises them - see graph.cpp).
 */
struct Graph
{
	int n;
	std::vector<std::vector<unsigned int>> adjmat;
	std::vector<unsigned int> label;

	unsigned int *degree;
	unsigned int **adjlist;

	std::vector<std::vector<std::pair<std::pair<unsigned int, unsigned int>, std::vector<int>>>> leaves;

	Graph(unsigned int n);
	~Graph();
	Graph(const Graph&) = delete;            // owns raw pointers; copying after
	Graph& operator=(const Graph&) = delete; // set_adjlist() would double-free
	Graph(Graph&&) = default;                // moves are fine (pointers transfer, not duplicate)
	Graph& operator=(Graph&&) = default;
};

// Constructs the induced subgraph of g on the vertex subset vv, reindexed 0..vv.size()-1
Graph induced_subgraph(Graph &g, std::vector<int> vv);

// Reads a graph from file; format is 'B' (MIVIA binary), 'L' (LAD), or 'D' (DIMACS)
Graph readGraph(char *filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

// Builds adjacency lists (g.degree, g.adjlist) from g.adjmat; required before GetEqClass
void set_adjlist(Graph &g);

// Computes negative modular symmetry classes; EqClass[v] is the class label of vertex v
void GetEqClass(Graph &g, ui *&EqClass);

// Groups leaf neighbours of each vertex by (edge label, vertex label) pair;
// required before search for McSplit+LL and McSplit+DAL (LUM strategy)
void pack_leaves(Graph &g);