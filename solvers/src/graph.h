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
 * and freed by the destructor. They are nullptr until set_adjlist() is called.
 */
struct Graph
{
	int n;										   // number of vertices
	std::vector<std::vector<unsigned int>> adjmat; // n x n adjacency matrix; adjmat[u][v] = edge label (0 = no edge)
	std::vector<unsigned int> label;			   // vertex labels; label[v] = 0 for unlabelled graphs

	unsigned int *degree;	// degree[v] = number of neighbours of v
	unsigned int **adjlist; // adjlist[v] = array of neighbour indices, length degree[v]

	std::vector<std::vector<std::pair<std::pair<unsigned int, unsigned int>, std::vector<int>>>> leaves;
	// leaves[u] = list of { (edge_label, vertex_label), [leaf neighbour indices] }

	Graph(unsigned int n);
	~Graph();
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