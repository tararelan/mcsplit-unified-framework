// #define _GNU_SOURCE

#include "graph.h"
#include <stdexcept>
#include <climits>
#include <algorithm>
#include <string>

// Constructor
Graph::Graph(unsigned int n) : n(n), adjmat(n, std::vector<unsigned int>(n, 0)), label(n, 0), degree(nullptr), adjlist(nullptr), leaves(n) {}

// Destructor
// Needs to free degree and adjlist if they were allocated
Graph::~Graph()
{
	if (adjlist != nullptr)
	{
		for (int i = 0; i < n; i++)
		{
			delete[] adjlist[i];
		}
		delete[] adjlist;
	}
	if (degree != nullptr)
	{
		delete[] degree;
	}
}

/**
 * Adds an edge between vertices v and w in the adjacency matrix.
 *
 * For self-loops (v == w), sets the most significant bit of the vertex label
 * instead of modifying the adjacency matrix.
 *
 * For directed graphs, encodes edge direction asymmetrically:
 *   - Forward edge (v -> w): stored in adjmat[v][w]
 *   - Reverse edge (w -> v): stored in adjmat[w][v], shifted left by 16 bits
 *     to distinguish direction in the same cell
 *
 * For undirected graphs, both directions are set identically.
 *
 * Uses bitwise OR (|=) rather than assignment for directed edges to allow
 * multiple edge labels to coexist in the same cell (e.g. both directions).
 *
 * @param g             Graph to modify
 * @param v             Source vertex
 * @param w             Target vertex
 * @param edge_val      Edge label value (used only if edge_labelled is true)
 * @param directed      If true, encode edge direction asymmetrically
 * @param edge_labelled If true, use edge_val; otherwise use fixed values (1/2)
 */
static void add_edge(Graph &g, int v, int w, unsigned int edge_val, bool directed, bool edge_labelled)
{
	if (v != w)
	{
		if (edge_labelled)
		{
			if (directed)
			{
				g.adjmat[v][w] |= edge_val;
				g.adjmat[w][v] |= (edge_val << 16);
			}
			else
			{
				g.adjmat[v][w] = edge_val;
				g.adjmat[w][v] = edge_val;
			}
		}
		else
		{
			if (directed)
			{
				g.adjmat[v][w] |= 1;
				g.adjmat[w][v] |= 2;
			}
			else
			{
				g.adjmat[v][w] = 1;
				g.adjmat[w][v] = 1;
			}
		}
	}
	else
	{
		// loop: set most significant bit of vertex label
		g.label[v] |= (1u << (sizeof(unsigned int) * CHAR_BIT - 1));
	}
}

/**
 * Constructs the induced subgraph of g on the vertex subset vv.
 *
 * The induced subgraph contains exactly the vertices in vv and all edges
 * between them that exist in g. Vertices are reindexed 0..vv.size()-1
 * in the order they appear in vv, preserving both adjacency structure
 * and vertex labels.
 *
 * @param g   Source graph
 * @param vv  Vertex subset to induce on (indices into g)
 * @return    New Graph containing the induced subgraph
 */
Graph induced_subgraph(Graph &g, std::vector<int> vv)
{
	Graph subg(vv.size());
	for (int i = 0; i < subg.n; i++)
	{
		for (int j = 0; j < subg.n; j++)
		{
			subg.adjmat[i][j] = g.adjmat[vv[i]][vv[j]];
		}
	}

	for (int i = 0; i < subg.n; i++)
	{
		subg.label[i] = g.label[vv[i]];
	}

	return subg;
}

/**
 * Reads a 2-byte little-endian word from a binary file stream.
 *
 * @param fp  File pointer to read from
 * @return    Integer value of the 2-byte word
 * @throws    std::runtime_error if fewer than 2 bytes could be read
 */
static int read_word(FILE *fp)
{
	unsigned char a[2];
	if (fread(a, 1, 2, fp) != 2)
	{
		throw std::runtime_error("Error reading file.");
	}
	return (int)a[0] | (((int)a[1]) << 8);
}

/**
 * Reads a graph from a MIVIA/VF binary format file.
 *
 * File format:
 *   - 2 bytes: number of vertices n
 *   - n x 2 bytes: vertex labels (one per vertex)
 *   - For each vertex i:
 *       - 2 bytes: number of outgoing edges
 *       - For each edge: 2 bytes target vertex, 2 bytes edge label
 *
 * Vertex and edge labels are extracted by shifting the raw 16-bit value
 * right by (16 - k1) bits, where k1 is derived from n to produce
 * approximately 33% label diversity relative to graph size. This matches
 * the labelling scheme used in the McSplit paper [McCreesh et al. 2017]
 * where the number of distinct labels is approximately 33% of n.
 * If labelled is false, all labels are ignored and edges are treated
 * as unlabelled (label value fixed to 1).
 *
 * @param filename  Path to the binary graph file
 * @param directed  If true, encode edge directions asymmetrically
 * @param labelled  If true, read and apply vertex and edge labels
 * @return          Parsed Graph object
 * @throws          std::runtime_error if the file cannot be opened or read
 */
static Graph readBinaryGraph(char *filename, bool directed, bool labelled)
{
	FILE *f;
	if ((f = fopen(filename, "rb")) == NULL)
	{
		throw std::runtime_error("Cannot open file.");
	}

	int n = read_word(f);
	Graph g(n);

	int m = n * 33 / 100;
	int p = 1, k1 = 0, k2 = 0;
	while (p < m && k1 < 16)
	{
		p *= 2;
		k1 = k2;
		k2++;
	}

	for (int i = 0; i < n; i++)
	{
		int label = (read_word(f) >> (16 - k1));
		if (labelled)
		{
			g.label[i] |= label;
		}
	}

	for (int i = 0; i < n; i++)
	{
		int len = read_word(f);
		for (int j = 0; j < len; j++)
		{
			int target = read_word(f);
			int label = (read_word(f) >> (16 - k1)) + 1;
			add_edge(g, i, target, labelled ? label : 1, directed, labelled);
		}
	}

	fclose(f);
	return g;
}

/**
 * Reads an unlabelled graph from a LAD format text file.
 *
 * File format:
 *   - Line 1: number of vertices n
 *   - n subsequent lines, one per vertex i:
 *       - First value: number of outgoing edges
 *       - Remaining values: target vertex indices
 *
 * LAD format does not support vertex or edge labels; all edges are
 * treated as unlabelled (edge_val = 1).
 *
 * @param filename  Path to the LAD format text file
 * @param directed  If true, encode edge directions asymmetrically
 * @return          Parsed Graph object
 * @throws          std::runtime_error if the file cannot be opened or
 *                  any value cannot be read
 */
static Graph readLadGraph(char *filename, bool directed)
{
	FILE *f;
	if ((f = fopen(filename, "r")) == NULL)
	{
		throw std::runtime_error("Cannot open file.");
	}

	int n = 0;
	if (fscanf(f, "%d", &n) != 1)
	{
		throw std::runtime_error("Number of vertices not read correctly.");
	}

	Graph g(n);

	for (int i = 0; i < n; i++)
	{
		int edge_count;
		if (fscanf(f, "%d", &edge_count) != 1)
		{
			throw std::runtime_error("Number of edges not read correctly.");
		}
		for (int j = 0; j < edge_count; j++)
		{
			int w;
			if (fscanf(f, "%d", &w) != 1)
			{
				throw std::runtime_error("an edge was not read correctly.");
			}
			add_edge(g, i, w, 1, directed, false);
		}
	}
	fclose(f);
	return g;
}

/**
 * Reads a graph from a DIMACS format text file.
 *
 * Supported line types:
 *   - 'p edge <n> <m>': problem line declaring n vertices and m edges;
 *     must appear before any edge or label lines
 *   - 'e <v> <w>':      edge between vertices v and w (1-indexed);
 *     converted to 0-indexed internally
 *   - 'n <v> <label>':  vertex label assignment (1-indexed);
 *     applied only if labelled is true
 *   - All other lines (comments etc.) are silently ignored
 *
 * @param filename  Path to the DIMACS format text file
 * @param directed  If true, encode edge directions asymmetrically
 * @param labelled  If true, read and apply vertex labels from 'n' lines
 * @return          Parsed Graph object
 * @throws          std::runtime_error if the file cannot be opened,
 *                  the problem line is malformed or missing, any edge
 *                  or label line is malformed, or the file is empty
 */
static Graph readDimacsGraph(char *filename, bool directed, bool labelled)
{
	FILE *f;
	if ((f = fopen(filename, "r")) == NULL)
	{
		throw std::runtime_error("Cannot open file.");
	}

	char *line = NULL;
	size_t nchar = 0;
	int n = 0, m = 0, edges_read = 0;
	Graph *g = nullptr;

	while (getline(&line, &nchar, f) != -1)
	{
		if (nchar > 0)
		{
			switch (line[0])
			{
			case 'p':
			{
				if (sscanf(line, "p edge %d %d", &n, &m) != 2)
				{
					throw std::runtime_error("Error reading p line.");
				}
				g = new Graph(n);
				break;
			}
			case 'e':
			{
				if (g == nullptr)
				{
					throw std::runtime_error("Graph size must be specified first.");
				}
				int v, w;
				if (sscanf(line, "e %d %d", &v, &w) != 2)
				{
					throw std::runtime_error("Error reading e line.");
				}
				add_edge(*g, v - 1, w - 1, 1, directed, labelled);
				edges_read++;
				break;
			}
			case 'n':
			{
				if (g == nullptr)
				{
					throw std::runtime_error("Graph size must be specified first.");
				}
				int v, label;
				if (sscanf(line, "n %d %d", &v, &label) != 2)
				{
					throw std::runtime_error("Error reading n line.");
				}
				if (labelled)
				{
					g->label[v - 1] |= label;
				}
				break;
			}
			}
		}
	}
	free(line);
	fclose(f);

	if (g == nullptr)
	{
		throw std::runtime_error("Empty graph file.");
	}

	Graph result = *g;
	delete g;
	return result;
}

/**
 * Reads a graph from file in one of the supported formats.
 *
 * Dispatches to the appropriate reader based on the format character:
 *   - 'B': MIVIA/VF binary format
 *   - 'L': LAD text format (unlabelled only)
 *   - 'D': DIMACS text format
 *
 * Vertex and edge labelling are combined into a single labelled flag
 * passed to the underlying reader, as the binary and DIMACS formats
 * do not distinguish between the two at the file level.
 *
 * @param filename        Path to the graph file
 * @param format          Format character: 'B', 'L', or 'D'
 * @param directed        If true, encode edge directions asymmetrically
 * @param edge_labelled   If true, read edge labels
 * @param vertex_labelled If true, read vertex labels
 * @return                Parsed Graph object
 * @throws                std::runtime_error if the format is unrecognised
 *                        or the underlying reader fails
 */
Graph readGraph(char *filename, char format, bool directed, bool edge_labelled, bool vertex_labelled)
{
	bool labelled = edge_labelled || vertex_labelled;
	switch (format)
	{
	case 'B':
		return readBinaryGraph(filename, directed, labelled);
	case 'L':
		return readLadGraph(filename, directed);
	case 'D':
		return readDimacsGraph(filename, directed, labelled);
	default:
		throw std::runtime_error("Unknown graph format.");
	}
}

/**
 * Builds adjacency lists from the adjacency matrix of g.
 *
 * Allocates and populates g.degree and g.adjlist by scanning each row
 * of g.adjmat. A non-zero entry at adjmat[i][j] indicates an edge from
 * i to j, which is recorded in adjlist[i]. For undirected graphs this
 * means both adjlist[i] and adjlist[j] will contain each other.
 *
 * Caller is responsible for freeing g.degree, each g.adjlist[i], and
 * g.adjlist when no longer needed.
 *
 * @param g  Graph whose adjacency lists are to be built;
 *           g.adjmat and g.n must already be initialised
 */
void set_adjlist(Graph &g)
{
	g.degree = new unsigned int[g.n]();
	g.adjlist = new unsigned int *[g.n];
	for (int i = 0; i < g.n; i++)
	{
		for (int j = 0; j < g.n; j++)
		{
			if (g.adjmat[i][j])
			{
				g.degree[i]++;
			}
		}
		g.adjlist[i] = new unsigned int[g.degree[i]];
		int idx = 0;
		for (int j = 0; j < g.n; j++)
		{
			if (g.adjmat[i][j])
			{
				g.adjlist[i][idx++] = j;
			}
		}
	}
}

/**
 * Computes negative modular symmetry classes for the vertices of g.
 *
 * Two vertices u and v are negatively symmetric if they share the same
 * open neighbourhood N^-(u) = N^-(v), meaning they are non-adjacent to
 * each other and have identical neighbour sets. Such vertices are
 * structurally interchangeable in any induced subgraph isomorphism and
 * form an independent set within their class.
 *
 * The algorithm identifies negative symmetry by checking, for each
 * unclassified vertex i, whether any neighbour-of-neighbour node_neg
 * has the same degree as i and the same neighbourhood (excluding i
 * itself). Only pairs with node_neg > i are considered to avoid
 * duplicate assignments.
 *
 * Vertices with degree 0 are skipped and left with class label 0,
 * as isolated vertices cannot share a neighbourhood with others.
 *
 * @param g        Graph with adjmat, adjlist, and degree initialised
 * @param EqClass  Output array of length g.n; EqClass[v] holds the
 *                 symmetry class label of vertex v. Vertices in the
 *                 same non-zero class are negatively symmetric.
 *                 Caller is responsible for freeing this array.
 */
void GetEqClass(Graph &g, ui *&EqClass)
{
	ui graph_size = g.n, label = 1, node = 0, node_neg = 0, node_nneg;
	EqClass = new ui[graph_size];
	bool equiv = true;
	for (ui i = 0; i < graph_size; i++)
	{
		EqClass[i] = 0;
	}
	for (ui i = 0; i < graph_size; i++)
	{
		if (EqClass[i] != 0 || g.degree[i] == 0)
		{
			continue;
		}
		node = g.adjlist[i][0];
		for (ui j = 0; j < g.degree[node]; j++)
		{
			node_neg = g.adjlist[node][j];
			if (g.degree[i] == g.degree[node_neg] && node_neg > i)
			{
				equiv = true;
				for (ui k = 0; k < g.degree[node_neg]; k++)
				{
					node_nneg = g.adjlist[node_neg][k];
					if (node_nneg != i && g.adjmat[i][node_nneg] == 0)
					{
						equiv = false;
						break;
					}
				}
				if (equiv)
				{
					EqClass[node_neg] = label;
				}
			}
		}
		EqClass[i] = label;
		++label;
	}
}

/**
 * Groups the leaf neighbours of each vertex by their (edge label, vertex label) pair.
 *
 * A leaf is a vertex with degree exactly 1. For each vertex u, this
 * function finds all leaf neighbours v and groups them by the pair
 * (adjmat[u][v], label[v]) - that is, leaves connected by the same
 * edge label and carrying the same vertex label are placed in the same
 * group. Within each group, leaf indices are stored in insertion order.
 * Groups are sorted lexicographically by label pair after all leaves
 * of u have been processed.
 *
 * This grouping is used by the LUM (Leaf vertex Union Match) strategy
 * in McSplit+LL and McSplit+DAL, which matches as many leaf pairs as
 * possible in a single step when a vertex pair (u, w) is committed to
 * the mapping. Leaves in the same label group are always in the same
 * bidomain and can be matched freely without splitting any other domain.
 *
 * Degree is computed locally here from adjmat rather than using
 * g.degree, as pack_leaves may be called before set_adjlist.
 *
 * @param g  Graph to process; g.adjmat, g.label, and g.leaves must
 *           be initialised. g.leaves[u] is populated in place.
 */
void pack_leaves(Graph &g)
{
	std::vector<int> deg(g.n, 0);
	for (int i = 0; i < g.n; i++)
		for (int j = 0; j < g.n; j++)
			if (i != j && g.adjmat[i][j])
				deg[i]++;

	for (int u = 0; u < g.n; u++)
	{
		for (int v = 0; v < g.n; v++)
		{
			if (g.adjmat[u][v] && u != v && deg[v] == 1)
			{
				std::pair<unsigned int, unsigned int> labels(g.adjmat[u][v], g.label[v]);
				int pos = -1;
				for (int k = 0;; k++)
				{
					if (k == (int)g.leaves[u].size())
						g.leaves[u].push_back(std::make_pair(labels, std::vector<int>()));
					if (g.leaves[u][k].first == labels)
					{
						pos = k;
						break;
					}
				}
				g.leaves[u][pos].second.push_back(v);
			}
		}
		std::sort(g.leaves[u].begin(), g.leaves[u].end());
	}
}