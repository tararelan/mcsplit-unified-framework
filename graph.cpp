#define _GNU_SOURCE

#include "graph.h"
#include <stdexcept>
#include <climits>

// Constructor
Graph::Graph(unsigned int n): n(n), adjmat(n, std::vector<unsigned int>(n, 0)), label(n, 0), degree(nullptr), adjlist(nullptr) {}

// Destructor
// Needs to free degree and adjlist if they were allocated
Graph::~Graph() {
	if (adjlist != nullptr) {
		for (int i = 0; i < n; i ++) {
			delete[] adjlist[i];
		}
		delete[] adjlist;
	}
	if (degree != nullptr) {
		delete[] degree;
	}
}

static void add_edge(Graph& g, int v, int w, unsigned int edge_val, bool directed, bool edge_labelled) {
    if (v != w) {
        if (edge_labelled) {
            if (directed) {
                g.adjmat[v][w] |= edge_val;
                g.adjmat[w][v] |= (edge_val << 16);
            } else {
                g.adjmat[v][w] = edge_val;
                g.adjmat[w][v] = edge_val;
            }
        } else {
            if (directed) {
                g.adjmat[v][w] |= 1;
                g.adjmat[w][v] |= 2;
            } else {
                g.adjmat[v][w] = 1;
                g.adjmat[w][v] = 1;
            }
        }
    } else {
        // loop: set most significant bit of vertex label
        g.label[v] |= (1u << (sizeof(unsigned int) * CHAR_BIT - 1));
    }
}

Graph induced_subgraph(Graph& g, std::vector<int> vv) {
	Graph subg(vv.size());
	for (int i = 0; i < subg.n; i++) {
		for (int j = 0; j < subg.n; j ++) {
			subg.adjmat[i][j] = g.adjmat[vv[i]][vv[j]];
		}
	}
			
	for (int i = 0; i < subg.n; i ++) {
		subg.label[i] = g.label[vv[i]];
	}
		
	return subg;
}


static int read_word(FILE*fp) {
	unsigned char a[2];
	if (fread(a, 1, 2, fp) != 2) {
		throw std::runtime_error("Error reading file.");
	}
	return (int)a[0] | (((int)a[1]) << 8);
}

// MIVIA / VF reader
static Graph readBinaryGraph(char* filename, bool directed, bool labelled) {
	FILE* f;
	if ((f = fopen(filename, "rb")) == NULL) {
		throw std::runtime_error("Cannot open file.");
	}

	int n = read_word(f);
	Graph g(n);

	int m = n * 33 / 100;
	int p = 1, k1 = 0, k2 = 0;
	while (p < m && k1 < 16) {
		p *= 2;
		k1 = k2;
		k2++;
	}

	for (int i = 0; i < n; i++) {
		int label = (read_word(f) >> (16 - k1));
		if (labelled) {
			g.label[i] |= label;
		}
	}

	for (int i = 0; i < n; i++) {
		int len = read_word(f);
		for (int j = 0; j < len; j++) {
			int target = read_word(f);
			int label = (read_word(f) >> (16 - k1)) + 1;
			add_edge(g, i, target, labelled ? label : 1, directed, labelled);
		}
	}

	fclose(f);
	return g;
}

// LAD reader
static Graph readLadGraph(char* filename, bool directed) {
	FILE* f;
	if ((f = fopen(filename, "r")) == NULL) {
		throw std::runtime_error("Cannot open file.");
	}

	int n = 0;
	if (fscanf(f, "%d", &n) != 1) {
		throw std::runtime_error("Number of vertices not read correctly.");
	}

	Graph g(n);

	for (int i = 0; i < n; i++) {
		int edge_count;
		if (fscanf(f, "%d", &edge_count) != 1) {
			throw std::runtime_error("Number of edges not read correctly.");
		}	
		for (int j = 0; j < edge_count; j++) {
			int w;
			if (fscanf(f, "%d", &w) != 1) {
				throw std::runtime_error("an edge was not read correctly.");
			}
			add_edge(g, i, w, 1, directed, false);
		}
	}
	fclose(f);
	return g;
}

// DIMACS reader
static Graph readDimacsGraph(char* filename, bool directed, bool labelled) {
	FILE* f;
	if ((f = fopen(filename, "r")) == NULL) {
		throw std::runtime_error("Cannot open file.");
	}

	char* line = NULL;
	size_t nchar = 0;
	int n = 0, m = 0, edges_read = 0;
	Graph* g = nullptr;

	while (getline(&line, &nchar, f) != -1) {
		if (nchar > 0) {
			switch (line[0]) {
				case 'p': {
					if (sscanf(line, "p edge %d %d", &n, &m) != 2) {
						throw std::runtime_error("Error reading p line.");
					}
					g = new Graph(n);
					break;
				}
				case 'e': {
					if (g == nullptr) {
						throw std::runtime_error("Graph size must be specified first.");
					}
					int v, w;
					if (sscanf(line, "e %d %d", &v, &w) != 2) {
						throw std::runtime_error("Error reading e line.");
					}
					add_edge(*g, v-1, w-1, 1, directed, labelled);
					edges_read++;
					break;
				}
				case 'n': {
					if (g == nullptr) {
						throw std::runtime_error("Graph size must be specified first.");
					}
					int v, label;
					if (sscanf(line, "n %d %d", &v, &label) != 2) {
						throw std::runtime_error("Error reading n line.");
					}
					if (labelled) {
						g->label[v-1] |= label;
					}
					break;
				}
			}
		}
	}
	free(line);
	fclose(f);

	if (g == nullptr) {
		throw std::runtime_error("Empty graph file.");
	}

	Graph result = *g;
	delete g;
	return result;
}

Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled) {
	bool labelled = edge_labelled || vertex_labelled;
	switch (format) {
		case 'B': return readBinaryGraph(filename, directed, labelled);
		case 'L': return readLadGraph(filename, directed);
		case 'D': return readDimacsGraph(filename, directed, labelled);
		default: throw std::runtime_error("Unknown graph format.");
	}
}

void set_adjlist(Graph &g) {
	g.degree = new unsigned int [g.n]();
	g.adjlist = new unsigned int*[g.n];
	for (int i = 0; i < g.n; i++) {
		for (int j = 0; j < g.n; j++) {
			if (g.adjmat[i][j]) {
				g.degree[i]++;
			}
		}
		g.adjlist[i] = new unsigned int[g.degree[i]];
		int idx = 0;
		for (int j = 0; j < g.n; j++) {
			if (g.adjmat[i][j]) {
				g.adjlist[i][idx++] = j;
			}
		}
	}
}

void GetEqClass(Graph& g, ui*& EqClass) {
	EqClass = new ui[g.n];
	for (int i = 0; i < g.n; i++) {
		EqClass[i] = i;
	}
	for (int i = 0; i < g.n; i++) {
		for (int j = i + 1; j < g.n; j++) {
			if (EqClass[j] != (ui)j) {
				continue;
			}
			bool equiv = true;
			for (int k = 0; k < g.n; k++) {
				if (k == i || k == j) {
					continue;
				}
				if (g.adjmat[i][k] != g.adjmat[j][k]) {
					equiv = false;
					break;
				}
			}
			if (equiv) {
				EqClass[j] = i;
			}
		}
	}
}