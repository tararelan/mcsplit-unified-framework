#pragma once

#include <vector>
#include <chrono>

// A matched pair of vertices, one from each graph
struct VtxPair {
	int v; // vertex in G
	int w; // vertex in H
	VtxPair(int v, int w): v(v), w(w) {}
};

// A label class: paired sets of vertices from G and H sharing the same label
struct Bidomain {
	int l;	// start index of left set in left buffer
	int r;	// start index of right set in right buffer
	int left_len;
	int right_len;
	bool is_adjacent;
	Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
		l(l), r(r), left_len(left_len), right_len(right_len), is_adjacent(is_adjacent) {}
};

// Statistics collected during search
struct Stats {
    unsigned long long nodes = 0;           // total nodes explored
    unsigned long long cut_branches = 0;    // branches pruned by bound
    unsigned long long conflicts = 0;       // reward signal updates (RL algorithms)
    unsigned long long bound_pruned = 0;    // branches pruned by bound specifically
    unsigned long long sym_pruned = 0;      // branches pruned by symmetry/equivalence
    double time_elapsed = 0.0;             // wall clock time in seconds
    double time_to_best = 0.0;             // time to find incumbent solution
    unsigned long long nodes_to_best = 0;  // nodes when incumbent was found
    int incumbent_size = 0;                // size of best solution found
    double root_upper_bound = 0.0;         // upper bound at root node
    bool aborted = false;                  // true if timeout was reached
};