#pragma once

#include <vector>
#include <string>

using ui = unsigned int;

struct Graph {
    int n;
    std::vector<std::vector<unsigned int>> adjmat;
    std::vector<unsigned int> label;
    
    // For SymSplit/RRSplit - adjacency list representation
    unsigned int *degree;
    unsigned int **adjlist;
    
    // For McSplit+LL - leaf vertex precomputation
    std::vector<std::vector<std::pair<std::pair<unsigned int, unsigned int>, std::vector<int>>>> leaves;

    Graph(unsigned int n);
    ~Graph();
};

Graph induced_subgraph(Graph& g, std::vector<int> vv);
Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);
void set_adjlist(Graph& g);
void GetEqClass(Graph& g, ui*& EqClass);