#ifndef GRAPH_H
#define GRAPH_H

#include <limits.h>
#include <stdbool.h>
#include <vector>

struct Node {
    unsigned int id;
    unsigned int original_id;
    unsigned int label;
    std::vector<Node> adjNodes;

    Node(unsigned int id, unsigned int label);
};

struct Graph {
    int n, e;
    std::vector<Node> adjlist;
    std::vector<std::vector<std::pair<std::pair<unsigned int, unsigned int>, std::vector<int>>>> leaves;

    Graph(unsigned int n);

    unsigned int get(const int u, const int v) const;

    void pack_leaves();

    int computeNumEdges();

    float computeDensity();
};

Graph induced_subgraph(struct Graph &g, std::vector<int> vv);

Graph readGraph(char *filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

#endif
