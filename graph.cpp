#include "graph.h"

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <string>

constexpr int BITS_PER_UNSIGNED_INT(CHAR_BIT * sizeof(unsigned int));

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

Node::Node(unsigned int id, unsigned int label) {
    this->id = id;
    this->original_id = id;
    this->label = label; // bidomain label
}

Graph::Graph(unsigned int n) {
    this->n = n;
    this->e = -1;
    for (unsigned int i = 0; i < n; ++i)
        adjlist.emplace_back(i, 0);
    leaves = std::vector<std::vector<std::pair<std::pair<unsigned int, unsigned int>, std::vector<int>>>> (n);
}

unsigned int Graph::get(const int u, const int v) const {
    if (u < this->n) {
        for (auto &edge : this->adjlist[u].adjNodes) {
            if (static_cast<int>(edge.id) == v) {
                return 1;
            }
        }
    }
    return 0;
}

void Graph::pack_leaves() {
    std::vector<int> deg(this->n, 0);

    for (int i = 0; i < this->n; i++)
        deg[i] += (int) this->adjlist[i].adjNodes.size();

    for (int u = 0; u < this->n; u++) {
        for (auto v: this->adjlist[u].adjNodes)
            if (deg[v.id] == 1) {
                std::pair<unsigned int, unsigned int> labels(1, this->adjlist[v.id].label);
                int pos = -1;
                for (int k = 0;; k++) {
                    if (k == int(this->leaves[u].size())) {
                        this->leaves[u].push_back(std::make_pair(labels, std::vector<int>()));
                    }
                    if (this->leaves[u][k].first == labels) {
                        pos = k;
                        break;
                    }
                }
                //            assert(pos != -1);
                this->leaves[u][pos].second.push_back(v.id);
            }
        sort(this->leaves[u].begin(), this->leaves[u].end());
    }
}

Graph induced_subgraph(struct Graph &g, std::vector<int> vv) {
    Graph subg(g.n);

#pragma omp parallel for
    for (int i = 0; i < subg.n; ++i) {
        subg.adjlist[i] = g.adjlist[vv[i]];
        subg.adjlist[i].id = i;
        for (int j = 0; j < (int) subg.adjlist[i].adjNodes.size(); ++j) {
            subg.adjlist[i].adjNodes[j].id =
                    std::find(vv.begin(), vv.end(), subg.adjlist[i].adjNodes[j].id) - vv.begin();
        }

        std::stable_sort(std::begin(subg.adjlist[i].adjNodes), std::end(subg.adjlist[i].adjNodes),
                         [&](Node a, Node b) { return a.id < b.id; });
    }

    subg.e = g.e;
    return subg;
}

void add_edge(Graph &g, int v, int w, bool directed = false, unsigned int val = 1) {
    if (v != w) {
        if (directed || val != 1) {
            std::cerr << "Error: this McSplit only supports undirected graphs with val=1" << std::endl;
            exit(1);
        } else {
            g.adjlist[v].adjNodes.push_back(Node(w, 0));
            g.adjlist[w].adjNodes.push_back(Node(v, 0));
        }
    } else {
        // To indicate that a vertex has a loop, we set the most
        // significant bit of its label to 1
        g.adjlist[v].label |= (1u << (BITS_PER_UNSIGNED_INT - 1));
    }
}

int Graph::computeNumEdges(){
    int nedges = 0;
    for (int i=0; i<this->n; i++)
        nedges += this->adjlist[i].adjNodes.size();
    this->e = nedges;
    return nedges;
}

/**
 * @return density = 2E / n(n-1)
 */
float Graph::computeDensity(){
    if (this->e < 0)
        this->computeNumEdges();
    return 2*float(this->e)/float(this->n*(this->n-1));
}

struct Graph readDimacsGraph(char *filename, bool directed, bool vertex_labelled) {
    struct Graph g(0);

    FILE *f;

    if ((f = fopen(filename, "r")) == NULL)
        fail("Cannot open file");

    char *line = NULL;
    size_t nchar = 0;

    int nvertices = 0;
    int medges = 0;
    int v, w;
    int edges_read = 0;
    int label;

    while (getline(&line, &nchar, f) != -1) {
        if (nchar > 0) {
            switch (line[0]) {
                case 'p':
                    if (sscanf(line, "p edge %d %d", &nvertices, &medges) != 2)
                        fail("Error reading a line beginning with p.\n");
                    g = Graph(nvertices);
                    break;
                case 'e':
                    if (sscanf(line, "e %d %d", &v, &w) != 2)
                        fail("Error reading a line beginning with e.\n");
                    add_edge(g, v - 1, w - 1, directed);
                    edges_read++;
                    break;
                case 'n':
                    if (sscanf(line, "n %d %d", &v, &label) != 2)
                        fail("Error reading a line beginning with n.\n");
                    if (vertex_labelled)
                        g.adjlist[v - 1].label |= label;
                    break;
            }
        }
    }

    if (medges > 0 && edges_read != medges)
        fail("Unexpected number of edges.");

    fclose(f);
    return g;
}

struct Graph readLadGraph(char *filename, bool directed) {
    struct Graph g(0);
    FILE *f;

    if ((f = fopen(filename, "r")) == NULL)
        fail("Cannot open file");

    int nvertices = 0;
    int w;

    if (fscanf(f, "%d", &nvertices) != 1)
        fail("Number of vertices not read correctly.\n");
    g = Graph(nvertices);

    for (int i = 0; i < nvertices; i++) {
        int edge_count;
        if (fscanf(f, "%d", &edge_count) != 1)
            fail("Number of edges not read correctly.\n");
        for (int j = 0; j < edge_count; j++) {
            if (fscanf(f, "%d", &w) != 1)
                fail("An edge was not read correctly.\n");
            add_edge(g, i, w, directed);
        }
    }

    fclose(f);
    return g;
}

int read_word(FILE *fp) {
    unsigned char a[2];
    if (fread(a, 1, 2, fp) != 2)
        fail("Error reading file.\n");
    return (int) a[0] | (((int) a[1]) << 8);
}

struct Graph readBinaryGraph(char *filename, bool directed, bool edge_labelled,
                             bool vertex_labelled) {
    struct Graph g(0);
    FILE *f;

    if ((f = fopen(filename, "rb")) == NULL)
        fail("Cannot open file");

    int nvertices = read_word(f);
    std::cout << "Nvertices: " << nvertices << std::endl;
    g = Graph(nvertices);

    // Labelling scheme: see
    // https://github.com/ciaranm/cp2016-max-common-connected-subgraph-paper/blob/master/code/solve_max_common_subgraph.cc
    int m = g.n * 33 / 100;
    int p = 1;
    int k1 = 0;
    int k2 = 0;
    while (p < m && k1 < 16) {
        p *= 2;
        k1 = k2;
        k2++;
    }
    //std::cout << "labelled: " << vertex_labelled << std::endl;
    for (int i = 0; i < nvertices; i++) {
        int label = (read_word(f) >> (16 - k1));
        //std::cout << "label: " << label << std::endl;
        if (vertex_labelled)
            g.adjlist[i].label |= label;
    }
    //std::cout << "edge_labelled: " << edge_labelled << std::endl;
    for (int i = 0; i < nvertices; i++) {
        int len = read_word(f);
        //std::cout << "len: " << len << std::endl;
        for (int j = 0; j < len; j++) {
            int target = read_word(f);
            int label = (read_word(f) >> (16 - k1)) + 1;
            add_edge(g, i, target, directed, edge_labelled ? label : 1);
        }
    }
    fclose(f);
    return g;
}

struct Graph readASCIIGraph(char *filename) {
    struct Graph g(0);
    FILE *f;

    if ((f = fopen(filename, "r")) == NULL)
        fail("Cannot open file");

    int nvertices, nedges;
    if (fscanf(f, "%d %d", &nvertices, &nedges) != 2)
        fail("Number of nodes and edges not read correctly.\n");

    std::cout << "nvertices: " << nvertices << std::endl;
    g = Graph(nvertices);
    int v1, v2;
    for (int i = 0; i < nedges; i++) {
        if (fscanf(f, "%d %d", &v1, &v2) != 2)
            fail("Bad edge format.\n");

        add_edge(g, v1, v2, false, 1);
    }
    g.e = nedges;
    fclose(f);
    return g;
}

struct Graph readGraph(char *filename, char format, bool directed, bool edge_labelled, bool vertex_labelled) {
    struct Graph g(0);
    if (format == 'D')
        g = readDimacsGraph(filename, directed, vertex_labelled);
    else if (format == 'L')
        g = readLadGraph(filename, directed);
    else if (format == 'B')
        g = readBinaryGraph(filename, directed, edge_labelled, vertex_labelled);
    else if (format == 'A')
        g = readASCIIGraph(filename);
    else
        fail("Unknown graph format\n");
    return g;
}
