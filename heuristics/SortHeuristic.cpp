#include "SortHeuristic.h"
#include <algorithm>
#include <thread>
#include <stack>
#include <queue>
#include <set>

#define VERBOSE false

namespace SortHeuristic {
    vector<int> Degree::sort(const Graph &g) {
        if (VERBOSE) std::cout << "Sorting by degree" << std::endl;
        vector<int> degree(g.n, 0);
        for (int v = 0; v < g.n; v++) {
            degree[v] = g.adjlist[v].adjNodes.size() - 1;
        }
        return degree;
    }

    vector<int> PageRank::sort(const Graph &g) {
        if (VERBOSE) std::cout << "Sorting by PageRank" << std::endl;
        constexpr float damping_factor = 0.85f;
        constexpr float epsilon = 0.00001f;
        std::vector<int> out_links = std::vector(g.n, 0);
        for (int i = 0; i < g.n; i++) {
            out_links[i] = g.adjlist[i].adjNodes.size();
        }
        // create a stochastic matrix (inefficient for big/sparse graphs, could be transformed into adj list (or just use the Graph's internal adj_list?))
        std::vector<std::vector<float>> stochastic_g = std::vector(g.n, std::vector(g.n, 0.0f));
        for (int i = 0; i < g.n; i++) {
            if (!out_links[i]) {
                for (int j = 0; j < g.n; j++) {
                    stochastic_g[i][j] = 1.0f / (float) g.n;
                }
            } else {
                for (auto &w: g.adjlist[i].adjNodes) {
                    stochastic_g[i][w.id] = 1.0f / (float) out_links[i];
                }
            }
        }
        std::vector<int> result(g.n, 0);
        std::vector<float> ranks(g.n, 0);
        std::vector<float> p(g.n, 1.0 / g.n);
        std::vector<std::vector<float>> transposed = std::vector(g.n, std::vector(g.n, 0.0f));
        // transpose matrix
        for (int i = 0; i < g.n; i++) {
            for (int j = 0; j < g.n; j++) {
                transposed[i][j] = stochastic_g[j][i];
            }
        }
        while (true) {
            std::fill(ranks.begin(), ranks.end(), 0);
            for (int i = 0; i < g.n; i++) {
                for (int j = 0; j < g.n; j++) {
                    ranks[i] = ranks[i] + transposed[i][j] * p[j];
                }
            }
            for (int i = 0; i < g.n; i++) {
                ranks[i] = damping_factor * ranks[i] + (1.0 - damping_factor) / (float) g.n;
            }
            float error = 0.0f;
            for (int i = 0; i < g.n; i++) {
                error += std::abs(ranks[i] - p[i]);
            }
            if (error < epsilon) {
                break;
            }

            for (int i = 0; i < g.n; i++) {
                p[i] = ranks[i];
            }
        }
        for (int i = 0; i < ranks.size(); i++) {
            result[i] = ranks[i] / epsilon;
        }
        return result;
    }

    // https://en.wikipedia.org/wiki/Clustering_coefficient
    vector<int> LocalClusteringCoefficient::sort(const Graph &g) {
        if (VERBOSE) std::cout << "Sorting by Local Clustering Coefficient" << std::endl;
        vector<int> result(g.n, 0);
        for (int i = 0; i < g.n; i++) {
            int degree = g.adjlist[i].adjNodes.size();
            if (degree < 2) {
                result[i] = 0;
                continue;
            }
            int num_triangles = 0;
            for (int j = 0; j < degree; j++) {
                int v = g.adjlist[i].adjNodes[j].id;
                for (int k = j + 1; k < degree; k++) {
                    if (g.get(v, g.adjlist[i].adjNodes[k].id) == 1) {
                        num_triangles++;
                    }
                }
            }
            result[i] = 2 * num_triangles * 100 / (degree * (degree - 1));
        }
        return result;
    }


    vector<int> KatzCentrality::sort(const Graph &g) {
        if (VERBOSE) std::cout << "Sorting by Katz Centrality" << std::endl;
        constexpr float alpha = 0.5f;
        vector<int> result(g.n, 0);
        for (int i = 0; i < g.n; i++) {
            vector<int> visited(g.n, 0);
            vector<double> score(g.n, 0);
            set<int> next_layer, current_layer;
            score[i] = 1;
            next_layer.insert(i);
            // run BFS
            while (!next_layer.empty()) {
                current_layer = next_layer;
                next_layer.clear();
                for (auto &v: current_layer) {
                    for (auto &w: g.adjlist[v].adjNodes) {
                        if (visited[w.id] == 0) {
                            score[w.id] += score[v] * alpha;
                            next_layer.insert(w.id);
                        }
                    }
                }
                for (auto &w: next_layer) {
                    visited[w] = 1;
                }
            }
            double sum = 0;
            for (int j = 0; j < g.n; j++)
                sum += score[j] * 10;
            result[i] = (int) sum;
        }
        return result;
    }

    ///////////////////////////
    // PARALLEL IMPLEMENTATIONS
    ///////////////////////////


    vector<int> Parallel::sort(const Graph &g) {

        std::vector<std::thread> threads;
        std::atomic<size_t> index;
        BC_.resize(g.n);
        std::fill(begin(BC_), end(BC_), 0.0);

        index.store(g.n - 1);

        // run threads
        for (size_t i = 0; i < num_threads - 1; i++)
            threads.emplace_back(std::thread([this, &index, &g] { run_worker(&index, g); }));

        // start working
        run_worker(&index, g);

        // wait for others to finish
        for (auto &thread: threads)
            thread.join();
        vector<int> results = get_result_vector();
        return results;
    }

    // Adapted from https://github.com/chivay/betweenness-centrality
    void BetweennessCentrality::process(const Graph &g, const size_t &vertex_id, std::vector<double> *BC_local) {
        size_t vec_size = BC_local->size();
        std::stack<size_t> S;
        std::vector<std::vector<size_t> > P(vec_size);
        std::vector<int> sigma(vec_size);
        std::vector<int> d(vec_size);
        std::vector<double> delta(vec_size);

        for (size_t w = 0; w < g.n; w++) {
            sigma[w] = 0;
            d[w] = -1;
            delta[w] = 0;
        }

        sigma[vertex_id] = 1;
        d[vertex_id] = 0;

        std::queue<size_t> Q;
        Q.push(vertex_id);

        while (!Q.empty()) {
            size_t v = Q.front();
            Q.pop();
            S.push(v);

            for (const auto &node: g.adjlist[v].adjNodes) {
                size_t w = node.id;
                if (d[w] < 0) {
                    Q.push(w);
                    d[w] = d[v] + 1;
                }

                if (d[w] == d[v] + 1) {
                    sigma[w] += sigma[v];
                    P[w].emplace_back(v);
                }
            }
        }

        while (!S.empty()) {
            size_t v = S.top();
            S.pop();

            for (size_t p: P[v]) {
                double result = (double(sigma[p]) / sigma[v]) * (1.0 + delta[v]);
                delta[p] += result;
            }

            if (v != vertex_id) {
                (*BC_local)[v] += delta[v];
            }
        }
    }

    // Adapted from https://github.com/konstantinNovichenko/Modified-Dijkstra-Centrality-Closeness-Betweenness
    // Our version is multi threaded and does not support weighted graphs
    void ClosenessCentrality::process(const Graph &g, const size_t &vertex_id, std::vector<double> *BC_local) {
        int s, Minimum, u, infinity, temp[g.n], centralityScore = 0;;
        std::vector<bool> T = std::vector<bool>(g.n);//Nodes to be visted
        std::vector<int> L = std::vector<int>(g.n);//lambda
        std::vector<int> father = std::vector<int>(g.n);

        infinity = 1000;
        s = (int) vertex_id;

        //INITIALIZING DATA STRUCTURES
        for (int j = 0; j < g.n; j++) {
            T[j] = false;            //T=V; all the vertices are eligible to be visited
            L[j] = infinity;            // at the beginning every vertex is at distance ?, from s
            temp[j] = infinity;
            father[j] = -1;
        }
        //WE ASSUME THE SOURCE IS 0 FOR EVERY GRAPH
        L[s] = 0;                    // s is at distance 0 of itself
        temp[s] = 0;                // Temp has the same contents of L

        // Let u be the vertex in T that has minimum L clearly at the beginning u=s
        for (int j = 0; j < g.n; j++) {                //LOOP TROUGH ALL THE VERTICES OF THE GRAPH
            //cout<<endl<<"STEP "<<i<<":\n________ ";
            Minimum = infinity;
            for (int k = 0; k < g.n; k++) {
                if (T[k] == 0) {
                    if (Minimum > temp[k]) {
                        Minimum = temp[k];            //finding minimum L[s]
                        u = k;
                    }
                }
            }
            temp[u] = infinity;                //Assigning INFINITY to the data structure already visited to find the next minimum L
            for (int k = 0; k < g.adjlist[u].adjNodes.size(); k++) {
                const Node *w = &g.adjlist[u].adjNodes[k];
                if (!T[w->id]) {       // if w Exist in T, proceed
                    if (L[w->id] > L[u] + 1) {
                        L[w->id] = L[u] + 1; // w is closer to s by using u;
                        temp[w->id] = L[w->id];
                        father[w->id] = u;
                    }
                }
            }
            T[u] = true;                //Discard visited vertex u from T
        }

        for (int j = 0; j < g.n; j++) {
            centralityScore += L[j]; // get the total distance from the vertex to other vertices
        }

        (*BC_local)[vertex_id] = 1.0 / (double) centralityScore * 10 * g.n; // push to the list of all centrality scores
    }

    void Parallel::run_worker(std::atomic<size_t> *idx, const Graph &g) {
        std::vector<double> BC_local(g.n);
        std::fill(begin(BC_local), end(BC_local), 0.0);

        while (true) {
            int my_index = (int) (*idx)--;

            if (my_index < 0)
                break;

            process(g, my_index, &BC_local);
        }

        // Synchronized section
        {
            std::lock_guard<std::mutex> guard(bc_mutex_);
            for (size_t i = 0; i < BC_local.size(); i++) {
                BC_[i] += BC_local[i];
            }
        }
    }

    std::vector<int> Parallel::get_result_vector() const {
        std::vector<int> results;
        results.resize(BC_.size());

        for (size_t i = 0; i < BC_.size(); i++) {
            results[i] = static_cast<int>(BC_[i] * 100);
        }
        return results;
    }

}

