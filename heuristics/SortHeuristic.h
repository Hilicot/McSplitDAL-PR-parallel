#ifndef MCSPLITDAL_SORTHEURISTIC_H
#define MCSPLITDAL_SORTHEURISTIC_H

#include <atomic>
#include <mutex>
#include <vector>
#include <iostream>
#include "../graph.h"

using namespace std;

namespace SortHeuristic {
    class Base {
    protected:
        int num_threads = 1;
    public:
        virtual ~Base() {};
        virtual string name(){return "Base";};
        virtual vector<int> sort(const Graph &g){std::cout << "Warning: no sort is selected!" << std::endl; return vector<int>();}; 
        void set_num_threads(int num_threads_) { this->num_threads = num_threads_; }
    };

    class Degree : public Base {
    public:
        ~Degree() {};
        string name() override {return "Degree";};
        vector<int> sort(const Graph &g) override;
    };

    class PageRank : public Base {
    public:
        ~PageRank() {};
        string name() override {return "PageRank";};
        vector<int> sort(const Graph &g) override;
    };

    class LocalClusteringCoefficient : public Base {
    public:
        ~LocalClusteringCoefficient() {};
        string name() override {return "LocalClusteringCoefficient";};
        vector<int> sort(const Graph &g) override;
    };

    class KatzCentrality : public Base {
    public:
        ~KatzCentrality() {};
        string name() override {return "KatzCentrality";};
        vector<int> sort(const Graph &g) override;
    };

    /////////////// Parallel Heuristics ///////////////

    class Parallel : public Base {
    public:
        ~Parallel() {};
        vector<int> sort(const Graph &g) override;
        [[nodiscard]] vector<int> get_result_vector() const;
    private:
        virtual void process(const Graph &g, const size_t &vertex_id, std::vector<double> *BC_local) {cout << "Warning: no parallel sort is selected!" << endl;}
        void run_worker(std::atomic<size_t> *idx, const Graph &g);
        std::mutex bc_mutex_;
        std::vector<double> BC_;
        std::string type_;
    };

    class BetweennessCentrality : public Parallel {
    public:
        ~BetweennessCentrality() {};
        string name() override {return "BetweennessCentrality";};
        void process(const Graph &g, const size_t &vertex_id, std::vector<double> *BC_local) override;
    };

    class ClosenessCentrality : public Parallel {
    public:
        ~ClosenessCentrality() {};
        string name() override {return "ClosenessCentrality";};
        void process(const Graph &g, const size_t &vertex_id, std::vector<double> *BC_local) override;
    };

}
#endif //MCSPLITDAL_SORTHEURISTIC_H
