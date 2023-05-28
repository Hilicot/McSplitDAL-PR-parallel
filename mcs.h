#ifndef MCSPLIT_MCS_H
#define MCSPLIT_MCS_H

#include <vector>
#include "graph.h"
#include "args.h"
#include "stats.h"
#include <set>
#include <memory>
#include <list>
#include <condition_variable>
#include <atomic>


using namespace std;
using gtype = double;

struct VtxPair {
    int v;
    int w;

    VtxPair(int v, int w) : v(v), w(w) {}
};

struct Bidomain {
    list<int> left;
    list<int> right;
    bool is_adjacent;

    // Bidomain(list<int> left, list<int> right, bool is_adjacent) : left(left), right(right), is_adjacent(is_adjacent) {};
    Bidomain(const list<int>& left, const list<int>& right, bool is_adjacent): left(left), right(right), is_adjacent(is_adjacent) {}

    int get_max_len() const { return max(left.size(), right.size()); }
};

struct NewBidomainResult {
    vector<Bidomain> new_domains;
    int reward;

    NewBidomainResult(vector<Bidomain>& new_domains, int reward): new_domains(new_domains), reward(reward) {}

};

struct Step {
    vector<Bidomain> domains;
    set<int> wselected;
    int w_iter;
    Bidomain *bd;
    int v;
    vector<VtxPair> current;

    Step(vector<Bidomain> &domains, int w_iter, int v, vector<VtxPair> current) : domains(domains), w_iter(w_iter), bd(nullptr), v(v), current(current) {
        this->wselected = set<int>();
    };

};

vector<VtxPair> mcs(const Graph &g0, const Graph &g1, void *rewards_p, Stats *stats);

#endif