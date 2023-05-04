#ifndef MCSPLITDAL_STATS_H
#define MCSPLITDAL_STATS_H

#include <chrono>
#include <atomic>

typedef struct Stats {
    std::atomic<unsigned long long> nodes{0};
    unsigned long long cutbranches{0};
    unsigned long long conflicts = 0;
    clock_t bestfind;
    unsigned long long bestnodes = 0, bestcount = 0;
    int dl = 0;
    clock_t start;
    std::atomic<bool> abort_due_to_timeout;
    int sleeping_threads = 0;
    bool swapped_graphs = false;
} Stats;

#endif //MCSPLITDAL_STATS_H
