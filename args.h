#ifndef MCSPLITDAL_ARGS_H
#define MCSPLITDAL_ARGS_H

#include "heuristics/SortHeuristic.h"

#ifdef MCSPLITDAL_MCSPLIT_DAL_H
#define EXTERN
#else
#define EXTERN extern
#endif

enum SwapPolicy {
    NO_SWAP,
    McSPLIT_SD,
    McSPLIT_SO,
    ADAPTIVE
};

enum Heuristic {
    min_max,
    min_product,
    rewards_based,
    heuristic_based
};

enum RewardSwitchPolicy {
    NO_CHANGE,
    CHANGE,
    RESET,
    RANDOM,
    STEAL
};

enum DAL_RewardPolicy {
    DAL_REWARD_MAX_NUM_DOMAINS,
    DAL_REWARD_MIN_MAX_DOMAIN_SIZE,
    DAL_REWARD_MIN_AVG_DOMAIN_SIZE,
};

enum NeighborOverlap {
    NO_OVERLAP,
    DAL_OVERLAP,
    RL_DAL_OVERLAP
};

struct RewardPolicy {
    RewardSwitchPolicy switch_policy;
    float reward_coefficient;
    int reward_switch_policy_threshold;
    int reward_policies_num;
    int current_reward_policy;
    int policy_switch_counter;
    DAL_RewardPolicy dal_reward_policy;
    NeighborOverlap neighbor_overlap;

    RewardPolicy() : reward_coefficient(1.0), reward_switch_policy_threshold(0), reward_policies_num(2),
                     current_reward_policy(1), policy_switch_counter(0), neighbor_overlap(NO_OVERLAP) {}
};

enum MCS {
    RL_DAL, LL_DAL
};

EXTERN struct arguments {
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool ascii;
    bool connected;
    bool directed;
    bool edge_labelled;
    bool vertex_labelled;
    bool big_first;
    bool random_start;
    Heuristic heuristic;
    SortHeuristic::Base *sort_heuristic;
    bool initialize_rewards;
    MCS mcs_method;
    char *filename1;
    char *filename2;
    int timeout;
    int threads;
    int max_thread_blocks;
    int max_iter;
    int arg_num;
    SwapPolicy swap_policy;
    RewardPolicy reward_policy;
} arguments;

#endif // MCSPLITDAL_ARGS_H
