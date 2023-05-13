#include <iostream>
#include <algorithm>
#include "reward.h"

const int short_memory_threshold = 1e5;
const int long_memory_threshold = 1e9;

void Reward::update(gtype reward, gtype dal_reward) {
    this->rl_component += reward;
    this->ll_component += reward;
    this->dal_component += dal_reward;
}

void Reward::normalize(int reward_max, int dal_max, double factor) {
    if (arguments.mcs_method == RL_DAL)
        this->normalized_reward =
                factor * this->rl_component / reward_max + (1 - factor) * this->dal_component / dal_max;
    else
        this->normalized_reward =
                factor * this->ll_component / reward_max + (1 - factor) * this->dal_component / dal_max;
}

void Reward::reset(int value) {
    this->rl_component = value;
    this->ll_component = value;
    this->dal_component = 0;
}

gtype Reward::get_reward(int current_reward_policy, bool normalized) const {
    if (current_reward_policy == 0) // RL or LL policy
        if (arguments.mcs_method == RL_DAL)
            return this->rl_component;
        else // if (arguments.mcs_method == LL_DAL)
            return this->ll_component;
    else if (current_reward_policy == 1) // DAL policy
        if (normalized)
            return normalized_reward;
        else
            return ll_component + dal_component;
    else {
        std::cerr << "Error: unknown reward policy" << std::endl;
        exit(1);
    }
}

void Reward::decay() {
    this->ll_component /= 2;
    this->dal_component /= 2;
}

vector<Reward> DoubleQRewards::get_left_rewards() {
    return this->V;
}

vector<Reward> DoubleQRewards::get_right_rewards(int v) {
    if (arguments.mcs_method == RL_DAL && current_reward_policy == 0)
        return this->SingleQ;
    else // if (arguments.mcs_method == LL_DAL)
        return this->Q[v];
}

void DoubleQRewards::initialize(const vector<int> &left, const vector<int> &right) {
    left_initial_sort_order = left;
    right_initial_sort_order = right;
    for (int i = 0; i < (int) left.size(); i++) {
        V[i].rl_component = left[i];
        V[i].ll_component = left[i];
        for (int j = 0; j < (int) right.size(); j++) {
            Q[i][j].rl_component = right[j];
            Q[i][j].ll_component = right[j];
        }
    }

    if (arguments.mcs_method == RL_DAL) {
        for (int j = 0; j < (int) right.size(); j++) {
            SingleQ[j].rl_component = right[j];
            SingleQ[j].ll_component = right[j];
        }
    }
}

void DoubleQRewards::rotate_reward_policy() {
    current_reward_policy = (current_reward_policy + 1) %
                                                    arguments.reward_policy.reward_policies_num;
}

/**
 * reset the rewards to the initial sort order
 */
void DoubleQRewards::reset_rewards() {
    for (unsigned int i = 0; i < V.size(); i++) {
        V[i].reset(left_initial_sort_order[i]);
        for (unsigned int j = 0; j < Q[i].size(); j++) {
            Q[i][j].reset(right_initial_sort_order[j]);
        }
    }

    if (arguments.mcs_method == RL_DAL)
        for (int j = 0; j < (int) right_initial_sort_order.size(); j++)
            SingleQ[j].reset(right_initial_sort_order[j]);
}

void DoubleQRewards::randomize_rewards() {
    // TODO to verify if this is correct
    /*
    for (int i = 0; i < arguments.reward_policy.reward_policies_num; i++) {
        for (int j = 0; j < n; j++) {
            V[i][j] = (gtype) rand() / RAND_MAX;    // TODO what should be the upper bound? probably very low, so that it decays very fast
            for (int k = 0; k < m; k++) {
                Q[i][j][k] = (gtype) rand() / RAND_MAX;
            }
        }
    }*/
    std::cerr << "Randomize rewards not implemented yet" << std::endl;
    exit(1);
}

void DoubleQRewards::update_policy_counter(const bool restart_counter) {
    if (restart_counter) { // A better solution was found, reset the counter
        policy_switch_counter = 0;
    } else { // Increase the policy counter
        policy_switch_counter++;
        if (policy_switch_counter > arguments.reward_policy.reward_switch_policy_threshold) {
            policy_switch_counter = 0;
            switch (arguments.reward_policy.switch_policy) {
                case NO_CHANGE:
                    // Do nothing
                    break;
                case CHANGE:
                    rotate_reward_policy();
                    break;
                case RESET: {
                    rotate_reward_policy();
                    reset_rewards();
                    break;
                }
                case RANDOM: {
                    rotate_reward_policy();
                    randomize_rewards();
                    break;
                }
                case STEAL:
                    std::cerr << "Steal policy not implemented yet" << std::endl;
                    break;
            }
        }
    }
}

void DoubleQRewards::update_rewards(const NewBidomainResult &new_domains_result, int v, int w, Stats *stats) {
    gtype reward = new_domains_result.reward;
    const vector<Bidomain> &new_domains = new_domains_result.new_domains;

    // Compute DAL reward
    gtype dal_reward = 0;
    if (arguments.reward_policy.dal_reward_policy == DAL_REWARD_MAX_NUM_DOMAINS)
        dal_reward = new_domains.size();
    else if (arguments.reward_policy.dal_reward_policy == DAL_REWARD_MIN_MAX_DOMAIN_SIZE) {
        auto max_bidomain = std::max_element(new_domains.begin(), new_domains.end(),
                                             [](const Bidomain &bd1, const Bidomain &bd2) {
                                                 return bd1.get_max_len() < bd2.get_max_len();
                                             });
        dal_reward = -max_bidomain->get_max_len() / 100; // partial rounding + normalization (to remove?)
    } else if (arguments.reward_policy.dal_reward_policy == DAL_REWARD_MIN_AVG_DOMAIN_SIZE) {
        int total = 0;
        for (const Bidomain &bd: new_domains) {
            total += bd.get_max_len();
        }
        dal_reward = -total / new_domains.size() / 100; // partial rounding + normalization (to remove?)
    }

    // update rewards
    if (reward > 0) {
        stats->conflicts++;

        V[v].update(reward, dal_reward);
        SingleQ[w].update(reward, dal_reward);
        Q[v][w].update(reward, dal_reward);

        // Do not decay if current policy is RL!
        if (arguments.mcs_method != RL_DAL || current_reward_policy != 0) {
            // TODO if we normalize, we might have to adjust the thresholds
            if (get_vertex_reward(v, false) > short_memory_threshold)
                for (auto &r: V)
                    r.decay();
            if (get_pair_reward(v, w, false) > long_memory_threshold)
                for (auto &r: Q[v])
                    r.decay();
        }
    }
}

gtype DoubleQRewards::get_vertex_reward(int v, bool normalized) const {
    return V[v].get_reward(current_reward_policy, normalized);
}

gtype DoubleQRewards::get_pair_reward(int v, int w, bool normalized) const {
    if (arguments.mcs_method == RL_DAL && current_reward_policy == 0)
        return SingleQ[w].get_reward(current_reward_policy, normalized);
    else // if (arguments.mcs_method == LL_DAL)
        return Q[v][w].get_reward(current_reward_policy, normalized);
}

