#include <iostream>
#include <algorithm>
#include <stack>
#include "mcs.h"
#include "reward.h"
#include <thread>
#include <list>

#define DEBUG 0

using namespace std;

const int short_memory_threshold = 1e5;
const int long_memory_threshold = 1e9;
mutex steps_mutex;
mutex reward_mutex;
mutex incumbent_mutex;
condition_variable steps_cv;
int block_size = -1;

bool reached_max_iter(Stats *stats) {
    return 0 < arguments.max_iter && arguments.max_iter < (int) stats->nodes;
}

int calc_bound(const vector<Bidomain> &domains) {
    int bound = 0;
    for (const Bidomain &bd: domains) {
        bound += std::min((int)bd.left.size(), (int)bd.right.size());
    }
    return bound;
}

int selectV_index(const Bidomain *bd, const Rewards &rewards) {
    gtype max_g = -1;
    int best_vtx = INT_MAX;
    for (int vtx: bd->left) {
        //unique_lock rlk(reward_mutex);  // Optimistic version
        double vtx_reward = rewards.get_vertex_reward(vtx, false);
        //rlk.unlock();

        if (vtx_reward > max_g) {
            best_vtx = vtx;
            max_g = vtx_reward;
        } else if (vtx_reward == max_g) {
            if (vtx < best_vtx) {
                best_vtx = vtx;
            }
        }
    }
    return best_vtx;
}

int select_bidomain(const vector<Bidomain> &domains, const Rewards &rewards,
                    int current_matching_size) {
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    double max_reward = -1;
    int tie_breaker;
    int i;
    int current;
    int best = -1;

    for (i = 0; i < (int)domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size > 0 && !bd.is_adjacent)
            continue;
        if (arguments.heuristic == rewards_based) {
            current = 0;
            for (int vtx: bd.left) {
                // unique_lock rlk(reward_mutex);  // Optimistic version
                double vtx_reward = rewards.get_vertex_reward(vtx, false);
                // rlk.unlock();

                current += vtx_reward;
            }
            if (current < max_reward) {
                max_reward = current;
                best = i;
            }
        } else {
            if (arguments.heuristic == heuristic_based) {
                current = 0;
                for (int vtx: bd.left) {
                    current += vtx;
                }
            } else
                current = arguments.heuristic == min_max ? std::max((int) bd.left.size(), (int) bd.right.size()) :
                          (int) bd.left.size() * (int) bd.right.size();
            if (current < min_size) {
                min_size = current;
                min_tie_breaker = bd.left.front();
                best = i;
            } else if (current == min_size) {
                tie_breaker = bd.left.front();
                if (tie_breaker < min_tie_breaker) {
                    min_tie_breaker = tie_breaker;
                    best = i;
                }
            }
        }
    }
    return best;
}

// multiway is for directed and/or labelled graphs
NewBidomainResult generate_new_domains(const vector<Bidomain> &d, vector<VtxPair> &current, Bidomain &bd, const Graph &g0, const Graph &g1, const int v, const int w) {
    set<int> left_excluded, right_excluded;

    left_excluded.insert(v);
    right_excluded.insert(w);
    
    current.emplace_back(v, w);

    int v_leaf, w_leaf;
    for (unsigned int i = 0, j = 0; i < g0.leaves[v].size() && j < g1.leaves[w].size();) {
        if (g0.leaves[v][i].first < g1.leaves[w][j].first)
            i++;
        else if (g0.leaves[v][i].first > g1.leaves[w][j].first)
            j++;
        else {
            const vector<int> &leaf0 = g0.leaves[v][i].second;
            const vector<int> &leaf1 = g1.leaves[w][j].second;
            for (unsigned int p = 0, q = 0; p < leaf0.size() && q < leaf1.size();) {
                if (!std::binary_search(bd.left.begin(), bd.left.end(), leaf0[p]))
                    p++;
                else if (!std::binary_search(bd.right.begin(), bd.right.end(), leaf1[q]))
                    q++;
                else {
                    v_leaf = leaf0[p], w_leaf = leaf1[q];
                    p++, q++;
                    current.emplace_back(v_leaf, w_leaf);
                    left_excluded.insert(v_leaf);
                    right_excluded.insert(w_leaf);
                }
            }
            i++, j++;
        }
    }

    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    int j = -1;
    int temp, total = 0;
    for (const Bidomain &old_bd: d) {
        list<int> left_matched, left_unmatched, right_matched, right_unmatched;
        j++;

        for (auto node: old_bd.left) {
            if (left_excluded.find(node) != left_excluded.end())
                continue;
            if (g0.get(v, node))
                left_matched.emplace_back(node);
            else
                left_unmatched.emplace_back(node);
        }

        for (auto node: old_bd.right) {
            if (right_excluded.find(node) != right_excluded.end())
                continue;
            if (g1.get(w, node))
                right_matched.emplace_back(node);
            else
                right_unmatched.emplace_back(node);
        }

        // compute reward
        temp = std::min(old_bd.left.size() - 1, old_bd.right.size() - 1) - std::min(left_matched.size(), right_matched.size()) -
               std::min(left_unmatched.size(), right_unmatched.size());
        total += temp;

        if (!left_unmatched.empty() && !right_unmatched.empty()) {
            new_d.emplace_back(std::move(left_unmatched), std::move(right_unmatched), old_bd.is_adjacent);
        }
        if (!left_matched.empty() && !right_matched.empty()) {
            new_d.emplace_back(std::move(left_matched), std::move(right_matched), true);
        }
    }

    left_excluded.clear();
    right_excluded.clear();

    NewBidomainResult result = {new_d, total};
    return result;
}

int getNeighborOverlapScores(const Graph &g0, const Graph &g1, const vector<VtxPair> &current, int v, int w) {
    int overlap_v = 0;
    int overlap_w = 0;
    // get number of selected neighbors of v
    for (auto &neighbor: g0.adjlist[v].adjNodes)
        for (auto j: current)
            if (j.v == (int) neighbor.id) {
                overlap_v++;
                break;
            }
    // get number of selected neighbors of w
    for (auto &neighbor: g1.adjlist[w].adjNodes)
        for (auto j: current)
            if (j.w == (int) neighbor.id) {
                overlap_w++;
                break;
            }

    return overlap_v + overlap_w;
}

int selectW_index(const Graph &g0, const Graph &g1, const vector<VtxPair> &current, const Bidomain *bd,
                  const Rewards &rewards, const int v, const set<int> &wselected) {
    gtype max_g = -1;
    int best_vtx = INT_MAX;
    for (int vtx : bd->right) {
        if (wselected.find(vtx) == wselected.end()) {
            gtype pair_reward = 0;
            // Compute overlap scores
            if (arguments.reward_policy.neighbor_overlap != NO_OVERLAP) {
                int overlap_score = getNeighborOverlapScores(g0, g1, current, v, vtx);
                pair_reward += overlap_score * 100;
            }
            // Compute regular reward for pair
            // unique_lock rlk(reward_mutex);  // Optimistic version
            pair_reward += rewards.get_pair_reward(v, vtx, false);
            // rlk.unlock();

            // Check if this is the best pair so far
            if (pair_reward > max_g) {
                best_vtx = vtx;
                max_g = pair_reward;
            } else if (pair_reward == max_g) {
                if (vtx < best_vtx) {
                    best_vtx = vtx;
                }
            }
        }
    }
    return best_vtx;
}

vector<VtxPair>
solve(const Graph &g0, const Graph &g1, Rewards &rewards, vector<VtxPair> &incumbent, list<Step *> &global_steps,
      unsigned int matching_size_goal, Stats *stats) {

    list<Step *> steps;
    while (true) {
        // pop one step from the global stack
        unique_lock lk(steps_mutex);
        while (global_steps.empty() && !stats->abort_due_to_timeout && !reached_max_iter(stats)) {
            stats->sleeping_threads++;
            if (stats->sleeping_threads == arguments.threads) {
                steps_cv.notify_all();
                break;
            }
            steps_cv.wait(lk);
            stats->sleeping_threads--;
        }
        if (stats->abort_due_to_timeout || reached_max_iter(stats) || stats->sleeping_threads == arguments.threads) {
            lk.unlock();
            return incumbent;
        }
        steps.emplace_back(global_steps.back());
        global_steps.pop_back();
        lk.unlock();

        // end cycle when there are no more steps, or when we have a W step after a certain number of steps
        int local_iter_count = 0;
        while (!steps.empty() && (block_size < 0 || (int) steps.size() < block_size || steps.back()->w_iter == -1)) {
            Step *s = steps.back();

            // check timeout
            if (stats->abort_due_to_timeout) {
                steps_cv.notify_all();
                while(!steps.empty()) delete steps.front(), steps.pop_front(); // Dealloc memory
                return incumbent;
            }

            /* V-step */
            if (s->w_iter == -1) {
                stats->nodes++;
                local_iter_count++;

                // check max iterations
                if (reached_max_iter(stats)) {
                    cout << "Reached " << stats->nodes << " iterations" << endl;
                    steps_cv.notify_all();
                    while(!steps.empty()) delete steps.front(), steps.pop_front(); // Dealloc memory
                    return incumbent;
                }

                // If the current matching is larger than the incumbent matching, update the incumbent
                unique_lock ilk(incumbent_mutex);
                if (s->current.size() > incumbent.size()) {
                    incumbent.clear();
                    for (auto &pr: s->current) {
                        incumbent.push_back(pr);
                    }
                    if (!arguments.quiet) {
                        cout << "Incumbent size: " << incumbent.size() << " after " << stats->nodes << " iterations" << endl;
                    }

                    stats->bestcount = stats->cutbranches + 1;
                    stats->bestnodes = stats->nodes;
                    stats->bestfind = clock();

                    unique_lock rlk(reward_mutex); // NB Check possible deadlock/Starvation!
                    rewards.update_policy_counter(true);
                    rlk.unlock();
                }
                ilk.unlock();

                // Prune the branch if the upper bound is too small
                int bound = (int) s->current.size() + calc_bound(s->domains);
                if (bound <= (int) incumbent.size() || bound < (int) matching_size_goal) {
                    delete steps.back();
                    steps.pop_back();
                    // If I am the first thread, set the block_size
                    if(block_size < 0){
                        block_size = arguments.max_thread_blocks;
                    }
                    stats->cutbranches++;
                    continue;
                }

                // Select a bidomain based on the heuristic
                int bd_idx = select_bidomain(s->domains, rewards, (int) s->current.size());
                if (bd_idx == -1) {
                    // In the MCCS case, there may be nothing we can branch on
                    continue;
                }
                auto bd = &s->domains[bd_idx];

                // Select vertex v (vertex with max reward)
                int v = selectV_index(bd, rewards);
  
                unique_lock rlk(reward_mutex);
                rewards.update_policy_counter(false);
                rlk.unlock();

                // Next iteration try to select a vertex w to pair with v (convert this v step to a w step)
                s->bd = bd;
                s->w_iter = 0;
                s->v = v;
                continue;
            }

            /* W-step */
            if (s->w_iter < (int) s->bd->right.size()) {
                int w = selectW_index(g0, g1, s->current, s->bd, rewards, s->v, s->wselected);
                s->wselected.insert(w);
                
                unique_lock rlk(reward_mutex);
                rewards.update_policy_counter(false);
                rlk.unlock();

#if DEBUG
                if (stats->nodes % 1 == 0) {
                    cout << "nodes: " << stats->nodes << ", v: " << s->v << ", w: " << w << ", size: " << s->current.size()
                         << ", dom: " << s->bd->left.size() - 1 << " " << s->bd->right.size() - 1 << endl; // ", steps: " << steps.size()<< endl;
                }
#endif

                // TODO check these are deep copies
                vector<VtxPair> new_current = s->current;
                Bidomain new_bd = *s->bd;
                auto result = generate_new_domains(s->domains, new_current, new_bd, g0, g1, s->v, w);
                
                rlk.lock();
                rewards.update_rewards(result, s->v, w, stats);
                rlk.unlock();
                
                s->w_iter++;
                // if this is the last W vertex, remove current W step
                if(s->w_iter >= (int) s->bd->right.size()) {
                    delete steps.back();
                    steps.pop_back();
                }

                // next iterations select a new vertex v
                Step *s2 = new Step(result.new_domains, -1, -1, new_current);
                steps.emplace_back(s2);
                
                continue;
            }
        }
#if DEBUG
        cout << "local_iter_count: " << local_iter_count << endl;
#endif // DEBUG

        // If the stack is not empty, push all steps in global stack
        if (!steps.empty()) {
            unique_lock lk2(steps_mutex);
            // copy only W steps to global stack
            for (auto &step: steps) {
                if (step->w_iter > -1)
                    global_steps.push_back(step);
            }
            lk2.unlock();
            steps_cv.notify_all();
            steps.clear();
        }
    }
}

vector<VtxPair> mcs(const Graph &g0, const Graph &g1, void *rewards_p, Stats *stats) {
    Rewards &rewards = *(Rewards *) rewards_p;

    auto domains = vector<Bidomain>{};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (auto node: g0.adjlist)
        left_labels.insert(node.label);
    for (auto node: g1.adjlist)
        right_labels.insert(node.label);
    std::set<unsigned int> labels; // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    // Create a bidomain for each label that appears in both graphs (only one at the start)
    for (unsigned int label: labels) {
        list<int> left;
        list<int> right;

        for (int i = 0; i < g0.n; i++)
            if (g0.adjlist[i].label == label)
                left.push_back(i);
        for (int i = 0; i < g1.n; i++)
            if (g1.adjlist[i].label == label)
                right.push_back(i);

        domains.emplace_back(left, right, false);
    }

    // Start threads
    stats->nodes = 0;
    if(!arguments.first_thread_goes_until_pruning) // by default the first thread goes until pruning because block_size = -1
        block_size = arguments.max_thread_blocks;   // so we set block_size to max_thread_blocks to disable it
    list<Step *> steps;
    // TODO remove g0_matched
    Step *sp = new Step(std::move(domains), -1, -1, vector<VtxPair>());
    steps.emplace_back(sp);
    vector<VtxPair> incumbent;
    vector<thread> threads;
    for (int i = 0; i < arguments.threads; i++) {
        cout << "Starting thread " << i+1 << " out of " << arguments.threads << endl; 
        //stats->start = clock();
        stats->nodes = 0;
        threads.emplace_back(solve, g0, g1, ref(rewards), ref(incumbent), ref(steps), 1, stats);
    }

    for (std::thread &t: threads)
        if (t.joinable())
            t.join();

    while(!steps.empty()) delete steps.front(), steps.pop_front(); // Dealloc memory

    if (arguments.timeout && double(clock() - stats->start) / CLOCKS_PER_SEC > arguments.timeout) {
        cout << "time out" << endl;
    }

    return incumbent;
}
