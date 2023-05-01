#include <iostream>
#include <algorithm>
#include <stack>
#include "mcs.h"
#include "reward.h"
#include <thread>
#include <condition_variable>
#include <list>

#define DEBUG true

using namespace std;

const int short_memory_threshold = 1e5;
const int long_memory_threshold = 1e9;

/*
void show(const vector<VtxPair> &current, const vector<Bidomain> &domains,
          const vector<int> &left, const vector<int> &right, Stats *stats) {
    cout << "Nodes: " << stats->nodes << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i = 0; i < current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    for (unsigned int i = 0; i < domains.size(); i++) {
        struct Bidomain bd = domains[i];
        cout << "Left  ";
        for (int j = 0; j < bd.left_len; j++)
            cout << left[bd.l + j] << " ";
        cout << std::endl;
        cout << "Right  ";
        for (int j = 0; j < bd.right_len; j++)
            cout << right[bd.r + j] << " ";
        cout << std::endl;
    }
    cout << "\n"
         << std::endl;
}*/

void delete_first_value_from_list(list<int> list, int value){
    for (auto it = list.begin(); it != list.end(); ++it) {
        if (*it == value) {
            list.erase(it);
            break;
        }
    }
}

bool reached_max_iter(Stats *stats) {
    return 0 < arguments.max_iter && arguments.max_iter < stats->nodes;
}

int calc_bound(const vector<Bidomain> &domains) {
    int bound = 0;
    for (const Bidomain &bd: domains) {
        bound += std::min(bd.left.size(), bd.right.size());
    }
    return bound;
}

int selectV_index(const Bidomain *bd, const Rewards &rewards) {
    gtype max_g = -1;
    int best_vtx = INT_MAX;
    for (int vtx: bd->left) {
        double vtx_reward = rewards.get_vertex_reward(vtx, false);
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
                double vtx_reward = rewards.get_vertex_reward(vtx, false);
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

// Returns length of left half of array
int partition(vector<int> &all_vv, int start, int len, const Graph &g, int index) {
    int i = 0;
    for (int j = 0; j < len; j++) {
        if (g.get(index, all_vv[start + j])) {
            std::swap(all_vv[start + i], all_vv[start + j]);
            i++;
        }
    }
    return i;
}

int remove_matched_vertex(vector<int> &arr, int start, int len, const vector<int> &matched) {
    int p = 0;
    for (int i = 0; i < len; i++) {
        if (!matched[arr[start + i]]) {
            std::swap(arr[start + i], arr[start + p]);
            p++;
        }
    }
    return p;
}

// multiway is for directed and/or labelled graphs
NewBidomainResult
generate_new_domains(const vector<Bidomain> &d, vector<VtxPair> &current, vector<int> &g0_matched,
                     vector<int> &g1_matched,
                     Bidomain &bd,
                     const Graph &g0, const Graph &g1, const int v, const int w) {
    current.emplace_back(v, w);
    g0_matched[v] = 1;
    g1_matched[w] = 1;

    int leaves_match_size = 0, v_leaf, w_leaf;
    for (unsigned int i = 0, j = 0; i < g0.leaves[v].size() && j < g1.leaves[w].size();) {
        if (g0.leaves[v][i].first < g1.leaves[w][j].first)
            i++;
        else if (g0.leaves[v][i].first > g1.leaves[w][j].first)
            j++;
        else {
            const vector<int> &leaf0 = g0.leaves[v][i].second;
            const vector<int> &leaf1 = g1.leaves[w][j].second;
            for (unsigned int p = 0, q = 0; p < leaf0.size() && q < leaf1.size();) {
                if (g0_matched[leaf0[p]])
                    p++;
                else if (g1_matched[leaf1[q]])
                    q++;
                else {
                    v_leaf = leaf0[p], w_leaf = leaf1[q];
                    p++, q++;
                    current.push_back(VtxPair(v_leaf, w_leaf));
                    g0_matched[v_leaf] = 1;
                    g1_matched[w_leaf] = 1;
                    leaves_match_size++;
                }
            }
            i++, j++;
        }
    }

    vector<Bidomain> new_d;
    new_d.reserve(d.size());
    int l, r, j = -1;
    int temp, total = 0;
    int unmatched_left_len, unmatched_right_len;
    for (const Bidomain &old_bd: d) {
        j++;
        l = old_bd.l;
        r = old_bd.r;
        if (leaves_match_size > 0 && !old_bd.is_adjacent) {
            unmatched_left_len = remove_matched_vertex(left, l, old_bd.left_len, g0_matched);
            unmatched_right_len = remove_matched_vertex(right, r, old_bd.right_len, g1_matched);
        } else {
            unmatched_left_len = old_bd.left_len;
            unmatched_right_len = old_bd.right_len;
        }
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, unmatched_left_len, g0, v);
        int right_len = partition(right, r, unmatched_right_len, g1, w);
        int left_len_noedge = unmatched_left_len - left_len;
        int right_len_noedge = unmatched_right_len - right_len;

        // compute reward
        temp = std::min(old_bd.left_len, old_bd.right_len) - std::min(left_len, right_len) -
               std::min(left_len_noedge, right_len_noedge);
        total += temp;

        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }

    NewBidomainResult result = {new_d, total};
    return result;
}

int getNeighborOverlapScores(const Graph &g0, const Graph &g1, const vector<VtxPair> &current, int v, int w) {
    int overlap_v = 0;
    int overlap_w = 0;
    // get number of selected neighbors of v
    for (auto &neighbor: g0.adjlist[v].adjNodes)
        for (auto j: current)
            if (j.v == neighbor.id) {
                overlap_v++;
                break;
            }
    // get number of selected neighbors of w
    for (auto &neighbor: g1.adjlist[w].adjNodes)
        for (auto j: current)
            if (j.w == neighbor.id) {
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
            pair_reward += rewards.get_pair_reward(v, vtx, false);

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

void remove_bidomain(vector<Bidomain> &domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

mutex steps_mutex;
condition_variable steps_cv;

vector<VtxPair>
solve(const Graph &g0, const Graph &g1, Rewards &rewards, vector<VtxPair> &incumbent, list<Step *> &global_steps,
      unsigned int matching_size_goal, Stats *stats) {

    list<Step *> steps;

    while (true) {
        // pop one step from the global stack
        unique_lock lk(steps_mutex);
        while (global_steps.empty() && !stats->abort_due_to_timeout && !reached_max_iter(stats)) {
            steps_cv.wait(lk);
        }
        if (stats->abort_due_to_timeout || reached_max_iter(stats))
            return incumbent;
        steps.emplace_back(global_steps.back());
        global_steps.pop_back();
        lk.unlock();

        // end cycle when there are no more steps, or when we have a W step after a certain number of steps
        // TODO remove the threshold for the first thread.
        // TODO adjust the threshold based on graph size (possibly, based on the depth at which the first pruning occurs?)
        while (!steps.empty() && steps.size() < 1000) {
            Step *s = steps.back();

            // check timeout
            if (stats->abort_due_to_timeout) {
                steps_cv.notify_one();
                return incumbent;
            }
            // TODO: unused since g0, g1 and current local to step, to be checked
            // delete eventual extra vertices from previous iterations
            // while (current.size() > s->cur_len) {
            //     VtxPair pr = current.back();
            //     s->g0_matched[pr.v] = 0;
            //     s->g1_matched[pr.w] = 0;
            //     current.pop_back();
            // }

            /* V-step */
            if (s->w_iter == -1) {
                stats->nodes++;

                // check max iterations
                if (reached_max_iter(stats)) {
                    cout << "Reached " << stats->nodes << " iterations" << endl;
                    steps_cv.notify_all();
                    return incumbent;
                }

                // If the current matching is larger than the incumbent matching, update the incumbent
                if (s->current.size() > incumbent.size()) {
                    incumbent.clear();
                    for (auto &pr: s->current) {
                        incumbent.push_back(pr);
                    }
                    if (!arguments.quiet) {
                        cout << "Incumbent size: " << incumbent.size() << endl;
                    }

                    rewards.update_policy_counter(true);
                }

                // Prune the branch if the upper bound is too small
                int bound = (int) s->current.size() + calc_bound(s->domains);
                if (bound <= incumbent.size() || bound < matching_size_goal) {
                    delete steps.back();
                    steps.pop_back();
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
                int v = selectV_index(s->bd, rewards);
                delete_first_value_from_list(bd->left, v);
                rewards.update_policy_counter(false);

                // Next iteration try to select a vertex w to pair with v (convert this v step to a w step)
                s->setBd(bd, bd_idx);
                s->w_iter = 0;
                s->v = v;
                continue;
            }

            /* W-step */
            if (s->w_iter < s->bd->right.size()) {
                int w = selectW_index(g0, g1, s->current, s->bd, rewards, s->v, s->wselected);
                delete_first_value_from_list(s->bd->right, w);
                s->wselected.insert(w);
                rewards.update_policy_counter(false);

#if DEBUG
                if (stats->nodes % 1000 == 0 && stats->nodes > 00) {
                    cout << "nodes: " << stats->nodes << ", v: " << s->v << ", w: " << w << ", size: " << s->current.size()
                         << ", dom: " << s->bd->left.size() << " " << s->bd->right.size() << endl;
                }
#endif
                // TODO check these are deep copies
                vector<VtxPair> new_current = s->current;
                Bidomain new_bd = *s->bd;
                auto result = generate_new_domains(s->domains, new_current, s->g0_matched, s->g1_matched, new_bd, g0,
                                                   g1, s->v,
                                                   w);
                rewards.update_rewards(result, s->v, w, stats);
                s->w_iter++;
                // if this is the last W vertex, remove current W step
                if(s->w_iter >= s->bd->right.size())
                    steps.pop_back();

                // next iterations select a new vertex v
                Step *s2 = new Step(result.new_domains, -1, -1, new_current, s->g0_matched, s->g1_matched);
                steps.emplace_back(s2);
                continue;
            }
        }

        // If the stack is not empty, push all steps in global stack
        if (!steps.empty()) {
            unique_lock lk2(steps_mutex);
            // TODO copy local steps into global steps. We need to think about the order (depth first/priority first?)
            // copy only W steps to global stack
            for (auto &step: steps) {
                if (step->w_iter > -1)
                    global_steps.push_back(step);
                else
                    delete step;
            }
            lk2.unlock();
            steps_cv.notify_all();
        }
    }
}

vector<VtxPair> mcs(const Graph &g0, const Graph &g1, void *rewards_p, Stats *stats) {

    vector<int> g0_matched(g0.n, 0);
    vector<int> g1_matched(g1.n, 0);

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
    list<Step *> steps;
    // TODO remove g0_matched
    Step *sp = new Step(std::move(domains), -1, -1, vector<VtxPair>(), g0_matched,
                        g1_matched);
    steps.emplace_back(sp);
    vector<VtxPair> incumbent;
    vector<thread> threads;
    for (int i = 0; i < arguments.threads; i++) {
        stats[i].start = clock();
        stats[i].nodes = 0;
        threads.emplace_back(solve, g0, g1, ref(rewards), ref(incumbent), ref(steps), 1, stats);
    }

    for (std::thread &t: threads)
        t.join();

    if (arguments.timeout && double(clock() - stats->start) / CLOCKS_PER_SEC > arguments.timeout) {
        cout << "time out" << endl;
    }

    return incumbent;
}
