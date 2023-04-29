#include <iostream>
#include <algorithm>
#include <set>
#include "mcs.h"
#include "reward.h"

using namespace std;

const int short_memory_threshold = 1e5;
const int long_memory_threshold = 1e9;

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
}

int calc_bound(const vector<Bidomain> &domains) {
    int bound = 0;
    for (const Bidomain &bd: domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}

int selectV_index(const vector<int> &arr, const Rewards &rewards, int start_idx, int len) {
    int idx = -1;
    gtype max_g = -1;
    int vtx, best_vtx = INT_MAX;
    for (int i = 0; i < len; i++) {
        vtx = arr[start_idx + i];
        double vtx_reward = rewards.get_vertex_reward(vtx, false);
        if (vtx_reward > max_g) {
            idx = i;
            best_vtx = vtx;
            max_g = vtx_reward;
        } else if (vtx_reward == max_g) {
            if (vtx < best_vtx) {
                idx = i;
                best_vtx = vtx;
            }
        }
    }
    return idx;
}

int select_bidomain(const vector<Bidomain> &domains, const vector<int> &left, const Rewards &rewards,
                    int current_matching_size) {
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    double max_reward = -1;
    int tie_breaker;
    unsigned int i;
    int current;
    int best = -1;

    for (i = 0; i < domains.size(); i++) {
        const Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size > 0 && !bd.is_adjacent)
            continue;
        if (arguments.heuristic == rewards_based) {
            current = 0;
            for (int j = bd.l; j < bd.l + bd.left_len; j++) {
                int vtx = left[j];
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
                for (int j = bd.l; j < bd.l + bd.left_len; j++) {
                    current += left[j];
                }
            } else
                current = arguments.heuristic == min_max ? std::max(bd.left_len, bd.right_len) : bd.left_len *
                                                                                                 bd.right_len;
            if (current < min_size) {
                min_size = current;
                min_tie_breaker = left[bd.l + selectV_index(left, rewards, bd.l, bd.left_len)];
                best = i;
            } else if (current == min_size) {
                tie_breaker = left[bd.l + selectV_index(left, rewards, bd.l, bd.left_len)];
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
generate_new_domains(const vector<Bidomain> &d, int bd_idx, vector<VtxPair> &current, vector<int> &g0_matched,
                     vector<int> &g1_matched,
                     vector<int> &left, vector<int> &right,
                     const Graph &g0, const Graph &g1, int v, int w,
                     bool multiway, Stats *stats) {
    current.push_back(VtxPair(v, w));
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
        if (leaves_match_size > 0 && old_bd.is_adjacent == false) {
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

#ifdef DEBUG

        printf("adj=%d ,noadj=%d ,old=%d ,temp=%d \n", std::min(left_len, right_len),
               std::min(left_len_noedge, right_len_noedge), std::min(old_bd.left_len, old_bd.right_len), temp);
        cout << "j=" << j << "  idx=" << bd_idx << endl;
        cout << "gl=" << lgrade[v] << " gr=" << rgrade[w] << endl;
#endif
        if (left_len_noedge && right_len_noedge)
            new_d.push_back({l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin + left_len,
                      [&](int a, int b) { return g0.get(v, a) < g0.get(v, b); });
            std::sort(r_begin, r_begin + right_len,
                      [&](int a, int b) { return g1.get(w, a) < g1.get(w, b); });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l < l_top && r < r_top) {
                unsigned int left_label = g0.get(v, left[l]);
                unsigned int right_label = g1.get(w, right[r]);
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do {
                        l++;
                    } while (l < l_top && g0.get(v, left[l]) == left_label);
                    do {
                        r++;
                    } while (r < r_top && g1.get(w, right[r]) == left_label);
                    new_d.push_back({lmin, rmin, l - lmin, r - rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});
        }
    }

#ifdef DEBUG0
    cout << "new domains are " << endl;
    for (const Bidomain &testd : new_d)
    {
        l = testd.l;
        r = testd.r;
        for (j = 0; j < testd.left_len; j++)
            cout << left[l + j] << " ";
        cout << " ; ";
        for (j = 0; j < testd.right_len; j++)
            cout << right[r + j] << " ";
        cout << endl;
    }
#endif
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

int selectW_index(const Graph &g0, const Graph &g1, const vector<VtxPair> &current, const vector<int> &arr,
                  const Rewards &rewards, const int v, int start_idx, int len,
                  const vector<int> &wselected) {
    int idx = -1;
    gtype max_g = -1;
    int vtx, best_vtx = INT_MAX;
    for (int i = 0; i < len; i++) {
        vtx = arr[start_idx + i];
        if (wselected[vtx] == 0) {
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
                idx = i;
                best_vtx = vtx;
                max_g = pair_reward;
            } else if (pair_reward == max_g) {
                if (vtx < best_vtx) {
                    idx = i;
                    best_vtx = vtx;
                }
            }
        }
    }
    return idx;
}

void remove_vtx_from_array(vector<int> &arr, int start_idx, int &len, int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

void remove_bidomain(vector<Bidomain> &domains, int idx) {
    domains[idx] = domains[domains.size() - 1];
    domains.pop_back();
}

void solve(const Graph &g0, const Graph &g1, Rewards &rewards,
           vector<VtxPair> &incumbent,
           vector<VtxPair> &current, vector<int> &g0_matched, vector<int> &g1_matched,
           vector<Bidomain> &domains, vector<int> &left, vector<int> &right, unsigned int matching_size_goal,
           Stats *stats) {
    // FIXME we have 2 timeout systems, remove one of them (the first seems to not work...)
    /*if (arguments.timeout && double(clock() - stats->start) / CLOCKS_PER_SEC > arguments.timeout) {
        return;
    }*/
    if (stats->abort_due_to_timeout)
        return;
    stats->nodes++;
    if (arguments.max_iter > 0 && stats->nodes > arguments.max_iter) {
        cout << "max_iter" << endl;
        return;
    }

    if (current.size() > incumbent.size()) { // incumbent 现任的
        incumbent = current;
        stats->bestcount = stats->cutbranches + 1;
        stats->bestnodes = stats->nodes;
        if (!arguments.quiet)
            cout << "Incumbent size: " << incumbent.size() << endl;
        stats->bestfind = clock();

        rewards.update_policy_counter(true);
    }

    // prune branch if upper bound is too small
    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= incumbent.size() || bound < matching_size_goal) {
        stats->cutbranches++;
        // cout << "nodes: " << stats->nodes  << " pruned" << endl;
        return;
    }
    // exit branch if goal already reached in big_first policy
    if (arguments.big_first && incumbent.size() == matching_size_goal)
        return;

    // select bidomain based on heuristic
    int bd_idx = select_bidomain(domains, left, rewards, current.size());
    if (bd_idx == -1) { // In the MCCS case, there may be nothing we can branch on
        return;
    }
    Bidomain &bd = domains[bd_idx];

    int v, w;
    int tmp_idx;

    // select vertex v (vertex with max reward)
    if(arguments.random_start && incumbent.size() == 0) // First vertex can optionally be random
        tmp_idx = rand() % bd.left_len;
    else
        tmp_idx = selectV_index(left, rewards, bd.l, bd.left_len);
    v = left[bd.l + tmp_idx];
    remove_vtx_from_array(left, bd.l, bd.left_len, tmp_idx); // remove v from bidomain
    rewards.update_policy_counter(false);

    // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
    vector<int> wselected(g1.n, 0);
    bd.right_len--; //
    for (int i = 0; i <= bd.right_len; i++) {
        tmp_idx = selectW_index(g0, g1, current, right, rewards, v, bd.r, bd.right_len + 1, wselected);
        w = right[bd.r + tmp_idx];
        wselected[w] = 1;
        std::swap(right[bd.r + tmp_idx], right[bd.r + bd.right_len]);
        rewards.update_policy_counter(false);
#if (DEBUG)
        if(stats->nodes % 100000 == 0)
            std::cout << "nodes: " << stats->nodes << ", v: " << v << ", w: " << w << ", size: " << current.size() << ", dom: "<< bd.left_len << " " << bd.right_len << std::endl;
#endif
        unsigned int cur_len = current.size();
        auto result = generate_new_domains(domains, bd_idx, current, g0_matched, g1_matched, left, right, g0, g1,
                                           v, w,
                                           arguments.directed || arguments.edge_labelled, stats);
        auto new_domains = result.new_domains;
        rewards.update_rewards(result, v, w, stats);

        stats->dl++;
        solve(g0, g1, rewards, incumbent, current, g0_matched, g1_matched, new_domains, left, right,
              matching_size_goal,
              stats);
        if (stats->abort_due_to_timeout) // hard timeout (else it gets stuck when trying to end the program gracefully)
            return;
        while (current.size() > cur_len) {
            VtxPair pr = current.back();
            current.pop_back();
            g0_matched[pr.v] = 0;
            g1_matched[pr.w] = 0;
        }
    }
    bd.right_len++;
    if (bd.left_len == 0)
        remove_bidomain(domains, bd_idx);
    solve(g0, g1, rewards, incumbent, current, g0_matched, g1_matched, domains, left, right, matching_size_goal,
          stats);
}

vector<VtxPair> mcs(const Graph &g0, const Graph &g1, void *rewards_p, Stats *stats) {
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right; // the buffer of vertex indices for the right partitions

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
        int start_l = left.size();
        int start_r = right.size();

        for (int i = 0; i < g0.n; i++)
            if (g0.adjlist[i].label == label)
                left.push_back(i);
        for (int i = 0; i < g1.n; i++)
            if (g1.adjlist[i].label == label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }

    vector<VtxPair> incumbent;

    if (arguments.big_first) {
        for (int k = 0; k < g0.n; k++) {
            unsigned int goal = g0.n - k;
            auto left_copy = left;
            auto right_copy = right;
            auto domains_copy = domains;
            vector<VtxPair> current;
            solve(g0, g1, rewards, incumbent, current, g0_matched, g1_matched, domains_copy, left_copy, right_copy,
                  goal,
                  stats);
            if (incumbent.size() == goal || stats->abort_due_to_timeout)
                break;
            if (!arguments.quiet)
                cout << "Upper bound: " << goal - 1 << std::endl;
        }
    } else {
        vector<VtxPair> current;
        solve(g0, g1, rewards, incumbent, current, g0_matched, g1_matched, domains, left, right, 1, stats);
    }

    if (arguments.timeout && double(clock() - stats->start) / CLOCKS_PER_SEC > arguments.timeout) {
        cout << "time out" << endl;
    }

    return incumbent;
}
