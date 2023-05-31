#include "mcsplitDAL+PR-iter.h"

using namespace std;

#define VSCORE
// #define DEBUG
#define Best
using std::cout;
using std::endl;
using std::vector;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}


/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format\vHEURISTIC can be min_max or min_product or rewards_based or heuristic_based";
static char args_doc[] = "HEURISTIC FILENAME1 FILENAME2";
static struct argp_option options[] = {
        {"quiet",                'q', 0,                   0, "Quiet output"},
        {"verbose",              'v', 0,                   0, "Verbose output"},
        {"dimacs",               'd', 0,                   0, "Read DIMACS format"},
        {"lad",                  'l', 0,                   0, "Read LAD format"},
        {"ascii",                'A', 0,                   0, "Read ASCII format"},
        {"connected",            'c', 0,                   0, "Solve max common CONNECTED subgraph problem"},
        {"directed",             'i', 0,                   0, "Use directed graphs"},
        {"labelled",             'a', 0,                   0, "Use edge and vertex labels"},
        {"vertex-labelled-only", 'x', 0,                   0, "Use vertex labels, but not edge labels"},
        {"big-first",            'b', 0,                   0, "First try to find an induced subgraph isomorphism, then decrement the target size"},
        {"timeout",              't', "timeout",           0, "Specify a timeout (seconds)"},
        {"iterations",           'I', "iterations",        0, "Specify a maximum number of iterations"},
        {"threads",              'p', "threads",           0, "Specify the number of threads"},
        {"max_thread_blocks",    'B', "blocks",            0, "Specify the number of steps a thread can work on locally"},
        {"random_start",         'r', 0,                   0, "Set random start to true"},
        {"dal_reward_policy",    'D', "dal_reward_policy", 0, "Specify the dal reward policy (num, max, avg)"},
        {"sort_heuristic",       's', "sort_heuristic",    0, "Specify the sort heuristic (degree, pagerank, betweenness, closeness, clustering, katz)"},
        {"pruning",              'P', 0,                   0, "Specify if the first thread goes on until pruning or not, before pushing to global queue"},
        {0}};

void set_default_arguments() {
    arguments.quiet = false;
    arguments.verbose = false;
    arguments.dimacs = false;
    arguments.lad = false;
    arguments.ascii = false;
    arguments.connected = false;
    arguments.directed = false;
    arguments.edge_labelled = false;
    arguments.vertex_labelled = false;
    arguments.big_first = false;
    arguments.filename1 = NULL;
    arguments.filename2 = NULL;
    arguments.timeout = 0;
    arguments.threads = 1;
    arguments.max_iter = 0;
    arguments.max_thread_blocks = 512;
    arguments.first_thread_goes_until_pruning = false;
    arguments.random_start = false;
    arguments.arg_num = 0;
    arguments.sort_heuristic = new SortHeuristic::Degree();
    arguments.initialize_rewards = false; // if false, rewards are initialized to 0, else to sort_heuristic
    arguments.mcs_method = RL_DAL;
    arguments.swap_policy = McSPLIT_SD;
    arguments.reward_policy.current_reward_policy = 1; // set starting policy (0:RL/LL, 1:DAL)
    arguments.reward_policy.reward_policies_num = 2;
    arguments.reward_policy.switch_policy = CHANGE;
    arguments.reward_policy.dal_reward_policy = DAL_REWARD_MAX_NUM_DOMAINS;
    arguments.reward_policy.neighbor_overlap = NO_OVERLAP;    // use neighbor overlap to select W
}

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'd':
            if (arguments.lad)
                fail("The -d and -l options cannot be used together.\n");
            arguments.dimacs = true;
            break;
        case 'l':
            if (arguments.dimacs)
                fail("The -d and -l options cannot be used together.\n");
            arguments.lad = true;
            break;
        case 'A':
            if (arguments.dimacs || arguments.lad)
                fail("The -d or -l options cannot be used together with -as.\n");
            arguments.ascii = true;
            break;
        case 'q':
            arguments.quiet = true;
            break;
        case 'v':
            arguments.verbose = true;
            break;
        case 'c':
            if (arguments.directed)
                fail("The connected and directed options can't be used together.");
            arguments.connected = true;
            break;
        case 'i':
            if (arguments.connected)
                fail("The connected and directed options can't be used together.");
            arguments.directed = true;
            break;
        case 'a':
            if (arguments.vertex_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.edge_labelled = true;
            arguments.vertex_labelled = true;
            break;
        case 'x':
            if (arguments.edge_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.vertex_labelled = true;
            break;
        case 'b':
            arguments.big_first = true;
            break;
        case 'P':
            arguments.first_thread_goes_until_pruning = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case 'I':
            arguments.max_iter = std::stoi(arg);
            break;
        case 'B':
            arguments.max_thread_blocks = std::stoi(arg);
            break;
        case 'p':
            arguments.threads = std::stoi(arg);
            break;
        case 'r':
            arguments.random_start = true;
            break;
        case 'D':
            if (string(arg) == "num")
                arguments.reward_policy.dal_reward_policy = DAL_REWARD_MAX_NUM_DOMAINS;
            else if (string(arg) == "max")
                arguments.reward_policy.dal_reward_policy = DAL_REWARD_MIN_MAX_DOMAIN_SIZE;
            else if (string(arg) == "avg")
                arguments.reward_policy.dal_reward_policy = DAL_REWARD_MIN_AVG_DOMAIN_SIZE;
            else
                fail("Unknown dal reward policy (try num, max, avg)");
            break;
        case 's':
            if (string(arg) == "degree")
                arguments.sort_heuristic = new SortHeuristic::Degree();
            else if (string(arg) == "pagerank")
                arguments.sort_heuristic = new SortHeuristic::PageRank();
            else if (string(arg) == "betweenness")
                arguments.sort_heuristic = new SortHeuristic::BetweennessCentrality();
            else if (string(arg) == "closeness")
                arguments.sort_heuristic = new SortHeuristic::ClosenessCentrality();
            else if (string(arg) == "clustering")
                arguments.sort_heuristic = new SortHeuristic::LocalClusteringCoefficient();
            else if (string(arg) == "katz")
                arguments.sort_heuristic = new SortHeuristic::KatzCentrality();
            else
                fail("Unknown sort heuristic (try degree, pagerank, betweenness, closeness, clustering, katz)");
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                if (std::string(arg) == "min_max")
                    arguments.heuristic = min_max;
                else if (std::string(arg) == "min_product")
                    arguments.heuristic = min_product;
                else if (std::string(arg) == "rewards_based")
                    arguments.heuristic = rewards_based;
                else if (std::string(arg) == "heuristic_based")
                    arguments.heuristic = heuristic_based;
                else
                    fail("Unknown heuristic (try min_max or min_product or rewards_based)");
            } else if (arguments.arg_num == 1) {
                arguments.filename1 = arg;
            } else if (arguments.arg_num == 2) {
                arguments.filename2 = arg;
            } else {
                argp_usage(state);
            }
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

/*******************************************************************************
                                 Main
*******************************************************************************/

bool check_sol(const Graph &g0, const Graph &g1, const vector<VtxPair> &solution) {
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i = 0; i < solution.size(); i++)
    {
        struct VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
        {
            cout << "Used left (v= " << p0.v << "): " << used_left[p0.v] << endl;
            cout << "Used right (w= " << p0.w << "): " << used_right[p0.w] << endl;
            return false;
        }
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.adjlist[p0.v].label != g1.adjlist[p0.w].label)
        {
            cout << g0.adjlist[p0.v].label << " != " << g1.adjlist[p0.w].label << endl;
            return false;
        }
        for (unsigned int j = i + 1; j < solution.size(); j++)
        {
            struct VtxPair p1 = solution[j];
            if (g0.get(p0.v, p1.v) != g1.get(p0.w, p1.w))
            {
                cout << "Edge left (" << p0.v << " -> " << p1.v << ") = " << g0.get(p0.v, p1.v) << endl;
                cout << "Edge right (" << p0.w << " -> " << p1.w << ") = " << g1.get(p0.w, p1.w) << endl;
                return false;
            }
        }
    }
    return true;
}

/**
 * based on arguments.swap_policy, return true if the graphs needs to be swapped.
 * McSPLIT_SD and McSPLIT_SO are based on Trimble's PHD thesis https://theses.gla.ac.uk/83350/
 */
bool swap_graphs(Graph &g0, Graph &g1) {
    switch (arguments.swap_policy) {
        case McSPLIT_SD: { // swap if density extremeness of g1 is bigger than that of g0
            // get densities
            float d0 = g0.computeDensity();
            float d1 = g1.computeDensity();
            // compute density extremeness
            double de0 = abs(0.5 - d0);
            double de1 = abs(0.5 - d1);
            return de1 > de0;
        }
        case McSPLIT_SO:
            return g1.n > g0.n;
        default:
            cerr << "swap policy unknown" << endl;
        case NO_SWAP:
            return false;
    }
}

int main(int argc, char **argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L'
                                                         : arguments.ascii ? 'A'
                                                                           : 'B';
    struct Graph g0 = readGraph(arguments.filename1, format, arguments.directed,
                                arguments.edge_labelled, arguments.vertex_labelled);
    struct Graph g1 = readGraph(arguments.filename2, format, arguments.directed,
                                arguments.edge_labelled, arguments.vertex_labelled);

    arguments.reward_policy.reward_switch_policy_threshold = 2 * std::min(g0.n, g1.n);

    Stats stats_s;
    Stats *stats = &stats_s;
    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    stats->abort_due_to_timeout.store(false);

    bool aborted = false;
#if 1
    if (0 != arguments.timeout) {
        timeout_thread = std::thread([&] {
            auto abort_time = std::chrono::steady_clock::now() + std::chrono::seconds(arguments.timeout);
            {
                /* Sleep until either we've reached the time limit,
                 * or we've finished all the work. */
                std::unique_lock<std::mutex> guard(timeout_mutex);
                while (!stats->abort_due_to_timeout.load()) {
                    if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                        /* We've woken up, and it's due to a timeout. */
                        aborted = true;
                        break;
                    }
                }
            }
            stats->abort_due_to_timeout.store(true);
        });
    }
#endif

    // decide whether to swap the graphs based on swap_policy
    if (swap_graphs(g0, g1)) {
        swap(g0, g1);
        stats->swapped_graphs = true;
        cout << "Swapped graphs" << endl;
    }

    //  auto start = std::chrono::steady_clock::now();
    stats->start = clock();

    // static sort order
    arguments.sort_heuristic->set_num_threads(10);
    std::vector<int> g0_deg = arguments.sort_heuristic->sort(g0);
    std::vector<int> g1_deg = arguments.sort_heuristic->sort(g1);

    // As implemented here, g1_dense and g0_dense are false for all instances
    // in the Experimental Evaluation section of the paper.  Thus,
    // we always sort the vertices in descending order of degree (or total degree,
    // in the case of directed graphs.  Improvements could be made here: it would
    // be nice if the program explored exactly the same search tree if both
    // input graphs were complemented.
    //  #ifdef VSCORE
    //    vector<int> lgrade(g0.n,0);
    //    vector<int> rgrade(g1.n,0);
    //  #endif
    vector<int> vv0(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = false; //sum(g1_deg) > g1.n * (g1.n - 1);
    std::stable_sort(std::begin(vv0), std::end(vv0),
                     [&](int a, int b) { return g1_dense ? (g0_deg[a] < g0_deg[b]) : (g0_deg[a] > g0_deg[b]); });

    vector<int> vv1(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = false; //sum(g0_deg) > g0.n * (g0.n - 1);
    std::stable_sort(std::begin(vv1), std::end(vv1),
                     [&](int a, int b) { //????????????????????????????????????????????????????
                         return g0_dense ? (g1_deg[a] < g1_deg[b]) : (g1_deg[a] > g1_deg[b]);
                     });
    std::cout << "Sorting done" << std::endl;


    struct Graph g0_sorted = induced_subgraph(g0, vv0);
    struct Graph g1_sorted = induced_subgraph(g1, vv1);
    clock_t time_elapsed = clock() - stats->start;
    std::cout << "Induced subgraph calculated in " << time_elapsed * 1000 / CLOCKS_PER_SEC << "ms" << endl;

    g0_sorted.pack_leaves();
    g1_sorted.pack_leaves();

    DoubleQRewards rewards(g0.n, g1.n);
    if(arguments.initialize_rewards){
        rewards.initialize(g0_deg, g1_deg);
    }

    // start clock
    stats->start = clock();

    vector<VtxPair> solution = mcs(g0_sorted, g1_sorted, (void *) &rewards, stats);

    // Convert to indices from original, unsorted graphs
    for (auto &vtx_pair: solution) {
        vtx_pair.v = vv0[vtx_pair.v];
        vtx_pair.w = vv1[vtx_pair.w];
    }

    // auto stop = std::chrono::steady_clock::now();
    // auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    time_elapsed = clock() - stats->start;
    clock_t time_find = stats->bestfind - stats->start;
    /* Clean up the timeout thread */
#if 1
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            stats->abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }
#endif

    if (!check_sol(g0, g1, solution))
        cout << "*** Error: Invalid solution" << endl;

    cout << "Solution size " << solution.size() << std::endl;
    for (int i = 0; i < g0.n; i++)
        for (unsigned int j = 0; j < solution.size(); j++)
            if (solution[j].v == i)
                cout << "(" << solution[j].v << " -> " << solution[j].w << ") ";
    cout << std::endl;

    cout << "Arguments:" << endl;
    cout << "  -t:                      " << arguments.timeout << endl;
    cout << "  -random_start:           " << arguments.random_start << endl;
    cout << "  -sort_heuristic:         " << arguments.sort_heuristic->name() << endl;
    cout << "  -initialize_reward:      " << arguments.initialize_rewards << endl;
    cout << "  -mcs_method:             " << arguments.mcs_method << endl;
    cout << "  -swap_policy:            " << arguments.swap_policy << endl;
    cout << "  -current_reward_policy:  " << arguments.reward_policy.current_reward_policy << endl;
    cout << "  -reward_policies_num:    " << arguments.reward_policy.reward_policies_num << endl;
    cout << "  -switch_policy:          " << arguments.reward_policy.switch_policy << endl;
    cout << "  -dal_reward_policy:      " << arguments.reward_policy.dal_reward_policy << endl;
    cout << "  -neighbor_overlap:       " << arguments.reward_policy.neighbor_overlap << endl;
    cout << "  -until_pruning:          " << arguments.first_thread_goes_until_pruning << endl;
    cout << "  -threads:                " << arguments.threads << endl;
    cout << "  -max_thread_blocks:      " << arguments.max_thread_blocks << endl;
    cout << "  -max_iter:               " << arguments.max_iter << endl;

    cout << endl;

    cout << "Nodes:                      " << stats->nodes << endl;
    cout << "Cut branches:               " << stats->cutbranches << endl;
    cout << "Conflicts:                  " << stats->conflicts << endl;
    printf("CPU time (ms):               %15ld\n", time_elapsed * 1000 / CLOCKS_PER_SEC);
    printf("FindBest time (ms):          %15ld\n", time_find * 1000 / CLOCKS_PER_SEC);
#ifdef Best
    cout << "Best nodes:                 " << stats->bestnodes << endl;
    cout << "Best count:                 " << stats->bestcount << endl;
    cout << "Swapped:                    " << stats->swapped_graphs << endl;
#endif
    if (aborted)
        cout << "TIMEOUT" << endl;

    delete arguments.sort_heuristic;

    return 0;
}
