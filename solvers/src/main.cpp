// #define _GNU_SOURCE

#include "graph.h"
#include "common.h"
#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <argp.h>

// Forward declarations
std::vector<VtxPair> mcs(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_dal(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_rr(const Graph& g, const Graph& h,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_sym(const Graph& g, const Graph& h,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_ll(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_dsb(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_rl(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout);
std::vector<VtxPair> mcs_rllum(const Graph& g, const Graph& h, bool multiway,
        Stats& stats, std::atomic<bool>& abort_due_to_timeout);

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Unified MCS solver";
static char args_doc[] = "FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"quiet",           'q', 0,         0, "Quiet output"},
    {"dimacs",          'd', 0,         0, "Read DIMACS format"},
    {"lad",             'l', 0,         0, "Read LAD format"},
    {"directed",        'i', 0,         0, "Use directed graphs"},
    {"labelled",        'a', 0,         0, "Use edge and vertex labels"},
    {"vertex-labelled", 'x', 0,         0, "Use vertex labels only"},
    {"timeout",         't', "timeout", 0, "Timeout in seconds"},
    {"algorithm",       'A', "algo",    0, "Algorithm: mcsplit (default), dal"},
    { 0 }
};

static struct {
    bool quiet = false;
    bool dimacs = false;
    bool lad = false;
    bool directed = false;
    bool edge_labelled = false;
    bool vertex_labelled = false;
    int timeout = 0;
    std::string algorithm = "mcsplit";
    char* filename1 = nullptr;
    char* filename2 = nullptr;
    int arg_num = 0;
} arguments;

static error_t parse_opt(int key, char* arg, struct argp_state* state) {
    switch (key) {
        case 'q': arguments.quiet = true; break;
        case 'd': arguments.dimacs = true; break;
        case 'l': arguments.lad = true; break;
        case 'i': arguments.directed = true; break;
        case 'a': arguments.edge_labelled = true; arguments.vertex_labelled = true; break;
        case 'x': arguments.vertex_labelled = true; break;
        case 't': arguments.timeout = std::stoi(arg); break;
        case 'A': arguments.algorithm = arg; break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) arguments.filename1 = arg;
            else if (arguments.arg_num == 1) arguments.filename2 = arg;
            else argp_usage(state);
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num < 2) argp_usage(state);
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/*******************************************************************************
                                    Main
*******************************************************************************/

int main(int argc, char** argv) {
    argp_parse(&argp, argc, argv, 0, 0, 0);

    if (!arguments.quiet)
        std::cerr << "# solution_size solution_edges nodes time_elapsed aborted root_ub nodes_to_best time_to_best cut_branches bound_pruned sym_pruned conflicts" << std::endl;

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    Graph g = readGraph(arguments.filename1, format, arguments.directed,
            arguments.edge_labelled, arguments.vertex_labelled);
    Graph h = readGraph(arguments.filename2, format, arguments.directed,
            arguments.edge_labelled, arguments.vertex_labelled);

    if (g.n > h.n) std::swap(g, h);

    if (!arguments.quiet) {
        std::cout << g.n << " vertices" << std::endl;
        std::cout << h.n << " vertices" << std::endl;
    }

    std::atomic<bool> abort_due_to_timeout(false);
    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    bool aborted = false;

    if (arguments.timeout > 0) {
        timeout_thread = std::thread([&] {
            auto abort_time = std::chrono::steady_clock::now() +
                    std::chrono::seconds(arguments.timeout);
            std::unique_lock<std::mutex> guard(timeout_mutex);
            while (!abort_due_to_timeout.load()) {
                if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                    aborted = true;
                    break;
                }
            }
            abort_due_to_timeout.store(true);
        });
    }

    bool multiway = arguments.directed || arguments.edge_labelled;

    auto start = std::chrono::steady_clock::now();
    Stats stats;
    std::vector<VtxPair> solution;

    if (arguments.algorithm == "dal")
        solution = mcs_dal(g, h, multiway, stats, abort_due_to_timeout);
	else if (arguments.algorithm == "rrsplit")
		solution = mcs_rr(g, h, stats, abort_due_to_timeout);
	else if (arguments.algorithm == "symsplit")
		solution = mcs_sym(g, h, stats, abort_due_to_timeout);
	else if (arguments.algorithm == "ll")
		solution = mcs_ll(g, h, multiway, stats, abort_due_to_timeout);
	else if (arguments.algorithm == "dsb")
		solution = mcs_dsb(g, h, multiway, stats, abort_due_to_timeout);
	else if (arguments.algorithm == "rl")
		solution = mcs_rl(g, h, multiway, stats, abort_due_to_timeout);
	else if (arguments.algorithm == "rllum")
		solution = mcs_rllum(g, h, multiway, stats, abort_due_to_timeout);
    else
        solution = mcs(g, h, multiway, stats, abort_due_to_timeout);

    auto stop = std::chrono::steady_clock::now();
    stats.time_elapsed = std::chrono::duration<double>(stop - start).count();
    stats.aborted = aborted;

    // Count edges in solution
    int solution_edges = 0;
    for (int i = 0; i < (int)solution.size(); i++)
        for (int j = i + 1; j < (int)solution.size(); j++)
            if (g.adjmat[solution[i].v][solution[j].v])
                solution_edges++;

    // Clean up timeout thread
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }

    std::cout << solution.size() << " "
              << solution_edges << " "
              << stats.nodes << " "
              << stats.time_elapsed << " "
              << stats.aborted << " "
              << stats.root_upper_bound << " "
              << stats.nodes_to_best << " "
              << stats.time_to_best << " "
              << stats.cut_branches << " "
              << stats.bound_pruned << " "
              << stats.sym_pruned << " "
              << stats.conflicts << std::endl;

    return 0;
}