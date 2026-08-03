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

std::vector<VtxPair> mcs(const Graph &g, const Graph &h, bool multiway,
						 Stats &stats, std::atomic<bool> &abort_due_to_timeout);
std::vector<VtxPair> mcs_rl(const Graph &g, const Graph &h, bool multiway,
							Stats &stats, std::atomic<bool> &abort_due_to_timeout);
// use_lum: enable leaf union matching; use_lsm: enable long-short memory scoring
std::vector<VtxPair> mcs_ll(const Graph &g, const Graph &h, bool multiway, Stats &stats,
							std::atomic<bool> &abort_due_to_timeout, bool use_lum = true, bool use_lsm = true);
// policy_mode: 0 = hybrid (alternating), 1 = RL only, 2 = DAL only
std::vector<VtxPair> mcs_dal(const Graph &g, const Graph &h, bool multiway,
							 Stats &stats, std::atomic<bool> &abort_due_to_timeout,
							 int policy_mode = 0);
// bound_mode: 0 = gated (adaptive), 1 = always on, 2 = never on (plain McSplit bound)
std::vector<VtxPair> mcs_dsb(const Graph &g, const Graph &h, bool multiway,
							 Stats &stats, std::atomic<bool> &abort_due_to_timeout, int bound_mode = 0);
// red_mode: bitmask - bit 0 = vertex-equivalence, bit 1 = maximality, bit 2 = tighter bound
std::vector<VtxPair> mcs_rr(const Graph &g, const Graph &h, Stats &stats,
							std::atomic<bool> &abort_due_to_timeout, int red_mode = 7);
// sym_mode: bitmask - bit 0 = variable symmetry, bit 1 = value symmetry
std::vector<VtxPair> mcs_sym(const Graph &g, const Graph &h,
							 Stats &stats, std::atomic<bool> &abort_due_to_timeout, int sym_mode = 3);

static char doc[] = "Unified MCS solver";
static char args_doc[] = "FILENAME1 FILENAME2";

// Command-line options accepted by the solver.
// Graph format defaults to MIVIA binary ('B') if neither --dimacs nor --lad is specified.
// Timeout of 0 means no limit. Algorithm names correspond to the mcs_* function variants.
// -A accepts more values than shown below (ablation sub-variants like ll_lsm,
// dal_rl, dsb_always, rrsplit_noveq, symsplit_varonly, etc.) - see the dispatch
// logic later in this file, or the User Manual appendix, for the full list.
static struct argp_option options[] = {
	{"quiet", 'q', 0, 0, "Suppress per-instance output; print only the solution size"},
	{"dimacs", 'd', 0, 0, "Read graphs in DIMACS format (default: MIVIA binary)"},
	{"lad", 'l', 0, 0, "Read graphs in LAD format (default: MIVIA binary)"},
	{"directed", 'i', 0, 0, "Treat graphs as directed"},
	{"labelled", 'a', 0, 0, "Use both vertex and edge labels"},
	{"vertex-labelled", 'x', 0, 0, "Use vertex labels only (no edge labels)"},
	{"timeout", 't', "timeout", 0, "Abort search after this many seconds (0 = no limit)"},
	{"algorithm", 'A', "algo", 0, "Algorithm to run: mcsplit, rl, ll, dal, dsb, rrsplit, symsplit"},
	{0}};

// Parsed argument state, populated by the argp parser and read by main().
static struct
{
	bool quiet = false;				   // suppress verbose output
	bool dimacs = false;			   // use DIMACS reader
	bool lad = false;				   // use LAD reader
	bool directed = false;			   // treat input as directed graphs
	bool edge_labelled = false;		   // read and match edge labels
	bool vertex_labelled = false;	   // read and match vertex labels
	int timeout = 0;				   // search timeout in seconds (0 = unlimited)
	std::string algorithm = "mcsplit"; // algorithm identifier
	char *filename1 = nullptr;		   // path to first input graph
	char *filename2 = nullptr;		   // path to second input graph
	int arg_num = 0;				   // counter for positional arguments
} arguments;

/**
 * Argp callback that populates the arguments struct from parsed command-line input.
 *
 * Handles both option keys (single-character flags and their arguments) and
 * positional arguments (the two input graph filenames). Calls argp_usage() to
 * print usage and exit if more than two positional arguments are supplied or
 * if fewer than two are present at ARGP_KEY_END.
 *
 * @param key    Option character or ARGP_KEY_* control value
 * @param arg    Option argument string, or NULL if the option takes no argument
 * @param state  Argp parser state (used for argp_usage())
 * @return       0 on success, ARGP_ERR_UNKNOWN for unrecognised keys
 */
static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	switch (key)
	{
	case 'q':
		arguments.quiet = true;
		break;
	case 'd':
		arguments.dimacs = true;
		break;
	case 'l':
		arguments.lad = true;
		break;
	case 'i':
		arguments.directed = true;
		break;
	case 'a':
		// -a/--labelled sets both flags: "labelled" means both vertex
		// and edge labels are used together, unlike -x which is vertex-only
		arguments.edge_labelled = true;
		arguments.vertex_labelled = true;
		break;
	case 'x':
		arguments.vertex_labelled = true;
		break;
	case 't':
		arguments.timeout = std::stoi(arg);
		break;
	case 'A':
		arguments.algorithm = arg;
		break;
	case ARGP_KEY_ARG:
		if (arguments.arg_num == 0)
			arguments.filename1 = arg;
		else if (arguments.arg_num == 1)
			arguments.filename2 = arg;
		else
			argp_usage(state);
		arguments.arg_num++;
		break;
	case ARGP_KEY_END:
		if (arguments.arg_num < 2)
			argp_usage(state);
		break;
	default:
		return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

/**
 * Entry point. Parses command-line arguments, reads both input graphs,
 * runs the selected MCS algorithm, and prints results to stdout as a
 * single CSV line.
 *
 * Output format:
 *   instance_a, instance_b, algo, size, edges, nodes, time,
 *   aborted, nodes_to_best, time_to_best,
 *   cut_branches, bound_pruned, sym_pruned
 *
 * If --quiet is set, only the solution size is printed.
 * If the timeout fires, aborted=1 and the best incumbent found is reported.
 */
int main(int argc, char **argv)
{
	argp_parse(&argp, argc, argv, 0, 0, 0);

	if (!arguments.quiet)
		std::cerr << "# solution_size solution_edges nodes time_elapsed aborted nodes_to_best time_to_best cut_branches bound_pruned sym_pruned" << std::endl;

	char format = arguments.dimacs ? 'D' : arguments.lad ? 'L'
														 : 'B';
	Graph g = readGraph(arguments.filename1, format, arguments.directed,
						arguments.edge_labelled, arguments.vertex_labelled);
	Graph h = readGraph(arguments.filename2, format, arguments.directed,
						arguments.edge_labelled, arguments.vertex_labelled);

	// convention: g is always the smaller graph; several algorithms'
	// complexity bounds (see Section 2.13) assume |V(G)| <= |V(H)|
	if (g.n > h.n)
		std::swap(g, h);

	if (!arguments.quiet)
	{
		std::cout << g.n << " vertices" << std::endl;
		std::cout << h.n << " vertices" << std::endl;
	}

	std::atomic<bool> abort_due_to_timeout(false);
	std::thread timeout_thread;
	std::mutex timeout_mutex;
	std::condition_variable timeout_cv;
	bool aborted = false;

	// background watchdog thread: sleeps until either the timeout elapses
	// or the search finishes first (abort_due_to_timeout set elsewhere)
	if (arguments.timeout > 0)
	{
		timeout_thread = std::thread([&]
									 {
            auto abort_time = std::chrono::steady_clock::now() +
                    std::chrono::seconds(arguments.timeout);
            std::unique_lock<std::mutex> guard(timeout_mutex);
            while (!abort_due_to_timeout.load()) {
                if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                    aborted = true;
                    break;
                }
            }
            abort_due_to_timeout.store(true); });
	}

	bool multiway = arguments.directed || arguments.edge_labelled;

	auto start = std::chrono::steady_clock::now();
	Stats stats;
	std::vector<VtxPair> solution;

	// dispatch table: algorithm name -> mcs_* call with the matching
	// mode/bitmask (see the declarations and options[] table above for
	// what each mode value means)
	if (arguments.algorithm == "rl")
		solution = mcs_rl(g, h, multiway, stats, abort_due_to_timeout);

	else if (arguments.algorithm == "ll")
		solution = mcs_ll(g, h, multiway, stats, abort_due_to_timeout, true, true);
	else if (arguments.algorithm == "ll_lsm")
		solution = mcs_ll(g, h, multiway, stats, abort_due_to_timeout, false, true);
	else if (arguments.algorithm == "ll_lum")
		solution = mcs_ll(g, h, multiway, stats, abort_due_to_timeout, true, false);
	else if (arguments.algorithm == "dal")
		solution = mcs_dal(g, h, multiway, stats, abort_due_to_timeout, 0);
	else if (arguments.algorithm == "dal_rl")
		solution = mcs_dal(g, h, multiway, stats, abort_due_to_timeout, 1);
	else if (arguments.algorithm == "dal_dal")
		solution = mcs_dal(g, h, multiway, stats, abort_due_to_timeout, 2);

	else if (arguments.algorithm == "dsb")
		solution = mcs_dsb(g, h, multiway, stats, abort_due_to_timeout);
	else if (arguments.algorithm == "dsb_always")
		solution = mcs_dsb(g, h, multiway, stats, abort_due_to_timeout, 1);
	else if (arguments.algorithm == "dsb_never")
		solution = mcs_dsb(g, h, multiway, stats, abort_due_to_timeout, 2);

	else if (arguments.algorithm == "rrsplit")
		solution = mcs_rr(g, h, stats, abort_due_to_timeout, 7); // all on
	else if (arguments.algorithm == "rrsplit_noveq")
		solution = mcs_rr(g, h, stats, abort_due_to_timeout, 6); // max+bound
	else if (arguments.algorithm == "rrsplit_nomax")
		solution = mcs_rr(g, h, stats, abort_due_to_timeout, 5); // veq+bound
	else if (arguments.algorithm == "rrsplit_nobound")
		solution = mcs_rr(g, h, stats, abort_due_to_timeout, 3); // veq+max

	else if (arguments.algorithm == "symsplit")
		solution = mcs_sym(g, h, stats, abort_due_to_timeout, 3); // full
	else if (arguments.algorithm == "symsplit_valonly")
		solution = mcs_sym(g, h, stats, abort_due_to_timeout, 2); // value only
	else if (arguments.algorithm == "symsplit_varonly")
		solution = mcs_sym(g, h, stats, abort_due_to_timeout, 1); // variable only

	else
		solution = mcs(g, h, multiway, stats, abort_due_to_timeout);

	auto stop = std::chrono::steady_clock::now();
	stats.time_elapsed = std::chrono::duration<double>(stop - start).count();
	stats.aborted = aborted;

	// Count edges in solution (boolean adjacency check only - does not
	// distinguish direction on directed graphs, see Section 4.1.2.2)
	int solution_edges = 0;
	for (int i = 0; i < (int)solution.size(); i++)
		for (int j = i + 1; j < (int)solution.size(); j++)
			if (g.adjmat[solution[i].v][solution[j].v])
				solution_edges++;

	// Clean up timeout thread
	if (timeout_thread.joinable())
	{
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
			  << stats.nodes_to_best << " "
			  << stats.time_to_best << " "
			  << stats.cut_branches << " "
			  << stats.bound_pruned << " "
			  << stats.sym_pruned << " ";

	return 0;
}