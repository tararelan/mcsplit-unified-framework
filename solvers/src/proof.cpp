#include "proof.h"
#include <stdexcept>

static std::string match_var(int v, int w)
{
	return "x" + std::to_string(v) + "_" + std::to_string(w);
}

static std::string bottom_var(int v)
{
	return "x" + std::to_string(v) + "_b";
}

ProofLogger::ProofLogger(const std::string &prefix) : opb(prefix + ".opb"), pbp(prefix + ".pbp")
{
	if (!opb || !pbp)
		throw std::runtime_error("Cannot open proof files with prefix " + prefix);
}

void ProofLogger::write_model(const Graph &g, const Graph &h)
{
	g_n = g.n;
	h_n = h.n;
	compatible.assign(g_n, std::vector<bool>(h_n, false));
	for (int v = 0; v < g_n; v++)
		for (int w = 0; w < h_n; w++)
			compatible[v][w] = g.label[v] == h.label[w];

	// count variables and constraints up front for the OPB header
	long long n_vars = g_n, n_constraints = 2LL * g_n;
	std::vector<int> users(h_n, 0);
	for (int v = 0; v < g_n; v++)
		for (int w = 0; w < h_n; w++)
			if (compatible[v][w])
			{
				n_vars++;
				users[w]++;
			}
	for (int w = 0; w < h_n; w++)
		if (users[w] > 0)
			n_constraints++;

	// two assignments conflict iff their adjacency differs; adjmat values
	// encode presence, direction and edge label, so compare them exactly
	auto conflicts = [&](int v1, int w1, int v2, int w2) {
		return w1 != w2 && compatible[v1][w1] && compatible[v2][w2] && g.adjmat[v1][v2] != h.adjmat[w1][w2];
	};
	for (int v1 = 0; v1 < g_n; v1++)
		for (int v2 = v1 + 1; v2 < g_n; v2++)
			for (int w1 = 0; w1 < h_n; w1++)
				for (int w2 = 0; w2 < h_n; w2++)
					if (conflicts(v1, w1, v2, w2))
						n_constraints++;

	opb << "* #variable= " << n_vars << " #constraint= " << n_constraints << "\n";
	opb << "min:";
	for (int v = 0; v < g_n; v++)
		opb << " 1 " << bottom_var(v);
	opb << " ;\n";

	// IDs 1..g_n: each v takes at least one value
	for (int v = 0; v < g_n; v++)
	{
		opb << "1 " << bottom_var(v);
		for (int w = 0; w < h_n; w++)
			if (compatible[v][w])
				opb << " 1 " << match_var(v, w);
		opb << " >= 1 ;\n";
	}
	// IDs g_n+1..2g_n: each v takes at most one value
	for (int v = 0; v < g_n; v++)
	{
		opb << "1 " << bottom_var(v);
		for (int w = 0; w < h_n; w++)
			if (compatible[v][w])
				opb << " 1 " << match_var(v, w);
		opb << " <= 1 ;\n";
	}
	next_id = 2LL * g_n + 1;

	// each w is the image of at most one v
	inj_id.assign(h_n, -1);
	for (int w = 0; w < h_n; w++)
	{
		if (users[w] == 0)
			continue;
		bool first = true;
		for (int v = 0; v < g_n; v++)
			if (compatible[v][w])
			{
				opb << (first ? "" : " ") << "1 " << match_var(v, w);
				first = false;
			}
		opb << " <= 1 ;\n";
		inj_id[w] = next_id++;
	}

	// adjacency must be preserved in both directions (induced subgraph)
	for (int v1 = 0; v1 < g_n; v1++)
		for (int v2 = v1 + 1; v2 < g_n; v2++)
			for (int w1 = 0; w1 < h_n; w1++)
				for (int w2 = 0; w2 < h_n; w2++)
					if (conflicts(v1, w1, v2, w2))
					{
						opb << "1 ~" << match_var(v1, w1) << " 1 ~" << match_var(v2, w2) << " >= 1 ;\n";
						next_id++;
					}
	opb.close();

	pbp << "pseudo-Boolean proof version 3.0\n";
	pbp << "f " << n_constraints << " ;\n";

	// the solver's initial incumbent is the empty mapping
	pbp << "soli";
	for (int v = 0; v < g_n; v++)
		pbp << " " << bottom_var(v);
	pbp << " ;\n";
	obj_id = next_id++;
}

void ProofLogger::set_vertex_orders(const std::vector<int> &g_order, const std::vector<int> &h_order)
{
	this->g_order = g_order;
	this->h_order = h_order;
}

void ProofLogger::new_incumbent(const std::vector<VtxPair> &current)
{
	if (aborted)
		return;
	std::vector<int> image(g_n, -1);
	for (const VtxPair &p : current)
		image[g_order[p.v]] = h_order[p.w];
	pbp << "soli";
	for (int v = 0; v < g_n; v++)
		pbp << " " << (image[v] == -1 ? bottom_var(v) : match_var(v, image[v]));
	pbp << " ;\n";
	obj_id = next_id++;
}

void ProofLogger::push_match(int v, int w)
{
	trail.push_back(match_var(g_order[v], h_order[w]));
}

void ProofLogger::push_exclude(int v)
{
	trail.push_back(bottom_var(g_order[v]));
}

void ProofLogger::pop()
{
	trail.pop_back();
}

void ProofLogger::write_backtrack_clause()
{
	pbp << "rup";
	for (const std::string &lit : trail)
		pbp << " 1 ~" << lit;
	pbp << " >= 1 ;\n";
	long long id = next_id++;
	if (trail.empty())
		contra_id = id;

	// the clauses learned for this node's children are subsumed by this one
	size_t depth = trail.size();
	if (child_nogoods.size() < depth + 2)
		child_nogoods.resize(depth + 2);
	if (!child_nogoods[depth + 1].empty())
	{
		pbp << "del id";
		for (long long c : child_nogoods[depth + 1])
			pbp << " " << c;
		pbp << " ;\n";
		child_nogoods[depth + 1].clear();
	}
	child_nogoods[depth].push_back(id);
}

void ProofLogger::prune_by_bound(const std::vector<Bidomain> &domains,
								 const std::vector<int> &left, const std::vector<int> &right)
{
	if (aborted)
		return;

	// For a bidomain (L, R) with |L| > |R|, summing "v takes some value" over
	// L and "w is used at most once" over R gives sum_{v in L} x<v>_b >= |L| - |R|
	// once the trail falsifies every assignment leaving the bidomain. Adding
	// these to the objective bound leaves a constraint that unit propagation
	// refutes under the trail. Smaller bidomains, matched vertices and
	// vertices outside every bidomain are handled by propagation alone.
	bool any = false;
	for (const Bidomain &bd : domains)
	{
		if (bd.left_len <= bd.right_len)
			continue;
		if (!any)
		{
			pbp << "pol " << obj_id;
			any = true;
		}
		for (int i = 0; i < bd.left_len; i++)
			pbp << " " << alo_id(g_order[left[bd.l + i]]) << " +";
		for (int i = 0; i < bd.right_len; i++)
			pbp << " " << inj_id[h_order[right[bd.r + i]]] << " +";
	}
	long long pol_id = -1;
	if (any)
	{
		pbp << " ;\n";
		pol_id = next_id++;
	}

	write_backtrack_clause();

	if (any)
		pbp << "del id " << pol_id << " ;\n";
}

void ProofLogger::backtrack()
{
	if (aborted)
		return;
	write_backtrack_clause();
}

void ProofLogger::abort()
{
	aborted = true;
}

void ProofLogger::finish(int best_size)
{
	int best_obj = g_n - best_size;
	long long lb_id = contra_id;
	int lb = best_obj;
	if (aborted || contra_id == -1)
	{
		// search incomplete: only the best logged solution is certified, with
		// the trivial lower bound sum_v x<v>_b >= 0
		pbp << "rup";
		for (int v = 0; v < g_n; v++)
			pbp << " 1 " << bottom_var(v);
		pbp << " >= 0 ;\n";
		lb_id = next_id++;
		lb = 0;
	}
	pbp << "output NONE ;\n";
	pbp << "conclusion BOUNDS " << lb << " : " << lb_id << " " << best_obj << " ;\n";
	pbp << "end pseudo-Boolean proof ;\n";
	pbp.close();
}
