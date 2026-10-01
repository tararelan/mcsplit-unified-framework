#include "proof.h"
#include <algorithm>
#include <functional>
#include <stdexcept>

static std::string match_var(int v, int w)
{
	return "x" + std::to_string(v) + "_" + std::to_string(w);
}

static std::string bottom_var(int v)
{
	return "x" + std::to_string(v) + "_b";
}

static std::string unused_var(int w)
{
	return "y" + std::to_string(w);
}

// Whether x and y are related inside a bidomain side: adjacent, or with
// complement set, distinct and non-adjacent. Like the DSB bound, this
// ignores direction and edge labels.
static bool related(const std::vector<std::vector<unsigned int>> &adj, int x, int y, bool complement)
{
	return x != y && (complement ? adj[x][y] == 0 : adj[x][y] != 0);
}

// Degrees within vs of the relation, and the number of related pairs.
static std::vector<int> side_degrees(const std::vector<std::vector<unsigned int>> &adj,
									 const std::vector<int> &vs, bool complement, int &pairs)
{
	std::vector<int> deg(vs.size(), 0);
	int total = 0;
	for (size_t i = 0; i < vs.size(); i++)
		for (size_t j = 0; j < vs.size(); j++)
			if (related(adj, vs[i], vs[j], complement))
				deg[i]++;
	for (int d : deg)
		total += d;
	pairs = total / 2;
	return deg;
}

// The fewest vertices of S whose degrees sum to at least diff (0 if diff <= 0).
static int degree_cover(std::vector<int> deg, int diff)
{
	if (diff <= 0)
		return 0;
	std::sort(deg.begin(), deg.end(), std::greater<int>());
	int r = 0, acc = 0;
	while (acc < diff)
		acc += deg[r++];
	return r;
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
	g_adj = g.adjmat;
	h_adj = h.adjmat;
	users.assign(h_n, 0);
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
	h_rank.assign(h_n, 0);
	for (int i = 0; i < h_n; i++)
		h_rank[h_order[i]] = i;
	g_image.assign(g_n, -1);
	h_used.assign(h_n, false);
	g_excluded.assign(g_n, false);
}

void ProofLogger::enable_unused_vertex_vars()
{
	unused_alo_id.assign(h_n, -1);
	unused_amo_id.assign(h_n, -1);
	for (int w = 0; w < h_n; w++)
	{
		int k = users[w];
		if (k == 0)
			continue;
		// y<w> => no v is mapped to w
		pbp << "red " << k << " ~" << unused_var(w);
		for (int v = 0; v < g_n; v++)
			if (compatible[v][w])
				pbp << " 1 ~" << match_var(v, w);
		pbp << " >= " << k << " : " << unused_var(w) << " -> 0 ;\n";
		long long fwd = next_id++;
		// no v is mapped to w => y<w>
		pbp << "red 1 " << unused_var(w);
		for (int v = 0; v < g_n; v++)
			if (compatible[v][w])
				pbp << " 1 " << match_var(v, w);
		pbp << " >= 1 : " << unused_var(w) << " -> 1 ;\n";
		unused_alo_id[w] = next_id++;
		// with injectivity, y<w> + sum_v x<v>_<w> <= 1
		if (k == 1)
			unused_amo_id[w] = fwd;
		else
		{
			pbp << "pol " << fwd << " " << inj_id[w] << " " << k - 1 << " * + " << k << " d ;\n";
			unused_amo_id[w] = next_id++;
		}
	}
}

void ProofLogger::break_g_twins(const std::vector<std::vector<int>> &classes)
{
	// pos(a) is the sorted index of a's image plus one, or h_n + 1 if a is
	// unmapped. For twins a < b, swapping them is an automorphism of the
	// model, so any solution with pos(a) > pos(b) can be swapped into one
	// with pos(a) < pos(b): red O(a,b): pos(b) - pos(a) >= 0. Deriving
	// O(a,b) for later b creates goals for the earlier O(a,c), which become
	// pos(c) - pos(b) >= 0; they follow from O(a,c) and the negation of
	// O(a,b) by adding them up.
	auto write_pos = [&](int a, int sign) {
		for (int w = 0; w < h_n; w++)
			if (compatible[a][w])
				pbp << " " << sign * (h_rank[w] + 1) << " " << match_var(a, w);
		pbp << " " << sign * (h_n + 1) << " " << bottom_var(a);
	};

	g_twin_class.assign(g_n, -1);
	for (const std::vector<int> &sorted_class : classes)
	{
		std::vector<int> cls;
		for (int v : sorted_class)
			cls.push_back(g_order[v]);
		if (cls.size() > 1)
		{
			for (int v : cls)
				g_twin_class[v] = g_twin_classes.size();
			g_twin_classes.push_back(cls);
		}
		for (size_t i = 0; i + 1 < cls.size(); i++)
			for (size_t k = i + 1; k < cls.size(); k++)
				if (compatible[cls[i]] != compatible[cls[k]])
					throw std::logic_error("proof logging: twins with different vertex labels");

		std::vector<std::vector<long long>> order_id(cls.size(), std::vector<long long>(cls.size(), -1));
		for (size_t i = 0; i + 1 < cls.size(); i++)
			for (size_t k = i + 1; k < cls.size(); k++)
			{
				const int a = cls[i], b = cls[k];
				pbp << "red";
				write_pos(b, 1);
				write_pos(a, -1);
				pbp << " >= 0 :";
				for (int w = 0; w < h_n; w++)
					if (compatible[a][w])
						pbp << " " << match_var(a, w) << " -> " << match_var(b, w) << " "
							<< match_var(b, w) << " -> " << match_var(a, w);
				pbp << " " << bottom_var(a) << " -> " << bottom_var(b) << " "
					<< bottom_var(b) << " -> " << bottom_var(a);
				if (k == i + 1)
				{
					pbp << " ;\n";
					next_id++;
				}
				else
				{
					pbp << " : subproof\n";
					long long neg = next_id++;
					for (size_t m = i + 1; m < k; m++)
					{
						pbp << "proofgoal " << order_id[i][m] << "\n";
						long long goal = next_id++;
						pbp << "pol " << neg << " " << goal << " + " << order_id[i][m] << " + ;\n";
						next_id++;
						pbp << "qed : -1 ;\n";
					}
					pbp << "qed ;\n";
					next_id++;
				}
				order_id[i][k] = next_id - 1;
			}

		// For propagation, between consecutive twins a < b: if b takes w,
		// then a takes something earlier than w. With unit propagation, this
		// chains along the class, and an unmapped twin unmaps all later ones.
		for (size_t i = 0; i + 1 < cls.size(); i++)
		{
			const int a = cls[i], b = cls[i + 1];
			for (int w = 0; w < h_n; w++)
			{
				if (!compatible[b][w])
					continue;
				pbp << "rup 1 ~" << match_var(b, w);
				for (int u = 0; u < h_n; u++)
					if (compatible[a][u] && h_rank[u] < h_rank[w])
						pbp << " 1 " << match_var(a, u);
				pbp << " >= 1 : " << order_id[i][i + 1] << " " << amo_id(b) << " " << alo_id(a) << " " << inj_id[w] << " ;\n";
				next_id++;
			}
		}

		pbp << "del id";
		for (size_t i = 0; i + 1 < cls.size(); i++)
			for (size_t k = i + 1; k < cls.size(); k++)
				pbp << " " << order_id[i][k];
		pbp << " ;\n";
	}
}

void ProofLogger::prune_by_twin_bound(const std::vector<Bidomain> &domains, const std::vector<int> &left,
									  const std::vector<int> &right, size_t bd_idx,
									  const std::vector<int> &twins, const std::vector<int> &high)
{
	if (aborted)
		return;

	// The twins can only take values from high (unit propagation over the
	// twin ordering rules out the rest), so summing their at-least-one
	// constraints with high's injectivity constraints gives
	// sum_{twins} x<v>_b >= |twins| - |high|. The other bidomains are
	// handled as in prune_by_bound.
	pbp << "pol " << obj_id;
	for (int v : twins)
		pbp << " " << alo_id(g_order[v]) << " +";
	for (int w : high)
		pbp << " " << inj_id[h_order[w]] << " +";
	for (size_t i = 0; i < domains.size(); i++)
	{
		const Bidomain &bd = domains[i];
		if (i == bd_idx || bd.left_len <= bd.right_len)
			continue;
		for (int k = 0; k < bd.left_len; k++)
			pbp << " " << alo_id(g_order[left[bd.l + k]]) << " +";
		for (int k = 0; k < bd.right_len; k++)
			pbp << " " << inj_id[h_order[right[bd.r + k]]] << " +";
	}
	pbp << " ;\n";
	long long pol_id = next_id++;
	write_backtrack_clause();
	pbp << "del id " << pol_id << " ;\n";
}

void ProofLogger::new_incumbent(const std::vector<VtxPair> &current)
{
	if (aborted)
		return;
	incumbent_size = current.size();
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
	v = g_order[v];
	w = h_order[w];
	trail.push_back(match_var(v, w));
	trail_is_match.push_back(true);
	trail_vertex.push_back(v);
	trail_lum.emplace_back();
	trail_attached_ids.emplace_back();
	g_image[v] = w;
	h_used[w] = true;
	n_matched++;
}

void ProofLogger::push_exclude(int v)
{
	v = g_order[v];
	trail.push_back(bottom_var(v));
	trail_is_match.push_back(false);
	trail_vertex.push_back(v);
	trail_lum.emplace_back();
	trail_attached_ids.emplace_back();
	g_excluded[v] = true;
}

void ProofLogger::pop()
{
	for (int vl : trail_lum.back())
	{
		h_used[g_image[vl]] = false;
		g_image[vl] = -1;
		n_matched--;
	}
	if (!trail_attached_ids.back().empty())
	{
		pbp << "del id";
		for (long long id : trail_attached_ids.back())
			pbp << " " << id;
		pbp << " ;\n";
	}

	int v = trail_vertex.back();
	if (trail_is_match.back())
	{
		h_used[g_image[v]] = false;
		g_image[v] = -1;
		n_matched--;
	}
	else
		g_excluded[v] = false;

	trail.pop_back();
	trail_is_match.pop_back();
	trail_vertex.pop_back();
	trail_lum.pop_back();
	trail_attached_ids.pop_back();
}

long long ProofLogger::derive_forced_match(int v, int w)
{
	// Every vertex not on the trail is free. Under the trail, v can still
	// take any free z with its adjacency pattern (or bottom), and w can
	// still be taken by any free a with w's pattern.
	auto same_pattern = [&](int a, int b) {
		for (int u = 0; u < g_n; u++)
			if (g_image[u] != -1 && g_adj[a][u] != h_adj[b][g_image[u]])
				return false;
		return true;
	};
	std::vector<int> zs{-1}, as;
	for (int z = 0; z < h_n; z++)
		if (z != w && compatible[v][z] && !h_used[z] && same_pattern(v, z))
			zs.push_back(z);
	for (int a = 0; a < g_n; a++)
		if (a != v && compatible[a][w] && g_image[a] == -1 && !g_excluded[a] && same_pattern(a, w))
			as.push_back(a);

	// Rule out each z in turn, one holder a of w at a time: the witness
	// moves v onto w and a onto z (or just v, if w is unused). The caller
	// is responsible for this preserving adjacency.
	auto value = [&](int u, int z) { return z == -1 ? bottom_var(u) : match_var(u, z); };
	std::vector<long long> ruled_out;
	for (int z : zs)
	{
		std::vector<long long> cases;
		for (int a : as)
		{
			pbp << "red";
			for (const std::string &lit : trail)
				pbp << " 1 ~" << lit;
			pbp << " 1 ~" << value(v, z) << " 1 ~" << match_var(a, w) << " >= 1 : "
				<< value(v, z) << " -> 0 " << match_var(v, w) << " -> 1 "
				<< match_var(a, w) << " -> 0 " << value(a, z) << " -> 1 ;\n";
			cases.push_back(next_id++);
		}
		pbp << "red";
		for (const std::string &lit : trail)
			pbp << " 1 ~" << lit;
		pbp << " 1 ~" << value(v, z);
		for (int a : as)
			pbp << " 1 " << match_var(a, w);
		pbp << " >= 1 : " << value(v, z) << " -> 0 " << match_var(v, w) << " -> 1 ;\n";
		cases.push_back(next_id++);

		pbp << "rup";
		for (const std::string &lit : trail)
			pbp << " 1 ~" << lit;
		pbp << " 1 ~" << value(v, z) << " >= 1 ;\n";
		ruled_out.push_back(next_id++);
		pbp << "del id";
		for (long long id : cases)
			pbp << " " << id;
		pbp << " ;\n";
	}

	pbp << "rup";
	for (const std::string &lit : trail)
		pbp << " 1 ~" << lit;
	pbp << " 1 " << match_var(v, w) << " >= 1 ;\n";
	long long forced = next_id++;
	pbp << "del id";
	for (long long id : ruled_out)
		pbp << " " << id;
	pbp << " ;\n";
	return forced;
}

std::vector<long long> &ProofLogger::attached_ids()
{
	return trail.empty() ? root_attached_ids : trail_attached_ids.back();
}

void ProofLogger::lum_match(int vl, int wl)
{
	vl = g_order[vl];
	wl = h_order[wl];
	if (g_image[vl] != -1 || g_excluded[vl] || h_used[wl] || !trail_is_match.back())
		throw std::logic_error("proof logging: LUM pair is not between free leaves under a match");

	// vl is adjacent only to v, and wl only to w (the last decision), so in
	// any solution under the trail where vl takes z, moving vl to wl (and
	// whoever had wl to z) gives a solution of the same size, or a larger
	// one if z is bottom and wl was unused
	if (!aborted)
		attached_ids().push_back(derive_forced_match(vl, wl));

	g_image[vl] = wl;
	h_used[wl] = true;
	n_matched++;
	trail_lum.back().push_back(vl);
}

void ProofLogger::maximal_match(int v, int w)
{
	v = g_order[v];
	w = h_order[w];
	if (g_image[v] != -1 || g_excluded[v] || h_used[w])
		throw std::logic_error("proof logging: maximal match is not between free vertices");

	// v->w splits no bidomain, so v and w have the same adjacency to every
	// pair that could still be matched. In any solution under the trail
	// where v takes z and a takes w, swapping to v->w and a->z preserves
	// every adjacency.
	if (!aborted)
		attached_ids().push_back(derive_forced_match(v, w));
}

void ProofLogger::maximal_match_with_twins(int v, int w, bool and_bottom)
{
	v = g_order[v];
	w = h_order[w];
	if (g_image[v] != -1 || g_excluded[v] || h_used[w])
		throw std::logic_error("proof logging: maximal match is not between free vertices");
	if (aborted)
		return;

	auto same_pattern = [&](int a, int b) {
		for (int u = 0; u < g_n; u++)
			if (g_image[u] != -1 && g_adj[a][u] != h_adj[b][g_image[u]])
				return false;
		return true;
	};
	std::vector<int> zs, as;
	for (int z = 0; z < h_n; z++)
		if (h_rank[z] > h_rank[w] && compatible[v][z] && !h_used[z] && same_pattern(v, z))
			zs.push_back(z);
	// bottom comes after every value in the twin ordering
	if (and_bottom)
		zs.push_back(-1);
	auto value = [&](int u, int z) { return z == -1 ? bottom_var(u) : match_var(u, z); };
	// v's own twins can't hold w when v takes a later value: its earlier
	// twins are decided, and the twin ordering puts its later twins after
	// v's value (which unit propagation sees), so leave them out
	for (int a = 0; a < g_n; a++)
		if (a != v && compatible[a][w] && g_image[a] == -1 && !g_excluded[a] && same_pattern(a, w) &&
			(g_twin_class[a] == -1 || g_twin_class[a] != g_twin_class[v]))
			as.push_back(a);

	auto write_guard = [&](int coeff) {
		for (const std::string &lit : trail)
			pbp << " " << coeff << " ~" << lit;
	};

	// v->w splits no bidomain, so in a solution under the trail where v
	// takes some z after w and a takes w, v->w and a->z preserves every
	// adjacency (as for maximal_match). v's twin order survives too: its
	// earlier twins are below w, its later twins above z. But a's twins
	// a < s1 < s2 < ... whose images lie between w and z would now come
	// before a, so the rewrite also shifts that chain along (a takes s1's
	// image, s1 takes s2's, and the last one takes z), which is the same
	// solution up to a permutation of twins. The cases are by the length r
	// of the chain, longest first: case r rules out "s_r is between w and
	// z", so that unit propagation provides "s_{r+1} isn't" for case r.
	// When z is bottom, everything after w is between w and z, and v's
	// later twins are all unmapped.
	for (int z : zs)
	{
		std::vector<int> between; // values strictly between w and z
		for (int c = 0; c < h_n; c++)
			if (h_rank[c] > h_rank[w] && (z == -1 || h_rank[c] < h_rank[z]))
				between.push_back(c);

		std::vector<long long> cases;
		for (int a : as)
		{
			// a's later twins (classes are in increasing index order, and
			// twins are branched in that order, so these are all free)
			std::vector<int> chain{a};
			if (g_twin_class[a] != -1)
			{
				const std::vector<int> &cls = g_twin_classes[g_twin_class[a]];
				for (size_t i = std::find(cls.begin(), cls.end(), a) - cls.begin() + 1; i < cls.size(); i++)
					if (g_image[cls[i]] == -1 && !g_excluded[cls[i]])
						chain.push_back(cls[i]);
			}
			for (size_t r = chain.size() - 1; r >= 1; r--)
			{
				const int last = chain[r];
				std::vector<int> range;
				for (int c : between)
					if (compatible[last][c])
						range.push_back(c);
				// s1 < ... < s_r all strictly between w and z needs r values
				// there, and otherwise unit propagation along the twin
				// ordering already rules this case out
				if (range.size() < r)
					continue;
				const int d = range.size();
				pbp << "red";
				write_guard(d);
				pbp << " " << d << " ~" << value(v, z) << " " << d << " ~" << match_var(a, w);
				for (int c : range)
					pbp << " 1 ~" << match_var(last, c);
				pbp << " >= " << d << " : " << value(v, z) << " -> 0 " << match_var(v, w) << " -> 1";
				for (size_t i = 0; i < r; i++)
				{
					for (int c = 0; c < h_n; c++)
						if (compatible[chain[i]][c])
							pbp << " " << match_var(chain[i], c) << " -> " << match_var(chain[i + 1], c);
					pbp << " " << bottom_var(chain[i]) << " -> " << bottom_var(chain[i + 1]);
				}
				for (int c = 0; c < h_n; c++)
					if (compatible[last][c])
						pbp << " " << match_var(last, c) << " -> " << (c == z ? 1 : 0);
				pbp << " " << bottom_var(last) << " -> " << (z == -1 ? 1 : 0) << " ;\n";
				cases.push_back(next_id++);
			}
			pbp << "red";
			write_guard(1);
			pbp << " 1 ~" << value(v, z) << " 1 ~" << match_var(a, w) << " >= 1 : "
				<< value(v, z) << " -> 0 " << match_var(v, w) << " -> 1 "
				<< match_var(a, w) << " -> 0 " << value(a, z) << " -> 1 ;\n";
			cases.push_back(next_id++);
		}
		pbp << "red";
		write_guard(1);
		pbp << " 1 ~" << value(v, z);
		for (int a : as)
			pbp << " 1 " << match_var(a, w);
		pbp << " >= 1 : " << value(v, z) << " -> 0 " << match_var(v, w) << " -> 1 ;\n";
		cases.push_back(next_id++);

		pbp << "rup";
		write_guard(1);
		pbp << " 1 ~" << value(v, z) << " >= 1 ;\n";
		attached_ids().push_back(next_id++);
		pbp << "del id";
		for (long long id : cases)
			pbp << " " << id;
		pbp << " ;\n";
	}
}

void ProofLogger::break_value_symmetry(int v, const std::vector<std::pair<int, int>> &twins)
{
	if (aborted)
		return;
	v = g_order[v];
	for (auto [b, w] : twins)
	{
		b = h_order[b];
		w = h_order[w];
		// b and w are twins in H, so swapping them is an automorphism of
		// the model, and maps a solution with v->w to one with v->b; b and
		// w are both free, so the trail is unaffected
		pbp << "red";
		for (const std::string &lit : trail)
			pbp << " 1 ~" << lit;
		pbp << " 1 ~" << match_var(v, w) << " >= 1 :";
		for (int u = 0; u < g_n; u++)
		{
			if (compatible[u][b] != compatible[u][w])
				throw std::logic_error("proof logging: twins with different vertex labels");
			if (compatible[u][b])
				pbp << " " << match_var(u, b) << " -> " << match_var(u, w) << " "
					<< match_var(u, w) << " -> " << match_var(u, b);
		}
		pbp << " ;\n";
		attached_ids().push_back(next_id++);
	}
}

void ProofLogger::write_trail_clause()
{
	for (const std::string &lit : trail)
		pbp << " 1 ~" << lit;
	pbp << " >= 1";
}

void ProofLogger::write_backtrack_clause()
{
	pbp << "rup";
	write_trail_clause();
	pbp << " ;\n";
	record_backtrack_clause(next_id++);
}

void ProofLogger::record_backtrack_clause(long long id)
{
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

ProofLogger::BidomainPlan ProofLogger::plan_bidomain(const std::vector<int> &L, const std::vector<int> &R) const
{
	BidomainPlan best;
	best.bound = std::min(L.size(), R.size());

	if (L.size() == 2 && R.size() == 2 && g_adj[L[0]][L[1]] != h_adj[R[0]][R[1]] && best.bound > 1)
	{
		best.kind = BidomainPlan::Pair;
		best.bound = 1;
	}

	for (bool s_is_left : {true, false})
	{
		if (!s_is_left && unused_alo_id.empty())
			continue;
		const std::vector<int> &S = s_is_left ? L : R, &B = s_is_left ? R : L;
		for (bool complement : {false, true})
		{
			int s_pairs, b_pairs;
			std::vector<int> deg = side_degrees(s_is_left ? g_adj : h_adj, S, complement, s_pairs);
			side_degrees(s_is_left ? h_adj : g_adj, B, complement, b_pairs);
			int r = degree_cover(deg, s_pairs - b_pairs);
			int bound = (int)S.size() - r;
			if (r > 0 && bound < best.bound)
			{
				best.kind = BidomainPlan::Degree;
				best.bound = bound;
				best.s_is_left = s_is_left;
				best.complement = complement;
				best.r = r;
			}
		}
	}
	return best;
}

long long ProofLogger::derive_pair_bound(const std::vector<int> &L, const std::vector<int> &R)
{
	// L = {a, b}, R = {c, d}, and ab, cd disagree, so a->c and a->d each
	// leave b nowhere to go
	pbp << "pbc 1 " << bottom_var(L[0]) << " 1 " << bottom_var(L[1]) << " >= 1 : subproof\n";
	next_id++;
	for (int c : R)
	{
		pbp << "rup 1 ~" << match_var(L[0], c) << " >= 1 ;\n";
		next_id++;
	}
	pbp << "rup >= 1 ;\n";
	next_id++;
	pbp << "qed ;\n";
	return next_id++;
}

long long ProofLogger::derive_degree_bound(const std::vector<int> &L, const std::vector<int> &R, const BidomainPlan &plan)
{
	// Under the trail (we are inside the pbc for the backtrack clause), the
	// vertices of S can only be matched into B and vice versa. Write p(s,t)
	// for "s and t are matched", U(s) for "s is unmatched", and deg for
	// degrees of the relation (adjacency, or non-adjacency) inside S and B.
	//
	// If p(s,t), every related neighbour of s that is matched must go to a
	// related neighbour of t, so at least deg(s) - deg(t) of them are
	// unmatched. Case splitting on the partner of s gives
	//   V_s:  sum_t (deg(s) - deg(t)) p(s,t) <= sum_{s' ~ s} U(s').
	// Adding the V_s, deg(s) times "s is unmatched or matched into B", and
	// deg(t) times "t has at most one partner in S" cancels every p(s,t):
	//   2 sum_s deg(s) U(s) >= 2 (pairs(S) - pairs(B)),
	// and rounding this knapsack gives sum_s U(s) >= r.
	const std::vector<int> &S = plan.s_is_left ? L : R, &B = plan.s_is_left ? R : L;
	const auto &s_adj = plan.s_is_left ? g_adj : h_adj, &b_adj = plan.s_is_left ? h_adj : g_adj;
	const bool comp = plan.complement;
	auto p = [&](int s, int t) { return plan.s_is_left ? match_var(s, t) : match_var(t, s); };
	auto u = [&](int s) { return plan.s_is_left ? bottom_var(s) : unused_var(s); };

	int s_pairs, b_pairs;
	std::vector<int> s_deg = side_degrees(s_adj, S, comp, s_pairs);
	std::vector<int> b_deg = side_degrees(b_adj, B, comp, b_pairs);

	// each s is unmatched or matched into B
	std::vector<long long> alo(S.size());
	for (size_t i = 0; i < S.size(); i++)
	{
		pbp << "rup 1 " << u(S[i]);
		for (int t : B)
			pbp << " 1 " << p(S[i], t);
		pbp << " >= 1 ;\n";
		alo[i] = next_id++;
	}

	// each t has at most one partner in S
	std::vector<long long> amo(B.size());
	for (size_t j = 0; j < B.size(); j++)
	{
		pbp << "ia";
		for (int s : S)
			pbp << " 1 ~" << p(s, B[j]);
		pbp << " >= " << S.size() - 1 << " : " << (plan.s_is_left ? inj_id[B[j]] : amo_id(B[j])) << " ;\n";
		amo[j] = next_id++;
	}

	std::vector<long long> vs;
	for (size_t i = 0; i < S.size(); i++)
	{
		const int s = S[i];
		bool b_has_degree = std::any_of(b_deg.begin(), b_deg.end(), [](int d) { return d > 0; });
		if (s_deg[i] == 0 && !b_has_degree)
			continue; // V_s is empty

		auto write_v = [&]() {
			for (size_t i2 = 0; i2 < S.size(); i2++)
				if (related(s_adj, s, S[i2], comp))
					pbp << " 1 " << u(S[i2]);
			for (size_t j = 0; j < B.size(); j++)
				if (s_deg[i] != b_deg[j])
					pbp << " " << b_deg[j] - s_deg[i] << " " << p(s, B[j]);
			pbp << " >= 0";
		};

		if (s_deg[i] == 0)
		{
			// every coefficient is non-negative, so this is trivially true
			pbp << "rup";
			write_v();
			pbp << " ;\n";
			vs.push_back(next_id++);
			continue;
		}

		pbp << "pbc";
		write_v();
		pbp << " : subproof\n";
		long long neg = next_id++;
		for (size_t j = 0; j < B.size(); j++)
		{
			const int t = B[j];
			if (s_deg[i] > b_deg[j])
			{
				// p(s,t) => each related neighbour s' of s is unmatched or
				// goes to a related neighbour of t
				std::vector<long long> ts;
				for (int s2 : S)
				{
					if (!related(s_adj, s, s2, comp))
						continue;
					pbp << "rup 1 ~" << p(s, t) << " 1 " << u(s2);
					for (int t2 : B)
						if (related(b_adj, t, t2, comp))
							pbp << " 1 " << p(s2, t2);
					pbp << " >= 1 ;\n";
					ts.push_back(next_id++);
				}
				// with injectivity over t's related neighbours, p(s,t) forces
				// deg(s) - deg(t) of s's related neighbours to be unmatched
				pbp << "pol " << ts[0];
				for (size_t k = 1; k < ts.size(); k++)
					pbp << " " << ts[k] << " +";
				for (size_t j2 = 0; j2 < B.size(); j2++)
					if (related(b_adj, t, B[j2], comp))
						pbp << " " << amo[j2] << " +";
				pbp << " ;\n";
				long long forced = next_id++;
				pbp << "pol " << neg << " " << forced << " + ;\n";
				next_id++;
			}
			pbp << "rup 1 ~" << p(s, t) << " >= 1 ;\n";
			next_id++;
		}
		pbp << "rup >= 1 ;\n";
		next_id++;
		pbp << "qed ;\n";
		vs.push_back(next_id++);
	}

	// sum, then round: weaken the degrees above the r-th largest down to it
	// and divide, leaving sum_s U(s) >= r
	std::vector<int> sorted_deg = s_deg;
	std::sort(sorted_deg.begin(), sorted_deg.end(), std::greater<int>());
	const int cap = sorted_deg[plan.r - 1];
	pbp << "pol";
	bool first = true;
	auto operand = [&](const std::string &op) {
		pbp << " " << op << (first ? "" : " +");
		first = false;
	};
	for (long long v : vs)
		operand(std::to_string(v));
	for (size_t i = 0; i < S.size(); i++)
		if (s_deg[i] > 0)
			operand(std::to_string(alo[i]) + " " + std::to_string(s_deg[i]) + " *");
	for (size_t j = 0; j < B.size(); j++)
		if (b_deg[j] > 0)
			operand(std::to_string(amo[j]) + " " + std::to_string(b_deg[j]) + " *");
	for (size_t i = 0; i < S.size(); i++)
		if (s_deg[i] > cap)
			operand("~" + u(S[i]) + " " + std::to_string(2 * (s_deg[i] - cap)) + " *");
	pbp << " " << 2 * cap << " d";
	for (size_t i = 0; i < S.size(); i++)
		if (s_deg[i] == 0)
			pbp << " " << u(S[i]) << " +";
	pbp << " ;\n";
	long long counted = next_id++;

	if (plan.s_is_left)
		return counted;

	// S is the H side: at least r vertices of R are unused, so at least
	// r + |L| - |R| vertices of L are unmatched
	std::vector<long long> parts;
	for (int c : R)
	{
		// c is unused or has one partner in L, but not both
		pbp << "ia 1 ~" << unused_var(c);
		for (int a : L)
			pbp << " 1 ~" << match_var(a, c);
		pbp << " >= " << L.size() << " : " << unused_amo_id[c] << " ;\n";
		parts.push_back(next_id++);
	}
	for (int a : L)
	{
		// a is unmatched or matched into R
		pbp << "rup 1 " << bottom_var(a);
		for (int c : R)
			pbp << " 1 " << match_var(a, c);
		pbp << " >= 1 ;\n";
		parts.push_back(next_id++);
	}
	pbp << "pol " << counted;
	for (long long id : parts)
		pbp << " " << id << " +";
	pbp << " ;\n";
	return next_id++;
}

void ProofLogger::prune_by_bound(const std::vector<Bidomain> &domains,
								 const std::vector<int> &left, const std::vector<int> &right)
{
	if (aborted)
		return;

	std::vector<std::vector<int>> Ls, Rs;
	int mcsplit_bound = n_matched;
	for (const Bidomain &bd : domains)
	{
		Ls.emplace_back();
		Rs.emplace_back();
		for (int i = 0; i < bd.left_len; i++)
			Ls.back().push_back(g_order[left[bd.l + i]]);
		for (int i = 0; i < bd.right_len; i++)
			Rs.back().push_back(h_order[right[bd.r + i]]);
		mcsplit_bound += std::min(bd.left_len, bd.right_len);
	}

	// operands summing "v takes some value" over L and "w is used at most
	// once" over R: once the trail falsifies every assignment leaving the
	// bidomain, this gives sum_{v in L} x<v>_b >= |L| - |R|
	auto mcsplit_operands = [&](size_t i) {
		for (int v : Ls[i])
			pbp << " " << alo_id(v) << " +";
		for (int w : Rs[i])
			pbp << " " << inj_id[w] << " +";
	};

	if (mcsplit_bound <= incumbent_size)
	{
		// Adding the bidomains with |L| > |R| to the objective bound leaves a
		// constraint that unit propagation refutes under the trail. Smaller
		// bidomains, matched vertices and vertices outside every bidomain
		// are handled by propagation alone.
		bool any = false;
		for (size_t i = 0; i < domains.size(); i++)
		{
			if (Ls[i].size() <= Rs[i].size())
				continue;
			if (!any)
				pbp << "pol " << obj_id;
			any = true;
			mcsplit_operands(i);
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
		return;
	}

	// McSplit's bound is not enough: tighten the bidomains whose bound
	// improves most until it is
	std::vector<BidomainPlan> plans(domains.size());
	std::vector<size_t> order;
	for (size_t i = 0; i < domains.size(); i++)
	{
		plans[i] = plan_bidomain(Ls[i], Rs[i]);
		if (plans[i].kind != BidomainPlan::MinSize)
			order.push_back(i);
	}
	auto gain = [&](size_t i) { return (int)std::min(Ls[i].size(), Rs[i].size()) - plans[i].bound; };
	std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return gain(a) > gain(b); });
	int excess = mcsplit_bound - incumbent_size;
	std::vector<bool> tighten(domains.size(), false);
	for (size_t i : order)
	{
		if (excess <= 0)
			break;
		tighten[i] = true;
		excess -= gain(i);
	}
	if (excess > 0)
		throw std::logic_error("proof logging: cannot justify the bound used to prune this node");

	pbp << "pbc";
	write_trail_clause();
	pbp << " : subproof\n";
	next_id++;
	std::vector<long long> tightened;
	for (size_t i = 0; i < domains.size(); i++)
		if (tighten[i])
		{
			long long id = plans[i].kind == BidomainPlan::Pair
							   ? derive_pair_bound(Ls[i], Rs[i])
							   : derive_degree_bound(Ls[i], Rs[i], plans[i]);
			// pin the derivation to the bound it was meant to give
			pbp << "ia";
			for (int v : Ls[i])
				pbp << " 1 " << bottom_var(v);
			pbp << " >= " << Ls[i].size() - plans[i].bound << " : " << id << " ;\n";
			tightened.push_back(next_id++);
		}
	pbp << "pol " << obj_id;
	for (long long id : tightened)
		pbp << " " << id << " +";
	for (size_t i = 0; i < domains.size(); i++)
		if (!tighten[i] && Ls[i].size() > Rs[i].size())
			mcsplit_operands(i);
	pbp << " ;\n";
	next_id++;
	pbp << "rup >= 1 ;\n";
	next_id++;
	pbp << "qed ;\n";
	record_backtrack_clause(next_id++);
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
