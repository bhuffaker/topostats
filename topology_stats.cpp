/*
** Computes various graph statistics and metrics on the input graph.
**
** Usage: topology_stats <input_graph>
**
** ---------------------------------------------------------------------
** Copyright (C) 2010 The Regents of the University of California.
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <limits>
#include <fstream>
#include <assert.h>
#include <string>
#include <vector>
#include <set>
#include <ext/hash_map>
#include <algorithm>
#include <cmath> // For log(), pow(), and floor()
#include <gsl/gsl_fit.h> // For gsl_fit_linear()
#include <unistd.h> // For getopt()

using namespace std;
using namespace __gnu_cxx; // For hash_map

namespace __gnu_cxx {
template <>
struct hash<std::string> {
        size_t operator() (const std::string& x) const {
                return hash<const char*>()(x.c_str());
	// hash<const char*> already exists
        }
};
}

//typedef int node_id;
typedef string node_id;
hash_map<node_id, set<node_id> > g_graph;
typedef hash_map<node_id, set<node_id> >::iterator graph_it;


template <class typeA, class typeB>
inline typeA pair_first(const pair<typeA,typeB>& in) { return in.first; }

template <class typeA, class typeB>
inline typeB pair_second(const pair<typeA,typeB>& in) { return in.second; }

template <class typeA, class typeB>
inline bool pair_second_lt(const pair<typeA, typeB>& a,
						const pair<typeA, typeB>& b)
{
    return a.second < b.second;
}

inline bool node_cmp(graph_it a, graph_it b)
{
    if (a->second.size() == b->second.size()) {
	return a->first < b->first; // ID sort to break ties.
    } else {
	return a->second.size() > b->second.size(); // Descending sort
    }
}

template <class typeA, class typeB> // Both types must be convertable to double.
double calc_exponent(vector<typeA> xvals, vector<typeB> yvals)
{
    assert(xvals.size() == yvals.size());
    double x[xvals.size()];
    double y[yvals.size()];
    int orig_size = xvals.size();
    // XXX Kludgy way of stripping out invalid values.
    int j = 0;
    for (int i = 0; i < xvals.size(); i++) {
	double temp_x, temp_y;
	double bad_val = -numeric_limits<double>::infinity();
	temp_x = log(xvals[i]);
	temp_y = log(yvals[i]);
	// NB: NaN is only value that doesn't equal itself.
	if (temp_x == temp_x && temp_y == temp_y &&
	    temp_x != bad_val && temp_y != bad_val)
	{
	    x[j] = temp_x;
	    y[j] = temp_y;
	    j++;
	}
    }
    int array_size = j;
    if (array_size == 0) { // Unlikely, but possible
	return numeric_limits<double>::quiet_NaN();
    }

    double unused, slope;
    gsl_fit_linear(x, 1, y, 1, array_size, &unused, &slope,
		    &unused, &unused, &unused, &unused);

    return -slope;
}

// Takes a value and returns its logarithmic bin value, of size
// 10^(exp_bin).
inline double log10_bin(double value, double exp_bin) {
    return pow(10, exp_bin*floor(log10(value)/exp_bin));
}

// Requires that degrees be sorted in ascending order.
void print_log_bins(hash_map<int,double>& vals, vector<int>& degrees,
			hash_map<int, int>& deg_dist, ofstream& out_file)
{
    double exp_binsize = 0.25;
    double bin = 0;
    double val_sum = 0;
    int current_binsize = 0;
    for (int i = 0; i < degrees.size(); i++) {
	int degree = degrees[i];
	int deg_count = deg_dist[degree];
	double new_bin = log10_bin(degree, exp_binsize);
	if (bin && bin < new_bin) {
	    out_file << bin << "\t" << val_sum / current_binsize << endl;
	    current_binsize = 0;
	    val_sum = 0;
	}
	bin = new_bin;
	current_binsize += deg_count;
	val_sum += vals[degree] * deg_count;
    }
    out_file << bin << "\t" << val_sum / current_binsize << endl;
    // If using gnuplot's fstep, need to add a duplicate value at the end.
    double prev_exp = log10(bin);
    out_file << pow(10, prev_exp+exp_binsize) << "\t"
		<< val_sum / current_binsize << endl;
}

void print_voodoo_bins(hash_map<int,double>& vals,
			vector<int>& degrees,
			hash_map<int, int>& deg_dist,
			ofstream& out_file)
{
    int max_binsize = 50;
    int current_binsize = 0;
    double deg_sum = 0;
    int min_deg = 0;
    int max_deg = 0;
    int prev_min_deg = 0;
    double val_sum = 0;
    bool first_loop = true;
    // Runs backwards to accumulate less numerous values first.
    // Used with gnuplot's 'line' style to create plateaus.
    for (int i = degrees.size()-1; i >= 0; i--) {
	int degree = degrees[i];
	int deg_count = deg_dist[degree];
	if (first_loop) {
	    prev_min_deg = degree;
	    first_loop = false;
	}
	if (current_binsize + deg_count > max_binsize) {
//	    out_file << deg_sum / current_binsize << "\t" <<
//			val_sum / current_binsize << endl;
	    out_file << prev_min_deg << "\t" << val_sum / current_binsize << endl;
	    out_file << min_deg << "\t" << val_sum / current_binsize << endl;
	    current_binsize = 0;
	    deg_sum = 0;
	    val_sum = 0;
	    prev_min_deg = min_deg;
	    min_deg = max_deg = 0;
	}
	if (current_binsize == 0) {
	    min_deg = max_deg = degree;
	}
	min_deg = degree < min_deg ? degree : min_deg;
	max_deg = degree > max_deg ? degree : max_deg;
	current_binsize += deg_count;
	deg_sum += degree * deg_count;
	val_sum += vals[degree] * deg_count;
    }
    out_file << prev_min_deg << "\t" << val_sum / current_binsize << endl;
    out_file << min_deg << "\t" << val_sum / current_binsize << endl;
//    out_file << deg_sum / current_binsize << "\t" <<
//		val_sum / current_binsize << endl;
}

void print_node_stats(const char * filename, 
			hash_map<node_id, int> & node_deg,
			hash_map<node_id, double> & node_avg_neighbor_degree,
			hash_map<node_id, int> & node_coreness,
			hash_map<node_id, double> & node_clustering
			)
{
    ofstream ostream;
    ostream.open(filename);
    ostream << "# node_id degree coreness clustering avg_neighbor_degree\n";
    for (hash_map<node_id, int>::iterator node = node_deg.begin();
	node != node_deg.end(); node++) {
	ostream << node->first << " " << node->second << " " << node_coreness[node->first] 
	    << " " << node_clustering[node->first] 
	    << " " << node_avg_neighbor_degree[node->first] << endl;
    }
    ostream.close();
}


int main(int argc, char ** argv)
{
    long num_edges = 0;
    char opt;
    bool voodoo_bin = false;
    bool log_bin = false;
    bool dump_raw_data = false;
    string output_filename = "topostats_default";
    while ((opt = getopt(argc, argv, "dlvO:")) != -1) {
	switch (opt) {
	    case 'O':
		output_filename = optarg;
		break;
	    case 'l':
		log_bin = true;
		break;
	    case 'v':
		voodoo_bin = true;
		break;
	    case 'd':
		dump_raw_data = true;
		break;
	    default:
		cerr << "Unknown option given: " << opt << endl;
		cerr << "   -d dump datafiles\n";
		cerr << "   -l log bin\n";
		cerr << "   -v voodoo bin\n";
		cerr << "   -O output filename\n";
		exit(1);
	}
    }

    if (argc-optind != 1) {
	cerr << "Usage: " << argv[0] << " <input_graph>\n";
	exit(1);
    }
    ifstream links_file(argv[optind]);

    if (!links_file.is_open()) {
	cerr << "Couldn't open link file " << argv[optind] << endl;
	exit(1);
    }
    string link;
    getline(links_file, link);
    while (!links_file.eof()) {
        if (link[0] == '#') { // Comment line, skip.
            getline(links_file, link);
            continue;
        }
	// NB: Assumes input data is numeric.
	node_id node1, node2;
	char buf1[1000], buf2[1000];
//	sscanf(link.c_str(), "%u %u", &node1, &node2);
	sscanf(link.c_str(), "%s %s", buf1, buf2);
	node1 = buf1;
	node2 = buf2;
	if (node1 != node2) { // Cull obviously bad input data.
	    g_graph[node1].insert(node2);
	    g_graph[node2].insert(node1);
	}

	getline(links_file, link);
    }

    int num_nodes = g_graph.size();

    hash_map<node_id, int> node_degs; // Needed for coreness and rich club.
    hash_map<int, int> deg_dist; // Degree counts
    hash_map<int, double> ave_nbr_deg; // Degree -> average neighbor degree.
    hash_map<int, double> clustering; // Degree -> clustering.

    // Used to dump the individual hashs
    // They will only be filled if the -d flag is set
    hash_map<node_id, int> node_coreness;
    hash_map<node_id, double> node_clustering;
    hash_map<node_id, double> node_avg_neighbor_degree;
//    hash_map<int, double> coreness; // Degree -> coreness.

    double ave_deg = 0;
    double ave_ave_nbr_deg = 0;
    double ave_clustering = 0;
    double ave_node_coreness = 0;

    double assort_sum = 0;
    double assort_prod = 0;
    double assort_sq = 0;

    vector<graph_it> nodes;

    for (graph_it node = g_graph.begin(); node != g_graph.end(); node++) {
	vector<node_id> neighbors(node->second.size());
	copy(node->second.begin(), node->second.end(), neighbors.begin());

	int node_deg = neighbors.size();
	nodes.push_back(node);
	deg_dist[node_deg]++;
	node_degs[node->first] = node_deg;
	num_edges += node_deg; // Count here to avoid overcounting repeated links
	double node_ave_nbr_deg = 0;
	double node_clusteringing = 0;
	int neighbor_links = 0;

	for (int i = 0; i < node_deg; i++) {
	    set<node_id>& nbrs_ref = g_graph[neighbors[i]];
	    int nbr_deg = nbrs_ref.size();
	    node_ave_nbr_deg += nbr_deg;
	    assort_sum += (node_deg + nbr_deg)/2.0;
	    assort_prod += node_deg * nbr_deg;
	    assort_sq += (node_deg*node_deg + nbr_deg*nbr_deg)/2.0;
	    for (int j = i+1; j < node_deg; j++) {
		if (nbrs_ref.find(neighbors[j]) != nbrs_ref.end()) {
		    neighbor_links++;
		}
	    }
	}
	node_ave_nbr_deg /= node_deg;
	if (node_deg > 1) {
	    node_clusteringing = 2.0 * neighbor_links / (node_deg*(node_deg-1));
	}

	ave_nbr_deg[node_deg] += node_ave_nbr_deg;
	clustering[node_deg] += node_clusteringing;
	ave_ave_nbr_deg += node_ave_nbr_deg;
	ave_clustering += node_clusteringing;

	if (dump_raw_data) {
	    node_avg_neighbor_degree[node->first] = node_ave_nbr_deg;
	    node_clustering[node->first] = node_clusteringing;
	}
    }

    sort(nodes.begin(), nodes.end(), node_cmp);

    long rich_club_count = 0;
    typedef set<node_id>::const_iterator set_it;

    int top_clique_size = 0;
    vector<double> ranks;
    vector<double> rcc;
    for (int i = 0; i < nodes.size(); i++) {
	graph_it node = nodes[i];
	// Check if node links into existing rich club.
	for (set_it j = node->second.begin(); j != node->second.end(); j++) {
	    // If neighbor degree is greater than ours, it's in rich club.
	    if (node_degs[*j] > node->second.size()) {
		rich_club_count++;
	    // If neighbor degree is less than ours, ignore it.
	    // If neighbor degree is the same as ours, compare node IDs.
	    } else if (node_degs[*j] == node->second.size() &&
			*j < node->first)
	    {
		rich_club_count++;
	    }
	}
	// Size of club is i+1.
	int rho = i+1;
	double rich_club_conn = 1; // By definition, first node has conn of 1.
	if (rho > 1) {
	    rich_club_conn = 2.0 * rich_club_count / rho / (rho-1);
	}
	if (rich_club_conn == 1) {
	    top_clique_size = rho;
	}
//	double normalized_rank = static_cast<double>(rho) / num_nodes;
//	ranks.push_back(normalized_rank);
//	rcc.push_back(rich_club_conn);
    }
//    double rcc_exp = calc_exponent(ranks, rcc);

    // NB: g_graph will be emptied and invalidated in the following process.
    int k = 1;  // Calculate different levels of k-cores.
    vector<int> removed_degs;
    int max_fringe_deg = 0;
    int fringe_size = 0;
    int min_node_coreness = numeric_limits<int>::max();
    bool fringe_found = false;
    while (!g_graph.empty()) {
	bool graph_changed = false;
	removed_degs.clear();
	do {
	    graph_changed = false;
	    graph_it temp;
	    graph_it node = g_graph.begin();
	    while (node != g_graph.end()) {
		int node_deg = node->second.size();
		if (node_deg <= k) {
		    removed_degs.push_back(node_degs[node->first]);
		    // Remove from sets of neighbors
		    for (set_it nbrs = node->second.begin();
					    nbrs != node->second.end(); nbrs++)
		    {
			g_graph[*nbrs].erase(node->first);
		    }
		    if (dump_raw_data) {
			node_coreness[node->first] = k-1;
		    }
		    temp = node;
		    node++;
		    g_graph.erase(temp);
		    graph_changed = true;

		} else {
		    node++;
		}
	    }
	} while (graph_changed);
	if (!fringe_found && removed_degs.size() > 0) {
	    max_fringe_deg = *max_element(removed_degs.begin(),
							removed_degs.end());
	    min_node_coreness = k-1;
	    fringe_size = removed_degs.size();
	    fringe_found = true;
	}
/*
	for (int i = 0; i < removed_degs.size(); i++) {
	    coreness[removed_degs[i]] += k-1;
	}
*/
	ave_node_coreness += removed_degs.size()*(k-1);

	/*
	// XXX Hacked in for outputting k-cores.
	char core_filename[256];
	sprintf(core_filename, "%s.%03d-core", argv[1], k);
	ofstream core_file(core_filename);
	for (graph_it node = g_graph.begin(); node != g_graph.end(); node++) {
	    for (set_it nbr = node->second.begin(); nbr != node->second.end(); nbr++)
	    {
		core_file << node->first << "\t" << *nbr << endl;
	    }
	}
	*/

	k++;
    }
    int core_size = removed_degs.size();
    int min_core_deg = *min_element(removed_degs.begin(), removed_degs.end());
    int max_node_coreness = k-2; // Extra -1 for last k++ before loop end.

    // Double-counted edges due to graph bidirectionality.
    num_edges /= 2;
    assort_sum /= 2;
    assort_prod /= 2;
    assort_sq /= 2;
    double assort_coeff = (assort_prod - assort_sum*assort_sum/num_edges) / 
			    (assort_sq - assort_sum*assort_sum/num_edges);

    typedef vector<int>::const_iterator deg_it;
    vector<int> degrees(deg_dist.size());
    vector<double> deg_ccdf(deg_dist.size());
    transform(deg_dist.begin(), deg_dist.end(), degrees.begin(),
							pair_first<int,int>);

    sort(degrees.begin(), degrees.end());
    int max_deg = degrees[degrees.size()-1];
    int tot_degrees = 0;
    double mean_sq_degree = 0;
    double clustering_coeff = 0;
    hash_map<int, double> ave_nbr_deg_binned;
    hash_map<int, double> clustering_binned;
    ofstream nbr_out, clus_out, ccdf_out;
    if (dump_raw_data) {

	print_node_stats((output_filename + ".nodes").c_str(),
	    node_degs, node_avg_neighbor_degree, node_coreness, node_clustering);

	nbr_out.open((output_filename + ".nbr").c_str());
	clus_out.open((output_filename + ".clus").c_str());
	ccdf_out.open((output_filename + ".ccdf").c_str());
	nbr_out.precision(15);
	clus_out.precision(15);
	ccdf_out.precision(15);
    }
    for (int i = 0; i < degrees.size(); i++) {
	int degree = degrees[i];
	int deg_count = deg_dist[degree];
	deg_ccdf[i] = tot_degrees;
	tot_degrees += deg_count;
	ave_deg += degree*deg_count;
	clustering_coeff += clustering[degree]*degree*(degree-1);
	mean_sq_degree += degree*degree*deg_count;
	ave_nbr_deg[degree] /= deg_count;
	ave_nbr_deg[degree] /= num_nodes-1; // Size-normalizing for comparison
	clustering[degree] /= deg_count;
//	coreness[degree] /= deg_count;
    }
    for (int i = 0; i < degrees.size(); i++) {
	int degree = degrees[i];
	int deg_count = deg_dist[degree];
	deg_ccdf[i] = static_cast<double>(tot_degrees - deg_ccdf[i]) /
								tot_degrees;
	if (dump_raw_data) {
	    clus_out << degree << "\t" << clustering[degree] << endl;
	    nbr_out << degree << "\t" << ave_nbr_deg[degree] << endl;
	    ccdf_out << degree << "\t" << deg_ccdf[i] << endl;
	}
    }

/*
    for (int i = 0; i < deg_ccdf.size(); i++) {
	deg_ccdf[i] = static_cast<double>(tot_degrees - deg_ccdf[i]) /
								tot_degrees;
	ccdf_out << degrees[i] << "\t" << deg_ccdf[i] << endl;
    }
*/
    if (dump_raw_data) {
	nbr_out.close();
	clus_out.close();
	ccdf_out.close();
    }

    if (voodoo_bin) {
	ofstream nbr_bin_out, clus_bin_out;
	nbr_bin_out.open((output_filename + ".nbr_voodoo_bin").c_str());
	clus_bin_out.open((output_filename + ".clus_voodoo_bin").c_str());
	nbr_bin_out.precision(15);
	clus_bin_out.precision(15);
	print_voodoo_bins(ave_nbr_deg, degrees, deg_dist, nbr_bin_out);
	print_voodoo_bins(clustering, degrees, deg_dist, clus_bin_out);
    }
    if (log_bin) {
	ofstream nbr_bin_out, clus_bin_out;
	nbr_bin_out.open((output_filename + ".nbr_log_bin").c_str());
	clus_bin_out.open((output_filename + ".clus_log_bin").c_str());
	nbr_bin_out.precision(15);
	clus_bin_out.precision(15);
	print_log_bins(ave_nbr_deg, degrees, deg_dist, nbr_bin_out);
	print_log_bins(clustering, degrees, deg_dist, clus_bin_out);
    }
    double deg_dist_exp = calc_exponent(degrees, deg_ccdf)+1;
    double power_law_degree = round(pow(num_nodes, 1/(deg_dist_exp-1)));

    clustering_coeff /= mean_sq_degree - ave_deg;

    ave_deg /= num_nodes;
    ave_ave_nbr_deg /= num_nodes;
    ave_ave_nbr_deg /= num_nodes-1; // Size-normalized for comparison
    ave_clustering /= num_nodes;
    ave_node_coreness /= num_nodes;

    double max_ave_nbr_deg = max_element(ave_nbr_deg.begin(), ave_nbr_deg.end(),
					pair_second_lt<int,double>)->second;

/*
    vector<double> plot_vals(deg_dist.size());

    transform(ave_nbr_deg.begin(), ave_nbr_deg.end(), degrees.begin(),
					pair_first<int,double>);
    transform(ave_nbr_deg.begin(), ave_nbr_deg.end(), plot_vals.begin(),
					pair_second<int,double>);
    double ave_nbr_deg_exp = calc_exponent(degrees, plot_vals);

    transform(clustering.begin(), clustering.end(), degrees.begin(),
					pair_first<int,double>);
    transform(clustering.begin(), clustering.end(), plot_vals.begin(),
					pair_second<int,double>);
    double clustering_exp = calc_exponent(degrees, plot_vals);

    transform(coreness.begin(), coreness.end(),
			degrees.begin(), pair_first<int,double>);
    transform(coreness.begin(), coreness.end(),
			plot_vals.begin(), pair_second<int,double>);
    // We expect this correlation to cause a positive slope.
    double coreness_exp = -calc_exponent(degrees, plot_vals);
*/

    cout.precision(15);

    cout << "Number of nodes:\t" << num_nodes << endl;
    cout << "Number of edges:\t" << num_edges << endl;
    cout << "Avg node degree:\t" << ave_deg << endl;
    cout << "Max node degree:\t" << max_deg << endl;
    cout << "Degree dist exponent (via CCDF) [warning: can be inaccurate]:\t" << deg_dist_exp << endl;
    cout << "Power-law maximum degree [warning: can be inaccurate]:\t" << power_law_degree << endl;
    cout << "Normalized max avg neighbor degree:\t" << max_ave_nbr_deg << endl;
    cout << "Normalized avg avg neighbor degree:\t" << ave_ave_nbr_deg << endl;
    cout << "Avg avg neighbor degree:\t" << ave_ave_nbr_deg*(num_nodes-1) << endl;
    //cout << "Avg neighbor degree exponent:\t" << ave_nbr_deg_exp << endl;
    cout << "Assortative coefficient:\t" << assort_coeff << endl;
    cout << "Mean clustering:\t" << ave_clustering << endl;
    cout << "Clustering coefficient:\t" << clustering_coeff << endl;
    //cout << "Clustering exponent:\t" << clustering_exp << endl;
    cout << "Top clique size:\t" << top_clique_size << endl;
    //cout << "Rich club connectivity exponent:\t" << rcc_exp << endl;
    cout << "Min node coreness:\t" << min_node_coreness << endl;
    cout << "Avg node coreness:\t" << ave_node_coreness << endl;
    cout << "Max node coreness:\t" << max_node_coreness << endl;
    cout << "Core size:\t" << core_size << endl;
    cout << "Min degree in core:\t" << min_core_deg << endl;
    cout << "Fringe size:\t" << fringe_size << endl;
    cout << "Max degree in fringe:\t" << max_fringe_deg << endl;
    //cout << "Coreness exponent:\t" << coreness_exp << endl;

    return 0;
}
