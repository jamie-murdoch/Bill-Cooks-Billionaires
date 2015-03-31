#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "util.h"

#include "Graph.h"
#include "TSP_Solver.h"

#define USE_GRAPHICS false


#if USE_GRAPHICS
#include "graphics.h"
#endif

using namespace std;

static void usage (char *f);
static int load_graph (Graph &graph, int ac, char** av);

static int delete_edges(Graph &g);
static int delete_edges2(Graph &g);
bool are_compatible(int p, int q, int x, int y, Graph &g);
bool is_edge(int p, int q, Graph &g);
//bool set_contains(int s, vector<int> R1, vector<int> R2);
bool set_contains(int r, int s, int q, int p, double lq, double lp, Graph & g, double deltar);

bool compare_results = false;

int main(int argc, char* argv[]) {
    //Initialize the problem
    //Must always pass in a .tsp file in first arg for now
    char *input_file = argv[1];
    cout << "Loading graph..." << flush;
    Graph graph(input_file);
    cout << "done." << endl;

    #if USE_GRAPHICS
        setup_sdl(graph);
        draw_points(graph.points, 0, 0, 255, true);
        //draw_tree(graph.kd_tree->mRoot);
        SDL_RenderPresent(renderer);
    #endif

    // //Test code for Kd_tree
    // for(int j = 0; j < 1000; j++) {
    //     Point2D p(1000 + rand() %4000, 1000 + rand() % 4000);
    //     double dist;
    //     double min_dist = rand() % 1000;
    //     int close = graph.kd_tree->find_closest_point(p, dist, min_dist);

    //     vector<double> dists(graph.node_count());
    //     for(int i = 0; i < graph.node_count(); i++) {
    //         dists[i] = (p - graph.points[i]).length();
    //     }

    //     sort(dists.begin(), dists.end());
    //     int index = 0;
    //     while(index < graph.node_count()) {
    //         if(dists[index] >= min_dist) break;
    //         index++;
    //     }
    //     cout << close << endl;
    //     cout << "actual: " << dists[index] << " nn: " << dist << endl;
    //     if(dists[index] != dist) {
    //         cout << "FAILED! " << endl;
    //         break;
    //     }
    // }

    cout << "Entering delete graph" << endl;

    double running_time = CO759_zeit();
    delete_edges(graph);
    running_time = CO759_zeit() - running_time;

    cout << "Removed " << graph.count_useless() << " out of " << graph.edges.size() << " edges in " << running_time << "s" << endl;
    cout << "Now there are " << graph.edges.size() - graph.count_useless() << " edges left." << endl;

    cout << "Let's try step 2... WARNING: CURRENTLY REMOVING INCORRECT EDGES" << endl;
    running_time = CO759_zeit();
    delete_edges2(graph);
    running_time = CO759_zeit() - running_time;

    cout << "Removed " << graph.count_useless() << " out of " << graph.edges.size() << " edges in " << running_time << "s" << endl;
    cout << "Now there are " << graph.edges.size() - graph.count_useless() << " edges left." << endl;
    
    //Run a test by comparing the removed edges to the ones in the optimal tour
    if(compare_results) {
        cout << "Running the TSP solver to compare results." << endl;
        TSP_Solver solver(graph);
        vector<int> tour_indices;
        solver.find_min_tour(tour_indices);

        //N^2 check for now
        bool pass = true;
        for(int j = 0; j < (int)tour_indices.size(); j++) {
            Edge *tour_edge = &graph.edges[tour_indices[j]];

            if(tour_edge->useless) {
                cout << "***Failed test! Deleted edge " << tour_edge->end[0] << " " << tour_edge->end[1] << " that was in the optimal tour!" << endl;
                pass = false;
            }
        }

        if(pass) {
            cout << "Passed the test. Optimal tour contains none of the removed edges." << endl;
        }
    }

    //Save files
    string pruned_fnam = string(input_file) + "-pruned.edg";
    string orig_fnam = string(input_file) + "-orig.edg";
    //strip path
    pruned_fnam = pruned_fnam.substr(pruned_fnam.find_last_of("\\/") + 1, pruned_fnam.size());
    orig_fnam = orig_fnam.substr(orig_fnam.find_last_of("\\/") + 1, orig_fnam.size());

    graph.save_edges(pruned_fnam, false);
    graph.save_edges(orig_fnam, true);

    cout << "Saved edge files to " << pruned_fnam << " and " << orig_fnam << endl;

    #if USE_GRAPHICS
        getchar();
        SDL_Quit();
    #endif

    return 0;
}

void circle_proj(int r, int s, double deltar, Graph &g, vector<double>& result)
{
    /* This is the function for computing s_r, as described in the paper. That is, given vertices s and r, this function
    computes s_r on the line rs such that |rs_r| = delta_r. This can be solved by writing out the euclidean norm for lines between r and s,
    which is quadratic and can be solved using the quadratic formula. 
    */ 
    if(g.points[s].x() == g.points[r].x()){
        // We can't define a line y=mx + b in this case
        if(g.points[s].y() > g.points[r].y()){
            result.push_back(g.points[r].x()); result.push_back(deltar + g.points[r].y());
            return;
        }
        else{
            result.push_back(g.points[r].x()); result.push_back(g.points[r].y() - deltar);
            return;
        }
    }
    double m = (g.points[s].y() - g.points[r].y()) / (g.points[s].x() - g.points[r].x());
    double x1 = deltar / (m * m + 1) + g.points[r].x(), x2 = g.points[r].x() - deltar / (m * m + 1);

    if(g.points[s].x() > g.points[r].x()){
        result.push_back(x1); result.push_back(m * (x1 - g.points[r].x()) + g.points[r].y());
        return;
    }
    else {
        result.push_back(x2); result.push_back(m * (x2 - g.points[r].x()) + g.points[r].y());
        return;
    }
}

double compute_lemma_8(int p, int q, int r, double deltar, double lp, double lq, Graph &g)
{
  double cos_eps_q = (g.lengths[p][q] * g.lengths[p][q] + g.lengths[q][r] * g.lengths[q][r] - g.lengths[p][r] * g.lengths[p][r]) / (2 * g.lengths[p][q] * g.lengths[q][r]),\
    cos_theta_q = (lq * lq + g.lengths[q][r] * g.lengths[q][r] - deltar * deltar) / (2 * lq * g.lengths[q][r]);
  return deltar - 1 - sqrt(g.lengths[p][q] * g.lengths[p][q] + lq * lq - 2 * g.lengths[p][q]  * lq * (cos_eps_q * cos_theta_q - sqrt(1 - cos_eps_q * cos_eps_q) * sqrt(1 - cos_theta_q * cos_theta_q)));
}

int find_potential(Graph &g, int p, int q, const vector<double> &delta_r, int &i, double &last_dist, const int num_points) {
    const Point2D &pnt_p = g.points[p];
    const Point2D &pnt_q = g.points[q];
    Point2D midpoint((pnt_p + pnt_q) / 2.0);

    for(int k = i; k < num_points; i++, k++){    
        double dist_to_midpoint;

        int r = g.kd_tree->find_closest_point(midpoint, dist_to_midpoint, last_dist);
        last_dist = dist_to_midpoint;

        if(r != p && r != q){
            double l_p = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[q][r] - 1;
            double l_q = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[p][r] - 1;

            if(l_p + l_q >= g.int_lengths[p][q] - 0.5){
                double alpha_p = 2 * acos( (l_q * l_q - delta_r[r] * delta_r[r] - g.lengths[r][q] * g.lengths[r][q]) / (2 * delta_r[r] * g.lengths[r][q]));
                double alpha_q = 2 * acos( (l_p * l_p - delta_r[r] * delta_r[r] - g.lengths[r][p] * g.lengths[r][p]) / (2 * delta_r[r] * g.lengths[r][p]));
                double gamma_r = acos(1 - pow(l_p + l_q - g.int_lengths[p][q] + 0.5,2) / (2 * delta_r[r] * delta_r[r]));

                if(gamma_r > max(alpha_p,alpha_q)){
    		    return r;
                }
            }
        }
    }
    return -1;

}
static int delete_edges(Graph &g)
{
    vector<double> delta_r;

    for(int i = 0; i < g.node_count(); i++){
        delta_r.push_back(0.5 + *min_element(g.int_lengths[i].begin(), g.int_lengths[i].end()) - 1);
    }

    const int num_points = 10;
    //#pragma omp parallel for schedule(dynamic)
    for(int h = 0; h < g.edges.size(); h++){
        Edge *pq = &g.edges[h];
        int p = pq->end[0];
        int q = pq->end[1];
        
        vector<int> potential_points;
        potential_points.reserve(num_points);

	vector<double> eq_19, eq_20;
    	eq_19.reserve(num_points);
        eq_20.reserve(num_points);
	
	int i = -1;
	double last_dist = 0.0;
	for(int u = 0; u < num_points; u++){
	  int r = find_potential(g, p, q, delta_r, u, last_dist, num_points);
	  // cout << "i = " << i << " r = " << new_r << " last_dist = " << last_dist << " potential points size = " << potential_points.size() << endl;
	  if(r != -1) {potential_points.push_back(r); i++;}
	  else break; 
   
	  double l_p = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[q][r] - 1, l_q = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[p][r] - 1;
	  eq_19.push_back(compute_lemma_8(p, q, r, delta_r[r], l_p, l_q, g));
	  eq_20.push_back(compute_lemma_8(q, p, r, delta_r[r], l_q, l_p, g));
        

	  // check to see if we can eliminate the edge
	  bool all_break = false;

	  for(int j = 0; j < i; j++){
	    int s = potential_points[j];

	    // each of these clauses corresponds to the main EE conditions under some re-ordering of (S_1, S_2), (R_1, R_2). 
	    if((g.int_lengths[p][q] - g.int_lengths[r][s] + eq_19[j] + eq_20[i] > 0
		&& g.int_lengths[p][q] - g.int_lengths[r][s] + eq_19[i] + eq_20[j] > 0) ||
	       (g.int_lengths[p][q] - g.int_lengths[r][s] + eq_20[j] + eq_20[i] > 0
	       && g.int_lengths[p][q] - g.int_lengths[r][s] + eq_19[i] + eq_19[j] > 0) ||
	       (g.int_lengths[p][q] - g.int_lengths[r][s] + eq_20[j] + eq_19[i] > 0
	       && g.int_lengths[p][q] - g.int_lengths[r][s] + eq_20[i] + eq_19[j] > 0)){
	      double l_p_s = delta_r[s] + g.int_lengths[p][q] - g.int_lengths[q][s] - 1, l_q_s = delta_r[s] + g.int_lengths[p][q] - g.int_lengths[p][s] - 1; // we could probably cache these 
	      double l_p_r = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[q][r] - 1, l_q_r = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[p][r] - 1; // we could probably cache these
	      
	      if(!set_contains(r, s, q, p, l_q_r, l_p_r, g, delta_r[r]) && !set_contains(s, r, q, p, l_q_s, l_p_s, g, delta_r[s])){
                    
		pq->useless = true;
		g.useless[p][q] = true;
		g.useless[q][p] = true;

		all_break = true;
		break;
	      }
	    }
	  }
	  if(all_break) break;
        }
    }

    cout << "Done edge deletions!" << endl;
    return 0;
}



static int delete_edges2(Graph &g){
  int num_points = min(50, g.node_count() - 2);

#pragma omp parallel for schedule(dynamic)
  for(int h = 0; h < g.edges.size(); h++){
    Edge *pq = &g.edges[h];
    int p = pq->end[0]; int q = pq->end[1];
    vector<vector<Point2D> > edge_pairs; edge_pairs.resize(num_points);
    if (pq->useless) continue;

    const Point2D &pnt_p = g.points[p];
    const Point2D &pnt_q = g.points[q];
    Point2D midpoint((pnt_p + pnt_q) / 2.0);

    double last_dist = 0.0;

    vector<int> r_vec; r_vec.resize(num_points);
    for(int i = 0; i < num_points; i++){    
      double dist_to_midpoint;

      int r = g.kd_tree->find_closest_point(midpoint, dist_to_midpoint, last_dist);
      r_vec[i] = r;
      last_dist = dist_to_midpoint;

      if (r == p || r == q) {i--; continue;}

      vector<int> compatible_points;
      for(int x = 0; x < g.node_count(); x++)
	if(x != p && x != q && is_edge(r, x, g) && are_compatible(p, q, r, x, g)) compatible_points.push_back(x);
      
      for(vector<int>::iterator x = compatible_points.begin(); x != compatible_points.end(); ++x){
	for(vector<int>::iterator y = x; y != compatible_points.end(); ++y){
	  if ((p == *x && q == *y) || (p == *y && q == *x)) continue;
	  if (g.int_lengths[*x][*y] + g.int_lengths[p][r] + g.int_lengths[q][r] >= g.int_lengths[p][q] + g.int_lengths[*x][r] + g.int_lengths[*y][r]){
	    edge_pairs[i].push_back(Point2D(*x,*y));
	    //found_violated_r_inequality = true;
	  }
	}
      }
      if (edge_pairs[i].size() == 0){
	pq->useless = true;
	g.useless[p][q] = true;
	g.useless[q][p] = true;
	break;
      }
    

      bool abort = false;
      for(int j = i - 1; j >= 0; j--){
	int s = r_vec[j];

	//	if(are_compatible(p, q, r, s, g)) continue;
	
	bool all_break = false;
	for(int x = 0; x < edge_pairs[i].size(); x++){
	  for(int y = 0; y < edge_pairs[j].size(); y++){
	    if(edge_pairs[j][y][0] == r || edge_pairs[j][y][1] == r || edge_pairs[i][x][0] == s || edge_pairs[i][x][1] == s){
	      all_break = true; break; // main EE theorem can't be used here
	    }
	      
	    // Sorry team... but I swear its better this way!
	    if(g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[s][edge_pairs[j][y][0]] - g.int_lengths[p][edge_pairs[j][y][0]] + g.int_lengths[r][edge_pairs[i][x][1]] - g.int_lengths[q][edge_pairs[i][x][1]] <= 0 ||
	       g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[r][edge_pairs[i][x][0]] - g.int_lengths[p][edge_pairs[i][x][0]] + g.int_lengths[s][edge_pairs[j][y][1]] - g.int_lengths[q][edge_pairs[j][y][1]] <= 0){
	      if(g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[s][edge_pairs[j][y][1]] - g.int_lengths[p][edge_pairs[j][y][1]] + g.int_lengths[r][edge_pairs[i][x][1]] - g.int_lengths[q][edge_pairs[i][x][1]] <= 0 ||
	       g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[r][edge_pairs[i][x][0]] - g.int_lengths[p][edge_pairs[i][x][0]] + g.int_lengths[s][edge_pairs[j][y][0]] - g.int_lengths[q][edge_pairs[j][y][0]] <= 0){
		if(g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[s][edge_pairs[j][y][1]] - g.int_lengths[p][edge_pairs[j][y][1]] + g.int_lengths[r][edge_pairs[i][x][0]] - g.int_lengths[q][edge_pairs[i][x][0]] <= 0 ||
	       g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[r][edge_pairs[i][x][1]] - g.int_lengths[p][edge_pairs[i][x][1]] + g.int_lengths[s][edge_pairs[j][y][0]] - g.int_lengths[q][edge_pairs[j][y][0]] <= 0){
		  if(g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[s][edge_pairs[j][y][0]] - g.int_lengths[p][edge_pairs[j][y][0]] + g.int_lengths[r][edge_pairs[i][x][0]] - g.int_lengths[q][edge_pairs[i][x][0]] <= 0 ||
	       g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[r][edge_pairs[i][x][1]] - g.int_lengths[p][edge_pairs[i][x][1]] + g.int_lengths[s][edge_pairs[j][y][1]] - g.int_lengths[q][edge_pairs[j][y][1]] <= 0){
	      all_break = true; break;
		  }
		}
	      }
	    }
	  }
	  if(all_break) break;
	}
	if(!all_break){
	  pq->useless = true;
	  g.lengths[p][q] = numeric_limits<double>::infinity();
	  g.lengths[q][p] = numeric_limits<double>::infinity();
	  abort = true; break;
	}
      }
      if(abort) break;
    }
  }

  cout << "Done Step2 deletions!" << endl;
  return 0;
}

//tests whether pq~xy according to equation (1) in the paper
bool are_compatible(int p, int q, int x, int y, Graph &g){
  return max(g.int_lengths[p][x] + g.int_lengths[q][y], g.int_lengths[p][y] + g.int_lengths[q][x]) >= g.int_lengths[p][q] + g.int_lengths[x][y];
}

bool is_edge(int p, int q, Graph &g){
  return !g.useless[p][q];
}

bool set_contains(int r, int s, int q, int p, double lq, double lp, Graph & g, double deltar)
{
  // Check if s \in R_p \cup R_q
  vector<double> s_r; 
  s_r.reserve(2);
  circle_proj(r, s, deltar, g, s_r);
  return sqrt(pow(g.points[q][0] - s_r[0],2) + pow(g.points[q][1] - s_r[1],2)) >= lq || sqrt(pow(g.points[p][0] - s_r[0], 2) + pow(g.points[p][1] - s_r[1], 2)) >= lp;
}

