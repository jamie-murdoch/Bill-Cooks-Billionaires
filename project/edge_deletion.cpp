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


bool compare_results = true;

int main(int argc, char* argv[]) {
    //Initialize the problem
    //Must always pass in a .tsp file in first arg for now
    cout << "Loading graph..." << flush;
    Graph graph(argv[1]);
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
    It may not be a bad idea to write some tests for this bad boy :)
    */
    //  cout << "Entered circle proj " << endl;
    if(g.points[s].x() == g.points[r].x()){
        // We can't define a line y=mx + b in this case
        if(g.points[s].y() > g.points[r].y()){
            result.push_back(g.points[r].x()); result.push_back(deltar + g.points[r].y());
        //cout << "Exitted circle_proj" << endl;
            return;
        }
        else{
            result.push_back(g.points[r].x()); result.push_back(g.points[r].y() - deltar);
        //cout << "Exitted circle_proj" << endl;
            return;
        }
    }
    double m = (g.points[s].y() - g.points[r].y()) / (g.points[s].x() - g.points[r].x());
    double x1 = deltar / (m * m + 1) + g.points[r].x(), x2 = g.points[r].x() - deltar / (m * m + 1);
    //cout << "r = " << r << " " << g.points[r].x() << " " << g.points[r].y() << " s = " << s << " " << g.points[s].x() << " " << g.points[s].y() << " x1 " << x1 << " x2 " << x2 << " deltar " << deltar << " distance between r and s: " << sqrt(pow(g.points[s].y() - g.points[r].y(), 2.0) + pow(g.points[s].x() - g.points[r].x(), 2.0)) << endl;
    //cout << "X1 " << g.points[1].x() << endl;
    if(g.points[s].x() > g.points[r].x()){
        result.push_back(x1); result.push_back(m * (x1 - g.points[r].x()) + g.points[r].y());
        //cout << "Exitted circle_proj" << endl;
        //cout << g.points[1].x() << endl;
        return;
    }
    else {
        result.push_back(x2); result.push_back(m * (x2 - g.points[r].x()) + g.points[r].y());
        //    cout << "Exited circle_proj" << endl;
        //cout << result[0] << " " << result[1] << endl;
        //    cout << "X1 " << g.points[1].x() << endl;
        return;
    }
}

double compute_lemma_8(int p, int q, int r, double deltar, double lp, double lq, Graph &g)
{
  double cos_eps_q = (g.lengths[p][q] * g.lengths[p][q] + g.lengths[q][r] * g.lengths[q][r] - g.lengths[p][r] * g.lengths[p][r]) / (2 * g.lengths[p][q] * g.lengths[q][r]),\
    cos_theta_q = (lq * lq + g.lengths[q][r] * g.lengths[q][r] - deltar * deltar) / (2 * lq * g.lengths[q][r]);
  return deltar - 1 - sqrt(g.lengths[p][q] * g.lengths[p][q] + lq * lq - 2 * g.lengths[p][q]  * lq * (cos_eps_q * cos_theta_q - sqrt(1 - cos_eps_q * cos_eps_q) * sqrt(1 - cos_theta_q * cos_theta_q)));
    vector<double> max_vec;
    vector<double> t_r;

    for(int t = 0; t < g.node_count(); t++){
        if(t != r || true ){
            circle_proj(r, t, deltar, g, t_r);
            //cout << "first out X1 " << g.points[1].x() << endl;
            if(sqrt(pow(g.points[q].x() - t_r[0],2) + pow(g.points[q].y() - t_r[1],2)) >= lq)
            {
                //cout << "Added q " << t << endl;
                max_vec.push_back(sqrt(pow(g.points[p].x() - t_r[0],2) + pow(g.points[p].y() - t_r[1],2)));
            }
            t_r.clear();
        }
    }


    if(max_vec.size() == 0) {
        return 0; // I'm not really sure if this is the right behaviour here, but it seems to make sense - if R_p is empty, then don't make any moves, and so there's no change in the tree's value
    }

    return deltar - 1 - *max_element(max_vec.begin(), max_vec.end());
}

void find_potential(Graph &g, int p, int q, const vector<double> &delta_r, vector<int> &potential_points) {
    const Point2D &pnt_p = g.points[p];
    const Point2D &pnt_q = g.points[q];
    Point2D midpoint((pnt_p + pnt_q) / 2.0);

    int num_potential = 0;
    double last_dist = 0.0;

    for(int i = 0; i < 10; i++){    
        double dist_to_midpoint;

        int r = g.kd_tree->find_closest_point(midpoint, dist_to_midpoint, last_dist);
        last_dist = dist_to_midpoint;
        
        const Point2D &pnt_r = g.points[r];

        if(r != p && r != q){
            double l_p = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[q][r] - 1;
            double l_q = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[p][r] - 1;

            if(l_p + l_q >= g.int_lengths[p][q] - 0.5){
                double alpha_p = 2 * acos( (l_q * l_q - delta_r[r] * delta_r[r] - g.lengths[r][q] * g.lengths[r][q]) / (2 * delta_r[r] * g.lengths[r][q]));
                double alpha_q = 2 * acos( (l_p * l_p - delta_r[r] * delta_r[r] - g.lengths[r][p] * g.lengths[r][p]) / (2 * delta_r[r] * g.lengths[r][p]));
                double gamma_r = acos(1 - pow(l_p + l_q - g.int_lengths[p][q] + 0.5,2) / (2 * delta_r[r] * delta_r[r]));

                if(gamma_r > max(alpha_p,alpha_q)){
                    potential_points.push_back(r);
                    num_potential++;

                    if(num_potential >= 10) {
                        break; 
                    }
                }
            }
        }
    }

}
static int delete_edges(Graph &g)
{
    vector<double> delta_r;

    for(int i = 0; i < g.node_count(); i++){
        delta_r.push_back(0.5 + *min_element(g.int_lengths[i].begin(), g.int_lengths[i].end()) - 1);
    }

    for(vector<Edge>::iterator pq = g.edges.begin(); pq != g.edges.end(); ++pq){
        int p = pq->end[0];
        int q = pq->end[1];
        
        vector<int> potential_points;
        potential_points.reserve(10);

        find_potential(g, p, q, delta_r, potential_points);

        // Compute eq_19, eq_20 for those chosen edges
        vector<double> eq_19, eq_20;
    	eq_19.resize(potential_points.size());
        eq_20.resize(potential_points.size()); 

        int i = 0;
        for(vector<int>::iterator it = potential_points.begin(); it != potential_points.end(); ++it, i++){
            double l_p = delta_r[*it] + g.int_lengths[p][q] - g.int_lengths[q][*it] - 1, l_q = delta_r[*it] + g.int_lengths[p][q] - g.int_lengths[p][*it] - 1;
            eq_19[i] = compute_lemma_8(p, q, *it, delta_r[*it], l_p, l_q, g);

            eq_20[i] = compute_lemma_8(q, p, *it, delta_r[*it], l_q, l_p, g);
        }

        // check to see if we can eliminate the edge
        bool all_break = false;
        int j;
        i = 0;

        for(vector<int>::iterator r = potential_points.begin(); r != potential_points.end(); ++r, i++){
            j = 0;

            for(vector<int>::iterator s = potential_points.begin(); s != potential_points.end(); ++s, j++){
                if(*s != *r && g.int_lengths[p][q] - g.int_lengths[*r][*s] + eq_19[j] + eq_20[i] > 0
                    && g.int_lengths[p][q] - g.int_lengths[*r][*s] + eq_19[i] + eq_20[j] > 0){
		  double l_p_s = delta_r[*s] + g.int_lengths[p][q] - g.int_lengths[q][*s] - 1, l_q_s = delta_r[*s] + g.int_lengths[p][q] - g.int_lengths[p][*s] - 1; // we could probably cache these... 
		  double l_p_r = delta_r[*r] + g.int_lengths[p][q] - g.int_lengths[q][*r] - 1, l_q_r = delta_r[*r] + g.int_lengths[p][q] - g.int_lengths[p][*r] - 1; // we could probably cache these...
		  if(!set_contains(*r, *s, q, p, l_q_r, l_p_r, g, delta_r[*r]) && !set_contains(*s, *r, q, p, l_q_s, l_p_s, g, delta_r[*s])){
                    
                    pq->useless = true;
		    /*g.int_lengths[p][q] = numeric_limits<int>::max();
		      g.int_lengths[q][p] = numeric_limits<int>::max();*/
		    g.lengths[p][q] = numeric_limits<double>::infinity();
		    g.lengths[q][p] = numeric_limits<double>::infinity();
                    //cout << "Deleted edge: " << pq->end[0] << " " << pq->end[1] <<endl;

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
  for(vector<Edge>::iterator pq = g.edges.begin(); pq != g.edges.end(); ++pq){
    int p = pq->end[0]; int q = pq->end[1];
    vector<vector<Point2D> > edge_pairs; edge_pairs.resize(g.node_count());
    if (pq->useless) continue;
    bool found_violated_r_inequality = false;
    for(int r = 0; r < g.node_count(); r++){
      if (r == p || r == q) continue;
      for(int x = 0; x < g.node_count(); x++){
	if (!is_edge(r, x, g) || !are_compatible(p, q, r, x, g)) continue;
	for(int y = x + 1; y < g.node_count(); y++){
	  if (!is_edge(r, y, g) || !are_compatible(p, q, r, y, g) || y == x || (p == x && q == y) || (p == y && q == x)) continue;
	  if (g.int_lengths[x][y] + g.int_lengths[p][r] + g.int_lengths[q][r] >= g.int_lengths[p][q] + g.int_lengths[x][r] + g.int_lengths[y][r]){
	    edge_pairs[r].push_back(Point2D(x,y));
	    //found_violated_r_inequality = true;
	  }
	}
      }
      if (edge_pairs[r].size() == 0){
	pq->useless = true;
	g.lengths[p][q] = numeric_limits<double>::infinity();
	g.lengths[q][p] = numeric_limits<double>::infinity();
	break;
      }
    }

    if(!pq->useless){
      for(int r = 0; r < g.node_count(); r++){
	bool abort = false;
	if(r == p || r == q) continue;
	for(int s = 0; s < g.node_count(); s++){
	  if(r == s || s == q || s == p) continue;
	  bool all_break = false;
	  for(int i = 0; i < edge_pairs[r].size(); i++){
	    for(int j = 0; j < edge_pairs[s].size(); j++){
	    
	      if(g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[s][edge_pairs[s][j][0]] - g.int_lengths[p][edge_pairs[s][j][0]] + g.int_lengths[r][edge_pairs[r][i][1]] - g.int_lengths[q][edge_pairs[r][i][1]] <= 0 ||
		 g.int_lengths[p][q] - g.int_lengths[r][s] + g.int_lengths[r][edge_pairs[r][i][0]] - g.int_lengths[p][edge_pairs[r][i][0]] + g.int_lengths[s][edge_pairs[s][j][1]] - g.int_lengths[q][edge_pairs[s][j][1]] <= 0){
		all_break = true; break;
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
  }

  cout << "Done Step2 deletions!" << endl;
  return 0;
}

//tests whether pq~xy according to equation (1) in the paper
bool are_compatible(int p, int q, int x, int y, Graph &g){
  return max(g.int_lengths[p][x] + g.int_lengths[q][y], g.int_lengths[p][y] + g.int_lengths[q][x]) >= g.int_lengths[p][q] + g.int_lengths[x][y];
}

bool is_edge(int p, int q, Graph &g){
  return g.lengths[p][q] != numeric_limits<double>::infinity();
}

bool set_contains(int r, int s, int q, int p, double lq, double lp, Graph & g, double deltar)
{
  // Check if s \in R_p \cup R_q
  vector<double> s_r; 
  s_r.reserve(2);
  circle_proj(r, s, deltar, g, s_r);
  return sqrt(pow(g.points[q][0] - s_r[0],2) + pow(g.points[q][1] - s_r[1],2)) >= lq || sqrt(pow(g.points[p][0] - s_r[0], 2) + pow(g.points[p][1] - s_r[1], 2)) >= lp;
  //    return find(R1.begin(), R1.end(), s) != R1.end() || find(R2.begin(), R2.end(), s) != R2.end();
/*  vector<int> intersect, singleton;
singleton.push_back(s);
intersect.resize(max(R1.size(), R2.size())); intersect[0] = -1;

set_intersection(singleton.begin(), singleton.end(), R1.begin(), R1.end(), intersect.begin());
if(intersect[0] != -1) return true;

set_intersection(singleton.begin(), singleton.end(), R2.begin(), R2.end(), intersect.begin());
if(intersect[0] != -1) return true;

return false;*/

}

