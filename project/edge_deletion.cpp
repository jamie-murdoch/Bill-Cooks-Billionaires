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

using namespace std;

static void usage (char *f);
static int load_graph (Graph &graph, int ac, char** av);
static int delete_edges(Graph &g);
bool set_contains(int s, vector<int> R1, vector<int> R2);

bool compare_results = false;

int main(int argc, char* argv[]) {
    //Initialize the problem
    //Must always pass in a .tsp file in first arg for now
    Graph graph(argv[1]);
    // cout << graph.points[0] << endl;
    // cout << graph.points[1] << endl;

    // int node = 5;
    // int close = graph.kd_tree->nn(node);
    // double length = graph.lengths[node][close];
    // cout << length<< endl;

    // sort(graph.lengths[node].begin(), graph.lengths[node].end());
    // cout << graph.lengths[node][0] << endl;
    // exit(1);
    vector<int> tour_indices;

    cout << "Entering delete graph" << endl;

    double running_time = CO759_zeit();
    delete_edges(graph);
    running_time = CO759_zeit() - running_time;

    cout << "Removed " << graph.count_useless() << " out of " << graph.edges.size() << " edges in " << running_time << "s" << endl;
    
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

    return 0;
}

void circle_proj(int r, int s, double deltar, Graph g, vector<double>& result)
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

double compute_lemma_8(int p, int q, int r, double deltar, double lp, double lq, Graph g, vector<int> &R_p)
{
    vector<double> max_vec;
    vector<double> t_r;

    for(int t = 0; t < g.node_count(); t++){
        if(t != r || true ){
            circle_proj(r, t, deltar, g, t_r);
            //cout << "first out X1 " << g.points[1].x() << endl;
            if(sqrt(pow(g.points[q].x() - t_r[0],2.0) + pow(g.points[q].y() - t_r[1],2.0)) >= lq)
            {
                //cout << "Added q " << t << endl;
                R_p.push_back(t);
                max_vec.push_back(sqrt(pow(g.points[p].x() - t_r[0],2.0) + pow(g.points[p].y() - t_r[1],2.0)));
            }
            t_r.clear();
        }
    }

    sort(R_p.begin(), R_p.end());

    if(max_vec.size() == 0) {
        return 0; // I'm not really sure if this is the right behaviour here, but it seems to make sense - if R_p is empty, then don't make any moves, and so there's no change in the tree's value
    }

    return deltar - 1 - *max_element(max_vec.begin(), max_vec.end());
}

static int delete_edges(Graph &g)
{
    vector<double> delta_r;

    for(int i = 0; i < g.node_count(); i++){
        delta_r.push_back(0.5 + *min_element(g.int_lengths[i].begin(), g.int_lengths[i].end()) - 1);
    }

    for(vector<Edge>::iterator pq = g.edges.begin(); pq != g.edges.end(); ++pq){
        int p = pq->end[0], q = pq->end[1];
        vector<double> mid_point;
        vector<double> dist_to_mid;
        vector<int> potential_points;

        mid_point.push_back(((double)g.points[p].x() + g.points[q].x()) / 2); mid_point.push_back(((double)g.points[p].y() + g.points[q].y()) / 2);

        for(int r = 0; r < g.node_count(); r++){
            if(r != p && r != q){
                double l_p = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[q][r] - 1, l_q = delta_r[r] + g.int_lengths[p][q] - g.int_lengths[p][r] - 1;

                if(l_p + l_q >= g.int_lengths[p][q] - 0.5){
                    double alpha_p = 2 * acos( (l_q * l_q - delta_r[r] * delta_r[r] - g.lengths[r][q] * g.lengths[r][q]) / (2 * delta_r[r] * g.lengths[r][q]));
                    double alpha_q = 2 * acos( (l_p * l_p - delta_r[r] * delta_r[r] - g.lengths[r][p] * g.lengths[r][p]) / (2 * delta_r[r] * g.lengths[r][p]));
                    double gamma_r = acos(1 - pow(l_p + l_q - g.int_lengths[p][q] + 0.5,2.0) / (2 * delta_r[r] * delta_r[r]));

                    // cout << "r " << r << " l_p " << l_p << " l_q " << l_q << " alpha_p " << alpha_p << " alpha_q " << alpha_q << " gamma_r " << gamma_r << " deltar " << delta_r[r] << " g.lengths[r][q] " << g.lengths[r][q] << " alpha_p arg " << (l_q * l_q - delta_r[r] * delta_r[r] - g.lengths[r][q] * g.lengths[r][q]) / (2 * delta_r[r] * g.lengths[r][q]) << endl;

                    if(gamma_r > max(alpha_p,alpha_q)){
                        potential_points.push_back(r);
                        // updates dist_to_mid
                        dist_to_mid.push_back(sqrt(pow(mid_point[0] - g.points[r].x(),2) + pow(mid_point[1] - g.points[r].y(),2.0)));
                    }
                }
            }
        }

        vector<int> points_to_check;
        if(potential_points.size() < 2) {
            continue;
        }
        else if(potential_points.size() < 10) {
            points_to_check = potential_points;
        }
        else {
            for(int i = 0; i < 10; i++){
                int new_r_ind = distance(dist_to_mid.begin(), min_element(dist_to_mid.begin(), dist_to_mid.end()));
                int new_r = potential_points[new_r_ind];

                potential_points.erase(potential_points.begin() + new_r_ind); dist_to_mid.erase(dist_to_mid.begin() + new_r_ind);
                points_to_check.push_back(new_r);
            }
        }

        // Compute eq_19, eq_20 for those chosen edges
        vector<double> eq_19, eq_20;
        vector<vector<int> > R_p, R_q; vector<int> temp;
        
        eq_19.resize(points_to_check.size()); eq_20.resize(points_to_check.size()); R_q.resize(points_to_check.size()); R_p.resize(points_to_check.size()); // does this do it?

        int i = 0;
        for(vector<int>::iterator it = points_to_check.begin(); it != points_to_check.end(); ++it, i++){
            double l_p = delta_r[*it] + g.int_lengths[p][q] - g.int_lengths[q][*it] - 1, l_q = delta_r[*it] + g.int_lengths[p][q] - g.int_lengths[p][*it] - 1;
            eq_19[i] = compute_lemma_8(p, q, *it, delta_r[*it], l_p, l_q, g, temp);
            R_p[i] = temp; temp.clear();

            eq_20[i] = compute_lemma_8(q, p, *it, delta_r[*it], l_q, l_p, g, temp);
            R_q[i] = temp; 
            temp.clear();
        }

        // check to see if we can eliminate the edge
        bool all_break = false;
        int j;
        i = 0;
        for(vector<int>::iterator r = points_to_check.begin(); r != points_to_check.end(); ++r, i++){
            j = 0;

            for(vector<int>::iterator s = points_to_check.begin(); s != points_to_check.end(); ++s, j++){
                if(*s != *r && g.int_lengths[p][q] - g.int_lengths[*r][*s] + eq_19[j] + eq_20[i] > 0
                    && g.int_lengths[p][q] - g.int_lengths[*r][*s] + eq_19[i] + eq_20[j] > 0 && \
                    !set_contains(*r, R_q[j], R_p[j]) && !set_contains(*s, R_q[i], R_p[i])){
                    
                    pq->useless = true;
                    cout << "Deleted edge: " << pq->end[0] << " " << pq->end[1] <<endl;

                    all_break = true;
                    break;
                }
            }

            if(all_break) break;
        }
    }

    cout << "Done edge deletions!" << endl;
    return 0;
}


bool set_contains(int s, vector<int> R1, vector<int> R2)
{
    return find(R1.begin(), R1.end(), s) != R1.end() || find(R2.begin(), R2.end(), s) != R2.end();
/*  vector<int> intersect, singleton;
singleton.push_back(s);
intersect.resize(max(R1.size(), R2.size())); intersect[0] = -1;

set_intersection(singleton.begin(), singleton.end(), R1.begin(), R1.end(), intersect.begin());
if(intersect[0] != -1) return true;

set_intersection(singleton.begin(), singleton.end(), R2.begin(), R2.end(), intersect.begin());
if(intersect[0] != -1) return true;

return false;*/

}

