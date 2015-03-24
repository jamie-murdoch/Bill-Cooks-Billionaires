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

using namespace std;

static void usage (char *f);
static int load_graph (Graph &graph, int ac, char** av);
static double euclid_edgelen (int i, int j, vector<double> x, vector<double> y);
static int rounded_euclid_edgelen(int i, int j, vector<double> x, vector<double> y);
static int delete_edges(Graph g);
bool set_contains(int s, vector<int> R1, vector<int> R2);

int main(int argc, char* argv[]) {
	//Initialize the problem
  Graph graph;
  vector<int> tour_indices;

  //Load/generate the problem
  if(load_graph(graph, argc, argv)) {
    cerr << "Problem creating graph." << endl;
    exit(1);
  }
  cout << "Entering delete graph" << endl;
  delete_edges(graph);
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
  if(g.x[s] == g.x[r]){
    // We can't define a line y=mx + b in this case
    if(g.y[s] > g.y[r]){
      result.push_back(g.x[r]); result.push_back(deltar + g.y[r]);
      //cout << "Exitted circle_proj" << endl;
      return;
    }
    else{
      result.push_back(g.x[r]); result.push_back(g.y[r] - deltar);
      //cout << "Exitted circle_proj" << endl;
      return;
    }
  }
  double m = (g.y[s] - g.y[r]) / (g.x[s] - g.x[r]);
  double x1 = deltar / (m * m + 1) + g.x[r], x2 = g.x[r] - deltar / (m * m + 1);
  //cout << "r = " << r << " " << g.x[r] << " " << g.y[r] << " s = " << s << " " << g.x[s] << " " << g.y[s] << " x1 " << x1 << " x2 " << x2 << " deltar " << deltar << " distance between r and s: " << sqrt(pow(g.y[s] - g.y[r], 2.0) + pow(g.x[s] - g.x[r], 2.0)) << endl;
  //cout << "X1 " << g.x[1] << endl;
  if(g.x[s] > g.x[r]){
    result.push_back(x1); result.push_back(m * (x1 - g.x[r]) + g.y[r]);
    //cout << "Exitted circle_proj" << endl;
    //cout << g.x[1] << endl;
    return;
  }
  else {
    result.push_back(x2); result.push_back(m * (x2 - g.x[r]) + g.y[r]);
    //    cout << "Exited circle_proj" << endl;
    //cout << result[0] << " " << result[1] << endl;
    //    cout << "X1 " << g.x[1] << endl;
    return;
  }
}

double compute_lemma_8(int p, int q, int r, double deltar, double lp, double lq, Graph g, vector<int> &R_p)
{
  vector<double> max_vec;
  vector<double> t_r;
  for(int t = 0; t < g.node_count; t++){
    if(t != r || true ){
      circle_proj(r, t, deltar, g, t_r);
      //cout << "first out X1 " << g.x[1] << endl;
      if(sqrt(pow(g.x[q] - t_r[0],2.0) + pow(g.y[q] - t_r[1],2.0)) >= lq)
	{
	  //	  cout << "Added q " << t << endl;
	  R_p.push_back(t);
	  max_vec.push_back(sqrt(pow(g.x[p] - t_r[0],2.0) + pow(g.y[p] - t_r[1],2.0)));
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

static int delete_edges(Graph g)
{
  vector<double> delta_r;
  vector<vector<double> > rounded_len, len; // Using 2d vectors allows constant access, as opposed to walking through an edge list
  rounded_len.resize(g.node_count); len.resize(g.node_count);
  for(int i = 0; i < g.node_count; i++){
    rounded_len[i].resize(g.node_count); len[i].resize(g.node_count);
    rounded_len[i][i] = numeric_limits<int>::max();
    len[i][i] = numeric_limits<double>::infinity();
  }
  for(vector<Edge>::iterator e = g.edges.begin(); e != g.edges.end(); ++e){
    len[e->end[0]][e->end[1]] = e->len; len[e->end[1]][e->end[0]] = e->len;
    rounded_len[e->end[0]][e->end[1]] = e->rounded_len; rounded_len[e->end[1]][e->end[0]] = e->rounded_len;
  }

  for(int i = 0; i < g.node_count; i++){
    delta_r.push_back(0.5 + *min_element(rounded_len[i].begin(), rounded_len[i].end()) - 1);
  }
  for(vector<Edge>::iterator pq = g.edges.begin(); pq != g.edges.end(); ++pq){
    int p = pq->end[0], q = pq->end[1];
    vector<double> mid_point;
    mid_point.push_back(((double)g.x[p] + g.x[q]) / 2); mid_point.push_back(((double)g.y[p] + g.y[q]) / 2);
    vector<int> potential_points;
    vector<double> dist_to_mid;
    for(int r = 0; r < g.node_count; r++){
      if(r != p && r != q){
	double l_p = delta_r[r] + rounded_len[p][q] - rounded_len[q][r] - 1, l_q = delta_r[r] + rounded_len[p][q] - rounded_len[p][r] - 1;
	if(l_p + l_q >= rounded_len[p][q] - 0.5){
	  double alpha_p = 2 * acos( (l_q * l_q - delta_r[r] * delta_r[r] - len[r][q] * len[r][q]) / (2 * delta_r[r] * len[r][q])), \
	    alpha_q = 2 * acos( (l_p * l_p - delta_r[r] * delta_r[r] - len[r][p] * len[r][p]) / (2 * delta_r[r] * len[r][p])), \
	    gamma_r = acos(1 - pow(l_p + l_q - rounded_len[p][q] + 0.5,2.0) / (2 * delta_r[r] * delta_r[r]));
	  //	  cout << "r " << r << " l_p " << l_p << " l_q " << l_q << " alpha_p " << alpha_p << " alpha_q " << alpha_q << " gamma_r " << gamma_r << " deltar " << delta_r[r] << " len[r][q] " << len[r][q] << " alpha_p arg " << (l_q * l_q - delta_r[r] * delta_r[r] - len[r][q] * len[r][q]) / (2 * delta_r[r] * len[r][q]) << endl;
	  if(gamma_r > max(alpha_p,alpha_q)){
	    potential_points.push_back(r);
	    // updates dist_to_mid
	    dist_to_mid.push_back(sqrt(pow(mid_point[0] - g.x[r],2) + pow(mid_point[1] - g.y[r],2.0)));
	}
	}
      }
    }
    vector<int> points_to_check;
    if(potential_points.size() < 2) continue;
    if(potential_points.size() < 10) points_to_check = potential_points;
    else
      {
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
      double l_p = delta_r[*it] + rounded_len[p][q] - rounded_len[q][*it] - 1, l_q = delta_r[*it] + rounded_len[p][q] - rounded_len[p][*it] - 1;
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

	if(*s != *r && rounded_len[p][q] - rounded_len[*r][*s] + eq_19[j] + eq_20[i] > 0
	   && rounded_len[p][q] - rounded_len[*r][*s] + eq_19[i] + eq_20[j] > 0 && \
	   !set_contains(*r, R_q[j], R_p[j]) && !set_contains(*s, R_q[i], R_p[i])){
	  cout << "Delete edge " << p << " " << q << endl;
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

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] [prob_file]\n", f);
    fprintf (stderr, "   -b d  gridsize d for random problems\n");
    fprintf (stderr, "   -g    prob_file has x-y coordinates\n");
    fprintf (stderr, "   -k d  generate problem with d cities\n");
    fprintf (stderr, "   -s d  random seed\n");
}

static int load_graph (Graph &graph, int ac, char** av) {
	char *fname = (char *) NULL;
	int geometric_data = 0;
	int random_city_count = 0;
	int gridsize_rand = 100;
	int rval = 0;

	//Parse args
	int c;
    int seed = (int) CO759_real_zeit ();

    if (ac == 1) {
        usage (av[0]);
        return 1;
    }

    while ((c = getopt (ac, av, "ab:gk:s:")) != EOF) {
        switch (c) {
        case 'b':
            gridsize_rand = atoi (optarg); 
            break;
        case 'g':
            geometric_data = 1;
            break;
        case 'k':
            random_city_count = atoi (optarg);
            break;
        case 's':
            seed = atoi (optarg);
            break;
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (optind < ac) fname = av[optind++];

    if (optind != ac) {
        usage (av[0]);
        return 1;
    }

    if (!fname && !random_city_count) {
        printf ("Must specify a problem file or use -k for random prob\n");
        return 1;
    }

    printf ("Seed = %d\n", seed);
    srandom (seed);
    //Done parsing args


    FILE *f = (FILE *) NULL;
    int i, j, end1, end2, w, node_count, edge_count;
    vector<double> x, y;

    if (fname) {
        if ((f = fopen (fname, "r")) == NULL) {
    	    fprintf (stderr, "Unable to open %s for input\n",fname);
            rval = 1;  goto CLEANUP;
        }
    }

    if (fname && geometric_data == 0) {
      cout << "ERROR! ERROR! DON'T GO HERE! BAD THINGS MAY HAPPEN" << endl;
        if (fscanf (f, "%d %d", &node_count, &edge_count) != 2) {
       	    fprintf (stderr, "Input file %s has invalid format\n",fname);
            rval = 1;  goto CLEANUP;
        }

        graph.edges.resize(edge_count);

        printf ("Nodes: %d  Edges: %d\n", node_count, edge_count);
        fflush (stdout);

        for (i = 0; i < edge_count; i++) {
    	    if (fscanf(f,"%d %d %d",&end1, &end2, &w) != 3) {
	        	fprintf (stderr, "%s has invalid input format\n", fname);
                rval = 1;  goto CLEANUP;
	   		}
	   		graph.edges[i].end[0] = end1;
	   		graph.edges[i].end[1] = end2;
	   		graph.edges[i].len = w;
        }
    } else {
        if (fname) {
            if (fscanf (f, "%d", &node_count) != 1) {
       	        fprintf (stderr, "Input file %s has invalid format\n",fname);
                rval = 1;  goto CLEANUP;
            }
        } else {
            node_count = random_city_count;
        }

        x.resize(node_count);
        y.resize(node_count);

        if (fname) {
            for (i = 0; i < node_count; i++) {
    	        if (fscanf(f,"%lf %lf",&x[i], &y[i]) != 2) {
	            fprintf (stderr, "%s has invalid input format\n", fname);
                    rval = 1;  goto CLEANUP;
	        }
            }
        } else {
            rval = CO759_build_xy (node_count, x, y, gridsize_rand);
            if (rval) {
                fprintf (stderr, "CO759_build_xy failed\n");
                goto CLEANUP;
            }
    
            printf ("%d\n", node_count);
            for (i = 0; i < node_count; i++) {
                printf ("%.0f %.0f\n", x[i], y[i]);
            }
            printf ("\n");
        }

        edge_count = (node_count * (node_count - 1)) / 2;
        printf ("Complete Graph: %d nodes, %d edges\n", node_count, edge_count);

        graph.edges.resize(edge_count);

        edge_count = 0;
        for (i = 0; i < node_count; i++) {
            for (j = i+1; j < node_count; j++) {
                graph.edges[edge_count].end[0] = i;
		graph.edges[edge_count].end[1] = j;
		graph.edges[edge_count].rounded_len = rounded_euclid_edgelen (i, j, x, y);
		graph.edges[edge_count].len = euclid_edgelen(i,j,x,y);
		edge_count++;
            }
        }
    }

    graph.node_count = node_count;
    graph.edge_count = edge_count;
    graph.x = x;
    graph.y = y;
CLEANUP:
    if (f) fclose (f);
    return rval;
}

static int rounded_euclid_edgelen (int i, int j, vector<double> x, vector<double> y)
{
    double t1 = x[i] - x[j], t2 = y[i] - y[j];
    return (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
}

static double euclid_edgelen(int i, int j, vector<double> x, vector<double> y)
{
  double t1 = x[i] - x[j], t2 = y[i] - y[j];
  return sqrt(t1 * t1 + t2 * t2);
}
