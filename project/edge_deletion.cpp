#include <iostream>
#include <vector>
#include <limits>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include "util.h"

#include "Graph.h"

using namespace std;

static void usage (char *f);
static int load_graph (Graph &graph, int ac, char** av);
static double euclid_edgelen (int i, int j, double *x, double *y);
static int rounded_euclid_edgelen(int i, int j, double *x, double *y);
static int delete_edges(Graph g);

int main(int argc, char* argv[]) {
	//Initialize the problem
	Graph graph;
    vector<int> tour_indices;

    //Load/generate the problem
	if(load_graph(graph, argc, argv)) {
        cerr << "Problem creating graph." << endl;
        exit(1);
    }


	return 0;
}

static vector<double> circle_proj(int r, int s, double deltar, Graph g)
{
  /* This is the function for computing s_r, as described in the paper. That is, given vertices s and r, this function
     computes s_r on the line rs such that |rs_r| = delta_r. This can be solved by writing out the euclidean norm for lines between r and s,
     which is quadratic and can be solved using the quadratic formula. 
     It may not be a bad idea to write some tests for this bad boy :)
   */
  double m = (g.edges[s].end[1] - g.edges[r].end[1]) / (g.edges[s].end[0] - g.edges[r].end[0]);
  double b = g.edges[r].end[1] + m * g.edges[r].end[0];
  double x1 = (-m * b + sqrt(m * m * b * b - (m * m + 1) * (b * b - deltar * deltar))) / (m * m + 1); // quadratic formula
  double x2 = (-m * b - sqrt(m * m * b * b - (m * m + 1) * (b * b - deltar * deltar))) / (m * m + 1);
  vector<double> result;
  if(x1 > g.edges[r].end[0] && x1 < g.edges[s].end[0]){
    result.push_back(x1); result.push_back(m * (x1 - g.edges[r].end[0]) + g.edges[r].end[1]);
    return result;
  }
  else if(x2 > g.edges[r].end[0] && x2 < g.edges[s].end[0]){
    result.push_back(x2); result.push_back(m * (x2 - g.edges[r].end[0]) + g.edges[r].end[1]);
    return result;
  }
  else {
    cout << "Something broke in circle_proj =(" << endl;
    result.push_back(-1); result.push_back(-1); return result;
  }
}

double compute_lemma_8(int p, double deltar, Graph g)
{
  vector<double> max_vec;
  for(int t = 0; t < g.node_count; t++){
    vector<double> t_r = circle_proj(r, t, delta_r[r], g);
    if(sqrt(pow(g.x[q] - t_r[0],2.0) + pow(g.y[q] - t_r[1],2.0)) >= lp)
      {
	max_vec.push_back(sqrt(pow(g.x[p] - t_r[0],2.0) + pow(g.y[p] - t_r[1],2.0)));
      }
  }
  return deltar - 1 - *max_element(max_vec);
}

static int delete_edges(Graph g)
{
  vector<double> delta_r;
  vector<vector<double>> rounded_len, len; // Using 2d vectors allows constant access, as opposed to walking through an edge list
  rounded_len.resize(g.node_count); len.resize(g.node_count);
  for(int i = 0; i < g.node_count; i++){
    rounded_len[i].resize(g.node_count); len[i].resize(g.node_count);
    rounded_len[i][i] = numeric_limits<int>.max();
    len[i][i] = numeric_limits<double>::infinity();
  }
  for(vector<Edge>::iterator e = g.edges.begin(); e != g.edges.end(); ++e){
    len[e.end[0]][e.end[1]] = e.len; len[e.end[1]][e.end[0]] = e.len;
    rounded_len[e.end[0]][e.end[1]] = e.rounded_len; len[e.end[1]][e.end[0]] = e.rounded_len;
  }

  vector<double> delta_r;
  for(int i = 0; i < g.node_count; i++){
    delta_r.push_back(0.5 + *min_element(rounded_len[i]) - 1);
  }
  
  for(vector<Edge>::iterator pq = g.edges.begin(); pq != g.edges.end(); ++pq){
    int p = pq->end[0], q = pq->end[1];
    vector<int> potential_points;
    vector<double> dist_to_mid;
    for(int r = 0; r < g.node_count; r++){
      if(r != p && r != q){
	double l_p = delta_r[r] + rounded_len[p][q] - rounded_len[q][r] - 1, l_q = delta_r[r] + rounded_len[p][q] - rounded_len[p][r] - 1;
	double alpha_p = 2 * acos( (l_q * l_q - delta_r[r] * delta_r[r] - len[r][q] * len[r][q]) / (2 * delta_r[r] * len[r][q])), \
	  alpha_q = 2 * acos( (l_p * l_p - delta_r[r] * delta_r[r] - len[r][p] * len[r][p]) / (2 * delta_r[r] * len[r][p])), \
	  gamma_r = acos(1 - pow(l_p + l_q - rounded_len[p][q] + 0.5,2.0) / (2 * delta_r[r] * delta_r[r]));
	if(gamma_r > max(alpha_p,alpha_q)){
	  potential_points.push_back(r);
	  // updates dist_to_mid
	}
      }
    }

    // Compute a smart order to look through potential_points - distance from the midpoint of pq?

    // Compute eq_19, eq_20 for those chosen edges
    vector<double> eq_19, eq_20;
    
    // check to see if we can eliminate the edge
    
	  }
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
    double *x = (double *) NULL, *y = (double *) NULL;

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

        x = (double *) malloc (node_count * sizeof (double));
        y = (double *) malloc (node_count * sizeof (double));
        if (!x || !y) {
            fprintf (stdout, "out of memory for x or y\n");
            rval = 1; goto CLEANUP;
        }

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

static int rounded_euclid_edgelen (int i, int j, double *x, double *y)
{
    double t1 = x[i] - x[j], t2 = y[i] - y[j];
    return (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
}

static double euclid_edgelen(int i, int j, double *x, double *y)
{
  double t1 = x[i] - x[j], t2 = y[i] - y[j];
  return sqrt(t1 * t1 + t2 * t2);
}
