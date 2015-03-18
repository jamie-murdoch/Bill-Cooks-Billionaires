#include <iostream>
#include <vector>

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

static int delete_edges(Graph g)
{
  int p, q;
  vector<double> delta_r;

  for(int i = 0; i < g.node_count; i++){
    int d = 0;
  }
  
  for(vector<Edge>::iterator pq = g.edges.begin(); pq != g.edges.end(); ++pq){
    p = pq->end[0]; q = pq->end[1];
    
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

CLEANUP:
    if (f) fclose (f);
    if (x) free (x);
    if (y) free (y);
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
