#ifndef KDTREE_H
#define KDTREE_H

#include <vector>

#define BUCKETSIZE 5;

//used to initialize the nearest-neighbor search
#define HUGE 9999

using namespace std;

struct kdnode {
  char bucket;
  int cutdim;
  int cutval;
  kdnode *loson, *hison;
  int lopt, hipt;
};


//global permutation vector perm whose elements are 2-elem vectors
//corresponding to x and y coordinate of a pt
vector<vector<int> > perm;

// returns the coordinate (0 for x, 1 for y) in which the entries of
// perm[l...u] have the largest difference
int findmaxspread(int l, int u);

// permutes perm[l...u] so that perm[(l+u)/2] has dim-coordinate no greater
//than any entry to its left and no less than any entry to its right
// this is accomplished (probably inefficiently) by just sorting that range
// in perm based on x or y coord. 
void select(int l, int u, int dim);


#endif
