#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include "algebra.h"

#define BUCKETSIZE 3

using namespace std;


struct KdNode {
  char bucket;
  int cutdim;
  double cutval;
  KdNode *loson, *hison;
  int lopt, hipt;
};

class KdTree {
 public:
  KdTree(const vector<Point2D> &points, const vector<vector<double> > &lengths);

  double dist(int i, int j);

  int find_closest_point(const Point2D &point, double &dist);
  int nearest_neighbor(int j);

  void out_of_tree_nn(KdNode *p, const Point2D &targetp, int &nnptnum, double &nndist);
  void in_tree_nn(KdNode *p, int &nntarget, int &nnptnum, double &nndist);

  double px(int i, int j); 
  int findmaxspread(int l, int u);
  bool pt_less(int i, int j, int dim);
  void select(int l, int u, int m, int dim);
  KdNode* build(int l, int u);
  void print_point(int i);
  void print_tree(KdNode *node);

 //private:
  vector<int> mPerm;
  const vector<Point2D> *mPoints;
  const vector<vector<double> > *mLengths;
  KdNode *mRoot;
};

#endif
