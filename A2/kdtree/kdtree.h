#ifndef KDTREE_H
#define KDTREE_H

#include <vector>

#define BUCKETSIZE 5

//used to initialize the nearest-neighbor search
//#define HUGE 9999

using namespace std;

struct Point {
  double x, y;
};

class PointSet {
 public:
  void addpt(Point newpt);
  int numpts(void);

  
  float x(int i, int j);
  float dist(int i, int j);

 private:
  vector<Point> ptvec;
};

struct KdNode {
  char bucket;
  int cutdim;
  float cutval;
  KdNode *loson, *hison;
  int lopt, hipt;
};

class KdTree {
 public:
  KdTree(PointSet *ptset);
  float px(int i, int j); 
  int findmaxspread(int l, int u);
  bool pt_less(int i, int j, int dim);
  void select(int l, int u, int m, int dim);
  KdNode* build(int l, int u);

 private:
  vector<int> perm;
  PointSet *source_ptset;
  KdNode *root;
};

#endif
