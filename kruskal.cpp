#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include "util.h"
#include "Graph.h"

using namespace std;

static int kruskal_tree(int ncount, int ecount, vector<Edge>& elist, vector<int>& tlist)
{
  int i, j = 0;
  vector<Node*> block;
  block.reserve(ncount);
  for (i = 0; i < ncount; i++){
    block[i] = new Node();
  }

  tlist.reserve(ncount - 1);
  for (i = 0; i < ncount - 1; i++) tlist[i] = -1;

  sort(elist.begin(), elist.end());
  Node* p1, *p2;
  for (i = 0; j < ncount - 1; i++){
    p1 = block[elist[i].end[0]]->find_canonical(); p2 = block[elist[i].end[1]]->find_canonical();
    if (p1 != p2){
      Node::link(p1, p2);
      tlist[j++] = i;
    }
  }

  for (i = 0; i < ncount; i++)
    delete block[i];
  return 0;
}



int main(int argc, char **argv){

  // char *problem_filename;
  // if (argc < 2) {
  //     cerr << "Usage: " << argv[0] << " edge_file" << endl;
  //     return 1;
  // } else {
  //     problem_filename = argv[1]; 
  // }


  // Graph base_graph(problem_filename);
  // base_graph.print_edges();   
  
  // return 0;

  int ncount; int ecount;
  vector<int> tlist;
  vector<Edge> edgelist;
  Edge myedge; 

  int i;
  int e0, e1, cost;

  ifstream fin;

  if(argc > 1) {
    fin.open(argv[1]);
  }
  else {
    cerr << "Pass input file as first argument." << endl;
    return 1;
  }

  if (fin >> ncount >> ecount){
    while(fin >> e0 >> e1 >> cost){
      myedge.end[0] = e0; myedge.end[1] = e1;
      myedge.len = cost;
      edgelist.push_back(myedge);
    }
  }
  fin.close();

  double szeit = CO759_zeit();
  kruskal_tree(ncount, ecount, edgelist, tlist);
  cout << "Running time " << CO759_zeit() - szeit << endl;

  unsigned long weight = 0;
  for (i = 0; i < ncount - 1; i++){
    weight += edgelist[tlist[i]].len;
  }
  cout << "mst has weight " << weight << endl;
  
}
