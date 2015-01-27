#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include "util.h"
#include "Graph.h"
#include "radix.h"

using namespace std;

vector<int> get_kruskal_mst_indices(Graph *graph)
{
  vector<Node> block(graph->node_count);
  vector<int> mst_indices(graph->node_count - 1, -1);

  sort(graph->edges.begin(), graph->edges.end());

  Node* p1, *p2;
  int j = 0;
  for (int i = 0; j < graph->node_count - 1; i++){
    p1 = block[graph->edges[i].end[0]].find_canonical(); //find_canonical_with_compression results in slowdown
    p2 = block[graph->edges[i].end[1]].find_canonical();

    if (p1 != p2){
      Node::link(p1, p2);
      mst_indices[j++] = i;
    }
  }

  return mst_indices;
}

int main(int argc, char **argv){
  double szeit;
  vector<int> mst_indices;
  Graph *graph;

  char *filename;
  if (argc < 2) {
      cerr << "Usage: " << argv[0] << " edge_file" << endl;
      return 1;
  } else {
      filename = argv[1]; 
  }

  cout << "Reading in the dataset..." << endl;
  szeit = CO759_zeit();
  graph = new Graph(filename);
  cout << "Dataset loading time: " << CO759_zeit() - szeit << endl;

  szeit = CO759_zeit();
  mst_indices = get_kruskal_mst_indices(graph);
  cout << "New Kruskal running time: " << CO759_zeit() - szeit << endl;

  unsigned long weight = graph->sum_edge_weights(mst_indices);

  cout << "Minimum spanning tree has weight " << weight << endl;

  delete graph;
}
