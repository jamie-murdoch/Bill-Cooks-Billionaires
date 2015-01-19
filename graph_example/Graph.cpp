#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

class Graph {
public:
    Graph (const char *filename);

    void print_edges();

private:
    int node_count;
    vector< pair<int,int> > edges;
    vector<int> edge_lengths;
};

//Load a graph from a file
Graph::Graph(const char *filename) {
    int edge_count;

    fstream fin;
    fin.open(filename);

    fin >> node_count >> edge_count;

    edges.resize(edge_count);
    edge_lengths.resize(edge_count);

    for(int i = 0; i < edge_count; i++) {
        fin >> edges[i].first >> edges[i].second >> edge_lengths[i];
    }
}

void Graph::print_edges() {
    for (int i = 0; i < edges.size(); ++i)
    {
        pair<int,int> edge = edges[i];
        cout << edge.first << ", " << edge.second << endl;
    }
}


int main (int argc, char **argv)
{
    char *problem_filename;
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " edge_file" << endl;
        return 1;
    } else {
        problem_filename = argv[1]; 
    }


    Graph base_graph(problem_filename);
    base_graph.print_edges();   
    
    return 0;
}


