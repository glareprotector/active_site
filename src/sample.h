// TEST
/*
This class is responsible for reading in 
1. node features (numNodes x numNodeFeatures)
2. edge features (numEdges x numEdgeFeatures)
3. nodeMap (numStates x numNodeFeatures
4. edgeMap (numStates x numStates x numEdgeFeatures)
5. info - numNodes, numEdges, numStates, numNodeFeatures, numEdgeFeatures
6. edgeNumber to (i,j) map

It represents one sample.  It supports the following functions (given theta)
1. calculate likelihood
2. calculate node marginals
3. calculate edge marginals
4. calculate gradient
*/

#ifndef SAMPLE_H
#define SAMPLE_H


#include <string>
#include <cmath>
#include <ext/hash_map>
#include <utility>
#include <tr1/memory>
#include <algorithm>
#include <iostream>
//#include "model.h"
#include "nums.h"


using namespace std;
using namespace __gnu_cxx;

// hash function for pair<int, int>
struct pair_hash{
  size_t operator() (const pair<int,int>& p) const{
    return p.first + p.second;
  }
 };

//int mean_field_max_iter = 100;

class model;

class sample{

 public:
  
  arbi_array<num> node_potentials;
  arbi_array<num> edge_potentials;
  arbi_array<num> node_features;
  arbi_array<num> edge_features;
  arbi_array<int> true_states;
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;

  num likelihood;

  int num_nodes;
  int num_edges;
    
  hash_map<int, arbi_array<int> > node_to_neighbors;
  hash_map< pair<int, int>, int, pair_hash> edge_to_pos;
  hash_map< int, pair<int, int> > pos_to_edge;

  int mean_field_max_iter;

  string folder;

  model* p_model;
  num get_node_potential(int, int);
  num get_edge_potential(int, int, int, int);
  sample(model*, arbi_array<num>, arbi_array<num>, arbi_array<int>, arbi_array<int>, string);
  sample();
  void set_node_potentials();
  void set_edge_potentials();
  void set_marginals();
  num get_likelihood();
  arbi_array<num> get_gradient();
  num get_log_Z();
  
  
};

#endif
