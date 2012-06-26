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
#include <math.h>
#include "globals.h"

// globals for MPI
int proc_id;
int num_procs;




using namespace std;
using namespace __gnu_cxx;

// hash function for pair<int, int>
struct pair_hash{
  size_t operator() (const pair<int,int>& p) const{
    return p.first + p.second;
  }
 };



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
    
  arbi_array< arbi_array<int> > node_to_neighbors;
  arbi_array< pair<int,int> > pos_to_edge;
  arbi_array<int> edge_to_pos;

  string folder;

  model* p_model;
  num get_node_potential(int, int);
  num get_edge_potential(int, int, int, int);
  sample(model*, arbi_array<num>, arbi_array<num>, arbi_array<int>, arbi_array<int>, string);
  sample();
  void set_node_potentials();
  void set_edge_potentials();
  void set_marginals();
  
  // functions related to likelihood
  num get_likelihood();
  arbi_array<num> get_likelihood_gradient();
  num get_log_Z();
  num get_config_likelihood();

  // separate set of functions for each obj function
  num get_exp_dist_obj();
  arbi_array<num> get_exp_dist_gradient();
  
  
  
};

#endif
