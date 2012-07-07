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
#include <Python/Python.h>
#include <sstream>

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
  string pdb_name;
  string chain_letter;

  model* p_model;
  num get_node_potential(arbi_array<num>, int, int);
  num get_edge_potential(arbi_array<num>, int, int, int, int);
  sample(model*, arbi_array<num>, arbi_array<num>, arbi_array<int>, arbi_array<int>, string, string, string);
  sample();

  arbi_array<num> get_node_potentials(arbi_array<num> theta);
  arbi_array<num> get_edge_potentials(arbi_array<num> theta);
  void get_marginals(arbi_array<num> theta, arbi_array<num>& node_marginals, arbi_array<num>& edge_marginals);
  void get_marginals(arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num>& node_marginals, arbi_array<num>& edge_marginals);
  
  num get_L(int which_obj, arbi_array<num>& theta);
  void get_dL_dMu(int which_obj, arbi_array<num> node_marginals, arbi_array<num> edge_marginals, arbi_array<num>& dL_dNode_Mu, arbi_array<num>& dL_dEdge_Mu);
  arbi_array<num> get_dL_dTheta(int which_obj, arbi_array<num> theta);
  arbi_array<num> get_dL_dTheta_Perturb(int which_obj, arbi_array<num> theta);

  void get_dPot_dTheta(arbi_array<num> theta, arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num>& dNode, arbi_array<num>& dEdge);

  // functions related to likelihood
  num get_data_likelihood(arbi_array<num> theta);
  arbi_array<num> get_data_likelihood_gradient(arbi_array<num> node_marginals, arbi_array<num> edge_marginals);
  arbi_array<num> get_data_likelihood_gradient(arbi_array<num> theta);
  num get_log_Z(arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num> node_marginals, arbi_array<num> edge_marginals);
  num get_data_potential(arbi_array<num> node_potentials, arbi_array<num> edge_potentials);

  // separate set of functions for each obj function
  num get_L_expected_distance(arbi_array<num> node_marginals, arbi_array<num> edge_marginals, arbi_array<num> dL_dNode_Mu, arbi_array<num> dL_dEdge_Mu);
  arbi_array<num> get_dL_dMu_expected_distance(arbi_array<num> theta);
  
  num smooth_f(num x);
  num d_smooth_f(num x);
  
  
  
};

#endif
