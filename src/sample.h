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
//#include "nums.h"
#include <math.h>
#include "globals.h"
#include "cpp_caller.h"
#include <sstream>
#include "helpers.h"
#include <cstdlib>

// globals for MPI
// int proc_id;
// int num_procs;

// want to replace array<T[1][1][1]> with array_3d[T]


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


  void register_pys(PyObject* pMaker, PyObject* pParams, bool recalculate);
  void unregister_pys();
  PyObject* get_pMaker();
  PyObject* get_pParams();
  bool get_recalculate();
  PyObject* pMaker_cur;
  PyObject* pParams_cur;
  int recalculate_cur;


  
  void simulate_states(arbi_array<num1d> f_theta);


  arbi_array<num2d> node_features;
  arbi_array<num2d> edge_features;
  arbi_array<int1d> true_states;

  num likelihood;

  int num_nodes;
  int num_edges;
    
  arbi_array< arbi_array<int1d>[1] > node_to_neighbors;
  arbi_array< pair<int,int>[1] > pos_to_edge;
  arbi_array<int2d> edge_to_pos;


  string pdb_name;
  string chain_letter;

  model* p_model;
  num get_node_potential(arbi_array<num2d>&, int, int);
  num get_edge_potential(arbi_array<num3d>&, int, int, int, int);
  sample(PyObject*, PyObject*, bool, model*, arbi_array<num2d>, arbi_array<num2d>, arbi_array<int2d>, arbi_array<int1d>, string, string);
  sample();

  arbi_array<num2d> get_node_potentials(arbi_array<num1d> theta);
  arbi_array<num3d> get_edge_potentials(arbi_array<num1d> theta);
  void get_marginals(arbi_array<num1d> theta, arbi_array<num2d>& node_marginals, arbi_array<num3d>& edge_marginals, int which_infer);
  void get_marginals(arbi_array<num2d> node_potentials, arbi_array<num3d> edge_potentials, arbi_array<num2d>& node_marginals, arbi_array<num3d>& edge_marginals, int which_infer);


  // these 3 are the functions that interface to outside.  they are same as the private ones except they call register and unregister.  so these 3 are wrappers
  void get_marginals(PyObject* pMaker, PyObject* pParams, bool recalculate, arbi_array<num1d> theta, arbi_array<num2d>& node_marginals, arbi_array<num3d>& edge_marginals, int which_infer);
  num get_L(PyObject* pMaker, PyObject* pParams, bool recalculate, int which_obj, arbi_array<num1d>& theta);
  arbi_array<num1d> get_dL_dTheta(PyObject* pMaker, PyObject* pParams, bool recalculate, int which_obj, arbi_array<num1d> theta);



  
  num get_L(int which_obj, arbi_array<num1d>& theta);
  void get_dL_dMu(int which_obj, arbi_array<num2d> node_marginals, arbi_array<num3d> edge_marginals, arbi_array<num2d>& dL_dNode_Mu, arbi_array<num3d>& dL_dEdge_Mu);
  arbi_array<num1d> get_dL_dTheta(int which_obj, arbi_array<num1d> theta);
  arbi_array<num1d> get_dL_dTheta_Perturb(int which_obj, arbi_array<num1d> theta, int which_infer);

  void pseudo_likelihood_helper(arbi_array<num1d> theta, arbi_array<num2d>& node_pseudos);
  num get_L_pseudo(arbi_array<num1d> theta);
  arbi_array<num1d> get_pseudo_likelihood_gradient(arbi_array<num1d> theta);

  arbi_array<num1d> get_feature_values(arbi_array<int1d> states);

  void get_marginals_logistic_regression(arbi_array<num2d> node_potentials, arbi_array<num2d>& node_marginals);


  void get_dPot_dTheta(arbi_array<num1d> theta, arbi_array<num2d> node_potentials, arbi_array<num3d> edge_potentials, arbi_array<num3d>& dNode, vector< vector< vector< vector< num> > > > & dEdge);

  // functions related to likelihood
  num get_data_likelihood(arbi_array<num1d> theta, int which_infer);
  arbi_array<num1d> get_data_likelihood_gradient(arbi_array<num2d> node_marginals, arbi_array<num3d> edge_marginals);
  arbi_array<num1d> get_data_likelihood_gradient(arbi_array<num1d> theta, int which_infer);
  num get_log_Z(arbi_array<num2d> node_potentials, arbi_array<num3d> edge_potentials, arbi_array<num2d> node_marginals, arbi_array<num3d> edge_marginals);
  num get_data_potential(arbi_array<num2d> node_potentials, arbi_array<num3d> edge_potentials);

  // separate set of functions for each obj function
  void get_dL_dMu_expected_distance(arbi_array<num2d> node_marginals, arbi_array<num3d> edge_marginals, arbi_array<num2d>& dL_dNode_Mu, arbi_array<num3d>& dL_dEdge_Mu);
  num get_L_expected_distance(arbi_array<num1d> theta, int which_infer);
  arbi_array<num1d> get_L_expected_distance_node_importance();

  
  void get_dL_dMu_nodewise(arbi_array<num2d> node_marginals, arbi_array<num3d> edge_marginals, arbi_array<num2d>& dL_dNode_Mu, arbi_array<num3d>& dL_dEdge_Mu);
  num get_L_nodewise(arbi_array<num1d> theta, int which_infer);

  num smooth_f(num x);
  num d_smooth_f(num x);

  void get_marginals_mean_field(arbi_array<num2d> node_potentials, arbi_array<num3d> edge_potentials, arbi_array<num2d>& node_marginals, arbi_array<num3d>& edge_marginals);
  void get_marginals_BP(arbi_array<num2d> node_potentials, arbi_array<num3d> edge_potentials, arbi_array<num2d>& node_marginals, arbi_array<num3d>& edge_marginals);
  
  num& get_message(arbi_array<num3d>& msgs, int& i, int& j, int& s);

  // junk that i should eventually get rid of
  arbi_array<num3d> msgs1;
  arbi_array<num3d> msgs2;
  arbi_array<num3d>* old_msgs;
  arbi_array<num3d>* new_msgs;
  int times_called;
};

#endif
