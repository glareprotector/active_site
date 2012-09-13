#ifndef model_h
#define model_h

#include "sample.h"
#include <exception>

//#include <string>

using namespace std;


class model{

 public:


  arbi_array<sample[1]> data;
  int num_samples;
  arbi_array<int1d> training_indicies;
  arbi_array<int1d> testing_indicies;
  int num_training;
  int num_testing;
  int num_folds;
  int which_fold;

  int num_states;
  int num_node_features;
  int num_edge_features;
  int theta_length;
  arbi_array<int2d> node_map;
  arbi_array<int3d> edge_map;

  arbi_array<num1d> gradient;
  num likelihood;


  sample read_sample();
  void load_data();
  void assign(int _num_folds, int _which_fold);
  void normalize();
  model();
  num get_L(int which_obj, arbi_array<num1d> theta);
  num get_reg(arbi_array<num1d> theta);
  arbi_array<num1d> get_dReg_dTheta(arbi_array<num1d> theta);
  arbi_array<num1d> get_dL_dTheta(int which_obj, arbi_array<num1d> theta);
  void report(arbi_array<num1d> theta, int iteration, num obj);
  
  int mean_field_max_iter;
  num reg_constant;
  int which_obj;
  int which_infer;
  int which_reg;

  string results_folder;

};

#endif
