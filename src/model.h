#ifndef model_h
#define model_h

#include "sample.h"
#include <exception>

//#include <string>

using namespace std;


class model{

 public:


  arbi_array<num> f_theta;


  arbi_array<sample> data;
  int num_samples;
  arbi_array<int> training_indicies;
  arbi_array<int> testing_indicies;
  int num_training;
  int num_testing;
  int num_folds;
  int which_fold;

  int num_states;
  int num_node_features;
  int num_edge_features;
  int theta_length;
  arbi_array<int> node_map;
  arbi_array<int> edge_map;

  arbi_array<num> gradient;
  num likelihood;
  //arbi_array<num> theta;

  sample read_sample(string folder_name);
  void load_data(arbi_array<string> folder_names);
  void assign(int _num_folds, int _which_fold);
  void normalize();
  model(int _num_states, int _num_node_features, int _num_edge_features, arbi_array<string> _folder_names, int _mean_field_max_iter, int _num_folds, int _which_fold, string results_folder, num _reg_constant, int _which_obj, int _which_infer);
  num get_L(int which_obj, arbi_array<num> theta);
  num get_reg(arbi_array<num> theta);
  arbi_array<num> get_dReg_dTheta(arbi_array<num> theta);
  arbi_array<num> get_dL_dTheta(int which_obj, arbi_array<num> theta);
  void report(arbi_array<num> theta, int iteration, num obj);
  
  int mean_field_max_iter;
  num reg_constant;
  int which_obj;
  int which_infer;
  int which_reg;

  string results_folder;

};

#endif
