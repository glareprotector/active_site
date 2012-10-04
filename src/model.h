#ifndef model_h
#define model_h

#include "sample.h"
#include <exception>

//#include <string>

using namespace std;


class model{

 public:


  arbi_array<sample[1]> data;
  arbi_array< pdb_name_struct[1] > master_pdb_list;

  arbi_array<int1d> all_training_indicies;
  arbi_array<int1d> all_testing_indicies;
  arbi_array<int1d> training_indicies;
  arbi_array<int1d> testing_indicies;
  arbi_array<int1d> do_i_care_idx;


  arbi_array<pdb_name_struct[1]> get_master_pdb_list(string pdb_list_file);
  arbi_array<pdb_name_struct[1]> get_master_pdb_list(arbi_array< pdb_name_struct[1]> training_pdb_list, arbi_array< pdb_name_struct[1]> testing_pdb_list);

  arbi_array<int1d> get_do_i_care_idx(arbi_array< pdb_name_struct[1] > master_pdb_list);

  void load_data(arbi_array< pdb_name_struct[1] > master_pdb_list, arbi_array<int1d> do_i_care_idx);
  
  void get_training_and_testing_indicies(arbi_array< pdb_name_struct[1] > master_pdb_list, int which_fold, int num_folds, arbi_array<int1d>& training_indicies, arbi_array<int1d>& testing_indicies);
  void get_training_and_testing_indicies(arbi_array< pdb_name_struct[1] > master_pdb_list, arbi_array< pdb_name_struct[1] > training_pdb_list, arbi_array< pdb_name_struct[1] > testing_pdb_list, arbi_array<int1d>& training_indicies, arbi_array<int1d>& testing_indicies);

  void get_own_training_and_testing_indicies(arbi_array<int1d> idx_i_care, arbi_array<int1d> training_indicies, arbi_array<int1d> testing_indicies, arbi_array<int1d>& own_training_indicies, arbi_array<int1d>& own_testing_indicies);

  void get_num_node_and_edge_features(this->data, int& num_node_features, int& num_edge_features);

  int get_maps(int num_states, int num_node_features, int num_edge_features, this->node_map, this->edge_map);


  int num_states;
  int num_node_features;
  int num_edge_features;
  int theta_length;
  arbi_array<int2d> node_map;
  arbi_array<int3d> edge_map;

  //arbi_array<num1d> gradient;
  // num likelihood;


  void normalize();
  
  num get_L(int which_obj, arbi_array<num1d> theta);
  num get_reg(arbi_array<num1d> theta);
  arbi_array<num1d> get_dReg_dTheta(arbi_array<num1d> theta);
  arbi_array<num1d> get_dL_dTheta(int which_obj, arbi_array<num1d> theta);
  void report(arbi_array<num1d> theta, int iteration, num obj);


  pdb_results_struct get_results(PyObject* pMaker, PyObject* pParams, bool recalculate);

  
  // int mean_field_max_iter;
  // num reg_constant;
  // int which_obj;
  // int which_obj2;
  // int which_infer;
  // int which_reg;

  //string results_folder;

  num prev_obj;

};

#endif
