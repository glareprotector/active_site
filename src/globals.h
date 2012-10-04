
#include <Python.h>

#include <string.h>

#ifndef GLOBALS_H
#define GLOBALS_H

//#define SERIAL
#define PARAM

#include "lite_fixed.hpp"




//typedef  vector<vector<vector<vector<double>>>> num4d_array;

#define arbi_array lite::array

using namespace std;

typedef double num;

typedef int int1d[1];
typedef int int2d[1][1];
typedef int int3d[1][1][1];
typedef num num1d[1];
typedef num num2d[1][1];
typedef num num3d[1][1][1];
//typedef num num4d[1][1][1][1];
typedef string string1d[1];
typedef string string2d[1][1];

struct pdb_name_struct{
  
  string pdb_name;
  string chain_letter;

};

struct results_struct{
  arbi_array<int1d> scores;
  arbi_array<int1d> true_classes;
  arbi_array<pdb_name_struct[1]> pdb_structs;
  arbi_array<int1d> sample_lengths;
}



namespace globals{
  
  string data_folder;
  string pdb_list_file;
  int mean_field_max_iter;
  int num_folds;
  int which_fold;
  string results_folder;
  num reg_constant;
  int which_obj;
  int which_infer;
  double eps;
  PyObject* pParams;
  bool recalculate;
  bool recalculate_nodewise_loss_f;

  int fdsa;
  


  const int INT_TYPE = 0;
  const int NUM_TYPE = 1;
  const int STRING_TYPE = 2;
  const int STRING_VECT = 3;
  const int INT_VECT = 4;
  const int NUM_VECT = 5;
  const int STRING_MAT = 6;

  void init(int argc, char** argv){
    data_folder = std::string("/home/fultonw/active_site/active_site/test/");
    pdb_list_file = std::string("catres_six.pdb_list");
    mean_field_max_iter = 100;
    //bp_max_iter = 100;
    num_folds = 2;
    which_fold = 0;
    results_folder = std::string("/home/fultonw/active_site/active_site/src/test_data/");
    reg_constant = 0;
    which_obj = 0;
    which_infer = 0;
    eps = 1.11e-16;
    fdsa = 160;
  }

}

#endif

