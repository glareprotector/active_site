#include <string>
#include <Python.h>



#ifndef GLOBALS_H
#define GLOBALS_H

//#define SERIAL
#define PARAM

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
  


  


  const int INT_TYPE = 0;
  const int NUM_TYPE = 1;
  const int STRING_TYPE = 2;
  const int STRING_VECT = 3;
  const int INT_VECT = 4;
  const int NUM_VECT = 5;
  const int STRING_MAT = 6;

  void init(int argc, char** argv){
    data_folder = string("/home/fultonw/active_site/active_site/test/");
    pdb_list_file = string("catres_six.pdb_list");
    mean_field_max_iter = 100;
    //bp_max_iter = 100;
    num_folds = 2;
    which_fold = 0;
    results_folder = string("/home/fultonw/active_site/active_site/src/test_data/");
    reg_constant = 0;
    which_obj = 0;
    which_infer = 0;
    eps = 1.11e-16;
  }

}

#endif

