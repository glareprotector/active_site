#include <string>

#ifndef GLOBALS_H
#define GLOBALS_H

//#define SERIAL

namespace globals{
  
  string data_folder;
  string pdb_list_file;
  int mean_field_max_iter;
  int num_folds;
  int which_fold;
  string results_folder;


  void init(int argc, char** argv){
    data_folder = string("/home/fultonw/active_site/active_site/test/");
    pdb_list_file = string("/home/fultonw/active_site/active_site/data/catres.pdb_list");
    mean_field_max_iter = 100;
    num_folds = 2;
    which_fold = 0;
    results_folder = string("/home/fultonw/active_site/active_site/src/test_data/");
  }

}

#endif

