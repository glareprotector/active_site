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


  void init(int argc, char** argv){
    data_folder = string("/home/fultonw/active_site/active_site/test/");
    pdb_list_file = string("/home/fultonw/active_site/active_site/data/catres_medium_present.pdb_list");
    mean_field_max_iter = 100;
    num_folds = 1;
    which_fold = 0;
  }

}

#endif

