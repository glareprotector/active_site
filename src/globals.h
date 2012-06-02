#include <string>

#ifndef GLOBALS_H
#define GLOBALS_H

namespace globals{
  
  string data_folder;
  string pdb_list_file;
  int mean_field_max_iter;


  void init(int argc, char** argv){
    data_folder = string("/home/fultonw/active_site/active_site/test/");
    pdb_list_file = string("/home/fultonw/active_site/active_site/data/catres_six.pdb_list");
    mean_field_max_iter = 100;
  }

}

#endif

