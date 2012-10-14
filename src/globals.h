
#include <Python.h>

#include <string.h>

#ifndef GLOBALS_H
#define GLOBALS_H

#define SERIAL
#define PARAM
//#define USINGMAIN

#include "lite_fixed.hpp"



#include "helpers.h"


#include <set>
#include <string>

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


// MPI
extern int num_procs;
extern int proc_id;
// end


struct pdb_name_struct{
  
  pdb_name_struct(){
    this->pdb_name = string("");
    this->chain_letter = string("");
  }


  pdb_name_struct(string pdb_name, string chain_letter){
    this->pdb_name = pdb_name;
    this->chain_letter = chain_letter;
  }

  bool operator==(pdb_name_struct other){
    return other.pdb_name==this->pdb_name && other.chain_letter == this->chain_letter;
  }

  string pdb_name;
  string chain_letter;

};

typedef arbi_array<pdb_name_struct[1]> pdb_struct_list;

struct pdb_results_struct{

  pdb_results_struct(arbi_array<num1d> scores, arbi_array<int1d> true_classes, arbi_array<pdb_name_struct[1]>pdb_structs, arbi_array<int1d> sample_lengths){
    this->scores = scores;
    this->true_classes = true_classes;
    this->pdb_structs = pdb_structs;
    this->sample_lengths = sample_lengths;
  }

  pdb_results_struct(){
    this->scores = arbi_array<num1d>();
    this->true_classes = arbi_array<int1d>();
    this->pdb_structs = arbi_array<pdb_name_struct[1]>();
    this->sample_lengths = arbi_array<int1d>();
  }

  arbi_array<num1d> scores;
  arbi_array<int1d> true_classes;
  arbi_array<pdb_name_struct[1]> pdb_structs;
  arbi_array<int1d> sample_lengths;
};



namespace globals{
  
  extern string data_folder;
  extern string pdb_list_file;
  extern int mean_field_max_iter;
  extern int num_folds;
  extern int which_fold;
  extern string results_folder;
  extern num reg_constant;
  extern int which_obj;
  extern int which_infer;
  extern double eps;
  extern PyObject* pParams;
  extern bool recalculate;
  extern bool recalculate_nodewise_loss_f;


  extern int fdsa;
  


  const int INT_TYPE = 0;
  const int NUM_TYPE = 1;
  const int STRING_TYPE = 2;
  const int STRING_VECT = 3;
  const int INT_VECT = 4;
  const int NUM_VECT = 5;
  const int STRING_MAT = 6;

  void init(int, char**);

  /*
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
    }*/

}

#endif

