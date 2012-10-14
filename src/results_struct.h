#ifndef RESULTS_STRUCT_H
#define RESULTS_STRUCT_H


#include "globals.h"

#define arbi_array lite::array

struct results_struct{
  arbi_array<num1d> scores;
  arbi_array<int1d> true_states;
  arbi_array<string1d> pdb_names;
  arbi_array<string1d> chain_letter;
  arbi_array<int1d> sizes;
}




#endif
