
#ifndef CRFLIB_H
#define CRFLIB_H

#include "cpp_caller.h"


#ifndef SERIAL
void init_crf(MPI_Comm comm);
#else
void init_crf();
#endif

pdb_results_struct get_results_given_testing_data_and_theta(PyObject* pMaker, PyObject* pParams, bool recalculate, arbi_array<num1d> theta, int which_infer,  pdb_struct_list testing_pdb_list, int num_states);

arbi_array<num1d> get_theta_given_training_data_and_hypers(PyObject* pMaker, PyObject* pParams, bool recalculate, int num_states, int max_iter, int which_obj, int which_reg, int which_infer, pdb_struct_list training_list);




#endif        
