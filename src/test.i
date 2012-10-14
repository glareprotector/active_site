%module test
%{
#define SWIT_FILE_WITH_INIT
#include "test.h"
#include "crf_lib.h"
  //#include <mpi.h>
  %}

// %include "mpi4py.i"
// %mpi4py_typemap(Comm, MPI_Comm);


#include "cpp_caller.h"










//#include "crf_lib.h"

%typemap(out) arbi_array<num1d> {

  $result = cpp_caller::cpp_num_vect_to_py_float_list($1);

}

%typemap(in) arbi_array<num1d>{

  $1 = cpp_caller::py_float_list_to_cpp_num_vect($input, false);

 }

%typemap(in) pdb_struct_list {
  $1 = cpp_caller::py_pdb_name_struct_list_to_cpp_pdb_name_struct_list($input, false);
 }



%typemap(out) pdb_results_struct {
  $result = cpp_caller::cpp_pdb_results_struct_to_py_pdb_results_struct($1);
 }




void init_crf();


pdb_results_struct get_results_given_testing_data_and_theta(PyObject* pMaker, PyObject* pParams, bool recalculate, arbi_array<num1d> theta, int which_infer,  pdb_struct_list testing_pdb_list, int num_states);

arbi_array<num1d> get_theta_given_training_data_and_hypers(PyObject* pMaker, PyObject* pParams, bool recalculate, int num_states, int max_iter, int which_obj, int which_reg, int which_infer, pdb_struct_list training_list);




int f(int x, arbi_array<num1d> ar);

