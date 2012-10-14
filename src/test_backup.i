%module test
%{
#define SWIT_FILE_WITH_INIT
#include "test.h"
#include <mpi.h>
  %}

%include "mpi4py.i"
%mpi4py_typemap(Comm, MPI_Comm);

// pdb_results_struct refers to the c++ version
%typemap(out) pdb_results_struct {

   PyObject* pScores = cpp_helper::cpp_num_vect_to_py_float_list($l.scores);
   PyObject* pTrue_classes = cpp_helper::cpp_int_vect_to_py_int_list($l.true_classes);
   PyObject* pPdb_name_structs = cpp_helper::cpp_pdb_name_struct_list_to_py_pdb_name_struct_list($l.pdb_structs);
   PyObject* pSample_lengths = cpp_helper::cpp_int_vect_to_py_int_list($l.sample_lengths);
   
   PyObject* pConstructor = cpp_helper::get_module_PyObject(string("cross_validation"), string("pdb_results_struct"));
   $result = PyObject_CallMethodObjArgs(pConstructor, pScores, pTrue_classes, pPdb_name_structs, pSample_lengths);

}

%typemap(out) arbi_array<num1d> {

   $result = cpp_caller::cpp_num_vect_to_py_float_list($l);

}




int f(int x, MPI_Comm comm);