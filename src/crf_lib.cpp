
#include "sample.cpp"


// this function will initialize things (set cpp_param paths static variable)
#ifndef SERIAL
void init_crf(MPI_Comm comm){
#else
void init_crf(){
#endif
  // set paths
  PyObject* pSysPath= cpp_caller::_get_module_PyObject(string("sys"), string("path"));
  for(int i = 0; i < PyList_Size(pSysPath); i++){
    cpp_caller::added_paths.insert(string(PyString_AsString(PyList_GetItem(pSysPath, i))));
  }
  // set proc_id
#ifndef SERIAL
  MPI_Comm_rank(comm, &proc_id);
  MPI_Comm_size(comm, &num_procs);
#else
  proc_id = 0;
  num_procs = 1;
#endif
}

// this function will be part of shared library that is turned into module
// params should contain testing list(actual list will not be in parameters since it was set earlier) and theta

 pdb_results_struct get_results_given_testing_data_and_theta(PyObject* pMaker, PyObject* pParams, bool recalculate, arbi_array<num1d> theta, int which_infer,  pdb_struct_list testing_pdb_list, int num_states){
  arbi_array<pdb_name_struct[1]> training_pdb_list;
  model m(pMaker, pParams, recalculate, training_pdb_list, testing_pdb_list, num_states);
  cout<<"ONE"<<endl;

  return m.get_results_struct(pMaker, pParams, recalculate, theta, which_infer);
 }
 
 arbi_array<num1d> get_theta_given_training_data_and_hypers(PyObject* pMaker, PyObject* pParams, bool recalculate, int num_states, int max_iter, int which_obj, int which_reg, int which_infer, pdb_struct_list training_pdb_list){


   // first, merge hyper params into params



   cout<<"ZERO"<<endl;
   arbi_array<pdb_name_struct[1]> testing_pdb_list;
   for(int i = 0; i < training_pdb_list.size().i0; i++){
     cout<<training_pdb_list(i).pdb_name<<" "<<i<<endl;
   }

   cout<<"ONE"<<endl;
   model m(pMaker, pParams, recalculate, training_pdb_list, testing_pdb_list, num_states);
   My_Minimizer* minner = new My_Minimizer(&m);
   vector<num> w0(m.theta_length, 0);
   return minner->LBFGS(pMaker, pParams, recalculate, w0, max_iter, which_obj, which_reg, which_infer);
 }
 


