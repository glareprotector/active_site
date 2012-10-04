#include "test.h"


set<string> cpp_caller::added_paths;

int f(int x, MPI_Comm comm){
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  std::cout<<rank<<" "<<size<<std::endl;
  
  PyObject* pSysPath= cpp_caller::_get_module_PyObject(string("sys"), string("path"));
  //for(int i = 0; i < PyList_Size(pSysPath); i++){
  //cpp_caller::added_paths.insert(string(PyString_AsString(PyList_GetItem(pSysPath, i))));
  //}

  PyObject* pASDF = cpp_caller::get_module_PyObject(string("math"), string("pow"));
  cpp_caller::py_print(pASDF);

  double y = 1.0;
  double x_sum = 0;
  MPI_Barrier(comm);
  MPI_Reduce(&y, &x_sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  cout<<rank<<" "<<x_sum<<endl;
  cout<<" END "<<endl;

  return x+1;
}
