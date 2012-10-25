#include <iostream>

#include "cpp_caller.h"
#include <Python.h>
#include "nums.h"
#include <string>

using namespace std;


int main(){

  Py_Initialize();
  arbi_array<int> indicies;
  arbi_array<num> distances;
  string pdb_name("2jcw");
  string chain_letter("A");
  int aa = 2;
  sorted_distances_getter::get(pdb_name, chain_letter, aa, indicies, distances);
  cout<<indicies<<endl;
  cout<<distances<<endl;

  //read_mat_to_num(string("test_mat.txt"),2,3);

  /*

  double asdf = 1.0d;
  do{
    asdf /= 2.0d;
  }
  while((double)1.0 + asdf != 1.0);
  cout<<asdf<<endl;
  

  return 0;
  */
}
