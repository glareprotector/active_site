#ifndef arbi_array_h
#define arbi_array_h

#include <stdarg.h>
#include <assert.h>
#include <fstream>
#include <string>

typedef double num;
using namespace std;


template <class T>
class arbi_array{

 public:

  T* m_data;
  int* dims;
  int* shift_lengths;
  int dim;
  int linear_length;
  
  arbi_array(int, ...);

  arbi_array(const arbi_array& x);

  arbi_array();

  arbi_array<T> operator+(const arbi_array<T>&);

  arbi_array<T>& operator= (const arbi_array<T>&);

  void append(T);

  int size(int);

  T& operator()(...);

  ~arbi_array();

  void fill(const T&);

  arbi_array<T> static transpose(arbi_array<T>);

  


};

template <class T>
ostream& operator<<(ostream& os, const arbi_array<T>& ar);


arbi_array<int> read_mat_to_int(string file, int num_row, int num_col);

  
arbi_array<num> read_mat_to_num(string file, int num_row, int num_col);

arbi_array<int>  read_vect_to_int(string file, int size);


//template <class T>
//ostream& operator<<(ostream& os, const arbi_array<T>& ar);



#include "nums.cpp"    

#endif


