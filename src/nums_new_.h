#ifndef arbi_array_h
#define arbi_array_h

#include <stdarg.h>
#include <assert.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <iomanip>

typedef double num;
using namespace std;


template <class T>
class arbi_array{
 private:
  int pos;
  int i;
  va_list iter;

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
  bool operator==(const arbi_array<T>&);

  arbi_array<T>& operator+= (const arbi_array<T>&);

  void append(T);

  int size(int);

  inline T& operator()(...);

  ~arbi_array();

  void fill(const T&);

  T max();

  arbi_array<T> static transpose(arbi_array<T>);

  void scale(num c);

  void write(string file_name, char sep);

};

template <class T>
ostream& operator<<(ostream& os, const arbi_array<T>& ar);


arbi_array<int> read_mat_to_int(string file, int num_row, int num_col, const char* sep = ",");

  
arbi_array<num> read_mat_to_num(string file, int num_row, int num_col, const char* sep = ",");

arbi_array<int>  read_vect_to_int(string file, int size, char sep);

arbi_array<num>  read_vect_to_num(string file, int size, char sep);

arbi_array<string> read_vect_to_string(string file); 

//template <class T>
//ostream& operator<<(ostream& os, const arbi_array<T>& ar);



#include "nums.cpp"    

#endif


