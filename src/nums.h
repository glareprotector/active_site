#ifndef arbi_array_h
#define arbi_array_h

#include <stdarg.h>
#include <assert.h>

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

  arbi_array(int _dim, ...){
    this->dim = _dim;
    dims = new int[this->dim];
    va_list iter;
    va_start(iter, _dim);
    this->linear_length = 1;
    this->shift_lengths = new int[this->dim];
    for(int i = 0; i < this->dim; i++){
      dims[i] = va_arg(iter, int);
      this->linear_length *= dims[i];
    }
    shift_lengths[this->dim-1] = 1;
    for(int i = this->dim-1-1; i >= 0; i--){
      shift_lengths[i] = shift_lengths[i+1] * this->dims[i+1];
    }
    this->m_data = new T[this->linear_length];
    //cout<<"linear_length in constructor: "<<this->linear_length<<endl;
    for(int i = 0; i < dim; i++){
      //cout<<dims[i]<<" ";
    }
    //cout<<endl;
    for(int i = 0; i < dim; i++){
      //cout<<shift_lengths[i]<<" ";
    }
    //cout<<endl;
    //cout<<endl;
  }

  arbi_array(const arbi_array& x){
    this->dim = x.dim;
    this->linear_length = x.linear_length;
    this->m_data = new T[x.linear_length];
    this->dims = new int[x.dim];
    this->shift_lengths = new int[x.dim];
    for(int i = 0; i < x.linear_length; i++){
      this->m_data[i] = x.m_data[i];
    }
    for(int i = 0; i < x.dim; i++){
      this->dims[i] = x.dims[i];
      this->shift_lengths[i] = x.shift_lengths[i];
    }
  }

  arbi_array(){
    this->dim = 0;
    this->dims = 0;
    this->m_data = 0;
    this->shift_lengths = 0;
    this->linear_length = 0;
  }

  arbi_array<T> operator+(const arbi_array<T>& ar){
    arbi_array<T> ans(*this);
    // check that dimensions are equal
    assert(ans.dim == ar.dim);
    for(int i = 0; i < ans.dim; i++){
      assert(ans.dims[i] == ar.dims[i]);
    }
    assert(this->linear_length == ar.linear_length);
    for(int i = 0; i < this->linear_length; i++){
      ans.m_data[i] = this->m_data[i] + ar.m_data[i];
    }
    return ans;
  }

  arbi_array<T>& operator= (const arbi_array<T>& ar){
    assert(this != &ar);
    if (m_data != 0) { delete[] m_data; }
    if (dims != 0) { delete[] dims; }
    if (shift_lengths != 0) { /*cout<<"shift_lengths: "<<shift_lengths<<endl;*/delete[] shift_lengths; }
    dim = ar.dim;
    linear_length = ar.linear_length;
    //cout<<"assignment operator linear_length: "<<linear_length<<endl;
    m_data = new T[linear_length];
    dims = new int[dim];
    shift_lengths = new int[dim];
    for(int i = 0; i < linear_length; i++){
      m_data[i] = ar.m_data[i];
    }
    for(int i = 0; i < dim; i++){
      dims[i] = ar.dims[i];
      shift_lengths[i] = ar.shift_lengths[i];
    }
    return *this;
  }

  void append(T x){
    assert(dim==1);
    T* new_m_data = new T[linear_length+1];
    for(int i = 0; i < linear_length; i++){
      new_m_data[i] = m_data[i];
    }
    linear_length++;
    dims[0] = dims[0] + 1;
    new_m_data[linear_length-1] = x;
    delete[] m_data;
    m_data = new_m_data;
  }

  int size(int which){
    return dims[which];
  }

  T& operator()(...){
    va_list iter;
    va_start(iter, NULL);
    int pos = 0;
    for(int i = 0; i < this->dim; i++){
      pos += this->shift_lengths[i] * va_arg(iter,int);
    }
    return this->m_data[pos];
  }

  ~arbi_array(){
    delete[] m_data;
  }

  void fill(const T& val){
    for(int i = 0; i < this->linear_length; i++){
      this->m_data[i] = val;
    }
  }

  

};

template <class T>
ostream& operator<<(ostream& os, const arbi_array<T>& ar){
  //cout<<"ar.linear_length: "<<ar.linear_length<<endl;
  for(int i = 0; i < ar.linear_length; i++){
    os<<ar.m_data[i];
  }
  return os;
}

#endif
