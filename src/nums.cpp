#include "nums.h"


template <class T>
arbi_array<T>::arbi_array(int _dim, ...){
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

template <class T>
arbi_array<T>::arbi_array(const arbi_array<T>& x){
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

template <class T>
arbi_array<T>::arbi_array(){
    this->dim = 0;
    this->dims = 0;
    this->m_data = 0;
    this->shift_lengths = 0;
    this->linear_length = 0;
  }

template <class T>
arbi_array<T> arbi_array<T>::operator+(const arbi_array<T>& ar){
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

template <class T>
arbi_array<T>& arbi_array<T>::operator= (const arbi_array<T>& ar){
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


template <class T>
  void arbi_array<T>::append(T x){
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


template<class T>
  int arbi_array<T>::size(int which){
    return dims[which];
  }

template<class T>
T& arbi_array<T>::operator()(...){
    va_list iter;
    va_start(iter, NULL);
    int pos = 0;
    for(int i = 0; i < this->dim; i++){
      pos += this->shift_lengths[i] * va_arg(iter,int);
    }
    return this->m_data[pos];
  }

template<class T>
arbi_array<T>::~arbi_array(){
    delete[] m_data;
  }

template<class T>
  void arbi_array<T>::fill(const T& val){
    for(int i = 0; i < this->linear_length; i++){
      this->m_data[i] = val;
    }
  }

template<class T>
  arbi_array<T> arbi_array<T>::transpose(arbi_array<T> x){
    assert(x.dim == 2);
    arbi_array<T> ans(2, x.dims[1], x.dims[0]);
    for(int i = 0; i < x.dims[0]; i++){
      for(int j = 0; j < x.dims[1]; j++){
	ans(j,i) = x(i,j);
      }
    }
    return ans;
  }

arbi_array<int>  read_vect_to_int(string file, int size){
    ifstream in(file.c_str());
    string elt;
    arbi_array<int> ans(1,size);
    for(int i = 0; i < size; i++){
      getline(in, elt, ',');
      ans(i) = atoi(elt.c_str());
    }
    return ans;
  }


arbi_array<num> read_mat_to_num(string file, int num_row, int num_col){
    arbi_array<num> ans(2, num_row, num_col);
    ifstream in(file.c_str());
    string line;
    char* line_cstr;
    char* elt;
    int i = 0;
    while(in.good()){
      getline(in, line);
      strcpy(line_cstr, line.c_str());
      for(int j = 0; j < num_col; j++){
	if(j == 0){
	  elt = strtok(line_cstr, ",");
	}
	else{
	  elt = strtok(NULL, ",");
	}
	ans(i,j) = atof(elt);
      }
      i++;
    }
    in.close();
    assert(i == num_row);
    return ans;
  }

 arbi_array<int> read_mat_to_int(string file, int num_row, int num_col){
    arbi_array<int> ans(2, num_row, num_col);
    ifstream in(file.c_str());
    string line;
    char* line_cstr;
    char* elt;
    int i = 0;
    while(in.good()){
      getline(in, line);
      strcpy(line_cstr, line.c_str());
      for(int j = 0; j < num_col; j++){
	if(j == 0){
	  elt = strtok(line_cstr, ",");
	}
	else{
	  elt = strtok(NULL, ",");
	}
	ans(i,j) = atoi(elt);
      }
      i++;
    }
    in.close();
    assert(i == num_row);
    return ans;
  }



template <class T>
ostream& operator<<(ostream& os, const arbi_array<T>& ar){
  for(int i = 0; i < ar.linear_length; i++){
    os<<ar.m_data[i];
  }
  return os;
}

