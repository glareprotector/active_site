#include "nums.h"
#include <string.h>

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
void arbi_array<T>::write(string file_name, char sep = ','){
  ofstream myfile;
  myfile.open(file_name.c_str());
  myfile<<setprecision(26);

  if(dim == 2){
    for(int i = 0; i < dims[0]; i++){
      for(int j = 0; j < dims[1]; j++){
	myfile<<operator()(i,j);
	if(j != dims[1]-1){
	  myfile<<sep;
	}
      }
      if(i != dims[0]-1){
	myfile<<endl;
      }
    }
  }
  else{
    for(int i = 0; i < linear_length; i++){
      myfile<<m_data[i];
      if(i != linear_length-1){
	myfile<<sep;
      }
    }
  }
}

template <class T>
arbi_array<T>& arbi_array<T>::operator+= (const arbi_array<T>& ar){
  assert(this.linear_length == ar.linear_length);
  for(int i = 0; i < this->linear_length; i++){
    this->m_data[i] += ar.m_data[i];
  }
  return *this;
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
  if(linear_length == 0){
    dim = 1;
    dims = new int[1];
    dims[0] = 1;
    shift_lengths = new int[1];
    shift_lengths[0] = 1;
  }
  else{
    dims[0] = dims[0] + 1;
  }
  assert(dim==1);
  T* new_m_data = new T[linear_length+1];
  for(int i = 0; i < linear_length; i++){
    new_m_data[i] = m_data[i];
  }
  linear_length++;

  new_m_data[linear_length-1] = x;
  delete[] m_data;
  m_data = new_m_data;
}

template<class T>
bool arbi_array<T>::operator==(const arbi_array<T>& ar){

  if(linear_length != ar.linear_length){
    return false;
  }
  if(dim != ar.dim){
    return false;
  }
  for(int i = 0; i < dim; i++){
    if(dims[i] != ar.dims[i]){
      return false;
    }
  }
  for(int i = 0; i < linear_length; i++){
    if(m_data[i] != ar.m_data[i]){
      return false;
    }
  }
  return true;
}


template<class T>
int arbi_array<T>::size(int which){
  if(dim == 0){
    return 0;
  }
  else{
    return dims[which];
  }
}

template<class T>
T arbi_array<T>::max(){
  assert(linear_length > 0);
  T ans = m_data[0];
  for(int i = 0; i < linear_length; i++){
    if(ans < m_data[i]){
      ans = m_data[i];
    }
  }
  return ans;
}

template<class T>
inline T& arbi_array<T>::operator()(...){
  
  va_start(iter, NULL);
  
  pos = 0;  

  if(dim == 1){
    int pos1 = va_arg(iter,int);
    assert(pos1 < size(0));
    return this->m_data[pos1];
  }
  if(dim == 2){
    int pos1 = va_arg(iter,int);
    int pos2 = va_arg(iter,int);
    assert(pos1 < size(0));
    assert(pos2 < size(1));
    pos += this->shift_lengths[0] * pos1;
    pos += this->shift_lengths[1] * pos2;
    return this->m_data[pos];
  }
  if(dim == 3){
    int pos1 = va_arg(iter,int);
    int pos2 = va_arg(iter,int);
    int pos3 = va_arg(iter,int);
    assert(pos1 < size(0));
    assert(pos2 < size(1));
    assert(pos3 < size(2));
    pos += this->shift_lengths[0] * pos1;
    pos += this->shift_lengths[1] * pos2;
    pos += this->shift_lengths[2] * pos3;
    return this->m_data[pos];
  }



  for(i = 0; i < this->dim; i++){
    pos += this->shift_lengths[i] * va_arg(iter,int);
  }
  return this->m_data[pos];
}

template<class T>
arbi_array<T>::~arbi_array(){
  if(shift_lengths != 0) delete[] shift_lengths;
  if(dims != 0) delete[] dims;
  if(m_data != 0) delete[] m_data;
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

template<class T>
void arbi_array<T>::scale(num c){
  for(int i = 0; i < linear_length; i++){
    m_data[i] *= c;
  }
}

class file_read_exception: public exception{
  virtual const char* what() const throw(){
    return "file read exception";
  }
};

void get_ifstream(const char* file, ifstream& in){
  in.open(file);
  if(in.fail()){
    cout<<"while reading: "<<file<<" , encountered a: "<<endl;
    throw file_read_exception();
  }
}

arbi_array<int> read_vect_to_int(string file, int size, char sep){
  ifstream in;
  get_ifstream(file.c_str(), in);
  string elt;
  arbi_array<int> ans(1,size);
  for(int i = 0; i < size; i++){
    getline(in, elt, sep);
    ans(i) = atoi(elt.c_str());
  }
  return ans;
}


arbi_array<num> read_vect_to_num(string file, int size, char sep){
  ifstream in;
  get_ifstream(file.c_str(), in);
  string elt;
  arbi_array<num> ans(1,size);
  for(int i = 0; i < size; i++){
    getline(in, elt, sep);
    ans(i) = atof(elt.c_str());
  }
  return ans;
}

arbi_array<num> read_mat_to_num(string file, int num_row, int num_col, const char* sep){
  arbi_array<num> ans(2, num_row, num_col);
  ifstream in;
  get_ifstream(file.c_str(), in);
  string line;
  char* line_cstr = new char[1000000];
  char* elt;
  int i = 0;
  while(in.good() && i < num_row){
    getline(in, line);
    strcpy(line_cstr, line.c_str());
    for(int j = 0; j < num_col; j++){
      if(j == 0){
	elt = strtok(line_cstr, sep);
      }
      else{
	elt = strtok(NULL, sep);
      }
      //cout<<i<<" "<<j<<endl;
      ans(i,j) = atof(elt);
    }

    i++;
  }
  in.close();
  //  delete[] line_cstr;
  //  delete[] elt;
  if(i != num_row){
    cout<<"i and num_row: "<<i<<" "<<num_row<<endl;
  }
  assert(i == num_row);
  return ans;
}

arbi_array<int> read_mat_to_int(string file, int num_row, int num_col, const char* sep){
  arbi_array<int> ans(2, num_row, num_col);
  ifstream in;
  get_ifstream(file.c_str(), in);
  string line;
  char* line_cstr = new char[1000000];
  char* elt;
  int i = 0;
  while(in.good() && i < num_row){
    getline(in, line);
    //cout<<line<<endl;
    strcpy(line_cstr, line.c_str());
    for(int j = 0; j < num_col; j++){
      if(j == 0){
	elt = strtok(line_cstr, sep);
      }
      else{
	elt = strtok(NULL, sep);
      }
      
      //cout<<i<<' '<<j<<' '<<num_col<<' '<<num_row<<endl;
      ans(i,j) = atoi(elt);
    }
    i++;
  }
  in.close();
  assert(i == num_row);
  return ans;
}

// assumes the vect is in a column, not a row
arbi_array<string> read_vect_to_string(string file){
  arbi_array<string> ans;
  ifstream in;
  get_ifstream(file.c_str(), in);
  string line;
  char* line_cstr;
  char* elt;
  int i = 0;
  while(in.good()){
    std::getline(in, line);
    ans.append(line);
  }
  in.close();
  return ans;
}



template <class T>
ostream& operator<<(ostream& os,  arbi_array<T>& ar){

  // if it is 2d mat
  if(ar.dim == 2){
    for(int i = 0; i < ar.dims[0]; i++){
      for(int j = 0; j < ar.dims[1]; j++){
	cout<<ar(i,j)<<' ';
      }
      cout<<endl;
    }
  }
  else{
    for(int i = 0; i < ar.linear_length; i++){
      os<<ar.m_data[i]<<" ";
    }
  }
  return os;
}

