/*
Implement helper array classes here
*/

#include <vector>

typedef double num;

using namespace std;

template <class T>
class array_1d{
 
 public:

  vector<T> m_data;
  int len1;

  array_1d(int len, int fill = 0){
    m_data = vector<T>(len);
    for(int i = 0; i < len; i++){
      m_data[i] = fill;
    }
    this->len1 = len;
  }

  array_1d()
    :m_data(vector<T>(0)){}

  T& operator()(int x){
    return m_data[x];
  }

};

template <class T>
class array_2d{

 public:
  
  vector< array_1d<T> > m_data;
  int len1, len2;

  array_2d(int len1, int len2, int fill = 0){
    m_data = vector< array_1d<T> > (len1);
    for(int i = 0; i < len1; i++){
      m_data[i] = array_1d<T>(len2, fill);
    }
    this->len1 = len1;
    this->len2 = len2;
  }

  array_2d()
    :m_data(array_2d(0,0)){}

  T& operator() (int i, int j){
    return m_data[i](j);
  }

};

template <class T>
class array_3d{

 public:

  vector< array_2d<T> > m_data;
  int len1, len2, len3;

  array_3d(int len1, int len2, int len3, int fill = 0){
    m_data = vector< array_2d<T> > (len1);
    for(int i = 0; i < len1; i++){
      m_data[i] = array_2d<T>(len2, len3, fill);
    }
    this->len1 = len1;
    this->len2 = len2;
    this->len3 = len3;
  }

  array_3d()
    :m_data(array_3d(0,0,0)){}

  T& operator() (int i, int j, int k){
    return m_data[i](j,k);
  }

};

template <class T>
class array_4d{

 public:

  vector< array_3d<T> > m_data;
  int len1, len2, len3, len4;

  array_4d(int len1, int len2, int len3, int len4, int fill = 0){
    m_data = vector< array_3d<T> > (len1);
    for(int i = 0; i < len1; i++){
      m_data[i] = array_3d<T>(len2, len3, len4, fill);
    }
    this->len1 = len1;
    this->len2 = len2;
    this->len3 = len3;
    this->len4 = len4;
  }

  array_4d()
    :m_data(array_4d(0,0,0,0)){}

  T& operator() (int i, int j, int k, int l){
    return m_data[i](j,k,l);
  }

};
