/*
Implement helper array classes here
*/

#include <vector>
<<<<<<< HEAD
#include <iostream>
=======
>>>>>>> 0151e131f606b906a4012a731159d741f266ca71

typedef double num;

using namespace std;

template <class T>
class array_1d{
 
 public:

  vector<T> m_data;
  int len1;

  array_1d(int len, int fill = 0){
<<<<<<< HEAD
    m_data = vector<T>(len, fill);
=======
    m_data = vector<T>(len);
    for(int i = 0; i < len; i++){
      m_data[i] = fill;
    }
>>>>>>> 0151e131f606b906a4012a731159d741f266ca71
    this->len1 = len;
  }

  array_1d()
    :m_data(vector<T>(0)){}

  T& operator()(int x){
    return m_data[x];
  }

<<<<<<< HEAD
  void append(T x){
    m_data.push_back(x);
  }

  //friend ostream& operator<<(ostream& os, const array_1d<T>& ar);

};

template <class T>
ostream& operator<<(ostream& os, const array_1d<T>& ar){
  //os<<"length: "<<ar.len1<<endl;
  for(int i = 0; i < ar.len1; i++){
    //os<<i<<endl;
    os<<ar.m_data[i]<<" ";
  }
  //cout<<"end1d"<<endl;
  return os;
}

template <class T>
=======
};

template <class T>
>>>>>>> 0151e131f606b906a4012a731159d741f266ca71
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

<<<<<<< HEAD
  array_2d(){
    m_data = vector< array_1d<T> >(0, array_1d<T>());
  }
=======
  array_2d()
    :m_data(array_2d(0,0)){}
>>>>>>> 0151e131f606b906a4012a731159d741f266ca71

  T& operator() (int i, int j){
    return m_data[i](j);
  }

};

template <class T>
<<<<<<< HEAD
ostream& operator<<(ostream& os, const array_2d<T>& ar){
  //cout<<"2dlen: "<<ar.len1<<endl;
  for(int i = 0; i < ar.len1; i++){
    //cout<<"QQ"<<i<<endl;
    //cout<<ar.m_data.size()<<endl;
    //    cout<<ar.m_data[i]<<endl;
     os<<ar.m_data[i]<<endl;
    //os<<endl;
    //cout<<"WW"<<endl;
  }
  return os;
}

template <class T>
=======
>>>>>>> 0151e131f606b906a4012a731159d741f266ca71
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

<<<<<<< HEAD
  array_3d(){
    m_data = vector< array_2d<T> >(0, array_2d<T>());
  }
    
=======
  array_3d()
    :m_data(array_3d(0,0,0)){}
>>>>>>> 0151e131f606b906a4012a731159d741f266ca71

  T& operator() (int i, int j, int k){
    return m_data[i](j,k);
  }

};

template <class T>
<<<<<<< HEAD
ostream& operator<<(ostream& os, const array_3d<T>& ar){
  for(int i = 0; i < ar.len1; i++){
    os<<"("<<i<<")"<<endl;
    os<<ar.m_data[i]<<endl;
  }
  return os;
}
=======
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
>>>>>>> 0151e131f606b906a4012a731159d741f266ca71
