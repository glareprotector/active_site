#ifndef HELPERS_H
#define HELPERS_H

#include <string>
#include <iostream>
#include <algorithm>
#include "lite_fixed.hpp"

namespace helpers{

  template<typename T>
    void append(lite::array<T[1]>& ar, T x){

    lite::array<T[1]> copy(ar);
    int the_size = ar.size().i0;
    ar.resize(the_size + 1);
    for(int i = 0; i < the_size; i++){
      ar(i) = copy(i);
    }
    ar(the_size) = x;
    

  }
  




}

#endif


