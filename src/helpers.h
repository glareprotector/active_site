#ifndef HELPERS_H
#define HELPERS_H

#include <string>
#include <iostream>
#include <algorithm>

namespace helpers{
  
  // a and b can be 0.  
  string split_string(string s, char delim, int a, int b){
    // find the position of nth occurent of delim from the end
    assert(a < b);
    int len = s.size();
    int count = 0;
    int a_pos = len, b_pos = 0;
    for(int i = len-1; i >= 0; i--){
      if(s[i] == delim){
	if(count == a){
	  a_pos = i;
	}
	if(count == b){
	  b_pos = i+1;
	}
	count++;
      }
    }
    a_pos = max(a_pos, 0);
    b_pos = min(b_pos, len-1);
    assert(b_pos < a_pos);
    return s.substr(b_pos, a_pos - b_pos);
  }

}

#endif


