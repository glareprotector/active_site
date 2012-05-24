#include "LBFGS.h"
#include <vector>
#include "sample.h"
#include "nums.h"

using namespace std;

class my_minimizer: public Minimizer{
 public:
  
  model* p_model;

  void ComputeGradient(vector<double>& gradient, const vector<double>& x){
    arbi_array<num> ans = p_model->get_gradient();
    assert(gradient.size() == ans.linear_length);
    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = ans(i);
    }
  }

  double ComputeFunction(const vector<double>& x){
    return p_model->get_likelihood();
  }
}
