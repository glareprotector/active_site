#ifndef CRF_LIB_H
#define CRF_LIB_H

// this function will be part of shared library that is turned into module
// params should contain testing list(actual list will not be in parameters since it was set earlier) and theta
results_struct get_results_given_testing_data_and_theta(PyObject* pMaker, PyObject* pParams, bool recalculate){
  model m(pMaker, pParams, recalculate);
  return m.get_results_struct(pMaker, pParams, recalculate);
}

arbi_array<int1d> get_theta_given_training_data_and_hypers(PyObject* pMaker, PyObject* pParams, bool recalculate){

  int max_iter = cpp_params::get_param(pMaker, pParams, string("maxiter"));
  model m(pMaker, pParams, recalculate);
  My_Minimizer* minner = new My_Minimizer(&m);
  vector<num> w0(m.theta_length, 0);
  return minner.LBFGS(w0, max_iter);
}


#endif
