#ifndef CPP_PARAM_H
#define CPP_PARAM_H

class cpp_param{
 public:
  
  static arbi_array<int1d> get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    return cpp_caller::py_int_list_to_cpp_int_vect(pResult, true);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, arbi_array<int1d> it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = cpp_int_vect_to_py_int_list(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static arbi_array<int1d>  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    return cpp_caller::py_int_list_to_cpp_int_vect(pResult, true);
  }



  static arbi_array<num1d> get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    return cpp_caller::py_float_list_to_cpp_num_vect(pResult, true);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, arbi_array<num1d> it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = cpp_num_vect_to_py_float_list(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static arbi_array<num1d>  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    return cpp_caller::py_float_list_to_cpp_num_vect(pResult, true);
  }



  static arbi_array<string1d> get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    return cpp_caller::py_string_list_to_cpp_string_vect(pResult, true);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, arbi_array<string1d> it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = cpp_string_vect_to_py_string_list(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static arbi_array<string1d>  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    return cpp_caller::py_string_list_to_cpp_string_vect(PyObject* pResult, true);
  }




  
  static int get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    return PyInt_AsLong(pResult);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, int it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = PyInt_FromLong(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static int  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    int ans = PyInt_AsLong(pResult);
    Py_DECREF(pResult);
    return ans;
  }



  static num get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    num ans = PyFloat_AsDouble(pResult);
    Py_DECREF(pResult);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, num it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = PyFloat_FromDouble(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static num  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    num ans = PyFloat_AsDouble(pResult);
    Py_DECREF(pResult);
    return ans;
  }



  static string get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    string ans = cpp_caller::CPPString_From_PyString(pResult);
    Py_DECREF(pResult);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, string it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = PyString_FromCPPString(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static string  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    string ans = cpp_caller::CPPString_From_PyString(pResult);
    Py_DECREF(pResult);
    return ans;
  }




  static arbi_array<int2d> get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    return cpp_caller::py_int_mat_to_cpp_int_mat(pResult, true);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, arbi_array<int2d> it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = cpp_int_mat_to_py_int_mat(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static arbi_array<int2d>  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    return cpp_caller::py_int_mat_to_cpp_int_mat(pResult, true);
  }



  static arbi_array<num2d> get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    return cpp_caller::py_float_mat_to_cpp_num_mat(pResult, true);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, arbi_array<num1d> it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = cpp_num_mat_to_py_float_mat(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static arbi_array<num2d>  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    return cpp_caller::py_float_mat_to_cpp_num_mat(pResult, true);
  }



  static arbi_array<string2d> get_param(PyObject* pMaker, PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    // convert pResult to the proper data type
    return cpp_caller::py_string_mat_to_cpp_string(pResult, true);
  }

  static void set_param(PyObject* pMaker, PyObject* pParams, string key, arbi_array<string2d> it){
    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pIt = CPPString_mat_to_py_string_mat(it);
    PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pKey, pIt, NULL);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    Py_DECREF(pIt);
  }

  static arbi_array<string2d>  get_var_or_file(PyObject* pMaker, PyObject* pParams, string module_name, string wrapper_name, bool recalculate, bool to_pickle = false, bool to_filelize = false, bool always_recalculate = false){
    PyObject* pMethodName = PyString_FromCPPString(string("get_var_or_file"));
    PyObject* pWrapper = cpp_caller::get_module_PyObject(module_name, wrapper_name);
    PyObject* pResult = PyObject_CallMethodObjArgs(pMaker, pMethodName, pParams, pWrapper, NULL);
    Py_DECREF(pMethodName);
    // convert pResult to the proper data type
    return cpp_caller::py_string_mat_to_cpp_string(PyObject* pResult, true);
  }


  static arbi_array<pdb_name_struct[1]> get_param



}

#endif
    
