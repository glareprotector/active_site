#ifndef CPP_CALLER_H
#define CPP_CALLER_H

#include <Python.h>
#include <cstring>
#include <string>

// this function will allow one to call module functions(not member functions of any object)
// will need to specify the module name, function name.

class cpp_caller{

 public:
  
  static PyObject* call_PyFunc(string module_name, string fxn_name, PyObject* pArgs){
    // load module
    PyObject* pFunc = get_module_PyObject(module_name_char, fxn_name_char);
    pValue = PyObject_CallObject(pFunc, pArgs);
    return pValue;
  }

  // this should be the only place a module is imported.  returns borrowed reference
  static PyObject* get_module_PyObject(string module_name, string obj_name){
    
    // convert arguments to char*
    char* module_name_char = new char[1000];
    char* obj_name_char = new char[1000];
    strcpy(module_name_char, module_name.c_str());
    strcpy(obj_name_char, obj_name.c_str());

    PyObject *pName, *pModule, *pDict, *pObj;
    pName = PyString_FromString(module_name_char);
    // if the module is already in sys.modules, don't import again
    PyObject* pModulesDict = PyImport_GetModuleDict();
    if(PyDict_Contains(pModulesDict, pName)){
      pModule = PyImport_Import(pName);
      pDict = PyModule_GetDict(pModule); // borrowed reference
      pObj = PyDict_GetItemString(pDict, obj_name_char); // borrowed reference
    }
    Py_DECREF(pName);
    delete[] module_name_char;
    delete[] obj_name_char;
    return pObj;
  }

  // converts python list of float to c++ array of nums
  static arbi_array<num> list_to_floats(PyObject* pList){
    int len = PyList_Size(pList);
    arbi_array<num> result(1, len);
    for(int i = 0; i < len; i++){
      result(i) = (num)PyFloat_AsDouble(PyList_GetItem(pList, i));
    }
    return result;
  }  

  // converts python list of ints to c++ array of ints
  static arbi_array<int> list_to_ints(PyObject* pList){
    int len = PyList_Size(pList);
    arbi_array<int> result(1, len);
    for(int i = 0; i < len; i++){
      result(i) = PyInt_AsLong(PyList_GetItem(pList, i));
    }
    return result;
  }
  
};

class cached_obj_getter: public cpp_caller{

  static void add_string_to_params(PyObject* pParams, string param_name, string param_val){
    PyObject *pParamName, *pParamVal;
    char* char_temp = new char[1000];
    str_cpy(char_temp, param_name.c_str());
    pParamName = PyString_FromString(char_temp);
    str_cpy(char_temp, param_val.c_str());
    pParamVal = PyString_FromString(char_temp);
    PyDict_SetItem(pParams, pParamName, pParamVal);
    Py_DECREF(pParamName);
    Py_DECREF(pParamVal);
  }

  static void add_int_to_params(PyObject* pParams, string param_name, int param_val){
    PyObject *pParamName, *pParamVal;
    char* char_temp = new char[1000];
    str_cpy(char_temp, param_name.c_str());
    pParamName = PyString_FromString(char_temp);
    str_cpy(char_temp, param_val.c_str());
    pParamVal = PyInt_FromLong(param_val);
    PyDict_SetItem(pParams, pParamName, pParamVal);
    Py_DECREF(pParamName);
    Py_DECREF(pParamVal);
  }

}

// for every site, need to retrieve list of other sites sorted by distance as well as the distances
// input would be pdb_name, chain_letter, site number.  input to python function is a dictionary of arguments
class sorted_distances_getter: public cached_obj_getter{

  static void get(string pdb_name, string chain_letter, int aa, arbi_array<int>& positions, arbi_array<num>& distances){
    
    string wrapper_name = string("pdb_chain_site_sorted_distances_obj_wrapper");
    PyObject* pWrapper = get_module_PyObject(string("objects"), wrapper_name); // borrowed reference
    // make params dict
    PyObject* pParams = PyDict_New();
    add_string_to_params(pParams, string("pdb_name"), pdb_name);
    add_string_to_params(pParams, string("chain_letter"), chain_letter);
    add_int_to_params(pParams, string("aa"), aa);
    // make tuple of arguments
    PyObject* pArgs = PyTuple_New(2);
    PyTuple_SetItem(pArgs, 0, pWrapper);
    PyTuple_SetItem(pArgs, 1, pParams);
    PyObject* pResult = call_PyFunc(string("global_stuff"), string("get_object"), pArgs);
    DECREF(pParams);
    DECREF(pArgs);
    
    // prase the returned python object into c++ format
    PyObject* pPositions = PyTuple_GetItem(pResult, 0); // borrowed reference
    PyObject* pDistances = PyTuple_GetItem(pResult, 1); // borrowed reference
    positions = list_as_ints(pPositions);
    distances = list_as_floats(pDistances);
  }  

};				    


#endif
