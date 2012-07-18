#ifndef CPP_CALLER_H
#define CPP_CALLER_H

#include <Python.h>
//#include <cstring>
#include <string>
//#include <string.h>
#include "nums.h"
#include <set>

typedef double num;
using namespace std;

// this function will allow one to call module functions(not member functions of any object)
// will need to specify the module name, function name.

class cpp_caller{
 
 public:

  static set<string> added_paths;

  /*static void init(){
    // add to added_paths everything initially in sys.path
    PyObject* pSysPath= _get_module_PyObject(string("sys"), string("path"));
    for(int i = 0; i < PyList_Size(pSysPath); i++){
      added_paths.insert(string(PyString_AsString(PyList_GetItem(pSysPath, i))));
    }
    }*/

  static PyObject* call_PyFunc(string module_name, string fxn_name, PyObject* pArgs){
    // load module
    PyObject* pFunc = get_module_PyObject(module_name, fxn_name);
    PyObject* pValue = PyObject_CallObject(pFunc, pArgs);
    return pValue;
    
  }

  // imports module, assuming it is in path.  returns new reference
  static PyObject* _get_module(string module_name){
    char* module_name_char = new char[1000];
    strcpy(module_name_char, module_name.c_str());
    PyObject* pName = PyString_FromString(module_name_char); 
    PyObject* pModule = PyImport_Import(pName);
    Py_DECREF(pName);
    delete[] module_name_char;
    return pModule;
  }

  // gets object from module, assuming path is correctly set.  returns borrowed reference
  static PyObject* _get_module_PyObject(string module_name, string obj_name){
    // convert arguments to char*
    char* obj_name_char = new char[1000];
    strcpy(obj_name_char, obj_name.c_str());
    PyObject *pModule, *pDict, *pObj;
    pModule = _get_module(module_name); // new
    pDict = PyModule_GetDict(pModule); // borrowed
    pObj = PyDict_GetItemString(pDict, obj_name_char); // borrowed
    Py_DECREF(pModule);
    delete[] obj_name_char;
    return pObj;
  }

  // adds module path to sys.path if needed.  returns borrowed reference
  static PyObject* get_module_PyObject(string module_name, string obj_name, string module_path = string(".")){
   
    // if module path is not in paths, add it
    if(added_paths.find(module_path) == added_paths.end()){
      // get reference to sys.path
      PyObject* pSysPath = _get_module_PyObject(string("sys"), string("path"));
      PyObject* pPathString = PyString_FromString(module_path.c_str());
      PyList_Insert(pSysPath, 0, pPathString);
      Py_DECREF(pPathString);
    }
    PyObject* pObj = _get_module_PyObject(module_name, obj_name);
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
 public:
  static void add_string_to_params(PyObject* pParams, string param_name, string param_val){
    PyObject *pParamName, *pParamVal;
    char* char_temp = new char[1000];
    strcpy(char_temp, param_name.c_str());
    pParamName = PyString_FromString(char_temp);
    strcpy(char_temp, param_val.c_str());
    pParamVal = PyString_FromString(char_temp);
    PyDict_SetItem(pParams, pParamName, pParamVal);
    Py_DECREF(pParamName);
    Py_DECREF(pParamVal);
    delete[] char_temp;
  }

  static void add_int_to_params(PyObject* pParams, string param_name, int param_val){
    PyObject *pParamName, *pParamVal;
    char* char_temp = new char[1000];
    strcpy(char_temp, param_name.c_str());
    pParamName = PyString_FromString(char_temp);
    pParamVal = PyInt_FromLong(param_val);
    PyDict_SetItem(pParams, pParamName, pParamVal);
    Py_DECREF(pParamName);
    Py_DECREF(pParamVal);
    delete[] char_temp;
  }

};

// for every site, need to retrieve list of other sites sorted by distance as well as the distances
// input would be pdb_name, chain_letter, site number.  input to python function is a dictionary of arguments
class sorted_distances_getter: public cached_obj_getter{
  
 public:

  static void get(string pdb_name, string chain_letter, int aa, arbi_array<int>& positions, arbi_array<num>& distances){
    
    string wrapper_name = string("pdb_chain_site_sorted_distances_obj_wrapper");
    PyObject* pWrapper = get_module_PyObject(string("objects"), wrapper_name); // borrowed reference
    // increment pWrapper ref count to make it seem like a new reference
    Py_INCREF(pWrapper);
    // make params dict
    PyObject* pParamDict = PyDict_New();
    add_string_to_params(pParamDict, string("pdb_name"), pdb_name);
    add_string_to_params(pParamDict, string("chain_letter"), chain_letter);
    add_int_to_params(pParamDict, string("aa"), aa);

    // make param object out of the dict.  first make argument list
    PyObject* pParamConstructorArgs = PyTuple_New(1); // new
    PyObject* pParamConstructor = get_module_PyObject(string("objects"), string("param")); // borrowed
    PyTuple_SetItem(pParamConstructorArgs, 0, pParamDict); // reference to pParamDict stolen by pParamConstructorArgs
    PyObject* pParams = call_PyFunc(string("objects"), string("param"), pParamConstructorArgs);
    Py_DECREF(pParamConstructorArgs);

    // make tuple of arguments
    PyObject* pArgs = PyTuple_New(2);
    PyTuple_SetItem(pArgs, 0, pWrapper);
    PyTuple_SetItem(pArgs, 1, pParams);
    PyObject* pResult = call_PyFunc(string("global_stuff"), string("get_object"), pArgs);
    Py_DECREF(pArgs);
    
    // prase the returned python object into c++ format
    PyObject* pDistances = PyTuple_GetItem(pResult, 0); // borrowed reference
    PyObject* pPositions = PyTuple_GetItem(pResult, 1); // borrowed reference
    Py_DECREF(pResult);
    positions = list_to_ints(pPositions);
    distances = list_to_floats(pDistances);
  }  

};				    


#endif
