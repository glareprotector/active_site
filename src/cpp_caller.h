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


  // returns new empty param object
  static PyObject* get_new_empty_param(){
    PyObject* pParamDict = PyDict_New();
    PyObject* pMakeNewParamArgs = PyTuple_New(1);
    PyTuple_SetItem(pMakeNewParamArgs, 0, pParamDict);
    PyObject* pParams = call_PyFunc(string("param"), string("param"), pMakeNewParamArgs);
    Py_DECREF(pMakeNewParamArgs);
    return pParams;
  }

  // sets param key to val, with the type(in python) of val indicated
  static void set_param(PyObject* pParams, string key, void* val_p, int type){
    PyObject* pKey = PyString_FromCPPString(key);
    PyObject* pObj;
    switch(type){
    case globals::INT_TYPE:
      pObj = PyInt_FromLong(*((int*)val_p));
      break;
    case globals::NUM_TYPE:
      pObj = PyFloat_FromDouble(*(num*)val_p);
      break;
    case globals::STRING_TYPE:
      pObj = PyString_FromCPPString(*(string*)val_p);
      break;
    case globals::STRING_VECT:
      pObj = cpp_caller::cpp_string_vect_to_py_string_list(*(arbi_array<string>*)val_p);
      break;
    case globals::NUM_VECT:
      pObj = cpp_caller::cpp_num_vect_to_py_float_list(*(arbi_array<num>*)val_p);
      break;
    case globals::STRING_MAT:
      pObj = cpp_caller::CPPString_mat_to_py_string_mat(*(arbi_array<string>*)val_p);
      break;
    }


    PyObject* pMethodName = PyString_FromCPPString(string("set_param"));
    PyObject* pResult = PyObject_CallMethodObjArgs(pParams, pMethodName, pKey, pObj, NULL);
    Py_DECREF(pKey);
    Py_DECREF(pObj);
    Py_DECREF(pMethodName);
    Py_DECREF(pResult);
  }

  // retrieves param from pParams of given key
  static PyObject* get_param(PyObject* pParams, string key){
    
    PyObject* pMethodName = PyString_FromCPPString(string("get_param"));
    PyObject* pKey = cpp_caller::PyString_FromCPPString(key);
    PyObject* pResult = PyObject_CallMethodObjArgs(pParams, pMethodName, pKey, NULL);
    PyObject_Print(pKey, stdout, 0);
    PyObject_Print(pResult, stdout, 0);
    py_print(pKey);
    py_print(pParams);
    py_print(pResult);
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    return pResult;
  }
    
  
  // converts python list of float to c++ array of nums
  static arbi_array<num> py_float_list_to_cpp_num_vect(PyObject* pList, bool decref = false){
    int len = PyList_Size(pList);
    arbi_array<num> result(1, len);
    for(int i = 0; i < len; i++){
      result(i) = (num)PyFloat_AsDouble(PyList_GetItem(pList, i));
    }
    if(decref) Py_DECREF(pList);
    return result;
  }

  // convert c++ array of nums to python list of floats
  static PyObject* cpp_num_vect_to_py_float_list(arbi_array<num> x){
    int size = x.size(0);
    PyObject* pX = PyList_New(size);
    for(int i = 0; i < size; i++){
      PyList_SetItem(pX, i, PyFloat_FromDouble(x(i)));
    }
    return pX;
  }


  // converts python list of ints to c++ array of ints
  static arbi_array<int> py_int_list_to_cpp_int_vect(PyObject* pList, bool decref = false){
    int len = PyList_Size(pList);
    arbi_array<int> result(1, len);
    for(int i = 0; i < len; i++){
      result(i) = PyInt_AsLong(PyList_GetItem(pList, i));
    }
    if(decref) Py_DECREF(pList);
    return result;
  }

  // convert c++ array of nums to python list of floats
  static PyObject* cpp_int_vect_to_py_int_list(arbi_array<int> x){
    int size = x.size(0);
    PyObject* pX = PyList_New(size);
    for(int i = 0; i < size; i++){
      PyList_SetItem(pX, i, PyInt_FromLong(x(i)));
    }
    return pX;
  }

  // convert python float list to c++ list of strings
  static arbi_array<string> py_string_list_to_cpp_string_vect(PyObject* pList, bool decref = false){
    int len = PyList_Size(pList);
    arbi_array<string> result(1, len);
    for(int i = 0; i < len; i++){
      result(i) = cpp_caller::CPPString_From_PyString(PyList_GetItem(pList, i));
    }
    if(decref) Py_DECREF(pList);
    return result;
  }

  // convert c++ list of strings to python float list
  static PyObject* cpp_string_vect_to_py_string_list(arbi_array<string> x){
    int size = x.size(0);
    PyObject* pList = PyList_New(size);
    for(int i = 0; i < size; i++){
      PyList_SetItem(pList, i, PyString_FromCPPString(x(i)));
    }
    return pList;
  }

  // prints a python object
  static void py_print(PyObject* pObj){
    PyObject* pFunc = cpp_caller::get_module_PyObject(string("new_new_objects"), string("print_stuff"));
    cout<<"                        OHOHOHOHO"<<endl;
    PyObject_Print(pObj, stdout, 0);
    cout<<"                        8989898989"<<endl;
    PyObject_Print(pFunc, stdout, 0);
    PyObject* pResult = PyObject_CallFunctionObjArgs(pFunc, pObj, NULL);
    Py_DECREF(pResult);
    //Py_DECREF(pFunc);
  }

  // convert python mat of strings to c++ mat of strings
  static arbi_array<string> py_string_mat_to_cpp_string(PyObject* pMat, bool decref = false){
    int size0 = PyList_Size(pMat);
    cout<<"inside string mat to cpp mat"<<endl;
    cpp_caller::py_print(pMat);
    PyObject* pRow = PyList_GetItem(pMat, 0);
    int size1 = PyList_Size(pRow);
    arbi_array<string> xx(2, size0, size1);
    for(int a = 0; a < size0; a++){
      PyObject* pRow = PyList_GetItem(pMat, a);
      for(int j = 0; j < size1; j++){
	xx(a,j) = CPPString_From_PyString(PyList_GetItem(pRow, j));
      }
    }
    if(decref) Py_DECREF(pMat);
    return xx;
  }

  // convert c++ mat of strings to python mat of strings
  static PyObject* CPPString_mat_to_py_string_mat(arbi_array<string> x){
    PyObject* pMat = PyList_New(x.size(0));
    for(int i = 0; i < x.size(0); i++){
      PyObject* pVect = PyList_New(x.size(1));
      for(int j = 0; j < x.size(1); j++){
	PyList_SetItem(pVect, j, PyString_FromCPPString(x(i,j)));
      }
      PyList_SetItem(pMat, i, pVect);
    }
    return pMat;
  }

  // convert python float matrix(list of lists) to c++ matrix of nums
  static arbi_array<num> py_float_mat_to_cpp_num_mat(PyObject* pMat, bool decref = false){
    // get dimensions of mat
    int height = PyList_Size(pMat);
    int width = PyList_Size(PyList_GetItem(pMat, 0));
    arbi_array<num> x(2, height, width);
    PyObject* pTemp;
    for(int i = 0; i < height; i++){
      for(int j = 0; j < width; j++){
	 pTemp = PyList_GetItem(PyList_GetItem(pMat, i), j);
	 x(i,j) = PyFloat_AsDouble(pTemp);
	 //Py_DECREF(pTemp);
      }
    }
    if(decref) Py_DECREF(pMat);
    return x;
  }


  // convert c++ matrix of nums to python matrix(list of lists) of floats
  static PyObject* cpp_int_mat_to_py_int_mat(arbi_array<num> x){
    PyObject* pX = PyList_New(x.size(0));
    for(int i = 0; i < x.size(0); i++){
      PyObject* pY = PyList_New(x.size(1));
      for(int j = 0; j < x.size(1); j++){
	PyList_SetItem(pY, j, PyInt_FromLong(x(i,j)));
      }
      PyList_SetItem(pX, i, pY);
    }
    return pX;
  }  


  static arbi_array<int> py_int_mat_to_cpp_int_mat(PyObject* pMat, bool decref = false){
    // get dimensions of mat
    int height = PyList_Size(pMat);
    int width = PyList_Size(PyList_GetItem(pMat, 0));
    arbi_array<int> x(2, height, width);
    PyObject* pTemp;
    for(int i = 0; i < height; i++){
      for(int j = 0; j < width; j++){
	 pTemp = PyList_GetItem(PyList_GetItem(pMat, i), j);
	 x(i,j) = PyInt_AsLong(pTemp);
	 //Py_DECREF(pTemp);
      }
    }
    if(decref) Py_DECREF(pMat);
    return x;
  }

  
  // returns new reference to python string
  static PyObject* PyString_FromCPPString(string s){
    PyObject* pString;
    char* char_temp = new char[1000];
    strcpy(char_temp, s.c_str());
    pString = PyString_FromString(char_temp);
    delete[] char_temp;
    return pString;
  }
  

  // converts python string to c++ string
  static string CPPString_From_PyString(PyObject* pStr){
    cout<<"2222"<<endl;
    PyObject_Print(pStr, stdout, 0);
    cout<<"3333"<<endl;
    return string(PyString_AsString(pStr));
  }
  

};

class cached_obj_getter: public cpp_caller{
 public:
  
  // takes in param object that is already set and calls the specified wrapper_instance 
    static PyObject* call_wrapper(string wrapper_module, string wrapper_instance_name, PyObject* pParams, bool recalculate, bool to_pickle, bool to_filelize){



    // convert the bools to python bools
    PyObject* pRecalculate;
    PyObject* pTo_Pickle;
    PyObject* pTo_Filelize;

    py_print(pParams);
    if(recalculate){ pRecalculate = Py_True;} else{ pRecalculate = Py_False;}
    if(to_pickle){ pTo_Pickle = Py_True;} else{ pTo_Pickle = Py_False;}
    if(to_filelize){ pTo_Filelize = Py_True;} else{ pTo_Filelize = Py_False;}

    // get the wrapper instance
    PyObject* pWrapperInstance = get_module_PyObject(wrapper_module, wrapper_instance_name);

    // get method name
    PyObject* pMethodName = PyString_FromCPPString(string("constructor"));
    py_print(pWrapperInstance);
    py_print(pMethodName);

    // call method
    PyObject* pResult = PyObject_CallMethodObjArgs(pWrapperInstance, pMethodName, pParams, pRecalculate, pTo_Pickle, pTo_Filelize, NULL);
    py_print(pResult);
    PyObject* pVal = PyTuple_GetItem(pResult, 2);
    Py_INCREF(pVal);
    Py_DECREF(pResult);
    cout<<"                inside call_wrapper"<<endl;
    py_print(pVal);
    return pVal;
  }

  static PyObject* _get(string wrapper_name, arbi_array<string> param_names, arbi_array<void*> values, arbi_array<int> types){
    string wrapper_instance_name = string("pdb_chain_site_sorted_distances_obj_wrapper");
    PyObject* pWrapper = get_module_PyObject(string("objects"), wrapper_instance_name); // borrowed reference
    // increment pWrapper ref count to make it seem like a new reference.  will feed
    Py_INCREF(pWrapper);
    
    //
    
    
    // make param object out of the dict.  first make argument list
    PyObject* pParamDict = PyDict_New();
    PyObject* pParamConstructorArgs = PyTuple_New(1); // new
    PyTuple_SetItem(pParamConstructorArgs, 0, pParamDict); // reference to pParamDict stolen by pParamConstructorArgs
    PyObject* pParams = call_PyFunc(string("objects"), string("param"), pParamConstructorArgs);
    Py_DECREF(pParamConstructorArgs);


    int num_params = param_names.size(0);
    for(int i = 0; i < num_params; i++){
      set_param(pParams, param_names(i), values(i), types(i));
    }
    
    // make tuple of arguments
    PyObject* pArgs = PyTuple_New(2);
    PyTuple_SetItem(pArgs, 0, pWrapper);
    PyTuple_SetItem(pArgs, 1, pParams);
    PyObject* pResult = call_PyFunc(string("global_stuff"), string("get_object"), pArgs);
    Py_DECREF(pArgs);

    // pResult is a tuple of 3 things, but only want the first
    PyObject* pVal = PyTuple_GetItem(pResult, 2);
    Py_INCREF(pVal);
    Py_DECREF(pResult);
    cout<<"                              OOOOOOOOOOOOOOOOOOOOOOOO"<<endl;
    py_print(pVal);
    return pVal;
  }
};

class crf_info_getter: public cached_obj_getter{

 public:

  static void get(string pdb_name, string chain_letter, arbi_array<num>& node_features, arbi_array<num>& edge_features, arbi_array<int>& edges, arbi_array<int>& true_states, int& num_nodes, int& num_edges){

    string wrapper_name = string("pdb_chain_data_obj_wrapper");
      
    int num_params = 2;

    arbi_array<string> param_names(1, num_params);
    arbi_array<void*> values(1, num_params);
    arbi_array<int> types(1, num_params);

    param_names(0) = string("pdb_name");
    values(0) = &pdb_name;
    types(0) = globals::STRING_TYPE;

    param_names(1) = string("chain_letter");
    values(1) = &chain_letter;
    types(0) = globals::STRING_TYPE;

    PyObject* pResult = _get(wrapper_name, param_names, values, types);

    // parse result
    PyObject* pNode_Features = PyTuple_GetItem(pResult, 0); // borrowed reference
    PyObject* pEdge_Features = PyTuple_GetItem(pResult, 1); // borrowed reference
    PyObject* pEdge_List = PyTuple_GetItem(pResult, 2); // borrowed reference
    PyObject* pTrue_States = PyTuple_GetItem(pResult, 3); // borrowed reference
    PyObject* pNum_Nodes = PyTuple_GetItem(pResult, 4); // borrowed reference
    PyObject* pNum_Edges = PyTuple_GetItem(pResult, 5); // borrowed reference
    node_features = py_float_list_to_cpp_num_vect(pNode_Features);
    edge_features = py_float_list_to_cpp_num_vect(pEdge_Features);
    edges = py_int_list_to_cpp_int_vect(pEdge_List);
    true_states = py_int_list_to_cpp_int_vect(pTrue_States);
    num_nodes = PyInt_AsLong(pNum_Nodes);
    num_edges = PyInt_AsLong(pNum_Edges);

    Py_DECREF(pResult);
    
  }
  };
  
// for every site, need to retrieve list of other sites sorted by distance as well as the distances
// input would be pdb_name, chain_letter, site number.  input to python function is a dictionary of arguments
/*
class sorted_distances_getter: public cached_obj_getter{
  
 public:
  
  static void get(string pdb_name, string chain_letter, int aa, arbi_array<int>& positions, arbi_array<num>& distances){
    
    string wrapper_name = string("pdb_chain_site_sorted_distances_obj_wrapper");
    
    int num_params = 3;
    arbi_array<string> pdb_names(1, num_params);
    arbi_array<<void*> values(1, num_params);
    arbi_array<int> types(1, num_params);
    
    param_names(0) = string("pdb_name");
    values(0) = pdb_name;
    types(0) = cpp_caller::STRING_TYPE;
    
    param_names(1) = string("chain_letter");
    values(1) = chain_letter;
    types(1) = cpp_caller::STRING_TYPE;
    
    param_names(2) = string("aa");
    values(2) = aa;
    types(2) = cpp_caller::INT_TYPE;
    
    PyObject* pResult = _get(wrapper_name, param_names, values, types);
    
    // convert results back to c++ objects
    PyObject* pDistances = PyTuple_GetItem(pResult, 0); // borrowed reference
    PyObject* pPositions = PyTuple_GetItem(pResult, 1); // borrowed reference
    positions = py_list_to_cpp_ints(pPositions);
    distances = py_list_to_cpp_floats(pDistances);
    
    Py_DECREF(pResult);
  }
  
};				    
*/
class experiment_results_file_getter: public cached_obj_getter{
    
 public:
  static void get(arbi_array<num> scores, arbi_array<int> sizes, arbi_array<string> pdb_names, arbi_array<string> chain_letters){
    // convert stuff to write to PyObjects and add to params
    cpp_caller::set_param(globals::pParams, string("scores"), &scores, globals::NUM_VECT);
    cpp_caller::set_param(globals::pParams, string("sizes"), &sizes, globals::INT_VECT);
    cpp_caller::set_param(globals::pParams, string("pdb_names"), &pdb_names, globals::STRING_VECT);
    cpp_caller::set_param(globals::pParams, string("chain_letters"), &chain_letters, globals::STRING_VECT);
    PyObject* pResults = cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_experiment_results_file_w"), globals::pParams, true, true, true);
    Py_DECREF(pResults);
  }
};

class roc_curve_plot_file_getter: public cached_obj_getter{
 public:
  static void get(){
    PyObject* pResults = cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_roc_curve_plot_file_w"), globals::pParams, true, false, false);
    Py_DECREF(pResults);
  }
};

class node_features_obj_getter: public cached_obj_getter{
 public:
  static arbi_array<num> get(string pdb_name, string chain_letter, bool recalculate=false){
    cpp_caller::set_param(globals::pParams, string("pdb_name"), &pdb_name, globals::STRING_TYPE);
    cpp_caller::set_param(globals::pParams, string("chain_letter"), &chain_letter, globals::STRING_TYPE);
    PyObject* pResults = cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_node_features_obj_w"), globals::pParams, recalculate, true, true);
    // convert results to arbi_array
    arbi_array<num> results = cpp_caller::py_float_mat_to_cpp_num_mat(pResults);
    Py_DECREF(pResults);
  }
};

class edge_features_obj_getter: public cached_obj_getter{
 public:
  static arbi_array<num> get(string pdb_name, string chain_letter, bool recalculate=false){
    cpp_caller::set_param(globals::pParams, string("pdb_name"), &pdb_name, globals::STRING_TYPE);
    cpp_caller::set_param(globals::pParams, string("chain_letter"), &chain_letter, globals::STRING_TYPE);
    PyObject* pResults = cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_edge_features_obj_w"), globals::pParams, recalculate, true, true);
    // convert results to arbi_array
    arbi_array<num> results = cpp_caller::py_float_mat_to_cpp_num_mat(pResults);
    Py_DECREF(pResults);
  }
};

  
  


#endif


/*

static void add_string_to_params(PyObject* pParams, string param_name, string param_val){
  PyObject *pParamName, *pParamVal;
  pParamName = PyString_FromCPPString(param_name);
  pParamVal = PyString_FromCPPString(char_val);
  PyDict_SetItem(pParams, pParamName, pParamVal);
  Py_DECREF(pParamName);
  Py_DECREF(pParamVal);
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

// ideally would then specify vector of param names, void pointers to objects, types, and for each type, need a way to convert (in c++) to python data struct
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
  positions = list_to_ints(pPositions);
  distances = list_to_floats(pDistances);
  Py_DECREF(pResult);
}  

*/
