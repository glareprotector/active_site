#ifndef CPP_CALLER_H
#define CPP_CALLER_H

#include <Python.h>
//#include <cstring>
#include <string>
//#include <string.h>
//#include "nums.h"
#include <set>
#include "globals.h"

typedef double num;
using namespace std;

// this function will allow one to call module functions(not member functions of any object)
// will need to specify the module name, function name.

class cpp_caller{
 
 public:

  static set<string> added_paths;

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
    if(pModule == NULL){
      cout<<"no module "<<module_name<<" "<<obj_name<<endl;
    }
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


  // returns attribute of object
  static PyObject* get_object_attr(PyObject* pObj, string attr_name){
    PyObject* pAttr_name = PyString_FromCPPString(attr_name);
    PyObject* pResult = PyObject_GetAttr(pObj, pAttr_name);
    Py_DECREF(pAttr_name);
    return pResult;
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
      pObj = cpp_caller::cpp_string_vect_to_py_string_list(*(arbi_array<string1d>*)val_p);
      break;
    case globals::NUM_VECT:
      pObj = cpp_caller::cpp_num_vect_to_py_float_list(*(arbi_array<num1d>*)val_p);
      break;
    case globals::STRING_MAT:
      pObj = cpp_caller::CPPString_mat_to_py_string_mat(*(arbi_array<string2d>*)val_p);
      break;
    case globals::INT_VECT:
      pObj = cpp_caller::cpp_int_vect_to_py_int_list(*(arbi_array<int1d>*)val_p);
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
    Py_DECREF(pMethodName);
    Py_DECREF(pKey);
    return pResult;
  }
    
  
  // converts python list of float to c++ array of nums
  static arbi_array<num1d> py_float_list_to_cpp_num_vect(PyObject* pList, bool decref = false){
    if(pList == NULL) throw 20;
    int len = PyList_Size(pList);
    arbi_array<num1d> result(len);
    for(int i = 0; i < len; i++){
      result(i) = (num)PyFloat_AsDouble(PyList_GetItem(pList, i));
    }
    if(decref) Py_DECREF(pList);
    return result;
  }

  // convert c++ array of nums to python list of floats
  static PyObject* cpp_num_vect_to_py_float_list(arbi_array<num1d> x){
    int size = x.size().i0;
    PyObject* pX = PyList_New(size);
    for(int i = 0; i < size; i++){
      PyList_SetItem(pX, i, PyFloat_FromDouble(x(i)));
    }
    return pX;
  }


  // converts python list of ints to c++ array of ints
  static arbi_array<int1d> py_int_list_to_cpp_int_vect(PyObject* pList, bool decref = false){
    if(pList == NULL) throw 20;
    int len = PyList_Size(pList);
    arbi_array<int1d> result; result.resize(len);
    for(int i = 0; i < len; i++){
      result(i) = PyInt_AsLong(PyList_GetItem(pList, i));
    }
    if(decref) Py_DECREF(pList);
    return result;
  }

  // convert c++ array of nums to python list of floats
  static PyObject* cpp_int_vect_to_py_int_list(arbi_array<int1d> x){
    int size = x.size().i0;
    PyObject* pX = PyList_New(size);
    for(int i = 0; i < size; i++){
      PyList_SetItem(pX, i, PyInt_FromLong(x(i)));
    }
    return pX;
  }

  // convert python float list to c++ list of strings
  static arbi_array<string1d> py_string_list_to_cpp_string_vect(PyObject* pList, bool decref = false){
    if(pList == NULL) throw 20;
    int len = PyList_Size(pList);
    arbi_array<string1d> result(len);
    for(int i = 0; i < len; i++){
      result(i) = cpp_caller::CPPString_From_PyString(PyList_GetItem(pList, i));
    }
    if(decref) Py_DECREF(pList);
    return result;
  }

  // convert c++ list of strings to python float list
  static PyObject* cpp_string_vect_to_py_string_list(arbi_array<string1d> x){
    int size = x.size().i0;
    PyObject* pList = PyList_New(size);
    for(int i = 0; i < size; i++){
      PyList_SetItem(pList, i, PyString_FromCPPString(x(i)));
    }
    return pList;
  }

  // prints a python object
  static void py_print(PyObject* pObj){
    PyObject* pFunc = cpp_caller::get_module_PyObject(string("new_new_objects"), string("print_stuff"));
    PyObject_Print(pObj, stdout, 0);
    PyObject_Print(pFunc, stdout, 0);
    PyObject* pResult = PyObject_CallFunctionObjArgs(pFunc, pObj, NULL);
    Py_DECREF(pResult);
  }

  // convert python mat of strings to c++ mat of strings
  static arbi_array<string2d> py_string_mat_to_cpp_string(PyObject* pMat, bool decref = false){
    if(pMat == NULL) throw 20;
    int size0 = PyList_Size(pMat);
    PyObject* pRow = PyList_GetItem(pMat, 0);
    int size1 = PyList_Size(pRow);
    arbi_array<string2d> xx(size0, size1);
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
  static PyObject* CPPString_mat_to_py_string_mat(arbi_array<string2d> x){
    PyObject* pMat = PyList_New(x.size().i0);
    for(int i = 0; i < x.size().i0; i++){
      PyObject* pVect = PyList_New(x.size().i1);
      for(int j = 0; j < x.size().i1; j++){
	PyList_SetItem(pVect, j, PyString_FromCPPString(x(i,j)));
      }
      PyList_SetItem(pMat, i, pVect);
    }
    return pMat;
  }


  static PyObject* cpp_num_mat_to_py_float_mat(arbi_array<num2d> x){
    PyObject* pX = PyList_New(x.size().i0);
    for(int i = 0; i < x.size().i0; i++){
      PyObject* pY = PyList_New(x.size().i1);
      for(int j = 0; j < x.size().i1; j++){
	PyList_SetItem(pY, j, PyFloat_FromDouble(x(i,j)));
      }
      PyList_SetItem(pX, i, pY);
    }
    return pX;
  }


  // convert python float matrix(list of lists) to c++ matrix of nums
  static arbi_array<num2d> py_float_mat_to_cpp_num_mat(PyObject* pMat, bool decref = false){
    if(pMat == NULL) {cout<<"NULL!";throw 20;}
    // get dimensions of mat
    int height = PyList_Size(pMat);
    int width = PyList_Size(PyList_GetItem(pMat, 0));
    arbi_array<num2d> x(height, width);
    PyObject* pTemp;
    for(int i = 0; i < height; i++){
      for(int j = 0; j < width; j++){
	 pTemp = PyList_GetItem(PyList_GetItem(pMat, i), j);
	 x(i,j) = PyFloat_AsDouble(pTemp);
      }
    }
    if(decref) Py_DECREF(pMat);
    return x;
  }


  // convert c++ matrix of nums to python matrix(list of lists) of floats
  static PyObject* cpp_int_mat_to_py_int_mat(arbi_array<int2d> x){
    PyObject* pX = PyList_New(x.size().i0);
    for(int i = 0; i < x.size().i0; i++){
      PyObject* pY = PyList_New(x.size().i1);
      for(int j = 0; j < x.size().i1; j++){
	PyList_SetItem(pY, j, PyInt_FromLong(x(i,j)));
      }
      PyList_SetItem(pX, i, pY);
    }
    return pX;
  }  


  static arbi_array<int2d> py_int_mat_to_cpp_int_mat(PyObject* pMat, bool decref = false){
    if(pMat == NULL) throw 20;
    // get dimensions of mat
    int height = PyList_Size(pMat);
    int width = PyList_Size(PyList_GetItem(pMat, 0));
    arbi_array<int2d> x(height, width);
    PyObject* pTemp;
    for(int i = 0; i < height; i++){
      for(int j = 0; j < width; j++){
	 pTemp = PyList_GetItem(PyList_GetItem(pMat, i), j);
	 x(i,j) = PyInt_AsLong(pTemp);
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
    //cout<<"2222"<<endl;
    PyObject_Print(pStr, stdout, 0);
    //cout<<"3333"<<endl;
    return string(PyString_AsString(pStr));
  }
  
  // returns new reference to python pdb_name_struct
  static PyObject* cpp_pdb_name_struct_to_py_pdb_name_struct(pdb_name_struct pdb){
    string pdb_name = pdb.pdb_name;
    string chain_letter = pdb.chain_letter;
    int start = pdb.start;
    int end = pdb.end;
    PyObject* pPdb_name = PyString_FromCPPString(pdb_name);
    PyObject* pChain_letter = PyString_FromCPPString(chain_letter);
    PyObject* pStart = PyInt_FromLong(start);
    PyObject* pEnd = PyInt_FromLong(end);
    // ERROR
    PyObject* pConstructor = cpp_caller::get_module_PyObject(string("cross_validation_pseudo"), string("pdb_name_struct"));
    PyObject* pResult = PyObject_CallFunctionObjArgs(pConstructor, pPdb_name, pChain_letter, pStart, pEnd, NULL);
    Py_DECREF(pConstructor);
    Py_DECREF(pPdb_name);
    Py_DECREF(pChain_letter);
    Py_DECREF(pStart);
    Py_DECREF(pEnd);
    return pResult;
  }

  // need to convert python pdb_name_struct list to c++ version.  don't need to do reverse but oh well
  static pdb_name_struct py_pdb_name_struct_to_cpp_pdb_name_struct(PyObject* pPdb){

    PyObject* pPdb_name = get_object_attr(pPdb, string("pdb_name"));
    PyObject* pChain_letter = get_object_attr(pPdb, string("chain_letter"));
    PyObject* pStart = get_object_attr(pPdb, string("start"));
    PyObject* pEnd = get_object_attr(pPdb, string("end"));
    pdb_name_struct ans = pdb_name_struct(CPPString_From_PyString(pPdb_name), CPPString_From_PyString(pChain_letter), PyInt_AsLong(pStart), PyInt_AsLong(pEnd));
    Py_DECREF(pPdb_name);
    Py_DECREF(pChain_letter);
    Py_DECREF(pStart);
    Py_DECREF(pEnd);
    return ans;
  }

  static arbi_array<pdb_name_struct[1]> py_pdb_name_struct_list_to_cpp_pdb_name_struct_list(PyObject* pList, bool decref=false){
    //cout<<"INSIDE!"<<endl;
    if(pList == NULL) throw 20;
    //cout<<"after"<<endl;
    int len = PyList_Size(pList);
    //cout<<"f"<<endl;
    arbi_array<pdb_name_struct[1]> result(len);
    //cout<<"g"<<endl;
    for(int i = 0; i < len; i++){
      result(i) = cpp_caller::py_pdb_name_struct_to_cpp_pdb_name_struct(PyList_GetItem(pList, i));
    }
    //cout<<"h"<<endl;
    if(decref) Py_DECREF(pList);
    //cout<<"i"<<endl;
    return result;
  }

  
      
  static PyObject* cpp_pdb_name_struct_list_to_py_pdb_name_struct_list(arbi_array<pdb_name_struct[1]> x){
    int size = x.size().i0;
    PyObject* pList = PyList_New(size);
    for(int i = 0; i < size; i++){
      PyList_SetItem(pList, i, cpp_pdb_name_struct_to_py_pdb_name_struct(x(i)));
    }
    return pList;
  }


  static PyObject* cpp_pdb_results_struct_to_py_pdb_results_struct(pdb_results_struct res){
    PyObject* pScores = cpp_num_vect_to_py_float_list(res.scores);
    PyObject* pTrue_classes = cpp_int_vect_to_py_int_list(res.true_classes);
    PyObject* pPdb_structs = cpp_pdb_name_struct_list_to_py_pdb_name_struct_list(res.pdb_structs);
    PyObject* pSample_lengths = cpp_int_vect_to_py_int_list(res.sample_lengths);
    PyObject* pConstructor = cpp_caller::get_module_PyObject(string("cross_validation_pseudo"), string("pdb_results_struct"));
    PyObject* pResult = PyObject_CallFunctionObjArgs(pConstructor, pScores, pTrue_classes, pPdb_structs, pSample_lengths, NULL);
    Py_DECREF(pConstructor);
    return pResult;
  }
};


//set<string> cpp_caller::added_paths;

class cached_obj_getter: public cpp_caller{
 public:
  
  // takes in param object that is already set and calls the specified wrapper_instance 
  static PyObject* call_wrapper(string wrapper_module, string wrapper_name, PyObject* pParams, bool recalculate, bool to_pickle, bool to_filelize, bool always_recalculate = false){

    // convert the bools to python bools
    PyObject* pRecalculate;
    PyObject* pTo_Pickle;
    PyObject* pTo_Filelize;
    PyObject* pAlways_Recalculate;

    if(recalculate){ pRecalculate = Py_True;} else{ pRecalculate = Py_False;}
    if(to_pickle){ pTo_Pickle = Py_True;} else{ pTo_Pickle = Py_False;}
    if(to_filelize){ pTo_Filelize = Py_True;} else{ pTo_Filelize = Py_False;}
    if(always_recalculate){ pAlways_Recalculate = Py_True;} else{ pAlways_Recalculate = Py_False;}

    // get the wrapper class reference(not an instance)
    PyObject* pWrapper = get_module_PyObject(wrapper_module, wrapper_name);
    // get wc.get_stuff reference
    PyObject* pGetStuff = get_module_PyObject(string("wc"), string("get_stuff")); 
    // call method
    PyObject* pResult = PyObject_CallFunctionObjArgs(pGetStuff, pWrapper, pParams, pRecalculate, pTo_Pickle, pTo_Filelize, pAlways_Recalculate, NULL);
    if(pResult == NULL){
      py_print(pParams);
      cout<<"BAD"<<wrapper_name<<endl;
    }
    return pResult;
  }

};

#endif

#include "cpp_param.h"
