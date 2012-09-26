#include "model.h"
#include "LBFGS.h"
#include <vector>
#include "helpers.h"
#include <time.h>
#include <cmath>


#ifndef SERIAL 
#include "mpi.h" 
#endif

sample model::read_sample(){

  // keep for now only bc sample for now takes these in as input
  string pdb_name = cpp_caller::CPPString_From_PyString(cpp_caller::get_param(globals::pParams, string("pdb_name")));
  string chain_letter = cpp_caller::CPPString_From_PyString(cpp_caller::get_param(globals::pParams, string("chain_letter")));

  cout<<"proc_id "<<proc_id<<" is loading "<<pdb_name<<" "<<chain_letter<<endl;
  
  PyObject* pNodeFeatures = cached_obj_getter::call_wrapper(string("new_new_objects"), string("jW"), globals::pParams, globals::recalculate, true, true);
  
  arbi_array<num2d> node_features = cpp_caller::py_float_mat_to_cpp_num_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("bmW"), globals::pParams, globals::recalculate, true, true), true);
  arbi_array<num2d> edge_features = cpp_caller::py_float_mat_to_cpp_num_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("bnW"), globals::pParams, globals::recalculate, true, true), true);  
  arbi_array<int1d> true_states = cpp_caller::py_int_list_to_cpp_int_vect(cached_obj_getter::call_wrapper(string("new_new_objects"), string("oW"), globals::pParams, globals::recalculate, true, true), true);
  arbi_array<int2d> edge_list = cpp_caller::py_int_mat_to_cpp_int_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("iW"), globals::pParams, globals::recalculate, true, true), true);

  this->num_node_features = node_features.size().i1;
  this->num_edge_features = edge_features.size().i1;

  return sample(this, node_features, edge_features, edge_list, true_states, pdb_name, chain_letter);
}

void model::normalize(){

  for(int i = 0; i < num_node_features; i++){

    bool to_normalize;
    //to_normalize = (i > 0) && (i < 4);
    //to_normalize = i >= num_node_features - 5;
    if(true){
      //if(to_normalize){

      // first need to calculate the mean and count
      num sum = 0;
      int count = 0;
      for(int j = 0; j < num_samples; j++){
	for(int k = 0; k < data(j).num_nodes; k++){
	  sum += data(j).node_features(k,i);
	  count++;
	}
      }
      
      #ifdef SERIAL
      
      num mean = sum / (num)count;

      #else

      num sum_sum;
      int count_sum;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(&sum, &sum_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&sum_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&count_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      num mean = sum_sum / (num)count_sum;
      
      #endif

      // now calculate sum of squared distance from mean
      num sqr_dist = 0;
      for(int j = 0; j < num_samples; j++){
	for(int k = 0; k < data(j).num_nodes; k++){
	  sqr_dist += pow( (data(j).node_features(k,i) - mean), 2);
	}
      }

      #ifdef SERIAL

      num variance = sqrt(sqr_dist / (num)count);

      #else

      num sqr_dist_sum;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(&sqr_dist, &sqr_dist_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&sqr_dist_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      num variance = sqrt(sqr_dist_sum / (num)count_sum);
  
     #endif

      // now actually normalize the features
      for(int j = 0; j < num_samples; j++){
	for(int k = 0; k < data(j).num_nodes; k++){
	  data(j).node_features(k,i) = data(j).node_features(k,i) - mean;
	  if(fabs(variance) > 1e-4){
	    data(j).node_features(k,i) /= variance;
	  }
	}
      }
    }
  }
}

// writes to specified folder the scores for each score as well as their true class
void model::report(arbi_array<num1d> theta, int iteration, num obj){



  #ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // arbi_array<int1d> true_classes;
  // arbi_array<num1d> scores;
  int* true_classes;
  num* scores;
  // sample_labels will be of length (4+1) * num_testing, with the 1 for the null characters, 4 for the length of sample_names
  char* sample_pdb_names;
  char* sample_chain_letters;
  int* sample_lengths;

  // need to retrieve the length of sample_pdb_names and sample_chain_letters later
  int num_testing_master;

  // first get total length of stuff to report
  int num_data = 0;
  cout<<num_testing<<" : "<<testing_indicies<<endl;
  for(int i = 0; i < num_testing; i++){
    num_data += data(testing_indicies(i)).num_nodes;
  }
  


  #ifndef SERIAL

  // give the root node the total number of sites
  MPI_Status status;
  int num_data_sum;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&num_data, &num_data_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&num_data_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // give the root node the total number of testing samples
  int num_testing_sum;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&num_testing, &num_testing_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&num_testing_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  num_testing_master = num_testing_sum;
  
  // root process gets larger true_classes and scores and sample name labels and letters
  if(proc_id == 0){
    // true_classes = arbi_array<int1d>; true_classes.resize(num_data_sum);
    // scores = arbi_array<num1d>(num_data_sum);
    true_classes = new int[num_data_sum];
    scores = new num[num_data_sum];
    int sample_pdb_names_length = num_testing_sum * 5;
    sample_pdb_names = new char[sample_pdb_names_length];
    int sample_chain_letters_length = num_testing_sum * 2;
    sample_chain_letters = new char[sample_chain_letters_length];
    sample_lengths = new int[num_testing_sum];
  }
  else{
    // true_classes = arbi_array<int1d>; true_classes.resize(num_data);
    // scores = arbi_array<num1d>(num_data);
    true_classes = new int[num_data];
    scores = new num[num_data];
    int sample_pdb_names_length = num_testing * 5;
    sample_pdb_names = new char[sample_pdb_names_length];
    int sample_chain_letters_length = num_testing * 2;
    sample_chain_letters = new char[sample_chain_letters_length];
    sample_lengths = new int[num_testing];
  }

  #else

  // true_classes.resize(num_data);
  // scores = arbi_array<num1d>(num_data);
  true_classes = new int[num_data];
  scores = new num[num_data];
  int sample_pdb_names_length = num_testing * 5;
  sample_pdb_names = new char[sample_pdb_names_length];
  int sample_chain_letters_length = num_testing * 2;
  sample_chain_letters = new char[sample_chain_letters_length];
  sample_lengths = new int[num_testing];
  num_testing_master = num_testing;

  #endif
  // each node gets list of sizes of its samples
  // put own scores/classes into the vector
  // score is the marginal of the marginal class
  int idx = 0;

  #ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  for(int i = 0; i < num_testing; i++){
    arbi_array<num2d> node_marginals;
    arbi_array<num3d> edge_marginals;
    data(testing_indicies(i)).get_marginals(theta, node_marginals, edge_marginals);
    for(int j = 0; j < data(testing_indicies(i)).num_nodes; j++){
      true_classes[idx] = data(testing_indicies(i)).true_states(j);
      scores[idx] = node_marginals(j,1);
      idx++;
    }
    sample_lengths[i] = data(testing_indicies(i)).num_nodes;
  }
  // likewise, each node gets a char of length 5 * num_testing for pdb_names, char of length 2 * num_testing for chain_letters
  for(int i = 0; i < num_testing; i++){
    strcpy(sample_pdb_names + (i*5), data(testing_indicies(i)).pdb_name.c_str());
    strcpy(sample_chain_letters + (i*2), data(testing_indicies(i)).chain_letter.c_str());
  }
    
     


  // if doing things in parallel, have to add more stuff to the vectors
  
  #ifndef SERIAL



  // each process takes turns sending messages to process 0
  int pos = idx;
  int sample_pos = num_testing;
  for(int i = 1 ; i < num_procs; i++){
    int* class_buf;
    num* score_buf;
    int next_score_size;
    int next_sample_size;
    if(proc_id == i){
      class_buf = true_classes;
      score_buf = scores;
      next_score_size = num_data;
      next_sample_size = num_testing;
      MPI_Send(&next_score_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&next_sample_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(class_buf, next_score_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(score_buf, next_score_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      
      MPI_Send(sample_pdb_names, 5 * next_sample_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
      MPI_Send(sample_chain_letters, 2 * next_sample_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
      MPI_Send(sample_lengths, num_testing, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    else if(proc_id == 0){
      class_buf = true_classes + pos;
      score_buf = scores + pos;
      MPI_Recv(&next_score_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&next_sample_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(class_buf, next_score_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(score_buf, next_score_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(sample_pdb_names + (5 * sample_pos), 5 * next_sample_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(sample_chain_letters + (2 * sample_pos), 2 * next_sample_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(sample_lengths + sample_pos, next_sample_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      pos += next_score_size;
      sample_pos += next_sample_size;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }



  #endif

  // now, parse the long strings into a name/letter for each sample.  only do this for proc_id 0
  // store each as arbi_array of strings
  arbi_array<string1d> parsed_pdb_names(num_testing_master);
  arbi_array<string1d> parsed_chain_letters(num_testing_master);
  arbi_array<int1d> sample_lengths_ar; sample_lengths_ar.resize(num_testing_master);
  if(proc_id == 0){
    for(int i = 0; i < num_testing_master; i++){
      parsed_pdb_names(i) = string(sample_pdb_names + (5*i));
      parsed_chain_letters(i) = string(sample_chain_letters + (2*i));
      sample_lengths_ar(i) = sample_lengths[i];
    }
  }




  cout<<"PROC_ID: "<<proc_id<<endl;

  if(proc_id == 0){


    // convert true_classes to array
    #ifndef SERIAL
    num_data = num_data_sum;
    #endif
    
    arbi_array<int1d> true_classes_ar; true_classes_ar.resize(num_data);
    arbi_array<num1d> scores_ar(num_data);
    for(int i = 0; i < num_data; i++){
      true_classes_ar(i) = true_classes[i];
      scores_ar(i) = scores[i];
    }


    

    cpp_caller::set_param(globals::pParams, string("scores"), &scores_ar, globals::NUM_VECT);
    cpp_caller::set_param(globals::pParams, string("true_states"), &true_classes_ar, globals::INT_VECT);
    cpp_caller::set_param(globals::pParams, string("pdb_names"), &parsed_pdb_names, globals::STRING_VECT);
    cpp_caller::set_param(globals::pParams, string("chain_letters"), &parsed_chain_letters, globals::STRING_VECT);
    cpp_caller::set_param(globals::pParams, string("sizes"), &sample_lengths_ar, globals::INT_VECT);
    // also send the current iteration
    cpp_caller::set_param(globals::pParams, string("iter"), &iteration, globals::INT_TYPE);
    cpp_caller::set_param(globals::pParams, string("obj_val"), &obj, globals::NUM_TYPE);
    // always recalculate this
    cached_obj_getter::call_wrapper(string("new_new_objects"), string("qW"), globals::pParams, true, false, false, true);
    cached_obj_getter::call_wrapper(string("new_new_objects"), string("blW"), globals::pParams, true, false, false, true);
    cached_obj_getter::call_wrapper(string("new_new_objects"), string("ahW"), globals::pParams, true, true, true, true);
  }


  

  #ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

}

void model::assign(int num_folds, int which_fold){
  


  training_indicies = arbi_array<int1d>();
  testing_indicies = arbi_array<int1d>();

  for(int i = 0; i < num_samples; i++){
    if((i % num_folds) == which_fold){
      //training_indicies.append(i);
      helpers::append(training_indicies,i);
    }
    else{
      //testing_indicies.append(i);
      helpers::append(testing_indicies,i);
    }
  }
  cout<<"testing: "<<testing_indicies<<endl;
  cout<<"training: "<<training_indicies<<endl;

  num_training = training_indicies.size().i0;
  num_testing = testing_indicies.size().i0;
}

// future usage: data_list is stored as param.  one by one, set pdb_name, call sample constructor with those params, which will be global

void model::load_data(){

  int num_samples;

  PyObject* pASDF;

  #ifndef SERIAL

  for(int i = 0; i < num_procs; i++){
    if(proc_id == i){
      MPI_Barrier(MPI_COMM_WORLD);
      pASDF = cached_obj_getter::call_wrapper(string("new_new_objects"), string("mW"), globals::pParams, globals::recalculate, false, false, false);
      cpp_caller::py_print(pASDF);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }


  #else
  pASDF = cached_obj_getter::call_wrapper(string("new_new_objects"), string("mW"), globals::pParams, globals::recalculate, true, true, true);

  #endif


  arbi_array<string2d> data_list = cpp_caller::py_string_mat_to_cpp_string(pASDF, true);
  num_samples = data_list.size().i0;

  #ifndef SERIAL

  // each node only cares about some of the samples.  store these in a index
  arbi_array<int1d> idx_i_care;
  for(int i = 0; i < num_samples; i++){
    if((i % num_procs) == proc_id){
      //idx_i_care.append(i);
      helpers::append(idx_i_care, i);
    }
  }

  #else

  num_samples = data_list.size().i0;
  // each node only cares about some of the samples.  store these in a index
  arbi_array<int1d> idx_i_care; idx_i_care.resize(0);
  for(int i = 0; i < num_samples; i++){
    //idx_i_care.append(i);
    helpers::append(idx_i_care,i);
  }

  
  #endif

  cout<<"building training samples"<<endl;
  //cout<<idx_i_care<<endl;

  // keep track of how many exceptions there were
  int num_exceptions = 0;

  this->training_indicies.resize(0);
  this->testing_indicies.resize(0);


  PyObject* pSelf = cached_obj_getter::get_param(globals::pParams, string("self"));
  bool self_test;
  if(pSelf == Py_True){
    self_test = true;
  }
  else{
    self_test = false;
  }
  //  Py_DECREF(pSelf);
  
  for(int i = 0; i < idx_i_care.size().i0; i++){
    //cout<<i<<endl;
    try{

      int j = idx_i_care(i);

      #ifdef PARAM
      cout<<data_list(j,0)<<endl;
      cout<<data_list(j,1)<<endl;
      string temp1 = data_list(j,0);
      string temp2 = data_list(j,1);
      cpp_caller::set_param(globals::pParams, string("pdb_name"), &temp1, globals::STRING_TYPE);
      cpp_caller::set_param(globals::pParams, string("chain_letter"), &temp2, globals::STRING_TYPE);
      #endif

      sample s = read_sample();


      cout<<"proc_id "<<proc_id<<" finished reading "<<temp1<<endl;

      //this->data.append(s);
      helpers::append(this->data, s);
      int next_sample_new_idx = this->data.size().i0 - 1;
      // assign testing/training based on original index
      if(self_test == false){
	if(j % this->num_folds == this->which_fold){
	//this->testing_indicies.append(next_sample_new_idx);
	  helpers::append(testing_indicies, next_sample_new_idx);
	}
	else{
	  //this->training_indicies.append(next_sample_new_idx);
	  helpers::append(training_indicies, next_sample_new_idx);
	}
      }
      else{
	helpers::append(testing_indicies, next_sample_new_idx);
	helpers::append(training_indicies, next_sample_new_idx);
      }
	
    }
    catch(...){
      num_exceptions++;
      cout<<"default exception while reading sample in folder: "<<endl;
    }
    
  }


  this->num_training = training_indicies.size().i0;
  this->num_testing = testing_indicies.size().i0;


  cout<<"training_indicies: "<<training_indicies<<endl;
  cout<<"testing_indicies: "<<testing_indicies<<endl;



  cout<<"NUM EXCEPTIONS: "<<proc_id<<" "<<num_exceptions<<endl;
  #ifndef SERIAL
  int exception_sum = 0;
  MPI_Barrier(MPI_COMM_WORLD);
  cout<<proc_id<<" num "<<num_exceptions<<endl;
  MPI_Reduce(&num_exceptions, &exception_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  cout<<proc_id<<" sum "<<exception_sum<<endl;
  MPI_Bcast(&exception_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
  cout<<proc_id<<" sumpost "<<exception_sum<<endl;
  MPI_Barrier(MPI_COMM_WORLD);
  num_exceptions = exception_sum;
  #endif

  if(proc_id == 0){
    cout<<"TOTAL NUMBER OF EXCEPTIONS: "<<num_exceptions<<endl;
  }
  #ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  this->num_samples = data.size().i0;

}


  
// params should have pdb list.  retrieve the list, set individual
model::model(){

  this->which_reg = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wreg")));
  this->num_states = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("ns")));
  this->num_folds = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("nfld")));
  this->which_fold = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wfld")));
  this->reg_constant = PyFloat_AsDouble(cpp_caller::get_param(globals::pParams, string("reg")));
  this->mean_field_max_iter = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("mfmi")));
  this->which_infer = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wif")));
  this->which_fold = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wfld")));
  this->num_folds = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("nfld")));
  this->which_obj = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wob")));
  this->which_obj2 = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wob2")));
  load_data();
  //assign(this->num_folds, this->which_fold);


  // set the node and edge map.  remember that edge map should be symmetric
  int idx = 0;
 
  this->node_map = arbi_array<int2d>(this->num_states, this->num_node_features);
  for(int i = 0; i < this->num_states; i++){
    for(int j = 0; j < this->num_node_features; j++){
      this->node_map(i,j) = idx;
      idx++;
    }
  }
  this->edge_map = arbi_array<int3d>(this->num_states, this->num_states, this->num_edge_features);
  for(int i = 0; i < this->num_states; i++){
    for(int j = 0; j <= i; j++){
      for(int k = 0; k < this->num_edge_features; k++){
	this->edge_map(i,j,k) = idx;
	this->edge_map(j,i,k) = idx;
	idx++;
      }
    }
  }
  
  this->theta_length = idx;

  cout<<node_map<<endl;
  //exit(1);

  
  normalize();

  this->prev_obj = -1;
  
}

arbi_array<num1d> model::get_dL_dTheta(int which_obj, arbi_array<num1d> theta){

  arbi_array<num1d> ans(theta_length);
  //ans.fill(0);
  ans = 0;
  for(int i = 0; i < num_training; i++){
      ans = ans + data(training_indicies(i)).get_dL_dTheta(which_obj, theta);
  }

  return ans;
}

num model::get_reg(arbi_array<num1d> theta){
  
  if(this->which_reg == 0){
    num reg = 0;
    for(int i = 0; i < theta.size().i0; i++){
      reg += theta(i) * theta(i);
    }
    return reg * reg_constant / 2.0;
  }
  else if(this->which_reg == 1){
    num reg = 0;
    for(int i = 0; i < theta.size().i0; i++){
      reg += abs(theta(i));
    }
    return reg * reg_constant;
  }
}

arbi_array<num1d> model::get_dReg_dTheta(arbi_array<num1d> theta){

  arbi_array<num1d> reg_grad(theta.size().i0);
  // have to add in term due to regularization
  if(this->which_reg == 0){
    for(int i = 0; i < theta.size().i0; i++){
      reg_grad(i) = theta(i) * reg_constant;
    }
    return reg_grad;
  }
  else if(this->which_reg == 1){
    for(int i = 0; i < theta.size().i0; i++){
      if(theta(i) > 0){
	reg_grad(i) = 1;
      }
      else{
	reg_grad(i) = -1.0;
      }
    }
    reg_grad *= this->reg_constant;
    return reg_grad;
  }
  else{
    assert(false);
  }
}
  
  

num model::get_L(int which_obj, arbi_array<num1d> theta){

  if(proc_id == 0){
    cout<<"which_obj: "<<which_obj<<endl;
  }

  num ans = 0, temp;
  for(int i = 0; i < num_training; i++){
    temp = data(training_indicies(i)).get_L(which_obj, theta);
    ans += temp;
  }

  return ans;
}

class My_Minimizer: public Minimizer{

 public:
  
  model* p_model;
  int which_obj;

  My_Minimizer(model* _p_model): Minimizer(false) {p_model = _p_model;}

  virtual void ComputeGradient(vector<double>& gradient, const vector<double>& x){
    
    //int which_obj = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wob")));
    int which_obj = p_model->which_obj;

    #ifdef SERIAL

    arbi_array<num1d> theta; theta.resize(x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    arbi_array<num1d> ans = p_model->get_dL_dTheta(which_obj, theta);
    
    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = ans(i);
    }
    #else

    arbi_array<num1d> theta(int(x.size()));
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }

    arbi_array<num1d> ans = p_model->get_dL_dTheta(which_obj, theta);
    num* ans_array = new num[p_model->theta_length];
    num* ans_sum_array = new num[p_model->theta_length];
    for(int i = 0; i < p_model->theta_length; i++){
      ans_array[i] = ans(i);
    }

   

    // send all ans to root for reducing.
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(ans_array, ans_sum_array, p_model->theta_length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    // broadcast sum back to everyone
    MPI_Bcast(ans_sum_array, p_model->theta_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //assert(gradient.size() == ans.linear_length);

    for(int i = 0; i < p_model->theta_length; i++){
      gradient[i] = ans_sum_array[i];
    }

    delete[] ans_array;
    delete[] ans_sum_array;
    
    #endif

    arbi_array<num1d> reg_grad = p_model->get_dReg_dTheta(theta);


    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = gradient[i] + reg_grad(i);
    }

    if(proc_id==0)cout<<"GRADIENT: "<<endl;
    for(int i = 0 ; i < gradient.size(); i++){
      if(proc_id==0)cout<<gradient[i]<<" ";
    }
    if(proc_id==0)cout<<endl;
  }

  virtual double ComputeFunction(const vector<double>& x){

    

    #ifndef SERIAL
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    arbi_array<num1d> theta(int(x.size()));
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }

    //int which_obj = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wob")));
    int which_obj = p_model->which_obj;

    double ans = p_model->get_L(which_obj, theta);
    assert(isfinite(ans));

    #ifndef SERIAL

    // now, send all messages to root for reducing.  then broadcast result back to everyone
    double ans_sum = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&ans, &ans_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ans_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ans = ans_sum;
    
    #endif


    // add in reg constant
    num reg_penalty = p_model->get_reg(theta);

    #ifndef SERIAL
    MPI_Barrier(MPI_COMM_WORLD);
    #endif

    return ans + reg_penalty;

  }

  virtual void Report (const vector<double> &theta, int iteration, double objective, double step_length){
    int s;

    if(p_model->prev_obj == -1){
      p_model->prev_obj = objective;
    }
    else{
      if(fabs(p_model->prev_obj-objective)< 0.00001){
	//exit(0);
      }
      p_model->prev_obj = objective;
    }
      

    if(iteration > 25000){
      p_model->which_obj = p_model->which_obj2;
    }

    if(iteration%5 == 1){
      arbi_array<num1d> theta_aa;  theta_aa.resize(theta.size());
      for(int i = 0; i < theta.size(); i++){
	theta_aa(i) = theta[i];
      }
      p_model->report(theta_aa, iteration, objective);
      
    }
  }

  virtual void Report (const string &s){
    int s1;
  }
};

set<string> cpp_caller::added_paths;

bool is_numeric(string s){
  for(int i = 0; i < s.length(); i++){
    if(isdigit(s[i])){
      return true;
    }
  }
  return false;
}

bool has_decimal(string s){
  for(int i = 0; i < s.length(); i++){
    if(s[i] == '.'){
      return true;
    }
  }
  return false;
}
    

// sets globals::pParams based on input arguments
// 
void set_input_params(PyObject* pParams, int argc, char* argv[]){
  // assume that arguments alternate 
  int i = 1;
  while(i < argc){
    string key = string(argv[i]);
    string val = string(argv[i+1]);
    cout<<i<<" "<<val<<" ok "<<key<<endl;
    // determine what type val is
    if(!is_numeric(val)){
      cout<<"string"<<endl;
      cpp_caller::set_param(pParams, key, &val, globals::STRING_TYPE);
    }
    else if(has_decimal(val)){
      cout<<"float"<<endl;
      num temp = atof(val.c_str());
      cpp_caller::set_param(pParams, key, &temp, globals::NUM_TYPE);
    }
    else{
      cout<<"int"<<endl;
      int temp = atoi(val.c_str());
      cpp_caller::set_param(pParams, key, &temp, globals::INT_TYPE);
    }
    i = i + 2;
   }
  cout<<i<<" whoz "<<argc<<endl;
  //assert(i == argc+1);
}




int main(int argc, char* argv[]){

  globals::init(argc, argv);
  Py_Initialize();



  #ifndef SERIAL
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  cout<<"PROC_ID_BEG: "<<proc_id<<endl;
  //exit(1);
  #else

  proc_id = 0;

  #endif

  cout<<"EE"<<endl;

  PyObject* pSysPath= cpp_caller::_get_module_PyObject(string("sys"), string("path"));
  
  for(int i = 0; i < PyList_Size(pSysPath); i++){
    cpp_caller::added_paths.insert(string(PyString_AsString(PyList_GetItem(pSysPath, i))));
  }

  cout<<"DD"<<endl;

  // for each node, set global_stuff.proc_id so that while writing files, can tell only root proc to do so

  PyObject* pProc_Id;
  PyObject* pBoolTemp;
  PyObject* pResult;

  cout<<"AA"<<endl;

  #ifdef SERIAL 

  pProc_Id = cpp_caller::get_module_PyObject(string("global_stuff"), string("proc_id"));
  pProc_Id = PyInt_FromLong(proc_id);
  globals::pParams = cpp_caller::get_module_PyObject(string("parameters"), string("the_params"));
  pBoolTemp = cpp_caller::get_module_PyObject(string("global_stuff"), string("recalculate"));
  if(pBoolTemp == Py_True){ 
    globals::recalculate = true;
  }
  else{
    globals::recalculate = false;
  }

  pBoolTemp = cpp_caller::get_module_PyObject(string("global_stuff"), string("recalculate_nodewise_loss_f"));
  if(pBoolTemp == Py_True){ 
    globals::recalculate_nodewise_loss_f = true;
  }
  else{
    globals::recalculate_nodewise_loss_f = false;
  }



  // create the experiment info file.  the experiment output results logically depend on this info file, since model reads info, creates model and then output.
  // result is a file handle which we don't do anything with, so decref it right away
  pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("lW"), globals::pParams, globals::recalculate, false, false);
  Py_DECREF(pResult);


  // test if param input is working properly
  set_input_params(globals::pParams, argc, argv);
  cpp_caller::py_print(globals::pParams);
  



  #else

  cout<<"BB"<<endl;

  for(int i = 0; i < num_procs; i++){
    MPI_Barrier(MPI_COMM_WORLD);
    if(proc_id == i){

      pProc_Id = cpp_caller::get_module_PyObject(string("global_stuff"), string("proc_id"));
      pProc_Id = PyInt_FromLong(proc_id);

      globals::pParams = cpp_caller::get_module_PyObject(string("parameters"), string("the_params"));
      set_input_params(globals::pParams, argc, argv);

      pBoolTemp = cpp_caller::get_module_PyObject(string("global_stuff"), string("recalculate"));
      if(pBoolTemp == Py_True){ 
	globals::recalculate = true;
      }
      else{
	globals::recalculate = false;
      }

      pBoolTemp = cpp_caller::get_module_PyObject(string("global_stuff"), string("recalculate_nodewise_loss_f"));
      if(pBoolTemp == Py_True){ 
	globals::recalculate_nodewise_loss_f = true;
      }
      else{
	globals::recalculate_nodewise_loss_f = false;
      }



      pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("lW"), globals::pParams, globals::recalculate, false, false);
      Py_DECREF(pResult);

    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  #endif

  model m;
  

  
  vector<num> w0(m.theta_length, 0);

  //srand(time(NULL));
  srand(0);
  for(int i = 0; i < m.theta_length; i++){
    w0[i] = ((num)(rand() % 100) / 100.0) - 0.5;
    //w0[i] = 1000;
    //w0[i] = 0;
    //w0[i]=w0ar(i);
    }
  
  My_Minimizer* minner = new My_Minimizer(&m);

  arbi_array<num1d> w0_aa(m.theta_length);

  ofstream asdf("theta_shorter.csv");
  for(int i = 0; i < m.theta_length; i++){
    asdf<<w0[i]<<',';
    w0_aa(i) = w0[i];
  }
  asdf.close();


  arbi_array<num3d> edge_potentials = m.data(0).get_edge_potentials(w0_aa);
  for(int i = 0; i < m.data(0).num_edges; i++){
    cout<<edge_potentials(i,0,0)<<" "<<edge_potentials(i,0,1)<<" "<<edge_potentials(i,1,0)<<" "<<edge_potentials(i,1,1)<<endl;
  }

  //cout<<m.data(0).edge_features<<endl;

  for(int i = 0; i < m.num_edge_features; i++){
    //cout<<m.edge_map(0,0,i)<<" "<<m.edge_map(0,1,i)<<" "<<m.edge_map(1,0,i)<<" "<<m.edge_map(1,1,i)<<endl;
  }

  
  cout<<w0_aa<<endl;
  //exit(1);


  arbi_array<num2d> node_marginals;
  arbi_array<num3d> edge_marginals;
  m.data(0).get_marginals(w0_aa, node_marginals, edge_marginals);
  cout<<node_marginals<<endl;
  cout<<edge_marginals(0,0,0)<<" "<<edge_marginals(0,1,0)<<" "<<edge_marginals(0,0,1)<<" "<<edge_marginals(0,1,1);
  cout<<w0_aa<<endl;


  cout<<m.data(0).get_edge_potential(edge_potentials,1,19,0,0)<<" "<<m.data(0).get_edge_potential(edge_potentials,1,19,0,1)<<" "<<m.data(0).get_edge_potential(edge_potentials,1,19,1,0)<<" "<<m.data(0).get_edge_potential(edge_potentials,1,19,1,1)<<endl;

  //exit(1);




  
  
  num val = minner->ComputeFunction(w0);  

  //exit(1);

  vector<num> grad(w0.size());
  minner->ComputeGradient(grad, w0);

  for(int i = 0; i < grad.size(); i++){
    cout<<grad[i]<<" ";
  }
  

  cout<<"theta: "<<endl;                                                                                                                                                           

  for(int i = 0; i < w0.size(); i++){                                                                                                                                              
    cout<<w0[i]<<" ";                                                                                                                                                              
  }
                                                                                                                                                                                 
 cout<<"\nfxn val: "<<val<<endl;                                                                                                                                                   
 //exit(1);





  
  //minner->LBFGS(w0,20000);

  minner->LBFGS(w0,5000);

  #ifndef SERIAL
  MPI_Finalize();
  #endif
  Py_Finalize();
  return 0;
}
