#include "model.h"
#include "LBFGS.h"
#include <vector>
#include "helpers.h"
#include <time.h>



#ifndef SERIAL 
#include <mpi.h> 
#endif

sample model::read_sample(string folder_name){

  // first read in info to figure out how many nodes, how many edges.  other params are read in already
  string info_file = folder_name + string("info.txt");
  string node_feature_file = folder_name + string("Xnode.csv");
  string edge_features_file = folder_name + string("Xedge.csv");
  string true_states_file = folder_name + string("true_y.csv");
  string edge_file = folder_name + string("edge_list.csv");

  string pdb_string = helpers::split_string(folder_name, '/', 0, 1);
  

  #ifdef PARAM
  // keep for now only bc sample for now takes these in as input
  string pdb_name = cpp_caller::CPPString_From_PyString(cpp_caller::get_param(globals::pParams, string("pdb_name")));
  string chain_letter = cpp_caller::CPPString_From_PyString(cpp_caller::get_param(globals::pParams, string("chain_letter")));

  PyObject* pNodeFeatures = cached_obj_getter::call_wrapper(string("new_new_objects"), string("jW"), globals::pParams, globals::recalculate, true, true);
  //cpp_caller::py_print(pNodeFeatures);
  arbi_array<num> node_features2 = cpp_caller::py_float_mat_to_cpp_num_mat(pNodeFeatures, true);


  #else
  
  string pdb_name = helpers::split_string(pdb_string, '_', 0, 1);
  string chain_letter = helpers::split_string(pdb_string, '_', -1, 0);
  #endif

  //arbi_array<int> info = read_vect_to_int(info_file, 2, ' ');
  //int num_nodes = info(0);
  //int num_edges = info(1);
  int xasdf;
  #ifdef PARAM

  //bool recalculate = false;
  arbi_array<num> node_features = cpp_caller::py_float_mat_to_cpp_num_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("jW"), globals::pParams, globals::recalculate, true, true), true);
  arbi_array<num> edge_features = cpp_caller::py_float_mat_to_cpp_num_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("kW"), globals::pParams, globals::recalculate, true, true), true);  
  arbi_array<int> true_states = cpp_caller::py_int_list_to_cpp_int_vect(cached_obj_getter::call_wrapper(string("new_new_objects"), string("oW"), globals::pParams, globals::recalculate, true, true), true);
  arbi_array<int> edge_list = cpp_caller::py_int_mat_to_cpp_int_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("iW"), globals::pParams, globals::recalculate, true, true), true);


  




  this->num_node_features = node_features.size(1);
  this->num_edge_features = edge_features.size(1);


  //cout<<node_features<<endl;
  //cout<<edge_list<<endl;
  //cout<<edge_features<<endl;
  //cout<<true_states<<endl;
 

  #else
  arbi_array<num> node_features_transposed = read_mat_to_num(node_feature_file, this->num_node_features, num_nodes);
  arbi_array<num> node_features = arbi_array<num>::transpose(node_features_transposed);

  arbi_array<num> edge_features_transposed = read_mat_to_num(edge_features_file, this->num_edge_features, num_edges);
  arbi_array<num> edge_features = arbi_array<num>::transpose(edge_features_transposed);
  arbi_array<int> true_states = read_vect_to_int(true_states_file, num_nodes, ',');

  // the true_states file has 1 and 2's, but this program wants 0 and 1's
  for(int i = 0; i < num_nodes; i++){
    true_states(i)--;
  }

  arbi_array<int> edges = read_mat_to_int(edge_file, num_edges, 2);

  #endif
  //cout<<node_features<<endl;

  //node_features.fill(0);
  //edge_features.fill(0);
  return sample(this, node_features, edge_features, edge_list, true_states, folder_name, pdb_name, chain_letter);
}

void model::normalize(){

  for(int i = 0; i < num_node_features; i++){

    bool to_normalize;
    //to_normalize = (i > 0) && (i < 4);
    to_normalize = i >= num_node_features - 4;
    //if(true){
     if(to_normalize){

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
	  data(j).node_features(k,i) /= variance;
	}
      }
    }
  }
}

// writes to specified folder the scores for each score as well as their true class
void model::report(arbi_array<num> theta, int iteration, num obj){

  /*if(proc_id == 0){
    cout<<"f_theta: "<<f_theta<<endl;
    cout<<"cur_theta: "<<theta<<endl;
    }*/

  arbi_array<int> true_classes;
  arbi_array<num> scores;
  // sample_labels will be of length (4+1) * num_testing, with the 1 for the null characters, 4 for the length of sample_names
  char* sample_pdb_names;
  char* sample_chain_letters;
  int* sample_lengths;

  // need to retrieve the length of sample_pdb_names and sample_chain_letters later
  int num_testing_master;

  // first get total length of stuff to report
  int num_data = 0;
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
    true_classes = arbi_array<int>(1, num_data_sum);
    scores = arbi_array<num>(1, num_data_sum);
    int sample_pdb_names_length = num_testing_sum * 5;
    sample_pdb_names = new char[sample_pdb_names_length];
    int sample_chain_letters_length = num_testing_sum * 2;
    sample_chain_letters = new char[sample_chain_letters_length];
    sample_lengths = new int[num_testing_sum];
  }
  else{
    true_classes = arbi_array<int>(1, num_data);
    scores = arbi_array<num>(1, num_data);
    int sample_pdb_names_length = num_testing * 5;
    sample_pdb_names = new char[sample_pdb_names_length];
    int sample_chain_letters_length = num_testing * 2;
    sample_chain_letters = new char[sample_chain_letters_length];
    sample_lengths = new int[num_testing];
  }

  #else

  true_classes = arbi_array<int>(1, num_data);
  scores = arbi_array<num>(1, num_data);
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
  for(int i = 0; i < num_procs; i++){
    #ifndef SERIAL
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
    //if(i == proc_id){
    //cout<<"FFF"<<endl<<"num_Testing"<<num_testing;
    //}
    #ifndef SERIAL
    MPI_Barrier(MPI_COMM_WORLD);
    #endif
  }
  //exit(1);
  for(int i = 0; i < num_testing; i++){
    arbi_array<num> node_marginals, edge_marginals;
    data(testing_indicies(i)).get_marginals(theta, node_marginals, edge_marginals);
    for(int j = 0; j < data(testing_indicies(i)).num_nodes; j++){
      true_classes(idx) = data(testing_indicies(i)).true_states(j);
      scores(idx) = node_marginals(j,1);
      //if(true_classes(idx) == 1){
      //cout<<scores(idx)<<" s ";
      //}
      idx++;
      //cout<<"scores1111"<<endl<<scores<<endl;
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
      class_buf = true_classes.m_data;
      score_buf = scores.m_data;
      next_score_size = true_classes.size(0);
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
      class_buf = true_classes.m_data + pos;
      score_buf = scores.m_data + pos;
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
  arbi_array<string> parsed_pdb_names(1, num_testing_master);
  arbi_array<string> parsed_chain_letters(1, num_testing_master);
  arbi_array<int> sample_lengths_ar(1, num_testing_master);
  if(proc_id == 0){
    for(int i = 0; i < num_testing_master; i++){
      parsed_pdb_names(i) = string(sample_pdb_names + (5*i));
      parsed_chain_letters(i) = string(sample_chain_letters + (2*i));
      sample_lengths_ar(i) = sample_lengths[i];
    }
  }

  //
  string class_file = results_folder + string("true_classes.csv");
  string score_file = results_folder + string("scores.csv");
  string theta_file = results_folder + string("theta.csv");
  //cout<<sample_lengths_ar<<endl;


  #ifdef PARAM

  //cout<<"scores"<<endl<<parsed_pdb_names<<endl;
  //exit(0);


  if(proc_id == 0){
    cpp_caller::set_param(globals::pParams, string("scores"), &scores, globals::NUM_VECT);
    cpp_caller::set_param(globals::pParams, string("true_states"), &true_classes, globals::INT_VECT);
    cpp_caller::set_param(globals::pParams, string("pdb_names"), &parsed_pdb_names, globals::STRING_VECT);
    cpp_caller::set_param(globals::pParams, string("chain_letters"), &parsed_chain_letters, globals::STRING_VECT);
    cpp_caller::set_param(globals::pParams, string("sizes"), &sample_lengths_ar, globals::INT_VECT);
    // also send the current iteration
    cpp_caller::set_param(globals::pParams, string("iter"), &iteration, globals::INT_TYPE);
    cpp_caller::set_param(globals::pParams, string("obj_val"), &obj, globals::NUM_TYPE);
    // always recalculate this
    cached_obj_getter::call_wrapper(string("new_new_objects"), string("qW"), globals::pParams, true, false, false, true);
  }
  #ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  #else
  
  if(proc_id == 0){
    true_classes.write(class_file, ',');
    scores.write(score_file, ',');
    theta.write(theta_file, ',');
  }

  #endif
  
  // now, feed this info to a cpp_caller.  first need to convert sample_lengths to arbi_array
  if(proc_id == 0){
    arbi_array<int> sample_lengths_ar(1, num_testing_master);
    for(int i = 0; i < num_testing_master; i++){
      sample_lengths_ar(i) = sample_lengths[i];
    }
    
  }
}

void model::assign(int num_folds, int which_fold){
  


  training_indicies = arbi_array<int>(1,0);
  testing_indicies = arbi_array<int>(1,0);

  for(int i = 0; i < num_samples; i++){
     if((i % num_folds) == which_fold){
      training_indicies.append(i);
        }
      else{
      testing_indicies.append(i);
      }
  }
  cout<<"testing: "<<testing_indicies<<endl;
  cout<<"training: "<<training_indicies<<endl;

  //assert(false);

  num_training = training_indicies.size(0);
  num_testing = testing_indicies.size(0);
  cout<<"8888888888888888"<<endl;
  //cout<<endl<<"proc_id "<<proc_id<<" training "<<training_indicies<<" testing "<<testing_indicies<<endl;
  #ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
}

// future usage: data_list is stored as param.  one by one, set pdb_name, call sample constructor with those params, which will be global

void model::load_data(arbi_array<string> folder_names){

  int num_samples;

  #ifdef PARAM
  //bool recalculate = false;


  // nodes have to take turns reading this in, since this object is not chain specific, and involves wile writing
  cout<<"00000000000000099999"<<endl;
  PyObject* pASDF;

  #ifndef SERIAL

  for(int i = 0; i < num_procs; i++){
    if(proc_id == i){
      MPI_Barrier(MPI_COMM_WORLD);
      pASDF = cached_obj_getter::call_wrapper(string("new_new_objects"), string("mW"), globals::pParams, globals::recalculate, true, true, true);
      cpp_caller::py_print(pASDF);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  cout<<"111111111222222222"<<endl;

  #else
  pASDF = cached_obj_getter::call_wrapper(string("new_new_objects"), string("mW"), globals::pParams, globals::recalculate, true, true, true);

  #endif




  cpp_caller::py_print(pASDF);
  arbi_array<string> data_list = cpp_caller::py_string_mat_to_cpp_string(pASDF, true);
  num_samples = data_list.size(0);
  //cpp_caller::py_print(pASDF);
  cout<<data_list<<endl;

  cout<<"                                      0000001111"<<endl;


  #endif

  #ifndef SERIAL

  // each node only cares about some of the samples.  store these in a index
  arbi_array<int> idx_i_care(1,0);
  for(int i = 0; i < num_samples; i++){
    if((i % num_procs) == proc_id){
      idx_i_care.append(i);
    }
  }

  cout<<idx_i_care<<endl<<"asdf"<<endl;
  //exit(1);

  #else

  num_samples = data_list.size(0);
    // each node only cares about some of the samples.  store these in a index
  arbi_array<int> idx_i_care(1,0);
  for(int i = 0; i < num_samples; i++){
    idx_i_care.append(i);
  }

  
  #endif

  #ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  cout<<proc_id<<endl<<idx_i_care<<endl;
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  //exit(0);

  cout<<"building training samples"<<endl;

  // keep track of how many exceptions there were
  int num_exceptions = 0;
  
  for(int i = 0; i < idx_i_care.size(0); i++){
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

      sample s = read_sample(folder_names(0));
      //cout<<s.node_features<<endl;
      this->data.append(s);
    }
    catch(...){
      num_exceptions++;
      cout<<"default exception while reading sample in folder: "<<endl;
    }

  }

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
    cout<<"NUMBER OF EXCEPTIONS: "<<num_exceptions<<endl;
  }
  #ifndef SERIAL
  MPI_Barrier(MPI_COMM_WORLD);
  #endif
  this->num_samples = data.size(0);
  
}


  
// params should have pdb list.  retrieve the list, set individual
model::model(int _num_states, int _num_node_features, int _num_edge_features, arbi_array<string> folder_names, int _mean_field_max_iter, int _num_folds, int _which_fold, string _results_folder, num _reg_constant, int _which_obj, int _which_infer){

  this->f_theta = arbi_array<num>(1, 57);
  
    
    for(int i = 0; i < 57; i++){
      f_theta(i) = ((num)(rand() % 100) / 100.0) - 0.5;
    }




  this->num_states = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("ns")));
  this->num_folds = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("nfld")));
  this->which_fold = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wfld")));
  load_data(folder_names);
  assign(this->num_folds, this->which_fold);
  cout<<"folds: "<<this->num_folds<<" "<<this->which_fold<<endl;
  //exit(1);
  //this->gradient = arbi_array<num>(1, this->theta_length);
  
  // set the node and edge map.  remember that edge map should be symmetric
  int idx = 0;
 
  this->node_map = arbi_array<int>(2, this->num_states, this->num_node_features);
  for(int i = 0; i < this->num_states; i++){
    for(int j = 0; j < this->num_node_features; j++){
      this->node_map(i,j) = idx;
      idx++;
    }
  }
  this->edge_map = arbi_array<int>(3, this->num_states, this->num_states, this->num_edge_features);
  for(int i = 0; i < this->num_states; i++){
    for(int j = 0; j <= i; j++){
      for(int k = 0; k < this->num_edge_features; k++){
	this->edge_map(i,j,k) = idx;
	this->edge_map(j,i,k) = idx;
	idx++;
      }
    }
  }
  cout<<node_map<<endl;
  this->theta_length = idx;
  //this->theta = arbi_array<num>(1, this->theta_length);


  #ifndef PARAM
  cout<<"ifnfdefPARAMS"<<endl;
  //exit(1);
  this->reg_constant = _reg_constant;
  this->mean_field_max_iter = _mean_field_max_iter;
  this->results_folder = _results_folder;
  this->which_obj = _which_obj;
  this->which_infer = _which_infer;

  #else

  this->reg_constant = PyFloat_AsDouble(cpp_caller::get_param(globals::pParams, string("reg")));
  this->mean_field_max_iter = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("mfmi")));
  this->which_infer = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wif")));
  this->which_fold = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wfld")));
  this->num_folds = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("nfld")));
  #endif

  
  normalize();
  //update_training();

  // results file wrapper should depend on something - the experiment info wrapper, which just prints out all params.  actually the model is the object created by the experiment wrapper
  // create params in c++, send to python, create experiment into wrapper using those params, which (fake calls model).  output spit to writer that depends on experiment wrapper
  // 

  /*for(int i = 0 ; i < num_samples; i++){
    data(i).simulate_states(f_theta);
    }*/


  



}

arbi_array<num> model::get_dL_dTheta(int which_obj, arbi_array<num> theta){

  arbi_array<num> ans(1, theta_length);
  ans.fill(0);
  for(int i = 0; i < num_training; i++){
      ans = ans + data(training_indicies(i)).get_dL_dTheta(which_obj, theta);
  }



  return ans;
}

num model::get_reg(arbi_array<num> theta){
  num reg = 0;
  for(int i = 0; i < theta.size(0); i++){
    reg += theta(i) * theta(i);
  }
  return reg * reg_constant / 2.0;
}

arbi_array<num> model::get_dReg_dTheta(arbi_array<num> theta){
  arbi_array<num> reg_grad(1, theta.size(0));
  // have to add in term due to regularization
  for(int i = 0; i < theta.size(0); i++){
    reg_grad(i) = theta(i) * reg_constant;
  }
  //  cout<<"reg_const: "<<reg_constant;
  return reg_grad;
}
  
  

num model::get_L(int which_obj, arbi_array<num> theta){

  //theta.write(string("theta_asdf"));


  //cout<<"num_training: "<<num_training<<" "<<training_indicies<<endl;

  //cout<<"which_obj: "<<which_obj<<endl;

  num ans = 0, temp;
  for(int i = 0; i < num_training; i++){
    //cout<<"a data: "<<data(training_indicies(i)).node_features<<endl;
    temp = data(training_indicies(i)).get_L(which_obj, theta);
    //cout<<endl<<proc_id<<" "<<i<<" "<<temp<<endl;
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
    
    int which_obj = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wob")));

    //cout<<"which_obj: "<<which_obj<<endl;


    //assert(false);

    //cout<<"gradient: "<<which_obj<<endl;

    #ifdef SERIAL

    arbi_array<num> theta(1, x.size());
    cout<<endl<<"THETA: ";
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
      cout<<x[i]<<' ';
    }
    cout<<endl;
    //p_model->set_theta(theta);
    arbi_array<num> ans = p_model->get_dL_dTheta(which_obj, theta);
    //cout<<"GRAD333"<<ans<<endl;
    assert(gradient.size() == ans.linear_length);
    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = ans(i);
    }
    //cout<<endl<<"GRADIENT"<<endl;
    cout<<"GRAD: "<<ans<<endl;
    //exit(1);
    #else

    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    //p_model->set_theta(theta);
    arbi_array<num> ans = p_model->get_dL_dTheta(which_obj, theta);
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

    assert(gradient.size() == ans.linear_length);

    for(int i = 0; i < p_model->theta_length; i++){
      gradient[i] = ans_sum_array[i];
    }

    delete[] ans_array;
    delete[] ans_sum_array;
    
    #endif

    arbi_array<num> reg_grad = p_model->get_dReg_dTheta(theta);
    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = gradient[i] + reg_grad(i);
    }
    if(proc_id == 0){
      cout<<"THETA: "<<theta<<endl;
      cout<<"GRADIENT: ";
    }
    if(proc_id==0){
      for(int i = 0; i < gradient.size(); i++){
	cout<<gradient[i]<<" ";
      }
    }
    if(proc_id==0)cout<<endl;
  }

  virtual double ComputeFunction(const vector<double>& x){
    //cout<<"COMPUTE FUNCTION"<<endl;
    //cout<<"fxn: "<<which_obj<<endl;

    //cout<<"THETAXXX"<<endl;
    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    //p_model->set_theta(theta);
    
    //cout<<"222222222"<<theta<<endl;
    int which_obj = PyInt_AsLong(cpp_caller::get_param(globals::pParams, string("wob")));

    double ans = p_model->get_L(which_obj, theta);
    //cout<<"ans: "<<ans<<endl;
    //exit(1);
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
    //cout<<"reg penalty: "<<reg_penalty;

    return ans + reg_penalty;

  }

  virtual void Report (const vector<double> &theta, int iteration, double objective, double step_length){
    int s;
    if(iteration%2 == 1){
      //cout<<"THETASIZE: "<<theta.size()<<endl;
      arbi_array<num> theta_aa(1, theta.size());
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

int main(int argc, char** argv){

  globals::init(argc, argv);
  Py_Initialize();



  #ifndef SERIAL
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  #else

  proc_id = 0;

  #endif



  PyObject* pSysPath= cpp_caller::_get_module_PyObject(string("sys"), string("path"));
  
  for(int i = 0; i < PyList_Size(pSysPath); i++){
    cpp_caller::added_paths.insert(string(PyString_AsString(PyList_GetItem(pSysPath, i))));
  }



  // for each node, set global_stuff.proc_id so that while writing files, can tell only root proc to do so

  PyObject* pProc_Id;
  PyObject* pBoolTemp;
  PyObject* pResult;

  #ifdef SERIAL 

  pProc_Id = cpp_caller::get_module_PyObject(string("global_stuff"), string("proc_id"));
  pProc_Id = PyInt_FromLong(proc_id);
  //exit(1);
  globals::pParams = cpp_caller::get_module_PyObject(string("parameters"), string("the_params"));
  pBoolTemp = cpp_caller::get_module_PyObject(string("global_stuff"), string("recalculate"));
  if(pBoolTemp == Py_True){ 
    globals::recalculate = true;
  }
  else{
    globals::recalculate = false;
  }

  // create the experiment info file.  the experiment output results logically depend on this info file, since model reads info, creates model and then output.
  // result is a file handle which we don't do anything with, so decref it right away
  pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("lW"), globals::pParams, globals::recalculate, false, false);
  Py_DECREF(pResult);

  #else
  for(int i = 0; i < num_procs; i++){
    MPI_Barrier(MPI_COMM_WORLD);
    if(proc_id == i){

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
      pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("lW"), globals::pParams, globals::recalculate, false, false);
      Py_DECREF(pResult);

    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  #endif
  cout<<"recalculate: "<<globals::recalculate<<endl;
  //exit(1);

  // insert params into the python environmenta
  /*globals::pParams = cpp_caller::get_new_empty_param();
  int dist_cut_off = 5;
  cpp_caller::set_param(globals::pParams, string("dist_cut_off"), &dist_cut_off, globals::INT_TYPE);
  num evalue = 1e-10;
  cpp_caller::set_param(globals::pParams, string("evalue"), &evalue, globals::NUM_TYPE);
  string data_list_file("catres_six.pdb_list");
  cpp_caller::set_param(globals::pParams, string("data_list_file"), &data_list_file, globals::STRING_TYPE);
  */

  // get globals::pParams from new_new_objects
  



  



    

  arbi_array<string> pdb_folders = read_vect_to_string(globals::pdb_list_file);
  
  for(int i = 0; i < pdb_folders.size(0); i++){
    pdb_folders(i) = globals::data_folder + pdb_folders(i) + '/';
  }

  cout<<proc_id<<": "<<pdb_folders<<endl;

  int num_states = 2;
  int num_node_features = 27;
  int num_edge_features = 1;

  model m(num_states, num_node_features, num_edge_features, pdb_folders, globals::mean_field_max_iter, globals::num_folds, globals::which_fold, globals::results_folder, globals::reg_constant, globals::which_obj, globals::which_infer);
  
  
  #ifndef SERIAL 
  MPI_Barrier(MPI_COMM_WORLD); 
  #endif

  

  #ifndef SERIAL 
  MPI_Barrier(MPI_COMM_WORLD); 
  #endif

  arbi_array<num> w0ar = read_vect_to_num(string("theta_asdf"), 57, ',');


  vector<num> w0(m.theta_length, 0);

  /* for(int i = 0; i < 20; i++){
    w0[i] = (num)(rand() % 40) / 40;
    w0[i] = .01;
    }*/
  //srand(time(NULL));
  srand(0);
  for(int i = 0; i < m.theta_length; i++){
    w0[i] = ((num)(rand() % 100) / 100.0) - 0.5;
    //w0[i] = 1000;
    //w0[i] = 0;
    w0[i]=w0ar(i);
    }
  
  My_Minimizer* minner = new My_Minimizer(&m);
  minner->LBFGS(w0,20000);

  #ifndef SERIAL
  MPI_Finalize();
  #endif
  Py_Finalize();
  return 0;
}
