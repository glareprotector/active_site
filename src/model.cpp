#include "model.h"
#include "LBFGS.h"
#include <vector>
#include "helpers.h"

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

  PyObject* pNodeFeatures = cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_node_features_obj_w"), globals::pParams, false, true, true);
  cpp_caller::py_print(pNodeFeatures);
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

  bool recalculate = false;
  arbi_array<num> node_features = cpp_caller::py_float_mat_to_cpp_num_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_node_features_obj_w"), globals::pParams, recalculate, true, true), true);
  arbi_array<num> edge_features = cpp_caller::py_float_mat_to_cpp_num_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_edge_features_obj_w"), globals::pParams, recalculate, true, true), true);  
  arbi_array<int> true_states = cpp_caller::py_int_list_to_cpp_int_vect(cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_true_states_obj_w"), globals::pParams, recalculate, true, true), true);
  arbi_array<int> edge_list = cpp_caller::py_int_mat_to_cpp_int_mat(cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_pdb_chain_edge_list_obj_w"), globals::pParams, recalculate, true, true), true);

  cout<<node_features<<endl;
  cout<<edge_list<<endl;
  cout<<edge_features<<endl;
  cout<<true_states<<endl;
 

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
  
  return sample(this, node_features, edge_features, edge_list, true_states, folder_name, pdb_name, chain_letter);
}

void model::normalize(){

  for(int i = 0; i < num_node_features; i++){

    bool to_normalize;
    to_normalize = (i > 0) && (i < 4);

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
void model::report(arbi_array<num> theta){

  
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
    sample_labels = new char[sample_pdb_names_length];
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
  for(int i = 0; i < num_testing; i++){
    arbi_array<num> node_marginals, edge_marginals;
    data(testing_indicies(i)).get_marginals(theta, node_marginals, edge_marginals);
    for(int j = 0; j < data(testing_indicies(i)).num_nodes; j++){
      true_classes(idx) = data(testing_indicies(i)).true_states(j);
      scores(idx) = node_marginals(j,1);
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
      MPI_Send(sample_lengths, num_testing, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    else if(proc_id == 0){
      class_buf = true_classes.m_data + pos;
      score_buf = scores.m_data + pos;
      MPI_Recv(&next_score_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&next_sample_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(class_buf, next_score_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(score_buf, next_score_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(sample_pdb_names[5 * sample_pos], 5 * next_sample_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(sample_chain_letters[5*sample_pos], 2 * next_sample_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(sample_lengths[sample_pos], next_sample_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      pos += next_size;
      sample_pos += next_sample_size;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  #endif

  // now, parse the long strings into a name/letter for each sample.  only do this for proc_id 0
  // store each as arbi_array of strings
  arbi_array<string> parsed_pdb_names(1, num_testing_master);
  arbi_array<string> parsed_chain_letters(1, num_testing_master);
  if(proc_id == 0){
    for(int i = 0; i < num_testing_master; i++){
      parsed_pdb_names(i) = string(sample_pdb_names + (5*i));
      parsed_chain_letters(i) = string(sample_chain_letters + (2*i));
    }
  }

  //
  string class_file = results_folder + string("true_classes.csv");
  string score_file = results_folder + string("scores.csv");
  string theta_file = results_folder + string("theta.csv");
  #ifdef SERIAL

  true_classes.write(class_file, ',');
  scores.write(score_file, ',');

  // call 

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

void model::assign(int _num_folds, int _which_fold){
  
  this->num_folds = _num_folds;
  this->which_fold = _which_fold;

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

}

// future usage: data_list is stored as param.  one by one, set pdb_name, call sample constructor with those params, which will be global

void model::load_data(arbi_array<string> folder_names){

  #ifdef SERIAL
  
  #ifdef PARAM
  bool recalculate = true;

  PyObject* pASDF = cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_formatted_data_list_obj_w"), globals::pParams, recalculate, false, true);
  cpp_caller::py_print(pASDF);
  arbi_array<string> data_list = cpp_caller::py_string_mat_to_cpp_string(pASDF, true);
  int num_samples = data_list.size(0);
  cpp_caller::py_print(pASDF);
  cout<<data_list<<endl;
  cout<<"                                      0000001111"<<endl;

  #else

  int num_samples = folder_names.size(0);
  
  #endif

  cout<<"building training samples"<<endl;

  for(int i = 0; i < num_samples; i++){
    try{

      #ifdef PARAM
      cout<<data_list(i,0)<<endl;
      cout<<data_list(i,1)<<endl;
      string temp1 = data_list(i,0);
      string temp2 = data_list(i,1);
      cpp_caller::set_param(globals::pParams, string("pdb_name"), &temp1, globals::STRING_TYPE);
      cpp_caller::set_param(globals::pParams, string("chain_letter"), &temp2, globals::STRING_TYPE);
      #endif

      sample s = read_sample(folder_names(i));
      this->data.append(s);
    }
    catch(exception& e){
      cout<<e.what()<<endl;
    }
    catch(...){
      cout<<"default exception while reading sample in folder: "<<folder_names(i)<<endl;
    }
  }
  this->num_samples = data.size(0);
  
  #else

  cout<<"building training samples"<<endl;
  // each process gets training data based on it's process id
  for(int i = 0; i < num_samples; i++){
    try{
      sample s = read_sample(folder_names(i));
      bool do_i_care = (i % num_procs) == proc_id;
      
      if(do_i_care){
	this->data.append(s);
      }
    }
    catch(exception& e){
      cout<<e.what()<<" "<<folder_names(i)<<endl;
    }
    catch(...){
      cout<<"default exception while reading sample in folder: "<<folder_names(i)<<endl;
    }
  }
  this->num_samples = data.size(0);
  
  #endif

}


  
// params should have pdb list.  retrieve the list, set individual
model::model(int _num_states, int _num_node_features, int _num_edge_features, arbi_array<string> folder_names, int _mean_field_max_iter, int _num_folds, int _which_fold, string _results_folder, num _reg_constant, int _which_obj, int _which_infer){

  this->num_states = _num_states;
  this->num_node_features = _num_node_features;
  this->num_edge_features = _num_edge_features;
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

  this->reg_constant = _reg_constant;
  this->mean_field_max_iter = _mean_field_max_iter;
  this->results_folder = _results_folder;
  this->which_obj = _which_obj;
  this->which_infer = _which_infer;

  load_data(folder_names);
  assign(_num_folds, _which_fold);
  //normalize();
  //update_training();

  // results file wrapper should depend on something - the experiment info wrapper, which just prints out all params.  actually the model is the object created by the experiment wrapper
  // create params in c++, send to python, create experiment into wrapper using those params, which (fake calls model).  output spit to writer that depends on experiment wrapper
  // 

}

arbi_array<num> model::get_dL_dTheta(int which_obj, arbi_array<num> theta){

  arbi_array<num> ans(1, theta_length);
  ans.fill(0);
  for(int i = 0; i < num_training; i++){
      ans = ans + data(training_indicies(i)).get_dL_dTheta(which_obj, theta);
  }

  // have to add in term due to regularization
  for(int i = 0; i < theta_length; i++){
    ans(i) += theta(i) * reg_constant;
  }

  return ans;
}

num model::get_L(int which_obj, arbi_array<num> theta){

  num ans = 0;
  for(int i = 0; i < num_training; i++){
    ans += data(training_indicies(i)).get_L(which_obj, theta);
  }

  // add in l2 penalty
  num penalty = 0;
  for(int i = 0; i < theta_length; i++){
    penalty += reg_constant * theta(i) * theta(i);
  }
  penalty /= 2.0;

  return ans + penalty;
}

class My_Minimizer: public Minimizer{

 public:
  
  model* p_model;
  int which_obj;

  My_Minimizer(model* _p_model): Minimizer(false) {p_model = _p_model;}

  virtual void ComputeGradient(vector<double>& gradient, const vector<double>& x){
    
    int which_obj = globals::which_obj;

    //assert(false);

    //cout<<"gradient: "<<which_obj<<endl;

    #ifdef SERIAL

    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    //p_model->set_theta(theta);
    arbi_array<num> ans = p_model->get_dL_dTheta(which_obj, theta);
    //cout<<"GRAD333"<<ans<<endl;
    assert(gradient.size() == ans.linear_length);
    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = ans(i);
    }
    //cout<<endl<<"GRADIENT"<<endl;
    //cout<<ans<<endl;
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

    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = ans_sum_array[i];
    }

    delete[] ans_array;
    delete[] ans_sum_array;
    
    #endif

  }

  virtual double ComputeFunction(const vector<double>& x){
    //cout<<"COMPUTE FUNCTION"<<endl;
    //cout<<"fxn: "<<which_obj<<endl;
    #ifdef SERIAL
    //cout<<"THETAXXX"<<endl;
    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
      //cout<<x[i]<<" ";
    }
    //p_model->set_theta(theta);
    
    //cout<<"222222222"<<theta<<endl;
    int which_obj = globals::which_obj;

    double ans = p_model->get_L(which_obj, theta);
    assert(isfinite(ans));
    return ans;
    
    #else
    
    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }

    //p_model->set_theta(theta);
    double ans = p_model->get_L(which_obj, theta);

    // now, send all messages to root for reducing.  then broadcast result back to everyone
    double ans_sum = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&ans, &ans_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ans_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    return ans_sum;
    
    #endif

  }

  virtual void Report (const vector<double> &theta, int iteration, double objective, double step_length){
    int s;
    if(iteration%2 == 0){
      arbi_array<num> theta_aa(1, theta.size());
      for(int i = 0; i < theta.size(); i++){
	theta_aa(i) = theta[i];
      }
      p_model->report(theta_aa);
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

  PyObject* pSysPath= cpp_caller::_get_module_PyObject(string("sys"), string("path"));
  
  for(int i = 0; i < PyList_Size(pSysPath); i++){
    cpp_caller::added_paths.insert(string(PyString_AsString(PyList_GetItem(pSysPath, i))));
  }


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
  globals::pParams = cpp_caller::get_module_PyObject(string("new_new_objects"), string("the_params"));
  



  // create the experiment info file.  the experiment output results logically depend on this info file, since model reads info, creates model and then output.
  // result is a file handle which we don't do anything with, so decref it right away
  PyObject* pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("the_experiment_info_file_w"), globals::pParams, true, false, false);
  Py_DECREF(pResult);

  
  
    

  #ifndef SERIAL
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  #else

  proc_id = 0;

  #endif

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

  cout<<"ALL FOLDERS:"<<endl;
  for(int i = 0; i < m.num_samples; i++){
    cout<<proc_id<<": "<<m.data(i).folder<<endl;
  }

  #ifndef SERIAL 
  MPI_Barrier(MPI_COMM_WORLD); 
  #endif

  vector<num> w0(m.theta_length, 0);

  /* for(int i = 0; i < 20; i++){
    w0[i] = (num)(rand() % 40) / 40;
    w0[i] = .01;
    }*/
  
  for(int i = 0; i < m.theta_length; i++){
    w0[i] = (num)(rand() % 40) / 40;
    //w0[i] = 1000;
    }

  My_Minimizer* minner = new My_Minimizer(&m);
  minner->LBFGS(w0,20000);

  #ifndef SERIAL
  MPI_Finalize();
  #endif
  Py_Finalize();
  return 0;
}
