#include "model.h"
#include "LBFGS.h"
#include <vector>
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

  string pdb_name;
  string chain_letter;
  istringstream split(folder_name);
  getline(split, pdb_name, '_');
  getline(split, chain_letter);

  arbi_array<int> info = read_vect_to_int(info_file, 2, ' ');
  int num_nodes = info(0);
  int num_edges = info(1);
  
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
  return sample(this, node_features, edge_features, edges, true_states, folder_name, pdb_name, chain_letter);
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

// calculates scores for every test based on the current value of theta
void model::test(){
  
  update_testing();
  report(results_folder);

}
    

// writes to specified folder the scores for each score as well as their true class
void model::report(arbi_array<num> theta, string report_folder){

  arbi_array<int> true_classes;
  arbi_array<num> scores;

  // first get total length of stuff to report
  int num_data = 0;
  for(int i = 0; i < num_testing; i++){
    num_data += data(testing_indicies(i)).num_nodes;
  }

  #ifndef SERIAL

  MPI_Status status;
  int num_data_sum;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&num_data, &num_data_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&num_data_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  
  // root process gets larger true_classes and scores
  if(proc_id == 0){
    true_classes = arbi_array<int>(1, num_data_sum);
    scores = arbi_array<num>(1, num_data_sum);
  }
  else{
    true_classes = arbi_array<int>(1, num_data);
    scores = arbi_array<num>(1, num_data);
  }

  #else

  true_classes = arbi_array<int>(1, num_data);
  scores = arbi_array<num>(1, num_data);

  #endif

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
  }

  // if doing things in parallel, have to add more stuff to the vectors
  
  #ifndef SERIAL

  // each process takes turns sending messages to process 0
  int pos = idx;
  for(int i = 1 ; i < num_procs; i++){
    int* class_buf;
    num* score_buf;
    int next_size;
    if(proc_id == i){
      class_buf = true_classes.m_data;
      score_buf = scores.m_data;
      next_size = true_classes.size(0);
      MPI_Send(&next_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(class_buf, next_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(score_buf, next_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else if(proc_id == 0){
      class_buf = true_classes.m_data + pos;
      score_buf = scores.m_data + pos;
      MPI_Recv(&next_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(class_buf, next_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(score_buf, next_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      pos += next_size;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  #endif

  //
  string class_file = report_folder + string("true_classes.csv");
  string score_file = report_folder + string("scores.csv");
  string theta_file = report_folder + string("theta.csv");
  #ifdef SERIAL

  true_classes.write(class_file, ',');
  scores.write(score_file, ',');

  #else
  
  if(proc_id == 0){
    true_classes.write(class_file, ',');
    scores.write(score_file, ',');
    theta.write(theta_file, ',');
  }

  #endif

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


void model::load_data(arbi_array<string> folder_names){

  #ifdef SERIAL
  
  int num_samples = folder_names.size(0);
  cout<<"building training samples"<<endl;

  for(int i = 0; i < num_samples; i++){
    try{
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
  for(int i = 0; i < folder_names.size(0); i++){
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


  

model::model(int _num_states, int _num_node_features, int _num_edge_features, arbi_array<string> folder_names, int _mean_field_max_iter, int _num_folds, int _which_fold, string _results_folder, num _reg_constant, int _which_obj){

  this->num_states = _num_states;
  this->num_node_features = _num_node_features;
  this->num_edge_features = _num_edge_features;
  this->gradient = arbi_array<num>(1, this->theta_length);
  
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
  
  this->theta_length = idx;
  this->theta = arbi_array<num>(1, this->theta_length);

  this->reg_constant = _reg_constant;
  this->mean_field_max_iter = _mean_field_max_iter;
  this->results_folder = _results_folder;
  this->which_obj = _which_obj;

  load_data(folder_names);
  assign(_num_folds, _which_fold);
  normalize();
  //update_training();
}

arbi_array<num> model::get_dL_dTheta(arbi_array<num> theta){

  arbi_array<num> ans(1, theta_length);
  ans.fill(0);
  for(int i = 0; i < num_training; i++){
      ans = ans + data(training_indicies(i)).get_dL_dTheta(which_obj, theta);
  }

  // have to add in term due to regularization
  for(int i = 0; i < theta_length; i++){
    ans(i) += theta(i) / reg_constant;
  }

  return ans;
}

num model::get_L(arbi_array<num> theta){

  num ans = 0;
  for(int i = 0; i < num_training; i++){
    ans += data(training_indicies(i)).get_L(which_obj, theta);
    }
  }

  // add in l2 penalty
  num penalty = 0;
  for(int i = 0; i < theta_length; i++){
    penalty += theta(i) * theta(i);
  }
  penalty /= (2.0 * reg_constant);

  return ans + penalty;
}

class My_Minimizer: public Minimizer{

 public:
  
  model* p_model;

  My_Minimizer(model* _p_model): Minimizer(false) {p_model = _p_model;};

  void ComputeGradient(vector<double>& gradient, const vector<double>& x){
    
    #ifdef SERIAL

    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    //p_model->set_theta(theta);
    arbi_array<num> ans = p_model->get_dL_dTheta(theta);
    assert(gradient.size() == ans.linear_length);
    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = ans(i);
    }
    
    #else

    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    //p_model->set_theta(theta);
    arbi_array<num> ans = p_model->get_dL_dTheta(theta);
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

  double ComputeFunction(const vector<double>& x){

    #ifdef SERIAL
      
    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    //p_model->set_theta(theta);
    double ans = p_model->get_L(theta);
    return ans;
    
    #else
    
    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    //p_model->set_theta(theta);
    double ans = p_model->get_L(theta);

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

  virtual void ReportAUC(){
    p_model->test();
  }

  virtual void Report (const vector<double> &theta, int iteration, double objective, double step_length){
    int s;
  }

  virtual void Report (const string &s){
    int s1;
  }
};


int main(int argc, char** argv){

  globals::init(argc, argv);
  
  PyObject* pName, pModule, pDict, pResult;

  pName = PyString_FromString("get_stuff");
  pModule = PyImport_Import(pName);
  pDict = PyModule_GetDict(pModule);
  

  #ifndef SERIAL
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  #endif

  arbi_array<string> pdb_folders = read_vect_to_string(globals::pdb_list_file);
  
  for(int i = 0; i < pdb_folders.size(0); i++){
    pdb_folders(i) = globals::data_folder + pdb_folders(i) + '/';
  }

  cout<<proc_id<<": "<<pdb_folders<<endl;

  int num_states = 2;
  int num_node_features = 27;
  int num_edge_features = 1;

  model m(num_states, num_node_features, num_edge_features, pdb_folders, globals::mean_field_max_iter, globals::num_folds, globals::which_fold, globals::results_folder, globals::reg_constant, globals::which_obj);
  
  
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
  My_Minimizer* minner = new My_Minimizer(&m);
  minner->LBFGS(w0,20000);

  #ifndef SERIAL
  MPI_Finalize();
  #endif

  return 0;
}
