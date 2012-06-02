#include "model.h"
#include "LBFGS.h"
#include <vector>
#ifndef SERIAL 
#include <mpi.h> 
#endif

// if you change theta, then need to set potentials/marginals
void model::set_theta(arbi_array<num> _theta){
  // only have to do something if this is a new value of theta
  if((_theta == theta) == false){
  
    this->theta = _theta;

    for(int i = 0; i < num_training; i++){
      training_data(i).set_node_potentials();
      training_data(i).set_edge_potentials();
      training_data(i).set_marginals();
    }
  }
}


sample model::read_sample(string folder_name){

  // first read in info to figure out how many nodes, how many edges.  other params are read in already
  string info_file = folder_name + string("info.txt");
  string node_feature_file = folder_name + string("Xnode.csv");
  string edge_features_file = folder_name + string("Xedge.csv");
  string true_states_file = folder_name + string("true_y.csv");
  string edge_file = folder_name + string("edge_list.csv");



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
  return sample(this, node_features, edge_features, edges, true_states, folder_name);
}

void model::load_data(arbi_array<string> folder_names){

  #ifdef SERIAL

  int num_samples = folder_names.size(0);
  cout<<"building training samples"<<endl;

  for(int i = 0; i < num_samples; i++){
    try{
      sample s = read_sample(folder_names(i));
      // for now assigning everything to training, but have to change this

      this->training_data.append(s);
    }
    catch(...){
      cout<<"default exception"<<endl;
    }
  }
  this->num_training = training_data.size(0);
  this->num_testing = testing_data.size(0);
  
  #else

  int num_samples = folder_names.size(0);

  cout<<"building training samples"<<endl;
  // each process gets training data based on it's process id
  for(int i = 0; i < num_samples; i++){
    try{
      sample s = read_sample(folder_names(i));
      // for now assigning everything to training, but have to change this
    
      //bool possibly_training = (i % 2) == 0;
      bool possibly_training = true;
      bool do_i_care = (i % num_procs) == proc_id;
      
      if(do_i_care){
	if(possibly_training){
	  this->training_data.append(s);
	}
	else{
	  this->testing_data.append(s);
	}
      }
    }
    catch(...){
      cout<<"default exception"<<endl;
    }
  }
  this->num_training = training_data.size(0);
  this->num_testing = testing_data.size(0);
  
  #endif

}


  

model::model(int _num_states, int _num_node_features, int _num_edge_features, arbi_array<string> folder_names, int _mean_field_max_iter){

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

  this->mean_field_max_iter = _mean_field_max_iter;

  load_data(folder_names);
}

arbi_array<num> model::get_gradient(){

  arbi_array<num> ans(1, theta_length);
  ans.fill(0);
  for(int i = 0; i < num_training; i++){
    ans = ans + training_data(i).get_gradient();
  }

  return ans;
}

num model::get_likelihood(){

  num ans = 0;
  for(int i = 0; i < num_training; i++){
    ans += training_data(i).get_likelihood();
  }

  return ans;
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
    p_model->set_theta(theta);
    arbi_array<num> ans = p_model->get_gradient();
    assert(gradient.size() == ans.linear_length);
    for(int i = 0; i < gradient.size(); i++){
      gradient[i] = ans(i);
    }
    
    #else

    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    p_model->set_theta(theta);
    arbi_array<num> ans = p_model->get_gradient();
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
    p_model->set_theta(theta);
    double ans = p_model->get_likelihood();
    return ans;
    
    #else
    
    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    p_model->set_theta(theta);
    double ans = p_model->get_likelihood();

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
  }

  virtual void Report (const string &s){
    int s1;
  }
};



//#include "sample.cpp"

int main(int argc, char** argv){

  globals::init(argc, argv);

  #ifndef SERIAL
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  #endif

  arbi_array<string> pdb_folders = read_vect_to_string(globals::pdb_list_file);

  for(int i = 0; i < pdb_folders.size(0); i++){
    pdb_folders(i) = globals::data_folder + pdb_folders(i) + '/';
  }

  cout<<pdb_folders<<endl;


  int num_states = 2;
  int num_node_features = 27;
  int num_edge_features = 1;

  model m(num_states, num_node_features, num_edge_features, pdb_folders, globals::mean_field_max_iter);

  #ifndef SERIAL 
  MPI_Barrier(MPI_COMM_WORLD); 
  #endif

  cout<<"ALL FOLDERS:"<<endl;
  for(int i = 0; i < m.num_training; i++){
    cout<<proc_id<<": "<<m.training_data(i).folder<<endl;
  }

  #ifndef SERIAL 
  MPI_Barrier(MPI_COMM_WORLD); 
  #endif

  vector<num> w0(m.theta_length, 0);
  My_Minimizer* minner = new My_Minimizer(&m);
  minner->LBFGS(w0,10);

  #ifndef SERIAL
  MPI_Finalize();
  #endif

  return 0;
}
