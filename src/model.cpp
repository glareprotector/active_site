#include "model.h"
//#include "nums.h"


void model::set_theta(arbi_array<num> _theta){
  this->theta = _theta;
}

void model::assign(int num_samples, arbi_array<int>& training_idx, arbi_array<int>& testing_idx, int& num_training, int& num_testing){
  // hardcode these in for now
  training_idx = arbi_array<int>(1,1);
  training_idx(0) = 0;
  num_training = 1;
  testing_idx = arbi_array<int>(1,1);
  testing_idx(0) = 0;
  num_testing = 0;
}

sample model::read_sample(string folder_name){
  /*
  // for now, just hardcode in the sample
  int num_nodes = 2;
  int num_edges = 1;
  
  arbi_array<num> node_features(2, num_nodes, this->num_node_features);
  node_features(0,0) = 1;
  node_features(1,0) = 1;
  
  arbi_array<num> edge_features(2, num_edges, this->num_edge_features);
  edge_features(0,0) = 1;
  
  arbi_array<int> edges(2, num_edges, 2);
  edges(0,0) = 0;
  edges(0,1) = 1;
  
  arbi_array<int> true_states(1, num_nodes);
  true_states(0) = 0;
  true_states(1) = 1;
  
  return sample(this, node_features, edge_features, edges, true_states);
  */

  // first read in info to figure out how many nodes, how many edges.  other params are read in already
  string info_file = folder_name + string("info.txt");
  string node_feature_file = folder_name + string("node_features.txt");
  string edge_features_file = folder_name + string("edge_features.txt");
  string true_states_file = folder_name + string("true_states.txt");
  string edge_file = folder_name + string("edges.txt");

  arbi_array<int> info = read_vect_to_int(info_file, 2);
  int num_nodes = info(0);
  int num_edges = info(1);
  
  arbi_array<num> node_features_transposed = read_mat_to_num(node_feature_file, this->num_node_features, num_nodes);
  arbi_array<num> node_features = arbi_array<num>::transpose(node_features_transposed);

  arbi_array<num> edge_features = read_mat_to_num(edge_features_file, num_edges, this->num_edge_features);
  arbi_array<int> true_states = read_vect_to_int(true_states_file, num_nodes);
  arbi_array<int> edges = read_mat_to_int(edge_file, num_edges, 2);
  return sample(this, node_features, edge_features, edges, true_states);
}

void model::load_data(arbi_array<string> folder_names){
  int num_samples = folder_names.size(0);
  cout<<"building training_idx"<<endl;
  arbi_array<int> training_idx;
  arbi_array<int> testing_idx;
  assign(num_samples, training_idx, testing_idx, this->num_training, this->num_testing);
  this->training_data = arbi_array<sample>(1,this->num_training);
  this->testing_data = arbi_array<sample>(1,this->num_testing);
  for(int i = 0; i < num_training; i++){
    this->training_data(i) = read_sample(folder_names(training_idx(i)));
  }
  for(int i = 0; i < num_testing; i++){
    this->testing_data(i) = read_sample(folder_names(testing_idx(i)));
  }
}

model::model(int _num_states, int _num_node_features, int _num_edge_features, arbi_array<int> _node_map, arbi_array<int> _edge_map, arbi_array<string> folder_names){
  this->theta_length = _num_states * _num_node_features + _num_states * _num_states * _num_edge_features;
  this->num_states = _num_states;
  this->num_node_features = _num_node_features;
  this->num_edge_features = _num_edge_features;
  this->node_map = _node_map;
  this->edge_map = _edge_map;
  this->gradient = arbi_array<num>(1, this->theta_length);
  this->theta = arbi_array<num>(1, this->theta_length);
  
  // would normally accept a list of folders to read from in load_data, but for now in load data just hardcode a sample
  
  // should read in num_states, num_node_features, num_edge_features first
  string model_info_file("asdf/");
  arbi_array<int> model_info = read_vect_to_int(model_info_file, 3);
  this->num_states = model_info(0);
  this->num_node_features = model_info(1);
  this->num_edge_features = model_info(2);


  load_data(folder_names);
}

arbi_array<num> model::get_gradient(){
  arbi_array<num> ans(1, theta_length);
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




//#include "sample.cpp"

int main(){

  arbi_array<num> tqt = read_mat_to_num(string("test_mat.txt"),2,3);
  cout<<tqt;
  return 0;
  
  cout<<"hello world"<<endl;

  int num_nodes = 2;
  int num_edges = 1;
  int num_states = 2;
  int num_node_features = 1;
  int num_edge_features = 1;
  int num_node_weights = num_states * num_node_features;
  int num_edge_weights = num_states * num_states * num_edge_features;
  
  arbi_array<int> node_map(2, num_states, num_node_features);
  arbi_array<int> edge_map(3, num_states, num_states, num_edge_features);
  int idx = 0;
  for(int i = 0; i < num_states; i++){
    for(int j = 0; j < num_node_features; j++){
      node_map(i,j) = idx;
      idx++;
    }
  }
  for(int i = 0; i < num_states; i++){
    for(int j = 0; j < num_states; j++){
      for(int k = 0; k < num_edge_features; k++){
	edge_map(i,j,k) = idx;
	idx++;
      }
    }
  }

  arbi_array<string> folder_names(1,1);
  folder_names(0) = string("asdf");


  model m(num_states, num_node_features, num_edge_features, node_map, edge_map, folder_names);


  int num_weights = idx;
  arbi_array<num> theta(1, num_weights);
  theta.fill(1);
  theta(0)=3;
  theta(1)=2;
  m.set_theta(theta);
  
  sample& s = m.training_data(0);

  cout<<"after set_theta"<<endl;
  s.set_node_potentials();
  cout<<"after set_node_potentials"<<endl;
  s.set_edge_potentials();
  cout<<"after set_edge_potentials"<<endl;
  
  


  //arbi_array<int> assignment(1,2);
  //assignment.fill(1);
  num likelihood = s.get_likelihood();
  cout<<"likelihood: "<<likelihood<<endl;
  
  cout<<"edge_potentials: "<<endl;
  cout<<s.edge_potentials;
  cout<<"node_potentials: "<<endl;
  cout<<s.node_potentials;
  
  s.set_node_marginals();
  cout<<"node_marginals"<<endl;
  cout<<s.node_marginals<<endl;
  cout<<"edge_marginals"<<endl;
  cout<<s.edge_marginals<<endl;

}
