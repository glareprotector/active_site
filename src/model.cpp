#include "model.h"

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
