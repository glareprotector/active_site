#include <iostream>
#include "model.h"


using namespace std;


int main(){

  //read_mat_to_num(string("test_mat.txt"),2,3);

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
