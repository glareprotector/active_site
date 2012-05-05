#include "sample.h"
//#include "arrays.h"
#include <iostream>



using namespace std;


int main(){
  
  cout<<"hello world"<<endl;

  int num_nodes = 2;
  int num_edges = 1;
  int num_states = 2;
  int num_node_features = 1;
  int num_edge_features = 1;
  int num_node_weights = num_states * num_node_features;
  int num_edge_weights = num_states * num_states * num_edge_features;
  
  array_1d<num> test(2);
  //cout<<endl;
  test(0)=1;
  test(1)=2;
  cout<<test<<endl;
  //return 0;

  
  array_2d<num> node_features(num_nodes, num_node_features);
  node_features(0,0) = 1;
  node_features(1,0) = 1;

  cout<<node_features<<endl;

  //  return 0;

  cout<<"FFFFFFFFFFFFFFFFFFFFFF"<<endl;

  array_2d<num> edge_features(num_edges, num_edge_features);
  edge_features(0,0) = 1;

  array_2d<int> node_map(num_states, num_node_features);
  array_3d<int> edge_map(num_states, num_states, num_edge_features);
  int idx = 0;
  for(int i = 0; i < num_states; i++){
    for(int j = 0; j < num_node_features; j++){
      node_map(i,j) = idx;
      cout<<"    "<<i<<" "<<j<<endl;
      cout<<"idx "<<idx<<endl; 
      idx++;
    }
  }
  for(int i = 0; i < num_states; i++){
    for(int j = 0; j < num_states; j++){
      for(int k = 0; k < num_edge_features; k++){
	edge_map(i,j,k) = idx;
	cout<<"    "<<i<<" "<<j<<" "<<k<<endl;
	cout<<"idx "<<idx<<endl;
	idx++;
      }
    }
  }

  array_2d<int> edge_to_ij(num_edges, 3);
  edge_to_ij(0,1) = 0;
  edge_to_ij(0,2) = 1;

  sample x(node_features, edge_features, node_map, edge_map, edge_to_ij, num_nodes, num_edges, num_states, num_node_features, num_edge_features);
  
  int num_weights = idx;
  array_1d<num> theta(num_weights,1);
  theta(0)=3;
  theta(1)=2;
  x.set_theta(theta);
  cout<<"after set_theta"<<endl;
  x.set_node_potentials();
  cout<<"after set_node_potentials"<<endl;
  x.set_edge_potentials();
  cout<<"after set_edge_potentials"<<endl;
  

  array_1d<int> assignment(2,1);
  cout<<"likelihood: "<<x.get_likelihood(assignment)<<endl;
  
  cout<<"edge_potentials: "<<endl;
  cout<<x.m_edge_potentials<<endl;
  cout<<"node_potentials: "<<endl;
  cout<<x.m_node_potentials<<endl;
  cout<<"theta: "<<endl;
  cout<<x.m_theta<<endl;
  cout<<"num_weights: "<<endl;
  cout<<idx<<endl;
  
  /*
  for(int i = 0; i < num_edges; i++){
    for(int j = 0; j < num_states; j++){
      for(int k = 0; k < num_states; k++){
	cout<<i<<j<<k<<"z"<<x.m_edge_potentials(i,j,k)<<" ";
      }
    }
  }

  
  cout<<x.m_node_potentials<<endl;
  cout<<x.m_edge_potentials<<endl;
  */

  x.set_node_marginals();
  cout<<"node_marginals"<<endl;
  cout<<x.m_node_marginals<<endl;
  cout<<"edge_marginals"<<endl;
  cout<<x.m_edge_marginals<<endl;

}
