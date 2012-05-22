// TEST
/*
This class is responsible for reading in 
1. node features (numNodes x numNodeFeatures)
2. edge features (numEdges x numEdgeFeatures)
3. nodeMap (numStates x numNodeFeatures
4. edgeMap (numStates x numStates x numEdgeFeatures)
5. info - numNodes, numEdges, numStates, numNodeFeatures, numEdgeFeatures
6. edgeNumber to (i,j) map

It represents one sample.  It supports the following functions (given theta)
1. calculate likelihood
2. calculate node marginals
3. calculate edge marginals
4. calculate gradient
*/

#ifndef SAMPLE_H
#define SAMPLE_H

#include "arrays.h"
#include <string>
#include <cmath>
#include <ext/hash_map>
#include <utility>
#include <tr1/memory>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace __gnu_cxx;

// hash function for pair<int, int>
struct pair_hash{
  size_t operator() (const pair<int,int>& p) const{
    return p.first + p.second;
  }
};

int mean_field_max_iter = 100;

class sample{

 public:
  
  array_2d<num> m_node_potentials;
  array_3d<num> m_edge_potentials;
  array_1d<num> m_theta;
  array_2d<num> m_node_features;
  array_2d<num> m_edge_features;
  array_2d<num> m_node_marginals;
  array_3d<num> m_edge_marginals;

  num likelihood;

  int num_states;
  int num_nodes;
  int num_edges;
  int num_node_features;
  int num_edge_features;



  array_2d<int> m_node_map;
  array_3d<int> m_edge_map;

  
  hash_map<int, array_1d<int> > node_to_neighbors;
  hash_map< pair<int, int>, int, pair_hash> edge_to_pos;
  hash_map< int, pair<int, int> > pos_to_edge;

  num get_node_potential(int node, int state){
    return m_node_potentials(node, state);
  }

  num get_edge_potential(int node1, int node2, int state1, int state2){
    if(node1 < node2){
      //cout<<edge_to_pos[pair<int,int>(node1,node2)]<<endl;
      //cout<<"GGGGG"<<endl;
      //cout<<m_edge_potentials<<endl;
      return m_edge_potentials(edge_to_pos[pair<int,int>(node1, node2)], state1, state2);
    }
    else{
      return get_edge_potential(node2, node1, state2, state1);
    }
  }

  sample(string file_path);

  sample(array_2d<num> node_features, array_2d<num> edge_features, array_2d<int> node_map, array_3d<int> edge_map, array_2d<int> edge_to_ij, int num_nodes, int num_edges, int num_states, int num_node_features, int num_edge_features){
    this->m_node_features = node_features;
    this->m_edge_features = edge_features;
    this->m_node_map = node_map;
    this->m_edge_map = edge_map;
    this->num_nodes = num_nodes;
    this->num_edges = num_edges;
    this->num_states = num_states;
    this->num_node_features = num_node_features;
    this->num_edge_features = num_edge_features;

    // initialize hash_maps
    node_to_neighbors = hash_map<int, array_1d<int> >();
    edge_to_pos = hash_map< pair<int, int>, int, pair_hash>();
    pos_to_edge = hash_map< int, pair<int,int> >();

    // assuming that i < j in edge_to_ij
    for(int i = 0; i < num_edges; i++){
      int node1 = edge_to_ij(i,1);
      int node2 = edge_to_ij(i,2);
      pos_to_edge[i] = pair<int,int>(node1, node2);
      edge_to_pos[pair<int,int>(node1, node2)] = i;
      if(node_to_neighbors.find(i) != node_to_neighbors.end()){
	node_to_neighbors[node1].append(node2);
      }
      else{
	node_to_neighbors[node1] = array_1d<int>(1, node2);
      }
    }
    
    // allocate node and edge potentials/marginals
    m_node_potentials = array_2d<num>(num_nodes, num_states);
    m_edge_potentials = array_3d<num>(num_edges, num_states, num_states);
    m_node_marginals = array_2d<num>(num_nodes, num_states);
    m_edge_marginals = array_3d<num>(num_edges, num_states, num_states);

  }
    

  void set_theta(array_1d<num> theta){
    m_theta = theta;
  }

  void set_node_potentials(){
    cout<<"num_nodes: "<<num_nodes<<endl;
    cout<<"num_states: "<<num_states<<endl;
    cout<<"num_node_features"<<num_node_features<<endl;
    for(int i = 0; i < num_nodes; i++){
      for(int j = 0; j < num_states; j++){
	num temp = 0;
	for(int k = 0; k < num_node_features; k++){
	  temp = temp + m_node_features(i,k) * m_theta(m_node_map(j,k));
	  m_node_potentials(i,j) = exp(temp);
	}
      }
    }
  }

  void set_edge_potentials(){
    for(int i = 0; i < num_edges; i++){
      for(int j = 0; j < num_states; j++){
	for(int k = 0; k < num_states; k++){
	  num temp = 0;
	  for(int l = 0; l < num_edge_features; l++){
	    temp = temp + m_edge_features(i,l) * m_theta(m_edge_map(j,k,l));
	  }
	  cout<<"z"<<i<<j<<k<<" "<<temp<<endl;
	  m_edge_potentials(i,j,k) = exp(temp);
	}
      }
    }
    //cout<<m_edge_potentials<<endl;
  }

  

  void set_node_marginals(){
    // make sure prev_marginals is normalized to begin with.  actually doesn't matter if all values are scaled up uniformly
    array_2d<num> prev_marginals(num_nodes, num_states, 1.0 / (double)num_states);
    array_2d<num> next_marginals(num_nodes, num_states);
    for(int i = 0; i < mean_field_max_iter; i++){
      for(int j = 0; j < num_nodes; j++){
	array_1d<int> neighbors = node_to_neighbors[j];
	for(int k = 0; k < num_states; k++){
	  num temp = get_node_potential(j,k);
	  for(int l = 0; l < neighbors.len1; l++){
	    int nbr = neighbors(l);
	    for(int m = 0; m < num_states; m++){
	      temp += get_node_potential(l,m) * get_edge_potential(j,nbr,k,m);
	    }
	  }
	  next_marginals(j,k) = exp(temp);
	}
	// normalize
	num sum = 0;
	for(int n = 0; n < num_states; n++){
	  sum += next_marginals(j,n);
	}
	for(int n = 0; n < num_states; n++){
	  next_marginals(j,n) /= sum;
	}
      }
      swap(prev_marginals, next_marginals);
    }
    m_node_marginals = next_marginals;

    // set edge marginals
    cout<<m_edge_marginals<<endl;
    for(int i = 0; i < num_edges; i++){
      int node1 = pos_to_edge[i].first;
      int node2 = pos_to_edge[i].second;
      for(int j = 0; j < num_states; j++){
	for(int k = 0; k < num_states; k++){
	  cout<<"OOOO"<<endl;
	  cout<<node1<<node2<<j<<k<<num_states<<endl;
	  
	  m_edge_marginals(i,j,k) = m_node_marginals(node1,j) * m_node_marginals(node2,k);
	}
      }
    }
  }

  
  num get_likelihood(array_1d<int> states){
    num temp = 0;
    for(int i = 0; i < num_nodes; i++){
      cout<<i<<endl;
      temp *= get_node_potential(i,states(i));
    }
    cout<<"FFFFFFFF"<<endl;
    for(int i = 0; i < num_edges; i++){
      int node1 = pos_to_edge[i].first;
      int node2 = pos_to_edge[i].second;
      cout<<node1<<" "<<node2<<" "<<states(node1)<<" "<<states(node2)<<endl;
      temp *= get_edge_potential(node1, node2, states(node1), states(node2));
    }
    return temp;
  }

};

#endif
