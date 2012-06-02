#include "sample.h"
// #include "model.h"
#include "Score.h"
#include "model.cpp"

num sample::get_node_potential(int node, int state){
  return node_potentials(node, state);
}

num sample::get_edge_potential(int node1, int node2, int state1, int state2){
  return edge_potentials(edge_to_pos(node1, node2), state1, state2);
}

sample::sample(){
  int x;
}

sample::sample(model* _p_model, arbi_array<num> _node_features, arbi_array<num> _edge_features, arbi_array<int> _edges, arbi_array<int> _true_states, string folder){
  this->p_model = _p_model;
  this->node_features = _node_features;
  this->edge_features = _edge_features;
  this->num_nodes = _node_features.size(0);
  this->num_edges = _edge_features.size(0);
  this->true_states = _true_states;
  this->folder = folder;

  // initialize hash_maps
  //this->node_to_neighbors = hash_map<int, arbi_array<int> >();
  //this->edge_to_pos = hash_map< pair<int, int>, int, pair_hash>();
  //this->pos_to_edge = hash_map< int, pair<int,int> >();

  this->node_to_neighbors = arbi_array< arbi_array<int> >(1, this->num_nodes);
  this->pos_to_edge = arbi_array< pair<int,int> >(1, this->num_edges);
  this->edge_to_pos = arbi_array<int>(2, this->num_nodes, this->num_nodes);

  for(int i = 0; i < this->num_edges; i++){
    int node1 = _edges(i,0); // should use hashmap version?
    int node2 = _edges(i,1); // ditto
    assert(node1 > node2);
    this->pos_to_edge(i) = pair<int,int>(node1, node2);

    this->edge_to_pos(node1,node2) = i;
    this->edge_to_pos(node2,node1) = i;
    this->node_to_neighbors(node1).append(node2);
    this->node_to_neighbors(node2).append(node1);
  }
    
  // allocate node and edge potentials/marginals
  this->node_potentials = arbi_array<num>(2, this->num_nodes, _p_model->num_states);
  this->edge_potentials = arbi_array<num>(3, this->num_edges, _p_model->num_states, _p_model->num_states);
  this->node_marginals = arbi_array<num>(2, this->num_nodes, _p_model->num_states);
  this->edge_marginals = arbi_array<num>(3, this->num_edges, _p_model->num_states, _p_model->num_states);

}

void sample::set_node_potentials(){
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      num temp = 0;
      for(int k = 0; k < p_model->num_node_features; k++){
	temp = temp + node_features(i,k) * p_model->theta(p_model->node_map(j,k));
      }
      node_potentials(i,j) = temp;
    }
  }

}

void sample::set_edge_potentials(){
  for(int i = 0; i < num_edges; i++){
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	num temp = 0;
	for(int l = 0; l < p_model->num_edge_features; l++){
	  temp = temp + edge_features(i,l) * p_model->theta(p_model->edge_map(j,k,l));
	}
	edge_potentials(i,j,k) = temp;
      }
    }
  }
}

// storing marginals in regular non-log form
void sample::set_marginals(){

  // set initial marginals to normalized potentials
  for(int i = 0; i < num_nodes; i++){
    num sum = 0;
    
    for(int j = 0; j < p_model->num_states; j++){
      (node_marginals)(i,j) = exp(get_node_potential(i,j));
      sum += (node_marginals)(i,j);
    }
    for(int j = 0; j < p_model->num_states; j++){
      (node_marginals)(i,j) /= sum;
    }
  }

  arbi_array<num> log_new_marginals(1,p_model->num_states);
  for(int i = 0; i < p_model->mean_field_max_iter; i++){
    for(int j = 0; j < num_nodes; j++){

      // calculate unnormalized log marginals
      log_new_marginals.fill(0);
      for(int k = 0; k < p_model->num_states; k++){
	num temp = 0;
	for(int l = 0; l < node_to_neighbors(j).size(0); l++){
	  int nbr = node_to_neighbors(j)(l);
	  for(int m = 0; m < p_model->num_states; m++){
	    log_new_marginals(k) += (node_marginals)(nbr,m) * (get_edge_potential(j,nbr,k,m));
	  }
	}
	log_new_marginals(k) += get_node_potential(j,k);
      }

      // normalize in log space
      num log_sum = 0;
      for(int n = 0; n < p_model->num_states; n++){
	LogScore_PLUS_EQUALS(log_sum, log_new_marginals(n));
      }
      for(int n = 0; n < p_model->num_states; n++){
	log_new_marginals(n) -= log_sum;
      }
      for(int n = 0; n < p_model->num_states; n++){
	(node_marginals)(j,n) = exp(log_new_marginals(n));
      }
      
      // normalize again due to rounding errors
      num sum = 0;
      for(int n = 0; n < p_model->num_states; n++){
	sum += (node_marginals)(j,n);
      }
      for(int n = 0; n < p_model->num_states; n++){
	(node_marginals)(j,n) /= sum;
      }
    }
  }
  
  // set edge marginals
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge(i).first;
    int node2 = pos_to_edge(i).second;
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	edge_marginals(i,j,k) = node_marginals(node1,j) * node_marginals(node2,k);
      }
    }
  }
}

// remember to take the negative of gradient
arbi_array<num> sample::get_gradient(){

  arbi_array<num> gradient = arbi_array<num>(1, p_model->theta_length);
  gradient.fill(0);

  // fill in gradient components for node parameters
  for(int i = 0; i < num_nodes; i++){
    for(int k = 0; k < p_model->num_node_features; k++){
      gradient(p_model->node_map(true_states(i),k)) += node_features(i,k);
      for(int j = 0; j < p_model->num_states; j++){
	gradient(p_model->node_map(j,k)) -= node_features(i,k) * node_marginals(i,j);
      }
    }
  }

  // fill in gradient components for edge parameters
  for(int i = 0; i < num_edges; i++){
    for(int k = 0; k < p_model->num_edge_features; k++){
      int node1 = pos_to_edge(i).first;
      int node2 = pos_to_edge(i).second;
      gradient(p_model->edge_map(true_states(node1),true_states(node2),k)) += edge_features(i,k);
      for(int j = 0; j < p_model->num_states; j++){
	for(int l = 0; l < p_model->num_states; l++){
	  gradient(p_model->edge_map(j,l,k)) -= edge_features(i,k) * edge_marginals(i,j,l);
	}
      }
    }
  }

  gradient.scale(-1.0);
  
  for(int i = 0; i < p_model->theta_length; i++){
    assert(isfinite(gradient(i)));
  }

  return gradient;
}

num sample::get_log_Z(){

  num energy = 0;
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      energy -= node_marginals(i,j) * get_node_potential(i,j);
    }
  }
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge(i).first;
    int node2 = pos_to_edge(i).second;
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	energy -= edge_marginals(i,j,k) * get_edge_potential(node1,node2,j,k);
      }
    }
  }
  num entropy = 0;
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      if(node_marginals(i,j) > 1e-10){
	entropy -= node_marginals(i,j) * log(node_marginals(i,j));
      }
    }
  }

  assert(isfinite(entropy));
  assert(isfinite(energy));

  return entropy - energy;
}
  

num sample::get_config_likelihood(){

  num temp = 0;
  for(int i = 0; i < num_nodes; i++){
    temp += get_node_potential(i,true_states(i));
  }
  
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge(i).first;
    int node2 = pos_to_edge(i).second;
    temp += get_edge_potential(node1, node2, true_states(node1), true_states(node2));
  }
  return temp;
}


// returns negative log likelihood
num sample::get_likelihood(){

  num log_z = get_log_Z();
  num config_likelihood = get_config_likelihood();

  assert(isfinite(log_z));
  assert(isfinite(config_likelihood));

  num ans = log_z - config_likelihood;
  assert(ans > 0);

  return ans;
}


