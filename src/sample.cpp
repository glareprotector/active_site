#include "sample.h"
// #include "model.h"
#include "Score.h"
#include "model.cpp"

num sample::get_node_potential(int node, int state){
  return node_potentials(node, state);
}

num sample::get_edge_potential(int node1, int node2, int state1, int state2){
  if(node1 < node2){
    return edge_potentials(edge_to_pos[pair<int,int>(node1, node2)], state1, state2);
  }
  else{
    return get_edge_potential(node2, node1, state2, state1);
  }
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
  this->node_to_neighbors = hash_map<int, arbi_array<int> >();
  this->edge_to_pos = hash_map< pair<int, int>, int, pair_hash>();
  this->pos_to_edge = hash_map< int, pair<int,int> >();

  // assuming that i < j in edge_to_ij.  but this is irrelevant since want each edge to be present twice
  for(int i = 0; i < this->num_edges; i++){
    int node1 = _edges(i,0); // should use hashmap version?
    int node2 = _edges(i,1); // ditto
    this->pos_to_edge[i] = pair<int,int>(node1, node2);
    this->edge_to_pos[pair<int,int>(node1, node2)] = i;
    if(this->node_to_neighbors.find(node1) != this->node_to_neighbors.end()){
      this->node_to_neighbors[node1].append(node2);
    }
    else{
      this->node_to_neighbors[node1] = arbi_array<int>(1,1);
      this->node_to_neighbors[node1](0) = node2;
    }
    if(this->node_to_neighbors.find(node2) != this->node_to_neighbors.end()){
      this->node_to_neighbors[node2].append(node1);
    }
    else{
      this->node_to_neighbors[node2] = arbi_array<int>(1,1);
      this->node_to_neighbors[node2](0) = node1;
    }
  }
    
  // allocate node and edge potentials/marginals
  this->node_potentials = arbi_array<num>(2, this->num_nodes, _p_model->num_states);
  this->edge_potentials = arbi_array<num>(3, this->num_edges, _p_model->num_states, _p_model->num_states);
  this->node_marginals = arbi_array<num>(2, this->num_nodes, _p_model->num_states);
  this->edge_marginals = arbi_array<num>(3, this->num_edges, _p_model->num_states, _p_model->num_states);

  // put this here for now
  this->mean_field_max_iter = 100;
}

void sample::set_node_potentials(){
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      num temp = 0;
      for(int k = 0; k < p_model->num_node_features; k++){
	temp = temp + node_features(i,k) * p_model->theta(p_model->node_map(j,k));
      }
      // node_potentials(i,j) = exp(temp); old
      node_potentials(i,j) = temp;
    }
  }
  //cout<<node_potentials<<endl;

}

void sample::set_edge_potentials(){
  for(int i = 0; i < num_edges; i++){
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	num temp = 0;
	for(int l = 0; l < p_model->num_edge_features; l++){
	  temp = temp + edge_features(i,l) * p_model->theta(p_model->edge_map(j,k,l));
	}
	// edge_potentials(i,j,k) = exp(temp); old
	edge_potentials(i,j,k) = temp;
      }
    }
  }
}

// storing marginals in regular non-log form
void sample::set_marginals(){
  // make sure prev_marginals is normalized to begin with.  actually doesn't matter if all values are scaled up uniformly
  arbi_array<num>* p_prev_marginals = new arbi_array<num>(2, num_nodes, p_model->num_states);
  p_prev_marginals->fill(1.0 / (double)p_model->num_states);
  arbi_array<num>* p_next_marginals = new arbi_array<num>(2, num_nodes, p_model->num_states);
  for(int i = 0; i < mean_field_max_iter; i++){
    for(int j = 0; j < num_nodes; j++){
      assert(node_to_neighbors.find(j) != node_to_neighbors.end()); // assuming graph is connected
      arbi_array<int> neighbors = node_to_neighbors[j];
      for(int k = 0; k < p_model->num_states; k++){
	num temp = 0;
	//num log_temp = LogScore_ZERO;
	for(int l = 0; l < neighbors.size(0); l++){
	  int nbr = neighbors(l);
	  for(int m = 0; m < p_model->num_states; m++){
	    temp += (*p_prev_marginals)(l,m) * (get_edge_potential(j,nbr,k,m)); // error here? should have log of that potential?
	    //LogScore_PLUS_EQUALS(log_temp, log( (*p_prev_marginals)(l,m) ) + get_edge_potential(j,nbr,k,m));
	  }
	}
	temp += get_node_potential(j,k);
	(*p_next_marginals)(j,k) = exp(temp);
	if((*p_next_marginals)(j,k) < 0){
	  cout<<"negative!: "<<(*p_next_marginals)(j,k)<<endl;
	}
	//(*p_next_marginals)(j,k) = exp(get_node_potential(j,k)) * exp(exp(log_temp));
      }
      // normalize
      num sum = 0;
      for(int n = 0; n < p_model->num_states; n++){
	sum += (*p_next_marginals)(j,n);
      }
      // don't normalize if sum is 0
      if(fabs(sum) > 1e-20){
	for(int n = 0; n < p_model->num_states; n++){
	  (*p_next_marginals)(j,n) /= sum;
	}
      }
      else{
	// marginals are uniform
	for(int k = 0; k < p_model->num_states; k++){
	  (*p_next_marginals)(j,k) = 1.0 / (num)p_model->num_states;
	}
      }
    }
    swap(p_prev_marginals, p_next_marginals);
  }
  node_marginals = *p_next_marginals;
  
  // set edge marginals
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge[i].first;
    int node2 = pos_to_edge[i].second;
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
  cout<<endl<<"BEFORE GRADIENT"<<endl;
  cout<<gradient<<endl;
  gradient.fill(0);
  for(int i = 0; i < num_nodes; i++){
    for(int k = 0; k < p_model->num_node_features; k++){
      gradient(p_model->node_map(true_states(i),k)) += node_features(i,k);
      for(int j = 0; j < p_model->num_states; j++){
	//cout<<node_features(i,k)<<" "<<node_marginals(i,j)<<endl;
	gradient(p_model->node_map(j,k)) -= node_features(i,k) * node_marginals(i,j);
      }
    }
  }

  for(int i = 0; i < num_edges; i++){
    for(int k = 0; k < p_model->num_edge_features; k++){
      int node1 = pos_to_edge[i].first;
      int node2 = pos_to_edge[i].second;
      gradient(p_model->edge_map(true_states(node1),true_states(node2),k)) += edge_features(i,k);
      for(int j = 0; j < p_model->num_states; j++){
	for(int l = 0; l < p_model->num_states; l++){
	  gradient(p_model->edge_map(j,l,k)) -= edge_features(i,k) * edge_marginals(i,j,l);
	}
      }
    }
  }

  gradient.scale(-1.0);

  return gradient;
}

num sample::get_log_Z(){
  cout<<endl<<"theta"<<endl;
  cout<<p_model->theta<<endl;

  //cout<<"node_marginals"<<endl;
  //cout<<node_marginals<<endl;
  num energy = 0;
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      energy -= node_marginals(i,j) * get_node_potential(i,j);
    }
  }
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge[i].first;
    int node2 = pos_to_edge[i].second;
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	energy -= edge_marginals(i,j,k) * get_edge_potential(node1,node2,j,k);
      }
    }
  }
  num entropy = 0;
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      entropy -= node_marginals(i,j) * log(node_marginals(i,j));
    }
  }
  cout<<"entropy: "<<entropy<<endl;
  cout<<"energy: "<<energy<<endl;
  cout<<endl<<"log Z: "<<(entropy - energy)<<endl;
  return entropy - energy;
}
  

// returns negative log likelihood
num sample::get_likelihood(){
  //num log_temp = LogScore_ZERO;
  //cout<<node_potentials<<endl;
  num temp = 0;
  for(int i = 0; i < num_nodes; i++){
    temp += get_node_potential(i,true_states(i));
    //cout<<temp<<endl;
    //LogScore_PLUS_EQUALS(log_temp, log( log(get_node_potential(i,true_states(i))));
  }
  
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge[i].first;
    int node2 = pos_to_edge[i].second;
    //cout<<node1<<node2<<"nodes!!"<<endl;
    temp += get_edge_potential(node1, node2, true_states(node1), true_states(node2));
    //LogScore_PLUS_EQUALS(log_temp, log(   log(get_edge_potential(node1, node2, true_states(node1), true_states(node2)))    )  );
    //cout<<temp<<endl;
  }
  //cout<<endl<<"config potential: "<<temp<<endl;
  //return temp - get_log_Z();

  //num log_NLL = log (get_log_Z()), log_temp
  cout<<get_log_Z()<<" asdf"<<temp;
  return get_log_Z() - temp;
}

//#include model.cpp
