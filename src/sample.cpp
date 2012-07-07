#include "sample.h"
#include "Score.h"
#include "model.cpp"

num sample::get_node_potential(arbi_array<num>& node_potentials, int node, int state){
  return node_potentials(node, state);
}

num sample::get_edge_potential(arbi_array<num>& edge_potentials, int node1, int node2, int state1, int state2){
  return edge_potentials(edge_to_pos(node1, node2), state1, state2);
}

sample::sample(){
  int x;
}

sample::sample(model* _p_model, arbi_array<num> _node_features, arbi_array<num> _edge_features, arbi_array<int> _edges, arbi_array<int> _true_states, string _folder, string _pdb_name, string _chain_letter){
  this->p_model = _p_model;
  this->node_features = _node_features;
  this->edge_features = _edge_features;
  this->num_nodes = _node_features.size(0);
  this->num_edges = _edge_features.size(0);
  this->true_states = _true_states;
  this->folder = _folder;
  this->pdb_name = _pdb_name;
  this->chain_letter = _chain_letter;

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

arbi_array<num> sample::get_node_potentials(arbi_array<num> theta){
  arbi_array<num> node_potentials(2, num_nodes, p_model->num_states);
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      num temp = 0;
      for(int k = 0; k < p_model->num_node_features; k++){
	temp = temp + node_features(i,k) * p_model->theta(p_model->node_map(j,k));
      }
      node_potentials(i,j) = temp;
    }
  }
  return node_potentials;
}


void sample::get_edge_potentials(arbi_array<num> theta){
  arbi_array<num> edge_potentials(3, num_edges, p_model->num_states, p_model->num_states);
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
  return edge_potentials;
}

void sample::get_marginals(arbi_array<num> theta, arbi_array<num>& node_marginals, arbi_array<num>& edge_marginals){

  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentials = get_edge_potentials(theta);
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;
  get_marginals(node_potentials, edge_potentials, node_marginals, edge_marginals);

}

// storing marginals in regular non-log form
// not assuming that node_marginals and edge_marginals are pre-allocated
void sample::get_marginals(arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num>& node_marginals, arbi_array<num>& edge_marginals){

  // allocate node_marginals and edge_marginals
  node_marginals = arbi_array<num>(2, num_nodes, p_model->num_states);
  edge_marginals = arbi_array<num>(3, num_edges, p_model->num_states, p_model->num_states);

  arbi_array<num> log_new_marginals(1,p_model->num_states);

  // set initial marginals to normalized potentials
  for(int i = 0; i < num_nodes; i++){
    num log_sum = 0;
    for(int j = 0; j < p_model->num_states; j++){      
      if(j == 1){
	log_sum = LogScore_ADD(get_node_potential(node_potentials,i,0), get_node_potential(node_potentials,i,1));
      }
      if(j > 1){
	LogScore_PLUS_EQUALS(log_sum, get_node_potential(node_potentials,i,j));
      }
    }

    num norm = 0;
    for(int j = 0; j < p_model->num_states; j++){
      node_marginals(i,j) = exp(get_node_potential(node_potentials,i,j) - log_sum);
      norm += node_marginals(i,j);
    }

    for(int j = 0; j < p_model->num_states; j++){
      node_marginals(i,j) /= norm;
    }

  }

  for(int i = 0; i < p_model->mean_field_max_iter; i++){
    for(int j = 0; j < num_nodes; j++){

      // calculate unnormalized log marginals
      log_new_marginals.fill(0);
      for(int k = 0; k < p_model->num_states; k++){
	num temp = 0;
	for(int l = 0; l < node_to_neighbors(j).size(0); l++){
	  int nbr = node_to_neighbors(j)(l);
	  for(int m = 0; m < p_model->num_states; m++){
	    log_new_marginals(k) += (node_marginals)(nbr,m) * (get_edge_potential(edge_potentials,j,nbr,k,m));
	  }
	}
	log_new_marginals(k) += get_node_potential(node_potentials,j,k);
      }

      // normalize in log space
      num log_sum = 0;
      for(int n = 0; n < p_model->num_states; n++){
	if(n == 1){
	  log_sum = LogScore_ADD(log_new_marginals(0), log_new_marginals(1));
	}
	if(n > 1){
	  LogScore_PLUS_EQUALS(log_sum, log_new_marginals(n));
	}
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
arbi_array<num> sample::get_data_likelihood_gradient(arbi_array<num> node_marginals, arbi_array<num> edge_marginals){

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

num sample::get_log_Z(arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num> node_marginals, arbi_array<num> edge_marginals){

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
      if(node_marginals(i,j) > 1e-80){
	entropy -= node_marginals(i,j) * log(node_marginals(i,j));
      }
    }
  }
  assert(isfinite(entropy));
  assert(isfinite(energy));

  return entropy - energy;
}
  

num sample::get_data_potential(arbi_array<num> node_potentials, arbi_array<num> edge_potentials){

  num temp = 0;
  for(int i = 0; i < num_nodes; i++){
    temp += get_node_potential(node_potentials, i,true_states(i));
  }
  
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge(i).first;
    int node2 = pos_to_edge(i).second;
    temp += get_edge_potential(edge_potentials, node1, node2, true_states(node1), true_states(node2));
  }
  return temp;
}


// returns negative log likelihood
num sample::get_data_likelihood(arbi_array<num> theta){

  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentials = get_edge_potentials(theta);
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;
  get_marginals(node_potentials, edge_potentials, node_marginals, edge_marginals);

  num log_z = get_log_Z(node_potentials, edge_potentials, node_marginals, edge_marginals);
  num data_potential = get_data_potential(node_potentials, edge_potentials);

  assert(isfinite(log_z));
  assert(isfinite(data_potential));

  num ans = log_z - data_potential;
  //if(ans <= 0){
  if(false){
    cout<<setprecision(32);
    cout<<this->folder<<" log_z: "<<log_z<<" config: "<<config_likelihood<<endl;
    cout<<"ans: "<<ans<<endl;
    cout<<p_model->theta<<endl;
    //cout<<endl<<"NODE FEATURES"<<endl;
    //cout<<node_features<<endl;
    //cout<<endl<<"NODE pOTENTIALS:"<<endl;
    //cout<<node_potentials<<endl;
    //cout<<endl<<"NODE MARGINALS:"<<endl;
    //cout<<node_marginals<<endl;
    cout<<this->folder<<" log_z: "<<log_z<<" config: "<<config_likelihood<<endl;
    p_model->theta.write(string("theta_shorter.csv"), ',');
    arbi_array<num> to_write = arbi_array<num>::transpose(node_features);
    to_write.write(this->folder + string("XnodeNormed.csv"), ',');
    cout<<"THETA: "<<p_model->theta<<endl;
    cout<<"FFFFFFFFFFFFF"<<endl;
    assert(false);
    }
  //assert(ans >= 0);

  return ans;
}

num sample::get_L(int which_obj, arbi_array<num>& theta){
  switch(which_obj){
  case 0:
    return get_data_likelihood(theta);
    break;
  case 1:
    return get_L_expected_distance(theta);
    break;
  }
}

num sample::smooth_f(num x){
  return exp(x*x);
}
 

num sample::get_L_expected_distance(arbi_array<num> theta){

  num loss = 0;
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;
  get_marginals(theta, node_marginals, edge_marginals);
  for(int i = 0 ; i < num_nodes; i++){
    arbi_array<num> sorted_distances, sorted_pos;
    sorted_distances_getter.get(pdb_name, chain_letter, i, sorted_pos, sorted_distances);
    // starting at position 0 to next to last, height from i(noninclusive) to i+1 is p(0 to i all not active).  so multiple the height by width which is d(i+1) - d(i)
    num cumulative = 1.0;
    num exp_val = 0;
    num closest_site_dist = -1;
    for(int j = 1; j < num_nodes; j++){
      cumulative *= node_marginals(sorted_pos(j-1), 0);
      exp_val += cumulative * (sorted_distances(j) - sorted_distances(j-1));
    }
    // find closest site
    for(int j = 0; j < num_nodes; j++){
      if(true_states(sorted_pos(j)) == 1){
	closest_site_dist = sorted_distances(sorted_pos(j));
	break;
      }
    }
    // is it safe to assume that every chain has at least one active site?  if not, set "true" closest distance to some big number
    assert(closest_site_dist != -1);
    loss += smooth_f(fabs(exp_val - closest_site_dist));
  }

  return loss;
}

num sample::d_smooth_f(num x){
  return exp(x*x) * 2x;
}
      
arbi_array<num> sample::get_dL_dMu_expected_distance(arbi_array<num> node_marginals, arbi_array<num> edge_marginals, arbi_array<num> dL_dNode_Mu, arbi_array<num> dL_dEdge_Mu){

  dL_dNode_Mu = arbi_array<num>(2, num_nodes, 2);
  assert(p_model->num_states == 2);
  dL_dNode_Mu.fill(0);
  for(int i = 0; i < num_nodes; i++){
    arbi_array<num> sorted_distances, sorted_pos;
    sorted_distances_getter.get(pdb_name, chain_letter, i, sorted_pos, sorted_distances);
    num cumulative = 1.0;
    num temp, first, der, closest_site_dist = -1;
    for(int j = 0; j < num_nodes; j++){
      if(true_states(sorted_pos(j)) == 1){
	closest_site_dist = sorted_distances(sorted_pos(j));
	break;
      }
    }
    assert(closest_site_dist != -1);
    for(int j = 1; j < num_nodes; j++){
      cumulative *= node_marginals(sorted_pos(j-1), 0);
      temp = (sorted_distances(j) - sorted_distances(j-1)) * cumulative;
      first = d_smooth_f(temp - closest_site_dist);
      for(int k = 0; k < j; k++){
	
	der = temp / node_marginals(sorted_pos(k));
	dL_dNode_Mu(sorted_pos(k),0) += first * der;
      }
      
    }
  }
  // gradient for 1 - marginal is negative of dL/dMu
  for(int i = 0; i < num_nodes; i++){
    dL_dNode_Mu(i,1) = -1.0 * dL_dNode_Mu(i,0);
  }
  // loss function doesn't depend on edge marginals, so dL_Edge_dMu should be zeros
  dL_dEdge_Mu = arbi_array<num>(num_edges, 2, 2);
  dL_dEdge_Mu.fill(0);

}
      
	


arbi_array<num> sample::get_data_likelihood_gradient(arbi_array<num> theta){
  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentails = get_edge_potentials(theta);
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;
  get_marginals(node_potentials, edge_potentials, node_marginals, edge_marginals);
  return get_data_likelihood_gradient(node_marginals, edge_marginals);
}

arbi_array<num> get_dL_dTheta(int which_obj, arbi_array<num> theta){
  
  if(which_obj == 0){
    return get_data_likelihood_gradient(theta);
  }
  else{
    return get_dL_dTheta_Perturb(which_obj, theta);
  }
}
    

void get_dPot_dTheta(arbi_array<num> theta, arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num>& dNode, arbi_array<num>& dEdge){
  
  dNode = arbi_array<num>(3, num_nodes, p-model->num_states, p_model->theta_length);
  dEdge = arbi_array<num>(4, num_edges, p_model->num_states, p_model->num_states, p_model->theta_length);

  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_node_features; k++){
	dNode(i,j,p_model->node_map(i,j,k)) = node_potentials(i,j) * node_features(i,k);
      }
    }
  }

  for(int i = 0; i < num_edges; i++){
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	for(int l = 0; l < p_model->num_edge_features; l++){
	  dEdge(i,j,k,p_model->edge_map(i,j,k,l)) = edge_potentials(i,j,k) * edge_features(i,l);
	}
      }
    }
  }

}

// implements the method of Domke
arbi_array<num> sample::get_dL_dTheta_Perturb(int which_obj, arbi_array<num> theta){
  
  arbi_array<num> dNode_Pot_dTheta;
  arbi_array<num> dEdge_Pot_dTheta;
  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentials = get_edge_potentials(theta);
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;
  get_marginals(node_potentials, edge_potentials, node_marginals, edge_marginals);
  arbi_array<num> dL_dNode_Mu;
  arbi_array<num> dL_dEdgeMu;

  get_dL_dMu(which_obj, node_marginals, edge_marginals, dL_dNode_Mu, dL_dEdgeMu);
  get_dPot_dTheta(theta, node_potentials, edge_potentials, dNode_Pot_dTheta, dEdge_Pot_dTheta);
  
  // calculate perturbed marginals
  num r = 0.01;
  dL_dNode_Mu.scale(r);
  dL_dEdge_Mu.scale(r);
  arbi_array<num> perturbed_node_potentials = node_potentials + dL_dNode_Mu;
  arbi_array<num> perturbed_edge_marginals = edge_potentials + dL_dEdge_Mu;
  arbi_array<num> perturbed_node_marginals;
  arbi_array<num> perturbed_edge_marginals;
  get_marginals(perturbed_node_potentials, perturbed_edge_potentials, perturbed_node_marginals, perturbed_edge_marginals);
  
  // calculate dL/dPot
  node_marginals.scale(-1.0);
  edge_marginals.scale(-1.0);
  perturbed_node_marginals = perturbed_node_marginals + node_marginals;
  perturbed_edge_marginals = perturbed_edge_marginals + edge_marginals;
  perturbed_node_marginals.scale(1.0 / r);
  perturbed_edge_marginals.scale(1.0 / r);
  arbi_array<num> dL_dNode_Pot = perturbed_node_marginals;
  arbi_array<num> dL_dEdge_Pot = perturbed_edge_marginals;

  // multiply dL/dPot and dPot/dTheta
  arbi_array<num> dL_dTheta(1, p_model->theta_length);
  for(int i = 0; i < p_model->theta_length; i++){
    num val = 0;
    for(int j = 0; j < num_nodes; j++){
      for(int k = 0; k < p_model->num_states; k++){
	val += dL_dNode_Pot(j,k) * dNode_Pot_dTheta(j,k,i);
      }
    }
    for(int j = 0; j < num_edges; j++){
      for(int k = 0; k < p_model->num_states; k++){
	for(int l = 0; l < p_model->num_states; l++){
	  val += dL_dEdge_Pot(j,k,l) * dEdge_Pot_dTheta(j,k,l,i);
	}
      }
    }
    dL_dTheta(i) = val;
  }

  return dL_dTheta;
}

 void sample::get_dL_dMu(int which_obj, arbi_array<num> node_marginals, arbi_array<num> edge_marginals, arbi_array<num> dL_dNode_Mu, arbi_array<num> dL_dEdge_Mu){
   switch(which_obj){
   case 1:
     get_dL_dMu_expected_distance(node_marginals, edge_marginals, dL_dNode_Mu, dL_dEdge_Mu);
     break;
   }
 }
