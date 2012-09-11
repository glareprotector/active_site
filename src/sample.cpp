#include "sample.h"
#include "Score.h"
#include "model.cpp"

num sample::get_node_potential(arbi_array<num>& node_potentials, int node, int state){
  return node_potentials(node, state);
}

num sample::get_edge_potential(arbi_array<num>& edge_potentials, int node1, int node2, int state1, int state2){
  if(node1 > node2){
    return edge_potentials(edge_to_pos(node1, node2), state1, state2);
  }
  else{
    return edge_potentials(edge_to_pos(node1, node2), state2, state1);
  }
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
  //this->node_potentials = arbi_array<num>(2, this->num_nodes, _p_model->num_states);
  //this->edge_potentials = arbi_array<num>(3, this->num_edges, _p_model->num_states, _p_model->num_states);
  //this->node_marginals = arbi_array<num>(2, this->num_nodes, _p_model->num_states);
  //this->edge_marginals = arbi_array<num>(3, this->num_edges, _p_model->num_states, _p_model->num_states);
  this->times_called = 0;

}


void sample::simulate_states(arbi_array<num> f_theta){
  // here fake modify the true_states.

  // this one puts residues close to true cat sites as 1
  for(int i = 0; i < true_states.size(0); i++){
    if(true_states(i) == 1){
      for(int j = i - 4; j< i + 4; j++){
	if(j >= 0 && j < true_states.size(0)-1){
	  if(rand()%3 == 0){
	    true_states(j) = 1;
	  }
	}
      }
    }
  }




  // this reads in fake weight vector, generates node potentials, decides states based on those
  
  arbi_array<num> f_node_potentials = get_node_potentials(f_theta);
  arbi_array<num> f_edge_potentials = get_edge_potentials(f_theta);
  arbi_array<num> f_node_marginals;
  arbi_array<num> f_edge_marginals;

  get_marginals_BP(f_node_potentials, f_edge_potentials, f_node_marginals, f_edge_marginals);
  for(int i = 0; i < num_nodes; i++){
    if(num(rand() % 1000) / 1000.0 > f_node_marginals(i,1)){
      true_states(i) = 1;
    }
    else{
      true_states(i) = 0;
    }
  }

}

arbi_array<num> sample::get_node_potentials(arbi_array<num> theta){
  arbi_array<num> node_potentials(2, num_nodes, p_model->num_states);

  node_potentials.fill(0);
  for(int i = 0; i < num_nodes; i++){
    //cout<<i<<" "<<node_features(i,3)<<endl;
    for(int j = 0; j < p_model->num_states; j++){
      num temp = 0;
      for(int k = 0; k < p_model->num_node_features; k++){
	temp = temp + node_features(i,k) * theta(p_model->node_map(j,k));
	//cout<<temp<<" ";
	//if(i == num_nodes-1 && j == 0){
	//if(k == p_model->num_node_features - 1 && j == 0){
	//	  cout<<k<<' '<<node_features(i,k)<<' '<<theta(p_model->node_map(j,k))<<' '<<temp<<" "<<num_nodes<<" "<<i<<" "<<node_features.size(0)<<" "<<node_features.size(1)<<endl;
	//}
      }
      //cout<<"nodepot: "<<proc_id<<" "<<i<<" "<<j<<temp<<endl;
      node_potentials(i,j) = temp;
      if(isfinite(temp) == false){
	for(int k = 0; k < p_model->num_node_features; k++){
	  //cout<<"POPOPOPOPOPOPOPOPPOPO: "<<k<<' '<<node_features(i,k)<<' '<<theta(p_model->node_map(j,k))<<' '<<temp<<" "<<num_nodes<<" "<<i<<" "<<node_features.size(0)<<" "<<node_features.size(1)<<endl;
	}
	//cout<<node_potentials<<endl;
	assert(false);
	//exit(1);
      }
      if(isfinite(temp) == false){

	for(int k = 0; k < p_model->num_node_features; k++){
	  cout<<k<<' '<<node_features(i,k)<<' '<<' '<<temp<<" "<<num_nodes<<" "<<i<<" "<<node_features.size(0)<<" "<<node_features.size(1)<<endl;
	}
	//cout<<node_potentials<<endl;
	assert(false);
	exit(1);
      }
    }
  }

  //cout<<"INSIDENODEPOTENTIALS INTEREST"<<proc_id<<" "<<pdb_name<<" "<<node_potentials(globals::fdsa,0)<<" "<<node_potentials(globals::fdsa,1);


  //cout<<endl<<"FFFFF "<<node_features(num_nodes-1, p_model->num_node_features-1)<<endl;
  //cout<<node_features<<endl;
  //cout<<node_potentials<<endl;
  return node_potentials;
}


arbi_array<num> sample::get_edge_potentials(arbi_array<num> theta){
  arbi_array<num> edge_potentials(3, num_edges, p_model->num_states, p_model->num_states);
  for(int i = 0; i < num_edges; i++){
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	num temp = 0;
	for(int l = 0; l < p_model->num_edge_features; l++){
	  temp = temp + edge_features(i,l) * theta(p_model->edge_map(j,k,l));
	  //cout<<edge_features(i,l)<<" "<<theta(p_model->edge_map(j,k,l))<<endl;
	}
	edge_potentials(i,j,k) = temp;
	if(isfinite(edge_potentials(i,j,k)) == false){
	  cout<<edge_features(i,0)<<" "<<theta(p_model->edge_map(j,k,0))<<" "<<p_model->edge_map(j,k,0)<<endl;
	  cout<<edge_features(i,1)<<" "<<theta(p_model->edge_map(j,k,1))<<" "<<p_model->edge_map(j,k,1)<<endl;
	  cout<<theta<<endl;
	  cout<<"quitting"<<endl;
	  assert(false);
	}
      }
    }
  }
  return edge_potentials;
}

void sample::get_marginals(arbi_array<num> theta, arbi_array<num>& node_marginals, arbi_array<num>& edge_marginals){
  //cout<<"\n\n\n\n\n\nthetehtehtehtehth: "<<endl<<theta<<endl;
  //  exit(1);



  //cout<<"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF "<<theta<<endl;
  #ifndef SERIAL
  
  #endif

  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentials = get_edge_potentials(theta);
  //cout<<"INTEREST NONE: "<<proc_id<<" "<<pdb_name<<" "<<node_potentials(globals::fdsa,0)<<" "<<node_potentials(globals::fdsa,1)<<endl;
  get_marginals(node_potentials, edge_potentials, node_marginals, edge_marginals);

  #ifndef SERIAL
  
  #endif

  //cout<<"GETMARGINALSWITHTHETA "<<proc_id<<" "<<pdb_name<<endl;

}

// which infer method to used is stored in model
void sample::get_marginals(arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num>& node_marginals, arbi_array<num>& edge_marginals){
  //cout<<"which_infer: "<<p_model->which_infer<<endl;


  //cout<<"\nGETMARGINALSWITHPOTPOTPOTPOTPOTPOTPOT "<<proc_id<<" "<<pdb_name<<endl;

  //cout<<"INTEREST11: "<<proc_id<<" "<<pdb_name<<" "<<node_potentials(globals::fdsa,0)<<" "<<node_potentials(globals::fdsa,1)<<endl;


  switch(p_model->which_infer){
  case 0:
    get_marginals_mean_field(node_potentials, edge_potentials, node_marginals, edge_marginals);
    break;
  case 1:
    get_marginals_BP(node_potentials, edge_potentials, node_marginals, edge_marginals);
    break;
  case 2:
    get_marginals_logistic_regression(node_potentials, node_marginals);
  }
  /*
  //cout<<node_marginals<<endl;
  arbi_array<num> node_marginals1;
  arbi_array<num> edge_marginals1;
  get_marginals_mean_field(node_potentials, edge_potentials, node_marginals1, edge_marginals1);
  //node_marginals = node_marginals1;
  //edge_marginals = edge_marginals1;
  for(int i = 0; i < num_nodes; i++){
    if(fabs(node_marginals(i,0)-node_marginals1(i,0)) > .1)cout<<node_marginals(i,0)<<" "<<node_marginals1(i,0)<<endl;
  }
  for(int i = 0; i < num_edges; i++){
    //cout<<edge_marginals(i,0,0)<<" "<<edge_marginals1(i,0,0)<<" "<<edge_marginals(i,0,1)<<" "<<edge_marginals1(i,0,1)<<edge_marginals(i,1,0)<<" "<<edge_marginals1(i,1,0)<<" "<<edge_marginals(i,1,1)<<" "<<edge_marginals1(i,1,1)<<endl;
  }
  //exit(1);
  */
}


void sample::get_marginals_logistic_regression(arbi_array<num> node_potentials, arbi_array<num>& node_marginals){
  // set initial marginals to normalized potentials
  node_marginals = arbi_array<num>(2, num_nodes, p_model->num_states);
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
      if(isfinite(node_marginals(i,j)) == false){
	cout<<"GGG "<<pdb_name<<endl;
	//cout<<node_marginals<<endl;
	cout<<"SEPEPEPEPEP"<<endl;
	cout<<node_potentials<<endl;
	assert(false);
	exit(1);
      }
    }
   
  }
  //cout<<node_marginals<<endl;
  //exit(1);
}


// storing marginals in regular non-log form
// not assuming that node_marginals and edge_marginals are pre-allocated
void sample::get_marginals_mean_field(arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num>& node_marginals, arbi_array<num>& edge_marginals){





  


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

    //cout<<proc_id<<" "<<i<<" norm "<<norm<<" "<<node_marginals(i,0)<<" "<<node_marginals(i,1)<<" "<<node_potentials(i,0)<<" "<<node_potentials(i,1)<<endl;
    
    for(int j = 0; j < p_model->num_states; j++){
      node_marginals(i,j) /= norm;
    }

  }



   for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      if(!isfinite(node_marginals(i,j))){
	cout<<node_marginals<<endl<<"before"<<endl;
	cout<<node_marginals(i,j)<<endl;
	for(int k = 0; k < p_model->num_node_features; k++){
	  cout<<node_features(i,k)<<" ";
	}
	cout<<endl;
	cout<<"proc_id: "<<proc_id<<" "<<i<<" "<<j<<endl;
	assert(false);
      }

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

  // check that all marginals are finite
  #ifndef SERIAL
  if(proc_id == 5) {
  #endif
    //cout<<"ZZZZZZZZZZZZZZ"<<endl<<node_potentials<<endl<<node_marginals<<endl;
  #ifndef SERIAL
  }
  #endif

  //cout<<node_potentials<<endl;

  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      if(!isfinite(node_marginals(i,j))){
	cout<<"AAAAAAAAAAAAAAAA"<<proc_id<<endl;
      }
    }
  }
  
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      if(!isfinite(node_marginals(i,j))){
	cout<<node_marginals<<endl;
	assert(false);
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

arbi_array<num> sample::get_feature_values(arbi_array<int> states){

  arbi_array<num> feature_values(1, p_model->theta_length);
  feature_values.fill(0);

  for(int i = 0; i < p_model->num_node_features; i++){
    for(int j = 0; j < num_nodes; j++){
      feature_values(p_model->node_map(states(j), i)) = feature_values(p_model->node_map(states(j), i)) + node_features(j,i);
    }
  }

  for(int i = 0; i < p_model->num_edge_features; i++){
    for(int j = 0; j < num_edges; j++){
      int u = pos_to_edge(j).first;
      int v = pos_to_edge(j).second;
      feature_values(p_model->edge_map(states(u),states(v),i)) += edge_features(j,i);
    }
  }

  return feature_values;
}


arbi_array<num> sample::get_pseudo_likelihood_gradient(arbi_array<num> theta){
  arbi_array<num> node_pseudos;
  pseudo_likelihood_helper(theta, node_pseudos);
  // have to compute feature value for each configuration that is 1 away from the true_states.  do this for each feature first
  // compute node_features for true_states
  arbi_array<num> grad(1, p_model->theta_length);
  grad.fill(0);
  
  assert(theta.size(0) == p_model->theta_length);
  
  for(int i = 0; i < num_nodes; i++){
    arbi_array<num> exp_grad(1, theta.size(0));
    exp_grad.fill(0);
    for(int k = 0; k < p_model->num_states; k++){
      arbi_array<int> states = true_states;
      states(i) = k;
      arbi_array<num> temp = get_feature_values(states);
      temp.scale(node_pseudos(i,k));
      exp_grad = exp_grad + temp;
    }
    exp_grad.scale(-1.0);
    grad = grad + exp_grad;
    grad = grad + get_feature_values(true_states);
  }
  grad.scale(-1.0);
  return grad;
}

// remember to take the negative of gradient
arbi_array<num> sample::get_data_likelihood_gradient(arbi_array<num> node_marginals, arbi_array<num> edge_marginals){

  arbi_array<num> gradient = arbi_array<num>(1, p_model->theta_length);
  gradient.fill(0);

  // cout<<node_marginals<<endl;

  // fill in gradient components for node parameters
  for(int i = 0; i < num_nodes; i++){
    for(int k = 0; k < p_model->num_node_features; k++){
      gradient(p_model->node_map(true_states(i),k)) += node_features(i,k);
      for(int j = 0; j < p_model->num_states; j++){
	gradient(p_model->node_map(j,k)) -= node_features(i,k) * node_marginals(i,j);
	//if(j == 0 && k == 0){
	//  cout<<j<<" "<<gradient(p_model->node_map(j,k))<<endl;
	//}
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
    //if(!isfinite(gradient(i))) {cout<<proc_id<<" "<<gradient(i)<<" ";}
    assert(isfinite(gradient(i)));
  }
  //assert(false);
  //cout<<"GRADIENT"<<gradient<<endl;
  //assert(false);
  return gradient;
}

num sample::get_log_Z(arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num> node_marginals, arbi_array<num> edge_marginals){

  num energy = 0;
  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      if(!isfinite(node_marginals(i,j))){
	cout<<i<<" "<<j<<" "<<node_marginals(i,j)<<endl;
	cout<<this->node_features<<endl;
	cout<<node_marginals<<endl;
	cout<<pdb_name<<" "<<pdb_name<<endl;
	assert(false);
      }
      energy -= node_marginals(i,j) * get_node_potential(node_potentials, i,j);
    }
  }
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge(i).first;
    int node2 = pos_to_edge(i).second;
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	assert(isfinite(edge_marginals(i,j,k)));

	energy -= edge_marginals(i,j,k) * get_edge_potential(edge_potentials, node1,node2,j,k);
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
  /*
  for(int i = 0; i < num_edges; i++){
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	if(edge_marginals(i,j,k) > 1e-80){
	  entropy -= edge_marginals(i,j,k) * log(edge_marginals(i,j,k));
	}
      }
    }
  }*/

  assert(isfinite(entropy));
  assert(isfinite(energy));
  num af = 0;
  for(int i = 0; i < node_features.size(0); i++){
    for(int j = 0; j < node_features.size(1); j++){
      af += node_features(i,j);
    }
  }
  //cout<<endl<<"Z: "<<proc_id<<" "<<entropy - energy<<" "<<af<<" GGGGGGGGG "<<node_features(num_nodes-1, p_model->num_node_features-1)<<endl;

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

  //cout<<"pot: "<<temp<<endl;
  return temp;
}


// returns negative log likelihood
num sample::get_data_likelihood(arbi_array<num> theta){




  //cout<<"333333333333"<<theta<<endl;
  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentials = get_edge_potentials(theta);
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;
  //cout<<theta<<endl;
  //cout<<node_potentials<<endl;

  get_marginals(node_potentials, edge_potentials, node_marginals, edge_marginals);
  //cout<<node_potentials<<endl;
  num log_z = get_log_Z(node_potentials, edge_potentials, node_marginals, edge_marginals);
  num data_potential = get_data_potential(node_potentials, edge_potentials);

  assert(isfinite(log_z));
  assert(isfinite(data_potential));
  //cout<<"log_z: "<<log_z<<endl;
  num ans = log_z - data_potential;
  //assert(ans > 0);
  if(ans <= 0){
    //if(false){
    cout<<setprecision(32);
    cout<<this->folder<<" log_z: "<<log_z<<" config: ";//<<config_likelihood<<endl;
    cout<<"ans: "<<ans<<endl;
    //cout<<theta<<endl;
    //cout<<endl<<"NODE FEATURES"<<endl;
    //cout<<node_features<<endl;
    //cout<<endl<<"NODE pOTENTIALS:"<<endl;
    //cout<<node_potentials<<endl;
    cout<<endl<<"NODE MARGINALS:"<<endl;
    cout<<node_marginals<<endl;
    cout<<this->folder<<" log_z: "<<log_z<<endl;//" config: "<<config_likelihood<<endl;
    //p_model->theta.write(string("theta_shorter.csv"), ',');
    arbi_array<num> to_write = arbi_array<num>::transpose(node_features);
    to_write.write(this->folder + string("XnodeNormed.csv"), ',');
    cout<<"THETA: "<<theta<<endl;
    cout<<"FFFFFFFFFFFFF"<<endl;
    //exit(1);
    //assert(false);
    }
  //assert(ans >= 0);
  /*
  arbi_array<num> asdf(1,3);
  asdf(0) = ans;
  asdf(1) = log_z;
  asdf(2) = data_potential;
  asdf.write(string("likelihood.txt"));
  theta.write(string("theta_shorter.csv"), ',');
  cout<<endl<<"THETA: "<<theta<<"LIKELIHOOD: "<<ans<<endl;
  node_potentials.write(string("node_potentials.csv"));
  arbi_array<num> exp_node_potentials(2,node_potentials.size(0), node_potentials.size(1));
  for(int i = 0; i < exp_node_potentials.size(0); i++){
    for(int j = 0; j < exp_node_potentials.size(1); j++){
      exp_node_potentials(i,j) = exp(node_potentials(i,j));
    }
  }
  cout<<exp_node_potentials<<endl;
  cout<<"edge_potentials"<<endl;
  cout<<exp(edge_potentials(0,0,0))<<" "<<exp(edge_potentials(0,0,1))<<" "<<exp(edge_potentials(0,1,0))<<" "<<exp(edge_potentials(0,1,1))<<endl;
  cout<<node_marginals<<endl;
  cout<<folder<<endl;
  if(fabs(theta(0)) > .000001){
    assert(false);
  }
  //  assert(false);
  */
  return ans;
}

num sample::get_L(int which_obj, arbi_array<num>& theta){

  num temp = 0;
  for(int i = 0; i < theta.size(0); i++){
    temp += theta(i);
  }
  //cout<<endl<<"E: "<<proc_id<<" "<<temp<<endl;
  switch(which_obj){
  case 0:
    return get_data_likelihood(theta);
    break;
  case 1:
    return get_L_expected_distance(theta);
    break;
  case 2:
    return get_L_nodewise(theta);
    break;
  case 3:
    return get_L_pseudo(theta);
    break;
  }
}

num sample::smooth_f(num x){
  //return exp(x);
  return -1.0 / (1.0 + x*x);
  return x*x;
  return exp(x*x);
}

num sample::get_L_pseudo(arbi_array<num> theta){

  arbi_array<num> node_pseudos;
  pseudo_likelihood_helper(theta, node_pseudos);
  num L = 0;
  for(int i = 0; i < num_nodes; i++){
    L += log(node_pseudos(i,true_states(i)));
  }
  if(isfinite(L) == false){
    cout<<node_pseudos<<endl;
    cout<<theta;
    cout<<endl<<L<<endl;
    exit(1);
  }
  return -1.0 * L;
}
    
 

num sample::get_L_nodewise(arbi_array<num> theta){

  // set param
  cpp_caller::set_param(globals::pParams, string("pdb_name"), &pdb_name, globals::STRING_TYPE);
  cpp_caller::set_param(globals::pParams, string("chain_letter"), &chain_letter, globals::STRING_TYPE);


  num loss = 0;
  arbi_array<num> node_marginals, edge_marginals;
  //cout<<"NODEWISE "<<proc_id<<" "<<pdb_name<<endl;
  get_marginals(theta, node_marginals, edge_marginals);

  //cout<<"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB "<<proc_id<<" theta: "<<theta<<endl;


  PyObject* pFakeTrueNum = cached_obj_getter::call_wrapper(string("new_new_objects"), string("bhW"), globals::pParams, globals::recalculate, true, true, false);
  arbi_array<num> fake_true_num = cpp_caller::py_float_list_to_cpp_num_vect(pFakeTrueNum);
  Py_DECREF(pFakeTrueNum);
  for(int i = 0; i < num_nodes; i++){
    arbi_array<int> sorted_pos;
    arbi_array<num> sorted_distances;
    //sorted_distances_getter::get(pdb_name, chain_letter, i, sorted_pos, sorted_distances);
    //cout<<i<<" "<<loss<<endl;
    loss += smooth_f(node_marginals(i,1) - fake_true_num(i));
  }
  if(isfinite(loss) == false){
    cout<<node_marginals<<endl;
    cout<<fake_true_num<<endl;
    cout<<theta<<endl;
    cout<<pdb_name<<endl;
    exit(1);
  }
  return loss;
}

void sample::get_dL_dMu_nodewise(arbi_array<num> node_marginals, arbi_array<num> edge_marginals, arbi_array<num>& dL_dNode_Mu, arbi_array<num>& dL_dEdge_Mu){


  // set param for pdb_name
  cpp_caller::set_param(globals::pParams, string("pdb_name"), &pdb_name, globals::STRING_TYPE);
  cpp_caller::set_param(globals::pParams, string("chain_letter"), &chain_letter, globals::STRING_TYPE);


  dL_dNode_Mu = arbi_array<num>(2, num_nodes, 2);
  dL_dNode_Mu.fill(0);
  PyObject* pFakeTrueNum = cached_obj_getter::call_wrapper(string("new_new_objects"), string("bhW"), globals::pParams, globals::recalculate, true, true, false);
  arbi_array<num> fake_true_num = cpp_caller::py_float_list_to_cpp_num_vect(pFakeTrueNum);
  Py_DECREF(pFakeTrueNum);
  for(int i = 0; i < num_nodes; i++){
    dL_dNode_Mu(i,1) = d_smooth_f(node_marginals(i,1) - fake_true_num(i));
    //dL_dNode_Mu(i,0) = -1.0 * dL_dNode_Mu(i,1);
  }
  //cout<<"\nDLDMU: "<<proc_id<<" "<<pdb_name<<dL_dNode_Mu(globals::fdsa,0)<<" "<<dL_dNode_Mu(globals::fdsa,1)<<" "<<node_marginals(globals::fdsa,1)<<" "<<fake_true_num(globals::fdsa)<<endl;
  //cout<<"\nfake true: "<<proc_id<<" faketruesize: "<<fake_true_num.size(0)<<endl;
  for(int i = 0; i < num_nodes; i++){
    //cout<<proc_id<<i<<" "<<fake_true_num(i)<<endl;
  }
  //cout<<endl<<"DLDMU"<<endl;
  //cout<<dL_dNode_Mu<<endl;
  dL_dEdge_Mu = arbi_array<num>(3, num_edges, 2, 2);
  dL_dEdge_Mu.fill(0);

}

num sample::get_L_expected_distance(arbi_array<num> theta){


  //exit(1);
  num loss = 0;
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;
  get_marginals(theta, node_marginals, edge_marginals);


  cpp_caller::set_param(globals::pParams, string("pdb_name"), &pdb_name, globals::STRING_TYPE);
  cpp_caller::set_param(globals::pParams, string("chain_letter"), &chain_letter, globals::STRING_TYPE);

  //  PyObject* pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("apW"), globals::pParams, globals::recalculate, false, false, false);
  PyObject* pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("apW"), globals::pParams, globals::recalculate, false, false, false);

  PyObject* pSorted_Dist = PyList_GetItem(pResult, 0);
  PyObject* pSorted_Pos = PyList_GetItem(pResult, 1);
  arbi_array<num> sorted_distances = cpp_caller::py_float_mat_to_cpp_num_mat(pSorted_Dist);
  arbi_array<int> sorted_pos = cpp_caller::py_int_mat_to_cpp_int_mat(pSorted_Pos);
  Py_DECREF(pResult);



  PyObject* pClosests = cached_obj_getter::call_wrapper(string("new_new_objects"), string("aqW"), globals::pParams, globals::recalculate, false, false, false);
  arbi_array<num> closest_dists = cpp_caller::py_float_list_to_cpp_num_vect(pClosests);
  


  Py_DECREF(pClosests);



  for(int i = 0 ; i < num_nodes; i++){
    
    //    cout<<"loss: "<<i<<endl;


    


    // starting at position 0 to next to last, height from i(noninclusive) to i+1 is p(0 to i all not active).  so multiple the height by width which is d(i+1) - d(i)
    num cumulative = 1.0;
    num exp_val = 0;
    num closest_site_dist = closest_dists(i);

    for(int j = 1; j < sorted_distances.size(1); j++){
      cumulative *= node_marginals(sorted_pos(i,j-1), 0);
      exp_val += cumulative * (sorted_distances(i,j) - sorted_distances(i,j-1));
    }

    loss += smooth_f(fabs(exp_val - closest_site_dist));

  }

  return loss;





}

num sample::d_smooth_f(num x){
  //return exp(x);
  return 2*x / ((x*x + 1) * (x*x + 1));
  return 2.0 * x;
  return exp(x*x) * 2.0 * x;
}
      
void sample::get_dL_dMu_expected_distance(arbi_array<num> node_marginals, arbi_array<num> edge_marginals, arbi_array<num>& dL_dNode_Mu, arbi_array<num>& dL_dEdge_Mu){




  assert(p_model->num_states == 2);
  dL_dNode_Mu = arbi_array<num>(2, num_nodes, 2);
  dL_dNode_Mu.fill(0);


  cpp_caller::set_param(globals::pParams, string("pdb_name"), &pdb_name, globals::STRING_TYPE);
  cpp_caller::set_param(globals::pParams, string("chain_letter"), &chain_letter, globals::STRING_TYPE);

  //PyObject* pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("apW"), globals::pParams, globals::recalculate, false, false, false);
  PyObject* pResult = cached_obj_getter::call_wrapper(string("new_new_objects"), string("apW"), globals::pParams, globals::recalculate, false, false, false);

  PyObject* pSorted_Dist = PyList_GetItem(pResult, 0);
  PyObject* pSorted_Pos = PyList_GetItem(pResult, 1);

  arbi_array<num> sorted_distances = cpp_caller::py_float_mat_to_cpp_num_mat(pSorted_Dist);
  arbi_array<int> sorted_pos = cpp_caller::py_int_mat_to_cpp_int_mat(pSorted_Pos);
  Py_DECREF(pResult);

  PyObject* pClosests = cached_obj_getter::call_wrapper(string("new_new_objects"), string("aqW"), globals::pParams, globals::recalculate, false, false, false);
  arbi_array<num> closest_dists = cpp_caller::py_float_list_to_cpp_num_vect(pClosests);
  Py_DECREF(pClosests);

  
  for(int i = 0; i < num_nodes; i++){
        
    //    cout<<"loss grad: "<<i<<endl;
    



    num cumulative = 1.0;
    num temp;

    num closest_site_dist = closest_dists(i);



    // for each term in summation for the site, calculate gradient (will be 0 for sites further away)
    num inside = 0;
    arbi_array<num> node_seconds(1, num_nodes);
    node_seconds.fill(0.0);



    for(int j = 1; j < sorted_pos.size(1); j++){
      cumulative *= node_marginals(sorted_pos(i,j-1), 0);
      temp = (sorted_distances(i,j) - sorted_distances(i,j-1)) * cumulative;
      inside += temp;
      assert(temp >= 0);
      for(int k = 0; k < j; k++){
	node_seconds(sorted_pos(i,k)) += temp / node_marginals(sorted_pos(i,k), 0);
      }      
    }
    //cout<<node_seconds<<endl;
    node_seconds.scale(d_smooth_f(inside - closest_site_dist));
    for(int j = 0; j < num_nodes; j++){
      dL_dNode_Mu(j,0) += node_seconds(j);
    }


  }

  // loss function doesn't depend on edge marginals, so dL_Edge_dMu should be zeros
  dL_dEdge_Mu = arbi_array<num>(3, num_edges, 2, 2);
  dL_dEdge_Mu.fill(0);
  

  

}
      
	
// function that takes in theta and computes node/edge potentials.  then for each node, for each node/each state, looks at neighbor true states to compute node's pseudolikelihood
void sample::pseudo_likelihood_helper(arbi_array<num> theta, arbi_array<num>& node_pseudos){

  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentials = get_edge_potentials(theta);
  
  //node_pseudos = arbi_array<num>(2, num_nodes, p_model->num_states);
  //node_pseudos.fill(0);
  arbi_array<num> node_pseudos1 = arbi_array<num>(2, num_nodes, p_model->num_states);
  node_pseudos1.fill(0);


  
  for(int i = 0; i < num_nodes; i++){
    for(int k = 0; k < p_model->num_states; k++){
      node_pseudos1(i,k) += node_potentials(i,k);
      arbi_array<int> nbrs = node_to_neighbors(i);
      for(int j = 0; j < nbrs.size(0); j++){
	int n = nbrs(j);
	node_pseudos1(i,k) += get_edge_potential(edge_potentials, i, n, k, true_states(n));
      }
    }
  }
 
  //cout<<node_pseudos1<<endl;


  // normalize proportional exp of their weights
  for(int i = 0; i < num_nodes; i++){
    num log_sum = 0;
    for(int j = 0; j < p_model->num_states; j++){      
      if(j == 1){
	log_sum = LogScore_ADD(node_pseudos1(i,0), node_pseudos1(i,1));
      }
      if(j > 1){
	LogScore_PLUS_EQUALS(log_sum, node_pseudos1(i,j));
      }
    }
    num norm = 0;
    for(int j = 0; j < p_model->num_states; j++){
      node_pseudos1(i,j) = exp(node_pseudos1(i,j) - log_sum);
      norm += node_pseudos1(i,j);
    }

    for(int j = 0; j < p_model->num_states; j++){
      node_pseudos1(i,j) /= norm;
    }
  }
  

  
  node_pseudos = node_pseudos1;

  //cout<<node_pseudos<<endl;
    //cout<<node_features<<endl;
    //cout<<"FSDFDFDF"<<endl;
  //exit(1);


  /*
  for(int i = 0; i < num_nodes; i++){
    arbi_array<int> states = true_states;
    for(int k = 0; k < p_model->num_states; k++){
      states(i) = k;
      arbi_array<num> feature_values = get_feature_values(states);
      for(int l = 0; l < p_model->theta_length; l++){
	node_pseudos(i,k) += feature_values(l) * theta(l);
      }
    }
    }

  //for(int i = 0; i < num_nodes; i++){
  //cout<<node_pseudos(i,1)-node_pseudos(i,0)<<" "<<node_pseudos1(i,1) - node_pseudos1(i,0)<<endl;
    //}


  

   // normalize proportional exp of their weights
  for(int i = 0; i < num_nodes; i++){
    num log_sum = 0;
    for(int j = 0; j < p_model->num_states; j++){      
      if(j == 1){
	log_sum = LogScore_ADD(node_pseudos(i,0), node_pseudos(i,1));
      }
      if(j > 1){
	LogScore_PLUS_EQUALS(log_sum, node_pseudos(i,j));
      }
    }
    num norm = 0;
    for(int j = 0; j < p_model->num_states; j++){
      node_pseudos(i,j) = exp(node_pseudos(i,j) - log_sum);
      norm += node_pseudos(i,j);
    }

    for(int j = 0; j < p_model->num_states; j++){
      node_pseudos(i,j) /= norm;
    }
    }


  
  //for(int i = 0; i < num_nodes; i++){
  //  cout<<node_pseudos(i,0)<<" "<<node_pseudos1(i,0)<<endl;
  //}
  */

}


arbi_array<num> sample::get_data_likelihood_gradient(arbi_array<num> theta){
  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentials = get_edge_potentials(theta);
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;
  //cout<<theta<<endl;
  //cout<<node_potentials<<endl;
  num af = 0;
  for(int i = 0; i < node_potentials.linear_length; i++){
    af += node_potentials.m_data[i];
  }
  //cout<<endl<<"G: "<<af<<endl<<theta<<endl;

  get_marginals(node_potentials, edge_potentials, node_marginals, edge_marginals);
  arbi_array<num> grad = get_data_likelihood_gradient(node_marginals, edge_marginals);
  /*if(fabs(theta(0)) > .00001){
    cout<<endl<<"GRADIENT"<<endl;
    cout<<grad<<endl;
    cout<<endl<<"THETA"<<endl;
    cout<<theta<<endl;
    theta.write(string("theta_shorter.csv"), ',');
    assert(false);
    }*/
  return grad;
}

arbi_array<num> sample::get_dL_dTheta(int which_obj, arbi_array<num> theta){
  
  if(which_obj == 0){
    return get_data_likelihood_gradient(theta);
  }
  else if(which_obj == 3){
    return get_pseudo_likelihood_gradient(theta);
  }
  else{
    return get_dL_dTheta_Perturb(which_obj, theta);
  }
}
    
// potential here refers to the theta in exp(theta * features)
void sample::get_dPot_dTheta(arbi_array<num> theta, arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num>& dNode, arbi_array<num>& dEdge){
  
  dNode = arbi_array<num>(3, num_nodes, p_model->num_states, p_model->theta_length);
  dEdge = arbi_array<num>(4, num_edges, p_model->num_states, p_model->num_states, p_model->theta_length);
  dNode.fill(0);
  dEdge.fill(0);


  for(int i = 0; i < num_nodes; i++){
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_node_features; k++){
	//cout<<i<<" "<<j<<" "<<k<<" "<<node_potentials(i,j)<<" "<<node_features(i,k)<<" "<<p_model->node_map(j,k)<<endl;
	//dNode(i,j,p_model->node_map(j,k)) = node_potentials(i,j) * node_features(i,k);
	assert(fabs(dNode(i,j,p_model->node_map(j,k))) < 0.000001);
	dNode(i,j,p_model->node_map(j,k)) = node_features(i,k);
      }
    }
  }

  for(int i = 0; i < num_edges; i++){
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	for(int l = 0; l < p_model->num_edge_features; l++){
	  dEdge(i,j,k,p_model->edge_map(j,k,l)) = edge_potentials(i,j,k) * edge_features(i,l);
	  dEdge(i,j,k,p_model->edge_map(j,k,l)) = edge_features(i,l);
	}
      }
    }
  }
  //cout<<dNode<<endl;
}

num& sample::get_message(arbi_array<num>& msgs, int i, int j, int s){
  // convention is that if i > j, store in row 0
  
  //cout<<proc_id<<" "<<i<<" "<<j<<endl;
  
  assert(i != j);
  if(i > j){
    return msgs(edge_to_pos(i,j), 0, s);
  }
  else{
    return msgs(edge_to_pos(i,j), 1, s);
  }
}
    

void sample::get_marginals_BP(arbi_array<num> node_potentials, arbi_array<num> edge_potentials, arbi_array<num>& node_marginals, arbi_array<num>& edge_marginals){
  // set aside data structure for old msgs and new msgs, and initialize all messages to 1
  /*static int times_called = 0;
  //static arbi_array<num> msgs1(3, num_edges, 2, p_model->num_states);
  //static arbi_array<num> msgs2(3, num_edges, 2, p_model->num_states);
  static arbi_array<num>* old_msgs;
  static arbi_array<num>* new_msgs;
  if(times_called == 0){
    msgs1 = arbi_array<num>(3, num_edges, 2, p_model->num_states);
    msgs2 = arbi_array<num>(3, num_edges, 2, p_model->num_states);
    old_msgs = &msgs1;
    new_msgs = &msgs2;
    (*new_msgs).fill(1.0);
  }
  (*new_msgs).fill(1.0);
  times_called++;
  */

  if(times_called == 0){
    msgs1 = arbi_array<num> (3, num_edges, 2, p_model->num_states);
    msgs2 = arbi_array<num> (3, num_edges, 2, p_model->num_states);
  }
  
  if(times_called == 0){
    old_msgs = &msgs1;
    new_msgs = &msgs2;
    (*new_msgs).fill(1.0);
  }
  //cout<<(*new_msgs)<<endl;
  times_called++;
  //(*old_msgs).fill(0.0);
  int bp_max_iter = 15;
  for(int i = 0; i < bp_max_iter; i++){
    swap(old_msgs, new_msgs);
    for(int j = 0; j < num_edges; j++){
      // edges go both ways
      int node1 = pos_to_edge(j).first;
      int node2 = pos_to_edge(j).second;
      num sum = 0;
      arbi_array<int> node1_nbrs = node_to_neighbors(node1);
      for(int k = 0; k < p_model->num_states; k++){
	get_message(*new_msgs, node1, node2, k) = 0;
	for(int l = 0; l < p_model->num_states; l++){
	  //	  cout<<get_edge_potential(edge_potentials, node1,node2,l,k);
	  num temp = exp(get_edge_potential(edge_potentials, node1,node2,l,k) + node_potentials(node1,l));
	  for(int m = 0; m < node1_nbrs.size(0); m++){
	    int nbr = node1_nbrs(m);
	    if(nbr != node2){
	      temp *= get_message(*old_msgs, nbr, node1, l);
	    }
	  }
	  //cout<<temp<<" ";
	  get_message(*new_msgs, node1, node2, k) += temp;
	  //cout<<get_message(new_msgs, node1, node2, k)<<" ";
	  sum += temp;
	}
      }
      for(int k = 0; k < p_model->num_states; k++){
	get_message(*new_msgs, node1, node2, k) /= sum;
	
      }
      
      //cout<<sum<<endl;
      // do the same for reverse edge.  everything the same except node1 and node2 switched
      swap(node1, node2);
      sum = 0;
      node1_nbrs = node_to_neighbors(node1);
      for(int k = 0; k < p_model->num_states; k++){
	get_message(*new_msgs, node1, node2, k) = 0;
	for(int l = 0; l < p_model->num_states; l++){
	  num temp = exp(get_edge_potential(edge_potentials, node1,node2,l,k) + node_potentials(node1,l));
	  for(int m = 0; m < node1_nbrs.size(0); m++){
	    int nbr = node1_nbrs(m);
	    if(nbr != node2){
	      temp *= get_message(*old_msgs, nbr, node1, l);
	    }
	  }
	  get_message(*new_msgs, node1, node2, k) += temp;
	  sum += temp;
	}
      }
      for(int k = 0; k < p_model->num_states; k++){
	get_message(*new_msgs, node1, node2, k) /= sum;
	//cout<<(*new_msgs);
	//cout<<get_message(*new_msgs, node1, node2, k)<<" "<<(*new_msgs).linear_length<<endl;
      }
      //cout<<*new_msgs<<endl;
    }
    //cout<<old_msgs<<endl;
    /*for(int a = 0; a < num_edges; a++){
      cout<<get_message(new_msgs, pos_to_edge(a).first, pos_to_edge(a).second, 0)<<" ";
      cout<<get_message(new_msgs, pos_to_edge(a).first, pos_to_edge(a).second, 0)<<" ";
      cout<<get_message(new_msgs, pos_to_edge(a).second, pos_to_edge(a).first, 1)<<" ";
      cout<<get_message(new_msgs, pos_to_edge(a).second, pos_to_edge(a).first, 1)<<" ";
    //cout<<new_msgs<<endl;
    }*/
  }
  // now use the messages to obtain marginals.  these messages are messages in original clique graph
  // only difference is that we went 2 steps in recursion before going to previous round of messages
  node_marginals = arbi_array<num>(2, num_nodes, p_model->num_states);
  node_marginals.fill(0);
  edge_marginals = arbi_array<num>(3, num_edges, p_model->num_states, p_model->num_states);
  edge_marginals.fill(0);
  //cout<<old_msgs<<endl;
  for(int i = 0; i < num_nodes; i++){
    num sum = 0;
    for(int j = 0; j < p_model->num_states; j++){
      num temp = exp(node_potentials(i,j));
      arbi_array<int> neighbors = node_to_neighbors(i);
      for(int k = 0; k < neighbors.size(0); k++){
	//cout<<temp<<" ";
	temp *= get_message(*new_msgs, neighbors(k), i, j);
      }
      node_marginals(i,j) = temp;
      sum += temp;
      //cout<<temp<<endl;
    }
    for(int j = 0; j < p_model->num_states; j++){
      node_marginals(i,j) /= sum;
    }
  }
  for(int i = 0; i < num_edges; i++){
    int node1 = pos_to_edge(i).first;
    int node2 = pos_to_edge(i).second;
    num sum = 0;
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	num temp = exp(node_potentials(node1,j) + edge_potentials(i,j,k) + node_potentials(node2,k));
	arbi_array<int> node1_neighbors = node_to_neighbors(node1);
	arbi_array<int> node2_neighbors = node_to_neighbors(node2);
	for(int l = 0; l < node1_neighbors.size(0); l++){
	  if(l != node2){
	    temp *= get_message(*new_msgs, node1_neighbors(l), node1, j);
	     }
	}
	for(int l = 0; l < node2_neighbors.size(0); l++){
	  if(l != node1){
	    temp *= get_message(*new_msgs, node2_neighbors(l), node2, k);
	     }
	}
	edge_marginals(i,j,k) = temp;
	sum += temp;
      }
    }
    // normalize edge marginals
    for(int j = 0; j < p_model->num_states; j++){
      for(int k = 0; k < p_model->num_states; k++){
	edge_marginals(i,j,k) /= sum;
      }
    }
  }

  /*for(int i = 0; i < num_nodes; i++){
    cout<<node_marginals(i,0)+node_marginals(i,1)<<endl;
    }*/

  //cout<<edge_marginals<<endl;
  //cout<<edge_marginals<<endl;
  //assert(false);
  //exit(1);
}

  


// implements the method of Domke
arbi_array<num> sample::get_dL_dTheta_Perturb(int which_obj, arbi_array<num> theta){
  
  arbi_array<num> dNode_Pot_dTheta;
  arbi_array<num> dEdge_Pot_dTheta;
  arbi_array<num> node_potentials = get_node_potentials(theta);
  arbi_array<num> edge_potentials = get_edge_potentials(theta);
  arbi_array<num> node_marginals;
  arbi_array<num> edge_marginals;



  //cout<<"\nINSIDENODEPOTENTIALSPERTURB INTEREST "<<proc_id<<" "<<pdb_name<<" "<<node_potentials(globals::fdsa,0)<<" "<<node_potentials(globals::fdsa,1);


  get_marginals(node_potentials, edge_potentials, node_marginals, edge_marginals);
  arbi_array<num> dL_dNode_Mu;
  arbi_array<num> dL_dEdge_Mu;

  get_dL_dMu(which_obj, node_marginals, edge_marginals, dL_dNode_Mu, dL_dEdge_Mu);
  get_dPot_dTheta(theta, node_potentials, edge_potentials, dNode_Pot_dTheta, dEdge_Pot_dTheta);

  // figure out what r should be
  num theta_max = max(node_potentials.max(), edge_potentials.max());
  num dL_dMu_max = max(dL_dNode_Mu.max(), dL_dEdge_Mu.max());
  num r = globals::eps * (1 + theta_max) / dL_dMu_max;
  //  cout<<"r: "<<r<<endl;
  
  /*if(!isfinite(r)){
    //cout<<node_potentials<<endl;
    cout<<dL_dNode_Mu<<endl;
    cout<<theta<<endl;
    cout<<"proc_id: "<<proc_id<<" theta_max:"<<theta_max<<" dl_dmu_max:"<<dL_dMu_max<<endl;
    assert(false);
    r = 1e-12;
    }*/

  r = 1e-8;
  // calculate perturbed marginals
  dL_dNode_Mu.scale(r);
  dL_dEdge_Mu.scale(r);
  //cout<<"NODEPOTENTIALS"<<endl;
  //cout<<node_potentials<<endl;
  //cout<<endl<<"dldnodemu"<<endl;
  //cout<<dL_dNode_Mu<<endl;
  arbi_array<num> perturbed_node_potentials = node_potentials + dL_dNode_Mu;
  arbi_array<num> perturbed_edge_potentials = edge_potentials + dL_dEdge_Mu;
  arbi_array<num> perturbed_node_marginals;
  arbi_array<num> perturbed_edge_marginals;

  //cout<<"\nPERTURB POTENTIAL SICK"<<proc_id<<" "<<pdb_name<<" "<<perturbed_node_potentials(globals::fdsa,0)<<" "<<perturbed_node_potentials(globals::fdsa,1)<<" "<<dL_dNode_Mu(globals::fdsa,0)<<" "<<dL_dNode_Mu(globals::fdsa,1)<<endl;
  


  get_marginals(perturbed_node_potentials, perturbed_edge_potentials, perturbed_node_marginals, perturbed_edge_marginals);
  
  //cout<<endl<<"perturbed node marginals"<<endl;
  //cout<<perturbed_node_marginals<<endl;

  // calculate dL/dPot
  node_marginals.scale(-1.0);
  edge_marginals.scale(-1.0);
  perturbed_node_marginals = perturbed_node_marginals + node_marginals;
  perturbed_edge_marginals = perturbed_edge_marginals + edge_marginals;
  perturbed_node_marginals.scale(1.0 / r);
  perturbed_edge_marginals.scale(1.0 / r);
  //cout<<endl<<"dl_dlnodepot"<<endl;
  //cout<<perturbed_node_marginals<<endl;
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
  //cout<<dL_dTheta<<endl;
  return dL_dTheta;
}

void sample::get_dL_dMu(int which_obj, arbi_array<num> node_marginals, arbi_array<num> edge_marginals, arbi_array<num>& dL_dNode_Mu, arbi_array<num>& dL_dEdge_Mu){
  switch(which_obj){
  case 1:
    get_dL_dMu_expected_distance(node_marginals, edge_marginals, dL_dNode_Mu, dL_dEdge_Mu);
    break;
  case 2:
    get_dL_dMu_nodewise(node_marginals, edge_marginals, dL_dNode_Mu, dL_dEdge_Mu);
  }
}


