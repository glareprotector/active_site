#include "model.h"
#include "LBFGS.h"
//#include "nums.h"

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
  arbi_array<int> true_states = read_vect_to_int(true_states_file, num_nodes, ' ');
  // the true_states file has 1 and 2's, but this program wants 0 and 1's
  for(int i = 0; i < num_nodes; i++){
    true_states(i)--;
  }
  arbi_array<int> edges = read_mat_to_int(edge_file, num_edges, 2);
  return sample(this, node_features, edge_features, edges, true_states, folder_name);
}

void model::load_data(arbi_array<string> folder_names){
  int num_samples = folder_names.size(0);
  cout<<"building training_idx"<<endl;

  for(int i = 0; i < num_samples; i++){
    try{
      sample s = read_sample(folder_names(i));
      // for now assigning everything to training, but have to change this
      /*if(i > num_samples / 2){
	this->training_data.append(s);
      }
      else{
	this->testing_data.append(s);
	}*/
      this->training_data.append(s);
    }
    catch(...){
      cout<<"default exception"<<endl;
    }
  }
  this->num_training = training_data.size(0);
  this->num_testing = testing_data.size(0);
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
  // string model_info_file("asdf/");
  // arbi_array<int> model_info = read_vect_to_int(model_info_file, 3);
  // this->num_states = model_info(0);
  // this->num_node_features = model_info(1);
  // this->num_edge_features = model_info(2);


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

class my_minimizer: public Minimizer{
 public:
  
  model* p_model;

  void ComputeGradient(vector<double>& gradient, const vector<double>& x){
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
  }

  double ComputeFunction(const vector<double>& x){
    arbi_array<num> theta(1, x.size());
    for(int i = 0; i < x.size(); i++){
      theta(i) = x[i];
    }
    p_model->set_theta(theta);
    return p_model->get_likelihood();
  }
};



//#include "sample.cpp"

int main(){

  //arbi_array<num> tqt = read_mat_to_num(string("test_mat.txt"),2,3);
  //cout<<tqt;
  //return 0;
  
  cout<<"hello world"<<endl;

  string pdb_list_file("/home/fultonw/active_site/active_site/data/catres_test.pdb_list");

  arbi_array<string> pdb_list = read_vect_to_string(pdb_list_file);
  cout<<pdb_list;

  string data_folder("/home/fultonw/active_site/active_site/test/");
  
  for(int i = 0; i < pdb_list.size(0); i++){
    pdb_list(i) = data_folder + pdb_list(i) + '/';
  }
  cout<<"pdb_list:"<<endl;
  cout<<pdb_list<<endl;

  int num_nodes = 2;
  int num_edges = 1;
  int num_states = 2;
  int num_node_features = 27;
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



  model m(num_states, num_node_features, num_edge_features, node_map, edge_map, pdb_list);


  int num_weights = idx;
  arbi_array<num> theta(1, num_weights);
  theta.fill(0);
  theta(0) = 0;
  theta(27) = 0;
  theta(1) = 7.3;
  theta(28) = 7;
  theta(2) = 1;
  theta(29) = 1.3;
  theta(3) = .08;
  theta(30) = .1;

  //theta.fill(1);
  //theta(0)=3;
  //theta(1)=2;
  m.set_theta(theta);

  for(int i = 0; i < m.num_training; i++){
    m.training_data(i).set_node_potentials();
    m.training_data(i).set_edge_potentials();
    m.training_data(i).set_marginals();
  }


  
  sample& s = m.training_data(0);
  cout<<"NODE FEATURRES: "<<endl;
  cout<<s.node_features<<endl;
  cout<<"NODE POTENTIALS:"<<endl;
  cout<<s.node_potentials<<endl;
  cout<<"edge potentials:"<<endl;
  cout<<s.edge_potentials<<endl;
  cout<<"node_marginals:"<<endl;
  cout<<s.node_marginals<<endl;
  cout<<"likelihood"<<endl;
  cout<<s.get_likelihood()<<endl;
  cout<<"gradient"<<endl;
  //arbi_array<num> grad = s.get_gradient();
  //cout<<grad<<endl;
  //cout<<s.get_likelihood();
  arbi_array<num> grad = s.get_gradient();
  operator<<(cout,grad);
  cout<<"log Z"<<endl;
  cout<<s.get_log_Z()<<endl;
  //  (cout<<s.get_gradient());
  /*


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
  */

  //cout<<endl<<"likelihood: ";
  //cout<<m.get_likelihood()<<endl;
  //cout<<"marginals: "<<endl;
  //cout<<s.node_marginals<<endl;
  cout<<"FOLDER: "<<s.folder<<endl;
  cout<<"ALL FOLDERS:"<<endl;
  for(int i = 0; i < m.num_training; i++){
    cout<<m.training_data(i).folder<<endl;
  }
}
