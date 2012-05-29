 #ifndef model_h
#define model_h

 #include "sample.h"
#include <exception>

//#include <string>

using namespace std;


class model{

 public:

  arbi_array<sample> training_data;
  arbi_array<sample> testing_data;
  int num_training;
  int num_testing;

  int num_states;
  int num_node_features;
  int num_edge_features;
  int theta_length;
  arbi_array<int> node_map;
  arbi_array<int> edge_map;

  arbi_array<num> gradient;
  num likelihood;
  arbi_array<num> theta;

  void set_theta(arbi_array<num> _theta);
  sample read_sample(string folder_name);
  void load_data(arbi_array<string> folder_names);
  model(int _num_states, int _num_node_features, int _num_edge_features, arbi_array<int> _node_map, arbi_array<int> _edge_map, arbi_array<string> folder_names);
  arbi_array<num> get_gradient();
  num get_likelihood();


};

#endif
