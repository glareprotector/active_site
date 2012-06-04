 #ifndef model_h
#define model_h

 #include "sample.h"
#include <exception>

//#include <string>

using namespace std;


class model{

 public:

  arbi_array<sample> data;
  int num_samples;
  arbi_array<int> training_indicies;
  arbi_array<int> testing_indicies;
  int num_training;
  int num_testing;
  int num_folds;
  int which_fold;

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
  void assign(int _num_folds, int _which_fold);
  void normalize();
  model(int _num_states, int _num_node_features, int _num_edge_features, arbi_array<string> _folder_names, int _mean_field_max_iter, int _num_folds, int _which_fold);
  arbi_array<num> get_gradient();
  num get_likelihood();

  int mean_field_max_iter;

};

#endif
