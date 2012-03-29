"""
created 3-26-2012
CRF will consist of potentials(node_id's, list of feature functions), values of feature/class associated weights, values of node_id associated variables.  potentials are really just a functions(the node_id's associated with function determine its C, and the values of node_id associated variables determine the values of y_C).  potentials can evaluate themselves.  also need a way to evaluate.  need to create the CRF by taking in a pdb chain
"""

import constants
import global_stuff
from objects import *
from cache import cache
import operator

class res_id(object):

    def __init__(self, pdb_name, chain_letter, pos):
        self.pdb_name = pdb_name
        self.chain_letter = chain_letter
        self.pos = pos

    def __hash__(self):
        return self.pdb_name.__hash__() + self.chain_letter.__hash__() + self.pos

    def __eq__(self, other):
        return (self.pdb_name == other.pdb_name) and (self.chain_letter == other.chain_letter) and (self.pos == other.pos)

class res_potential(object):

    def get_var_ranges(self):
        return [[0],[1]]

    # returns list of features you actually compute(aka don't depend on y values)
    def get_sub_feature_wrappers(self):
        return [conservation_feature_function_wrapper(), inverse_average_distances_feature_function_wrapper()]

    # returns list of feature_wrappers
    def get_feature_wrappers(self):
        to_return = []
        temp_params = param({})
        for val in self.get_var_ranges():
            temp_params.set_param('indicator_value', val)
            temp_params.set_param('wrapper1', indicator_function_wrapper(temp_params))
            temp_params.set_param('wrapper2', y_value_feature_function_wrapper(temp_params))
            temp_params.set_param('wrapper1', function_composition_wrapper(temp_params))
            for f in get_sub_feature_wrappers(params):
                temp_params.set_param('wrapper2', f)
                to_return.append(function_multiplier_wrapper(temp_params))
        return to_return

    # based on the information in params(equivalent of 'c') looks in global variables and retrieves feature values
    # for those params, including values of y.  for res_potential, params will be pdb_name, chain_letter, pos
    def get_feature_value_vector(self):
        vect = []
        for feature_wrapper in self.feature_wrappers:
            f = global_stuff.the_obj_manager.get_variable(feature_wrapper)
            vect.append(f(self.params))
        return vect

    # takes in weights for all kinds of nodes, retrieves the weight vector corresponding to itself
    # and returns potential value. needs ordered list of feature
    # self.info would contain the {feature} kind, the clique template
    def get_value(self, all_weights):
        weights = all_weights.get_param(self.info)
        v1 = self.get_feature_value_vector(self)
        return global_stuff.dot_product(v1, weights)
        

    # params should contain pdb_name, chain_letter, pos
    def __init__(self, params):
        self.params = params
        self.feature_wrappers = self.get_feature_wrappers(params)
    

class potential(object):

    def get_feature_wrappers(self):
        pass

    # assuming that var_info is sorted according to potential_info for now
    def __init__(self, potential_info, var_info, added_feature_names):

        self.potential_info = potential_info
        self.var_info = var_info
        feature_wrappers = self.get_feature_wrappers()
        self.feature_name_to_feature_function = {}
        for wrapper in feature_wrappers():
            self.feature_name_to_feature_function[wrapper.get_name()] = global_stuff.the_obj_manager.get_variable(wrapper)
        added_feature_names = feature_name_to_feature_function.keys()
        

class crf(object):

    def 

    def __init__(self, params, recalculate):
        self.potentials = []
        self.params = params
        
        res_dists = global_stuff.the_obj_manager.get_variable(pdb_chain_pairwise_distance_obj_wrapper(self.params, recalculate))
        aa_to_pos = global_stuff.the_obj_manager.get_variable(pdb_chain_aa_to_pos_obj_wrapper(self.params, recalculate))
        
