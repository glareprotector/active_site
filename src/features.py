"""
created 3-27-2012
contains the function wrapper classes(made them wrappers so that i could retrieve their name without any new code) as well as the functions needed to implement those functions
"""

from objects import *

# returns function that returns 1 if the argument is in indicator_values
# params should contain g_wrapper - a function wrapper for how to get the value, and indicator_values = the values to match
class indicator_function_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return ['g_wrapper', 'indicator_values']

    def constructor(self, recalculate):
        #pdb.set_trace()
        g = global_stuff.the_obj_manager.get_variable(self.get_param('g_wrapper'))
        
        # params is whatever g needs to fetch the value
        def f(params):
            return int(g(params) in self.get_param('indicator_values'))

        return f

# returns the list of indicator function wrappers associated with a categorical variable
class categorical_function_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return ['g_wrapper', 'indicator_values_list']

    def constructor(self, recalculate):
        to_return = []
        for iv in self.get_param('indicator_values_list'):
            self.params.set_param('indicator_values', iv)
            to_return.append(indicator_function_wrapper(self.params))
        return to_return

class get_residue_function_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return []

    def constructor(self, recalculate):
        
        def f(params):
            seq = global_stuff.the_obj_manager.get_variable(pdb_chain_seq_obj_wrapper(params), recalculate)
            map = global_stuff.the_obj_manager.get_variable(pdb_chain_pos_to_aa_dict_obj_wrapper(params), recalculate)
            aa = seq[map[params.get_param('pos')]]
            return aa
        
        return f

class inverse_average_distance_feature_function_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return []

    def constructor(self, recalculate):
        
        # params should have pdb name, chain letter, pos
        def f(params):
            inv_dists = global_stuff.the_obj_manager.get_variable(pdb_chain_inverse_average_distances_obj_wrapper(params), recalculate)
            map = global_stuff.the_obj_manager.get_variable(pdb_chain_pos_to_aa_dict_obj_wrapper(params), recalculate)
            return inv_dists[map[params.get_param('pos')]]
        
        return f
        
class conservation_feature_function_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return []

    def constructor(self, recalculate):

        # params is params object that will specify pdb_name, chain_letter, pos
        def f(params):

            scores = global_stuff.the_obj_manager.get_variable(pdb_chain_conservation_score_obj_wrapper(params), recalculate)
            map = global_stuff.the_obj_manager.get_variable(pdb_chain_pos_to_aa_dict_obj_wrapper(params), recalculate)
            return scores[map[params.get_param('pos')]]

        return f

class y_value_feature_function_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return []

    def constructor(self, recalculate):
        
        # params will specify pdb_name, chain_letter, pos
        def f(params):
            return the_global_stuff.get_y_value(params)

class function_composition_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return ['wrapper1', 'wrapper2']

    def constructor(self):

        # g always takes in a real(which will be called x
        # h always returns a real
        def f(params):
            g = global_stuff.the_obj_manager.get_variable(self.get_param('wrapper1'))
            h = global_stuff.the_obj_manager.get_variable(self.get_param('wrapper2'))
            params.set_param('x', h(params))
            return g(params)


class function_multiplier_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return ['wrapper1', 'wrapper2']

    def constructor(self):

        # params is params object that will specify information needed to compute both functions
        def f(params):
            g = global_stuff.the_obj_manager.get_variable(self.get_param('wrapper1'))
            h = global_stuff.the_obj_manager.get_variable(self.get_param('wrapper2'))
            return g(params) * h(params)
