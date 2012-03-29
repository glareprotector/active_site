"""
created 3-27-2012
contains the function wrapper classes(made them wrappers so that i could retrieve their name without any new code) as well as the functions needed to implement those functions
"""

from objects import *

class indicator_function_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return ['indicator_value']

    def constructor(self, recalculate):
        
        def f(x):
            if x == self.get_param('indicator_value'):
                return 1
            else:
                return 0

        return f

class conservation_feature_function_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return []

    def constructor(self):

        # params is params object that will specify pdb_name, chain_letter, pos
        def f(params):

            scores = global_stuff.the_obj_manager.get_variable(pdb_chain_conservation_score_obj_wrapper(params))
            map = global_stuff.the_obj_manager.get_variable(pdb_chain_pos_to_aa_dict_obj_wrapper(params))
            return scores[map[params.get_param('pos')]]

        return f

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
