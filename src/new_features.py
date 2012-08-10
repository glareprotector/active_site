#from new_new_objects import *
from wrapper_decorator import dec
from wrapper import *
import new_new_objects as objects

class indicator_fxn_w(obj_wrapper):

    # params will contain g, x where you return whether g.constructor((params)) is in the set x.  g is a wrapper
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        g = self.get_param(params, 'g')
        val = self.get_var_or_file(g, params, recalculate, False)
        return [float(val in self.get_param(params, 'x'))]

the_indicator_fxn_w = indicator_fxn_w()

class categorical_fxn_w(obj_wrapper):

    # params will contain g - the fxn, and 'values' - list of values
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        ans = []
        values = self.get_param(params, 'values')
        for x in values:
            self.set_param(params, 'x', x)
            ans = ans + self.get_var_or_file(the_indicator_fxn_w, params, recalculate, False)
        return ans

the_categorical_fxn_w = categorical_fxn_w()

class get_residue_fxn_w(obj_wrapper):

    # params will contain pdb_name, chain_letter, and pos, which is location in actual chain, not in my data structure
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        seq = self.get_var_or_file(objects.the_pdb_chain_seq_obj_w, params, recalculate, True)
        map = self.get_var_or_file(objects.the_pdb_chain_pos_to_aa_dict_obj_w, params, recalculate, True)
        aa = seq[map[self.get_param(params, 'pos')]]
        assert len(seq) == len(map.keys())
        aa = seq[map[params.get_param('pos')]]
        return aa

the_get_residue_fxn_w = get_residue_fxn_w()

class residue_categorical_fxn_w(obj_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        self.set_param(params, 'g', the_get_residue_fxn_w)
        self.set_param(params, 'values', [ [x] for x in constants.AMINO_ACID_LIST])
        return self.get_var_or_file(the_categorical_fxn_w, params, recalculate, False)

the_residue_categorical_fxn_w = residue_categorical_fxn_w()

class residue_class_categorical_fxn_w(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        self.set_param(params, 'g', the_get_residue_fxn_w)
        self.set_param(params, 'values', constants.AMINO_ACID_CATEGORIES)
        return self.get_var_or_file(the_categorical_fxn_w, params, recalculate, False)

the_residue_class_categorical_fxn_w = residue_class_categorical_fxn_w()

class inverse_avg_dist_fxn_w(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        inv_dists = self.get_var_or_file(objects.the_pdb_chain_inv_avg_dist_obj_w, params, recalculate, True)
        map = self.get_var_or_file(objects.the_pdb_chain_pos_to_aa_dict_obj_w, params, recalculate, True)
        return [inv_dists[map[self.get_param(params,'pos')]]]

the_inverse_avg_dist_fxn_w = inverse_avg_dist_fxn_w()

class ones_fxn_w(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        return 1

the_ones_fxn_w = ones_fxn_w()


dict = {'pdb_name':'1asy', 'chain_letter':'A', 'edge_feature_list':[the_ones_fxn_w], 'node_feature_list':[the_residue_class_categorical_fxn_w, the_residue_categorical_fxn_w, the_inverse_avg_dist_fxn_w], 'dist_cut_off':5}

dict = {'pdb_name':'1asy', 'chain_letter':'A', 'edge_feature_list':[the_ones_fxn_w], 'node_feature_list':[the_residue_class_categorical_fxn_w, the_inverse_avg_dist_fxn_w], 'dist_cut_off':5}

params = param(dict)
#print 'about to call node_features constructor'
used_keys, all_keys, node_features  = objects.the_pdb_chain_pairwise_dist_obj_w.constructor(params, True, True, True)
