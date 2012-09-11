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

class categorical_fxn_w(obj_wrapper):

    # params will contain g - the fxn, and 'values' - list of values
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        ans = []
        values = self.get_param(params, 'values')
        for x in values:
            self.set_param(params, 'x', x)
            ans = ans + self.get_var_or_file(indicator_fxn_w, params, recalculate, False)
        return ans

class get_residue_fxn_w(obj_wrapper):

    # params will contain pdb_name, chain_letter, and pos, which is location in actual chain, not in my data structure
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        seq = self.get_var_or_file(objects.pdb_chain_seq_obj_w, params, recalculate, True)
        map = self.get_var_or_file(objects.pdb_chain_pos_to_aa_dict_obj_w, params, recalculate, True)
        aa = seq[map[self.get_param(params, 'pos')]]
        assert len(seq) == len(map.keys())
        aa = seq[map[params.get_param('pos')]]
        return aa

class residue_categorical_fxn_w(obj_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        self.set_param(params, 'g', get_residue_fxn_w)
        self.set_param(params, 'values', [ [x] for x in constants.AMINO_ACID_LIST])
        return self.get_var_or_file(categorical_fxn_w, params, recalculate, False)

class residue_class_categorical_fxn_w(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        self.set_param(params, 'g', get_residue_fxn_w)
        self.set_param(params, 'values', constants.AMINO_ACID_CATEGORIES)
        return self.get_var_or_file(categorical_fxn_w, params, recalculate, False)

class inverse_avg_dist_fxn_w(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        inv_dists = self.get_var_or_file(objects.pdb_chain_inv_avg_dist_obj_w, params, recalculate, True)
        map = self.get_var_or_file(objects.pdb_chain_pos_to_aa_dict_obj_w, params, recalculate, True)
        return [inv_dists[map[self.get_param(params,'pos')]]]

class ones_fxn_w(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        return [1]

