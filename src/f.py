#from new_new_objects import *
from wrapper_decorator import dec
from wrapper import *
import new_new_objects as objects

class rW(obj_wrapper):

    # params will contain g, x where you return whether g.constructor((params)) is in the set x.  g is a wrapper
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        g = self.get_param(params, 'g')
        val = self.get_var_or_file(g, params, recalculate, False)
        return [float(val in self.get_param(params, 'x'))]

class sW(obj_wrapper):

    # params will contain g - the fxn, and 'values' - list of values
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        ans = []
        values = self.get_param(params, 'values')
        for x in values:
            self.set_param(params, 'x', x)
            ans = ans + self.get_var_or_file(rW, params, recalculate, False)
        return ans

class tW(obj_wrapper):

    # params will contain pdb_name, chain_letter, and pos, which is location in actual chain, not in my data structure
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        seq = self.get_var_or_file(objects.dW, params, recalculate, True)
        map = self.get_var_or_file(objects.eW, params, recalculate, True)
        aa = seq[map[self.get_param(params, 'pos')]]
        assert len(seq) == len(map.keys())
        aa = seq[map[params.get_param('pos')]]
        return aa

class uW(obj_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        self.set_param(params, 'g', tW)
        self.set_param(params, 'values', [ [x] for x in constants.AMINO_ACID_LIST])
        return self.get_var_or_file(sW, params, recalculate, False)

class vW(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        self.set_param(params, 'g', tW)
        self.set_param(params, 'values', constants.AMINO_ACID_CATEGORIES)
        return self.get_var_or_file(sW, params, recalculate, False)

class wW(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        inv_dists = self.get_var_or_file(objects.hW, params, recalculate, True)
        map = self.get_var_or_file(objects.eW, params, recalculate, True)
        return [inv_dists[map[self.get_param(params,'pos')]]]

class xW(obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        return [1]

class zW(obj_wrapper):

    # intrepid features
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        data = self.get_var_or_file(objects.yW, params, recalculate, True, True, always_recalculate)
        pos = self.get_param(params, 'pos')
        # position in data can't be higher than its position
        start_search_pos = min(len(data)-1, pos)
        i = start_search_pos
        idx_to_use = [4,5,8]

        asdf = 0
        while i >= -1:
            if i == -1:
                print "didn't find pos in intrepid data"
                assert False
            if data[i][1] == pos:
                break
            i = i-1
            asdf = asdf + 1
        line = data[i]
        return [line[j] for j in idx_to_use]

class aaW(obj_wrapper):

    # gets b-factor of residue
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False):
        pos_to_aa = self.get_var_or_file(objects.eW, params, recalculate, True, False, always_recalculate)
        chain = self.get_var_or_file(objects.cW, params, recalculate, True, False, always_recalculate)
        res = chain[pos_to_aa[self.get_param(params, 'pos')]]
        return [global_stuff.get_representative_atom(res).get_bfactor()]
    
