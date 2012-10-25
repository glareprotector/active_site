#from new_new_objects import *
from wrapper_decorator import dec
from wrapper import *
import new_new_objects as objects
import re
import string
from Bio.PDB import Polypeptide
import numpy
import global_stuff
import helper

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
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

#        pdb.set_trace()
        all_bs = self.get_var_or_file(objects.bfW, params, recalculate, True, False, False)

        
        pos_to_aa = self.get_var_or_file(objects.eW, params, recalculate, True, False, always_recalculate)
        return [all_bs[pos_to_aa[self.get_param(params, 'pos')]]]
#        return [global_stuff.get_representative_atom(res).get_bfactor()]
    
class akW(obj_wrapper):

    # computes KL div between 2 columns of 2 chains
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        
        msa = self.get_var_or_file(objects.agW, params, recalculate, True, True, False)
        pos_to_aa = self.get_var_or_file(objects.eW, params, recalculate, True, False, False)
        aa1 = pos_to_aa[self.get_param(params, 'pos1')]
        aa2 = pos_to_aa[self.get_param(params, 'pos2')]
        col1 = msa[:,aa1]
        col2 = msa[:,aa2]
        d1 = re.sub(r'-','',col1)
        d2 = re.sub(r'-','',col2)
        return [helper.get_KL(d1,d2)]

class axW(obj_wrapper):

    # g fxn for DSSP categorical wrapper
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        dssp_dict = self.get_var_or_file(objects.awW, params, recalculate, True, True, False)
        return dssp_dict[(self.get_param(params, 'pos'),self.get_param(params, c))]


class ayW(obj_wrapper):

    #categorical feature fxn for DSSP
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        self.set_param(params, 'g', axW)
        self.set_param(params, 'values', [ [x] for x in constants.DSSP_CATEGORIES])
        return self.get_var_or_file(sW, params, recalculate, False)
    

class bbW(obj_wrapper):

    # NACCESS feature
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        the_dict = self.get_var_or_file(objects.baW, params, recalculate, True, True, False)
        chain_letter = self.get_param(params, c)
        pos = self.get_param(params, 'pos')
        key = (pos, chain_letter)
        vals = the_dict[key]
        res_three = vals[0]
        chain_seq_in_one = self.get_var_or_file(objects.dW, params, recalculate, True, False)
        pos_to_aa = self.get_var_or_file(objects.eW, params, recalculate, True, False)
        res_one = chain_seq_in_one[pos_to_aa[pos]]
        assert(Polypeptide.three_to_one(res_three) == res_one)
#        pdb.set_trace()
        return [vals[1], vals[3], vals[5], vals[7], vals[9]]


class beW(obj_wrapper):

    # ligsite feature
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        sites = self.get_var_or_file(objects.bdW, params, recalculate, True, True, False)
        num_pockets_to_consider = self.get_param(params, 'lgn')
#        cutoff = self.get_param(params, 'lgc')
        coords = self.get_var_or_file(objects.bgW, params, recalculate, True, False, False)
        pos_to_aa = self.get_var_or_file(objects.eW, params, recalculate, True, False, False)
        coord = coords[pos_to_aa[self.get_param(params, 'pos')]]
        yes = 999
        for i in range(num_pockets_to_consider):
            pocket_coord = numpy.array(sites[i])
            dist = helper.physical_distance(coord - pocket_coord)

            if dist  < yes:
                yes = dist
                
        return [yes]
        

class bjW(obj_wrapper):

    # mutual information networks feature
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        answers = self.get_var_or_file(objects.biW, params, recalculate, True, True, False)
        pos_to_aa = self.get_var_or_file(objects.eW, params, recalculate, True, False, False)
        return answers[pos_to_aa[self.get_param(params, 'pos')]]
