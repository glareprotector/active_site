#from fimport *
import f as features
from wrapper_decorator import dec
#from wrapper import *
import pdb
#pdb.set_trace()
import wrapper
import param



import Bio.PDB
import constants
import global_stuff

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.PDB import Polypeptide
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

import math
import subprocess
import string
import os
import random
import pdb

# doesn't matter if pdb_name and chain_letter are capitalized
class fW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        pdb_file_name = self.get_param(params, 'pdb_name')
        pdbl = Bio.PDB.PDBList()
        pdbl.retrieve_pdb_file(pdb_file_name, pdir=self.get_holding_folder())
        subprocess.call(['mv', self.get_holding_folder() + string.lower('pdb'+pdb_file_name+'.ent'), self.get_holding_location()])

        return open(self.get_holding_location(), 'r')

        
class cW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        f = self.get_var_or_file(fW, params, recalculate, to_pickle)
        structure = Bio.PDB.PDBParser().get_structure(self.get_param(params, 'pdb_name'), f)
        chain = Bio.PDB.PPBuilder().build_peptides(structure[0][self.get_param(params, 'chain_letter')])
        to_return = []
        for chain_frag in chain:
            to_return = to_return + chain_frag
        # should raise exception here if error
#        pdb.set_trace()
        return to_return

class aW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        chain_obj = self.get_var_or_file(cW, params, recalculate, True)
        chain_positions = [chain_obj[j].get_id()[1] for j in range(len(chain_obj))]
        return chain_positions

class bW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        chain_obj = self.get_var_or_file(cW, params, recalculate, True, False)
        # write the seq file at location + name
        raw_seq_string = ''.join([Polypeptide.three_to_one(res.resname) for res in chain_obj])
        seq = Bio.Seq.Seq(raw_seq_string)
        seq_record = Bio.SeqRecord.SeqRecord(seq)
        SeqIO.write(seq_record, self.get_holding_location(), 'fasta')
        return open(self.get_holding_location(),'r')

class dW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        seq_file_handle = self.get_var_or_file(bW, params, recalculate, True, False)
        asdf = SeqIO.read(seq_file_handle, 'fasta')
        return asdf

class eW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        aa_to_pos = self.get_var_or_file(aW, params, recalculate, True)

        pos_to_aa_dict = {}
        for i in range(len(aa_to_pos)):
            pos_to_aa_dict[aa_to_pos[i]] = i
#        pdb.set_trace()
        return pos_to_aa_dict

class gW(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        residues = self.get_var_or_file(cW, params, recalculate, True)
        rep_atoms = [global_stuff.get_representative_atom(res) for res in residues]
        num_res = len(residues)
        dists = [[-1 for i in range(num_res)] for j in range(num_res)]
        for i in range(num_res):
            for j in range(num_res):
                try:
                    dists[i][j] = rep_atoms[i] - rep_atoms[j]
                except Exception as e:
                    print 'ERROR: distance fail', self.params, i, j, residues[i].child_dict.keys(), residues[j].child_dict.keys()
                    dists[i][j] = -1
#        pdb.set_trace()
        return dists

class hW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        dists = self.get_var_or_file(gW, params, recalculate, True)
        inv_avg_dists = [-1 for i in range(len(dists))]
        for i in range(len(dists)):
            val = 0;
            count = 0
            for j in range(len(dists)):
                if dists[i][j] != -1:
                    val = val + dists[i][j]
                    count = count + 1
            inv_avg_dists[i] = 1.0 / (val / count)
        return inv_avg_dists

class iW(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        res_dists = self.get_var_or_file(gW, params, recalculate, True)
        aa_to_pos = self.get_var_or_file(aW, params, recalculate, True)
        chain_seq = self.get_var_or_file(dW, params, recalculate, True)

        edges = []

        for i in range(len(chain_seq)):
            for j in range(i):
                if math.fabs(i-j) == 1:
                    edges.append([int(i),int(j)])
                elif res_dists[i][j] < self.get_param(params, 'dist_cut_off') and math.fabs(i-j) > 5:
                    edges.append([int(i),int(j)])

        return edges

class jW(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    # params will contain node_feature_list
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = True, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        feature_list = self.get_param(params, 'node_feature_list')
        aa_to_pos = self.get_var_or_file(aW, params, recalculate, True)
        node_features = []
        for i in range(len(aa_to_pos)):
            pos = aa_to_pos[i]
            self.set_param(params, 'pos', pos)
            temp = []
            for fxn_w in feature_list:
                temp = temp + self.get_var_or_file(fxn_w, params, recalculate, False)
            node_features.append(temp)
        #pdb.set_trace()
        return node_features

class kW(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    # params will contain edge_feature_list
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        feature_list = self.get_param(params, 'edge_feature_list')
        aa_to_pos = self.get_var_or_file(aW, params, recalculate, True)
        edge_list = self.get_var_or_file(iW, params, recalculate, True)
        edge_features = []
        for i in range(len(edge_list)):
            edge = edge_list[i]
            pos1 = aa_to_pos[edge[0]]
            pos2 = aa_to_pos[edge[1]]
            self.set_param(params, 'pos1', pos1)
            self.set_param(params, 'pos2', pos2)
            temp = []
            for fxn_w in feature_list:
                temp = temp + self.get_var_or_file(fxn_w, params, recalculate, False)
            edge_features.append(temp)
        return edge_features

class lW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        # this simply extracts every single parameter, puts them into a dict, and prints the dict.  it also sets params if they weren't gotten here
        print params
        #pdb.set_trace()
        the_dict = {}
        the_dict["dist_cut_off"] = self.get_param(params, "dist_cut_off");
        # data_list will be list of tuples of (pdb_name, chain_letter)
        the_dict["data_list_file"] = self.get_param(params, "data_list_file")
        the_dict['node_feature_list'] = self.get_param(params, "node_feature_list")
        the_dict['edge_feature_list'] = self.get_param(params, "edge_feature_list")

        return the_dict

class mW(wrapper.mat_obj_wrapper):

    # params should include file name want to read - 'data_list_file'
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(self.get_param(params, 'data_list_file'), 'r')
        ans = []
        for line in f:
            ans.append(string.split(string.strip(line), sep='_'))
       # pdb.set_trace()
        print ans
        return ans

        
class nW(wrapper.mat_obj_wrapper):

    # params will be [(scores, sizes, pdb_name, chain_letter),   ]
    # params which which it is stored does NOT include these things...only data_list, params for getting features
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        # in perfect world, this would have been generated model file, which then generates results.  so simulate this by getting parameters from experiment info
        # this would have been the root of call tree
        experiment_info_pretending_to_be_model_params = self.get_var_or_file(lW, params, recalculate, False, True)

        # node_features can get passed params, because params would have been given to this wrapper which would create model and use them to get node features
        # after getting scores, have wrappers depending on scores that calculate roc, other measures
        # i can't call them arbitrary order, so have to call scores, then other wrappers
#        pdb.set_trace()
        scores = self.get_param(params, 'scores', False)
        true_states = self.get_param(params, 'true_states', False)
        sizes = self.get_param(params, 'sizes', False)
        pdb_names = self.get_param(params, 'pdb_names', False)
        chain_letters = self.get_param(params, 'chain_letters', False)
        num_samples = len(pdb_names)
        pos = 0
        mat = []
        # alternate between pdb_names, chain_letter and scores
        for i in range(num_samples):
            mat.append([pdb_names[i], chain_letters[i], sizes[i]])
            mat.append(scores[pos:pos+sizes[i]])
            mat.append(true_states[pos:pos+sizes[i]])
            pos = pos + sizes[i]
        return mat

# for features i can't get automatically, assume that there is folder for each chain containing the features.  wrappers will then read from those folders
class oW(wrapper.vect_obj_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        pos_to_aa = self.get_var_or_file(eW, params, recalculate, True, False)
        #pdb.set_trace()
        pdb_name = self.get_param(params, 'pdb_name')
        chain_letter = self.get_param(params, 'chain_letter')
        aux_folder = global_stuff.get_aux_folder(pdb_name, chain_letter)
        states_file = aux_folder + string.lower(pdb_name) + '.catres'
        f = open(states_file, 'r')
        true_states = [0 for i in range(len(pos_to_aa.keys()))]
        chain_seq_in_one = self.get_var_or_file(dW, params, recalculate, True, False)
        #pdb.set_trace()
        for line in f:
            processed = string.split(string.split(line, sep='\t')[0], sep=' ')
            pos = int(processed[2])
            aa = pos_to_aa[pos]
            true_states[pos_to_aa[pos]] = 1
            three = processed[1]
            one = Polypeptide.three_to_one(three)
            assert one == chain_seq_in_one[aa]
            true_states[aa] = 1
        return true_states


class pW(wrapper.mat_obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        results = self.get_var_or_file(nW, params, recalculate, to_pickle, False, always_recalculate)

        assert len(results) % 3 == 0
        num_samples = len(results) / 3
        roc_classes = []
        roc_scores = []
        for i in range(num_samples):
            pdb_name = results[3 * i][0]
            chain_letter = results[ 3 * i][1]
            scores = results[(3 * i) + 1]
            true_states_from_cpp = results[(3 * i) + 2]
            self.set_param(params, "pdb_name", pdb_name)
            self.set_param(params, "chain_letter", chain_letter)
            true_states = self.get_var_or_file(oW, params, recalculate, True, True)
            roc_classes = roc_classes + true_states
            roc_scores = roc_scores + scores
            assert(len(roc_classes) == len(roc_scores))
            
        ans = global_stuff.get_transpose([roc_classes, roc_scores])
        return ans


class qW(wrapper.file_wrapper):

    # this part does not work
    # params will be those required by experiment_results, which are those required by experiment_info.
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = True, old_obj = None):
#        pdb.set_trace()
        # first get instance of pW
        from wc import wc
        self.set_param(params, 'which_wrapper_class', pW)
        pW_instance = self.old_get_var_or_file(wc, params, True, False, False)
        f  = self.old_get_var_or_file(pW_instance.cache.file_dumper_wrapper, params, recalculate, False, False, always_recalculate)
        # call R roc curve script, specifying input location and where to write curve
        source = f.name
        destination = self.get_holding_location()
        old_destination = self.get_backup_location()
#        pdb.set_trace()
        assert os.path.isfile(source)
        assert os.path.isfile(constants.ROC_CURVE_SCRIPT)
        iteration = self.get_param(params, 'iter', record=False)
        obj_val = self.get_param(params, 'obj_val', record=False)
        subprocess.call(['Rscript', constants.ROC_CURVE_SCRIPT, source, destination, str(iteration), str(obj_val)])
        # merge current version with old version.
        if old_obj != None:
            subprocess.call(['convert', '-append', old_obj.name, destination, destination])
        return open(destination)

    def get_file_location(self, object_key):
        return self.get_folder(object_key) + self.get_name(object_key) + '.jpg'

class yW(wrapper.obj_wrapper):

    # reads in intrepid file, stores as processed mat.  assumes that the required files are already in the correct folders.
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        pos_to_aa = self.get_var_or_file(eW, params, recalculate, True, False, True)
        pdb_name = self.get_param(params, 'pdb_name')
        chain_letter = self.get_param(params, 'chain_letter')
        aux_folder = global_stuff.get_aux_folder(pdb_name, chain_letter)
        intrepid_file = aux_folder + 'intrepid.aux'
        f = open(intrepid_file, 'r')
        all_lines = f.readlines()
        all_lines.pop(0)
        int_idx = [0,1]
        float_idx = [4,5,6,7,8,9]
        ans = []
        for line in all_lines:
            a_line = string.strip(line)
            a_line = string.split(a_line, sep='|')
            a_line.pop(4)
            for i in int_idx:
                a_line[i] = int(a_line[i])
            for i in float_idx:
                a_line[i] = float(a_line[i])
            ans.append(a_line)
        return ans
    
            
class abW(wrapper.obj_wrapper):

    # this wrapper doesn't output anything.  sole purpose is to call node features for stuff in data_list and store the results
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        data_list = self.get_var_or_file(mW, params, recalculate, False, False)
        which_wrapper = self.get_param(params, 'which_wrapperq')
        for i in range(len(data_list)):
            try:
                pdb_name = data_list[i][0]
                chain_letter = data_list[i][1]
                self.set_param(params, 'pdb_name', pdb_name)
                self.get_var_or_file(which_wrapper, params, recalculate, False, False)
            except:
                print 'error while downloading ', data_list[i][0]
            
    

class rW(wrapper.obj_wrapper):

    # takes in proc_id, num_procs, data_list, and returns which ones it tries calculate features for a node depending on its proc_id
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = True, old_obj = None):
        data_list = self.get_var_or_file(mW, params, True, False, False)
    

def print_stuff(x):
    #pdb.set_trace()
    print 'printing!!!!!!!!'
    print x


