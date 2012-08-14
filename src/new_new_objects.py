#from new_features import *
import new_features as features
from wrapper_decorator import dec
#from wrapper import *
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



class pdb_file_w(wrapper.file_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        pdb_file_name = self.get_param(params, 'pdb_name')
        pdbl = Bio.PDB.PDBList()
        pdbl.retrieve_pdb_file(pdb_file_name, pdir=self.get_holding_folder())
        subprocess.call(['mv', self.get_holding_folder() + string.lower('pdb'+pdb_file_name+'.ent'), self.get_holding_location()])
        return open(self.get_holding_location(), 'r')

the_pdb_file_w = pdb_file_w()
        
class pdb_chain_w(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        f = self.get_var_or_file(the_pdb_file_w, params, recalculate, to_pickle)
        structure = Bio.PDB.PDBParser().get_structure(self.get_param(params, 'pdb_name'), f)
        chain = Bio.PDB.PPBuilder().build_peptides(structure[0][self.get_param(params, 'chain_letter')])
        to_return = []
        for chain_frag in chain:
            to_return = to_return + chain_frag
        # should raise exception here if error
        return to_return

the_pdb_chain_w = pdb_chain_w()

class pdb_chain_aa_to_pos_obj_w(wrapper.obj_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        chain_obj = self.get_var_or_file(the_pdb_chain_w, params, recalculate, True)
        chain_positions = [chain_obj[j].get_id()[1] for j in range(len(chain_obj))]
        return chain_positions

the_pdb_chain_aa_to_pos_obj_w = pdb_chain_aa_to_pos_obj_w()

class pdb_chain_seq_file_w(wrapper.file_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        chain_obj = self.get_var_or_file(the_pdb_chain_w, params, recalculate, True, False)
        # write the seq file at location + name
        raw_seq_string = ''.join([Polypeptide.three_to_one(res.resname) for res in chain_obj])
        seq = Bio.Seq.Seq(raw_seq_string)
        seq_record = Bio.SeqRecord.SeqRecord(seq)
        SeqIO.write(seq_record, self.get_holding_location(), 'fasta')
        return open(self.get_holding_location(),'r')

the_pdb_chain_seq_file_w = pdb_chain_seq_file_w()

class pdb_chain_seq_obj_w(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        seq_file_handle = self.get_var_or_file(the_pdb_chain_seq_file_w, params, recalculate, True, False)
        asdf = SeqIO.read(seq_file_handle, 'fasta')
        return asdf


the_pdb_chain_seq_obj_w = pdb_chain_seq_obj_w()

class pdb_chain_pos_to_aa_dict_obj_w(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        aa_to_pos = self.get_var_or_file(the_pdb_chain_aa_to_pos_obj_w, params, recalculate, True)
        pos_to_aa_dict = {}
        for i in range(len(aa_to_pos)):
            pos_to_aa_dict[aa_to_pos[i]] = i
        return pos_to_aa_dict

the_pdb_chain_pos_to_aa_dict_obj_w = pdb_chain_pos_to_aa_dict_obj_w()

class pdb_chain_pairwise_dist_obj_w(wrapper.mat_obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        residues = self.get_var_or_file(the_pdb_chain_w, params, recalculate, True)
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
        return dists

the_pdb_chain_pairwise_dist_obj_w = pdb_chain_pairwise_dist_obj_w()

class pdb_chain_inv_avg_dist_obj_w(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        dists = self.get_var_or_file(the_pdb_chain_pairwise_dist_obj_w, params, recalculate, True)
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

the_pdb_chain_inv_avg_dist_obj_w = pdb_chain_inv_avg_dist_obj_w()

class pdb_chain_edge_list_obj_w(wrapper.mat_obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        pdb.set_trace()
        res_dists = self.get_var_or_file(the_pdb_chain_pairwise_dist_obj_w, params, recalculate, True)
        aa_to_pos = self.get_var_or_file(the_pdb_chain_aa_to_pos_obj_w, params, recalculate, True)
        chain_seq = self.get_var_or_file(the_pdb_chain_seq_obj_w, params, recalculate, True)

        edges = []

        for i in range(len(chain_seq)):
            for j in range(i):
                if math.fabs(i-j) == 1:
                    edges.append([int(i),int(j)])
                elif res_dists[i][j] < self.get_param(params, 'dist_cut_off') and math.fabs(i-j) > 5:
                    edges.append([int(i),int(j)])

        return edges

the_pdb_chain_edge_list_obj_w = pdb_chain_edge_list_obj_w()

class node_features_obj_w(wrapper.mat_obj_wrapper):

    # params will contain node_feature_list
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = True):
        feature_list = self.get_param(params, 'node_feature_list')
        aa_to_pos = self.get_var_or_file(the_pdb_chain_aa_to_pos_obj_w, params, recalculate, True)
        node_features = []
        for i in range(len(aa_to_pos)):
            pos = aa_to_pos[i]
            self.set_param(params, 'pos', pos)
            temp = []
            for the_fxn_w in feature_list:
                temp = temp + self.get_var_or_file(the_fxn_w, params, recalculate, False)
            node_features.append(temp)
        return node_features

#pdb.set_trace()
the_node_features_obj_w = node_features_obj_w()

class edge_features_obj_w(wrapper.mat_obj_wrapper):

    # params will contain edge_feature_list
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        pdb.set_trace()
        feature_list = self.get_param(params, 'edge_feature_list')
        aa_to_pos = self.get_var_or_file(the_pdb_chain_aa_to_pos_obj_w, params, recalculate, True)
        edge_list = self.get_var_or_file(the_pdb_chain_edge_list_obj_w, params, recalculate, True)
        edge_features = []
        for i in range(len(edge_list)):
            edge = edge_list[i]
            pos1 = aa_to_pos[edge[0]]
            pos2 = aa_to_pos[edge[1]]
            self.set_param(params, 'pos1', pos1)
            self.set_param(params, 'pos2', pos2)
            temp = []
            for the_fxn_w in feature_list:
                temp = temp + self.get_var_or_file(the_fxn_w, params, recalculate, False)
            edge_features.append(temp)
        return edge_features

the_edge_features_obj_w = edge_features_obj_w()

class experiment_info_file_w(wrapper.file_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        # this simply extracts every single parameter, puts them into a dict, and prints the dict
        print params
#        pdb.set_trace()
        the_dict = {}
        the_dict["dist_cut_off"] = self.get_param(params, "dist_cut_off");
        # data_list will be list of tuples of (pdb_name, chain_letter)
        the_dict["data_list_file"] = self.get_param(params, "data_list_file")
        the_dict['node_feature_list'] = self.get_param(params, "node_feature_list")
        the_dict['edge_feature_list'] = self.get_param(params, "edge_feature_list")
        f = open(self.get_holding_location(), 'w')
        f.write(str(the_dict))
        f.close()
        return f

the_experiment_info_file_w = experiment_info_file_w()

class formatted_data_list_obj_w(wrapper.mat_obj_wrapper):

    # params should include file name want to read - 'data_list_file'
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        f = open(self.get_param(params, 'data_list_file'), 'r')
        ans = []
        for line in f:
            ans.append(string.split(string.strip(line), sep='_'))
       # pdb.set_trace()
        print ans
        return ans

the_formatted_data_list_obj_w = formatted_data_list_obj_w()
        
class experiment_results_obj_w(wrapper.mat_obj_wrapper):

    # params will be [(scores, sizes, pdb_name, chain_letter),   ]
    # params which which it is stored does NOT include these things...only data_list, params for getting features
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        # in perfect world, this would have been generated model file, which then generates results.  so simulate this by getting parameters from experiment info
        # this would have been the root of call tree
        experiment_info_pretending_to_be_model_params = self.get_var_or_file(the_experiment_info_file_w, params, recalculate, False, True)
        # node_features can get passed params, because params would have been given to this wrapper which would create model and use them to get node features
        # after getting scores, have wrappers depending on scores that calculate roc, other measures
        # i can't call them arbitrary order, so have to call scores, then other wrappers
        pdb.set_trace()
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
            mat.append([pdb_names[i], chain_letters[i]])
            mat.append(scores[pos:pos+sizes[i]])
        return mat

the_experiment_results_obj_w = experiment_results_obj_w();

class true_states_obj_w(wrapper.vect_obj_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        pdb.set_trace()
        map = self.get_var_or_file(the_pdb_chain_aa_to_pos_obj_w, params, recalculate, True, False)
        return [random.randrange(0,2) for x in range(len(map))]

the_true_states_obj_w = true_states_obj_w()

class roc_curve_input_obj_w(wrapper.mat_obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        results = self.get_var_or_file(the_experiment_results_obj_w, params, recalculate, to_pickle, False)
        assert len(results) % 2 == 0
        num_samples = len(results) / 2
        roc_classes = []
        roc_scores = []
        for i in range(num_samples):
            pdb_name = results[2 * i][0]
            chain_letter = results[ 2 * i][1]
            scores = results[(2 * i) + 1]
            self.set_param(params, "pdb_name", pdb_name)
            self.set_param(params, "chain_letter", chain_letter)
            true_states = self.get_var_or_file(the_true_states_obj_w, params, recalculate, True, True)
            roc_classes = roc_classes + true_states
            roc_scores = roc_scores + scores
            
        ans = global_stuff.transpose([roc_classes, roc_scores])
        return ans

the_roc_curve_input_obj_w = roc_curve_input_obj_w()

class roc_curve_plot_file_w(wrapper.file_wrapper):

    # params will be those required by experiment_results, which are those required by experiment_info.
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False):
        f  = self.get_var_or_file(the_roc_curve_input_obj_w.cache.file_dumper, params, recalculate, False, False)
        subprocess.call(['yard-plot', '-o', self.get_holding_location(), f.name])
        return f

def print_stuff(x):
    #pdb.set_trace()
    print 'printing!!!!!!!!'
    print x

the_roc_curve_plot_file_w = roc_curve_plot_file_w()


the_dict = {'data_list_file':'catres_two.pdb_list', 'edge_feature_list':[features.the_ones_fxn_w], 'node_feature_list':[features.the_residue_categorical_fxn_w, features.the_residue_class_categorical_fxn_w, features.the_inverse_avg_dist_fxn_w], 'dist_cut_off':5}
the_params = param.param(the_dict)
#pdb.set_trace()        
the_experiment_info_file_w.constructor(the_params, False, False, False)

# py_initialize
# from c++, get params as python object
# have wrapper 
# create model using python params.  store params in python. model will then read in params
# sample will have to set some parameters
