#from fimport *
import f as features
from wrapper_decorator import dec
#from wrapper import *
import pdb
#pdb.set_trace()
import wrapper
import param
import numpy




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
import re

# doesn't matter if pdb_name and chain_letter are capitalized
class fW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        pdb.set_trace()
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
        #raw_seq_string = ''.join([Polypeptide.three_to_one(res.resname) for res in chain_obj])
        raw_seq_list = []
        for res in chain_obj:
            try:
                raw_seq_list.append(Polypeptide.three_to_one(res.resname))
            except:
                raw_seq_list.append('0')
        raw_seq_string = ''.join(raw_seq_list)
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

        #temp = self.get_var_or_file(adW, params, recalculate, False, False, False)
        

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
#        pdb.set_trace()
        return edge_features

class lW(wrapper.obj_wrapper, wrapper.shorten_name_wrapper):

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
        the_dict['wif'] = self.get_param(params, "wif")
        the_dict['wob'] = self.get_param(params, "wob")
        the_dict['reg'] = self.get_param(params, "reg")
        the_dict['wreg'] = self.get_param(params, "wreg")
        the_dict['wob2'] = self.get_param(params, "wob2")
        the_dict['nwc'] = self.get_param(params, "nwc")
        the_dict['wtpr'] = self.get_param(params, "wtpr")
        the_dict['posw'] = self.get_param(params, "posw")
        the_dict['sfc'] = self.get_param(params, "sfc")
        
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

        
class nW(wrapper.mat_obj_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

    # params will be [(scores, sizes, pdb_name, chain_letter),   ]
    # params which which it is stored does NOT include these things...only data_list, params for getting features
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        # in perfect world, this would have been generated model file, which then generates results.  so simulate this by getting parameters from experiment info
        # this would have been the root of call tree
        #pdb.set_trace()
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
class oW(wrapper.vect_obj_wrapper, wrapper.by_pdb_folder_wrapper):
    
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
            #assert one == chain_seq_in_one[aa]
            true_states[aa] = 1
        return true_states


class pW(wrapper.mat_obj_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

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
            true_states = self.get_var_or_file(oW, params, global_stuff.recalculate, True, True)
            roc_classes = roc_classes + true_states
            roc_scores = roc_scores + scores
            assert(len(roc_classes) == len(roc_scores))

        print [(roc_classes[i],roc_scores[i]) for i in range(len(roc_scores)) if roc_classes[i] == 1 ]
            
        ans = global_stuff.get_transpose([roc_classes, roc_scores])
        return ans


class qW(wrapper.file_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

    # this part does not work
    # params will be those required by experiment_results, which are those required by experiment_info.
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = True, old_obj = None):
#        pdb.set_trace()
        # first get instance of pW
        from wc import wc
        # figure out which roc file input wrapper to use
        which_classifier_roc_input = self.get_param(params, 'wclf')
        # if using svm, treat 'iter' as the fold number just so i don't have to change roc curve code yet
        #import try_svm
        #if which_classifier_roc_input == try_svm.atW:
        #    self.set_param(params, 'iter', self.get_param(params, 'wfld'))
        self.set_param(params, 'which_wrapper_class', which_classifier_roc_input)
        pW_instance = self.old_get_var_or_file(wc, params, True, False, False)
        if always_recalculate:
            f  = self.old_get_var_or_file(pW_instance.cache.file_dumper_wrapper, params, recalculate, False, False, 2)
        else:
            f  = self.old_get_var_or_file(pW_instance.cache.file_dumper_wrapper, params, recalculate, False, False, False)
      #  pdb.set_trace()
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

class yW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

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
                print pdb_name
                self.set_param(params, 'pdb_name', pdb_name)
                self.get_var_or_file(which_wrapper, params, recalculate, False, False)
            except:
                print 'error while downloading ', data_list[i][0]
            
    
# blast results file wrapper(xml format)
class adW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        seq_records = []
        f = self.get_var_or_file(bW, params, recalculate, False, False, False)
        query = SeqIO.parse(f, 'fasta')
        seq_records.append(query)
        print 'RUNNING BLAST!!!!!!!'
        psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLAST_PATH, outfmt = 5, query = '\''+f.name+'\'', db = global_stuff.BLASTDB_PATH, out = self.get_holding_location())
        #pdb.set_trace()
        subprocess.call(str(psi_blast_cline), shell=True, executable='/bin/bash')
        return open(self.get_holding_location())


# processsed blast results format(bunch of sequences in fasta in 1 file)
class aeW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        # parse blast xml file, then do processing
        #pdb.set_trace()
        blast_xml_handle = self.get_var_or_file(adW, params, recalculate, False, False, False)
        record = NCBIXML.read(blast_xml_handle)
        seen = set()
        seq_records = []
        # add the query sequence, and have a set so that only add each sequence once
        query = self.get_var_or_file(dW, params, recalculate, True, False, False)
        query.id = 'QUERY'
        seq_records.append(query)
        seen.add(query.seq.tostring())
        # add high scoring pairs in alignments with sufficiently low evalue that have not been seen
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < self.get_param(params, 'evalue') and not hsp.sbjct in seen:
                    seq_records.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(hsp.sbjct), id = alignment.hit_id))
        # write hits to fasta file
        output_handle = open(self.get_holding_location(), 'w')
        SeqIO.write(seq_records, output_handle, 'fasta')
        print 'WROTE ', self.get_holding_location()
        return output_handle


# gets the result of muscle
class afW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        msa_input_handle = self.get_var_or_file(aeW, params, recalculate, to_pickle, to_filelize, always_recalculate)
        cline = MuscleCommandline(cmd = global_stuff.MUSCLE_PATH, input = '\''+msa_input_handle.name+'\'', out = self.get_holding_location(), clw = False, maxiters = 2)
        #cline()

        subprocess.call(str(cline), shell=True, executable='/bin/bash')
        return open(self.get_holding_location())


# processed msa output(columns with skips removed)
class agW(wrapper.msa_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        f = self.get_var_or_file(afW, params, recalculate, False, False, False)
        msa = AlignIO.read(f.name, 'fasta')
        # search for the query sequence
        idx = -1
        for i in range(len(msa)):
            if msa[i].id == 'QUERY':
                idx = i
                break
        # find the first non-insertion column
        i = 0

        while msa[idx,i] == '-':
            #print msa[idx,i]
            i = i + 1
            #pdb.set_trace()
            
        to_return = msa[:,i:(i+1)]
        # add in all the other columns
        for k in range(i+1, msa.get_alignment_length()):
            if msa[idx,k] != '-':
                #print k
                to_return = to_return + msa[:,k:(k+1)]
        return to_return

# if a chain has k true sites, it sees if any of the k most highly ranked sites are within a certain cut-off
class ahW(wrapper.mat_obj_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = True, old_obj = None):
#        pdb.set_trace()
        results = self.get_var_or_file(nW, params, recalculate, True, True, always_recalculate)
        num_samples = len(results) / 3
        cutoffs = global_stuff.metric_cutoffs
        these_results = [0 for i in range(len(cutoffs))]
        total_num_sites = 0
#        pdb.set_trace()
        for i in range(num_samples):
            pdb_name = results[3*i][0]
 #           pdb.set_trace()
            chain_letter = results[3*i][1]
            scores = results[(3*i)+1]
            true_states = results[(3*i)+2]
            # sort scores by marginals
            temp = [ [scores[i], i] for i in range(len(scores))]
#            pdb.set_trace()

            sorted_temp = sorted(temp, key = lambda x: x[0], reverse = True)
            self.set_param(params, 'pdb_name', pdb_name)
            self.set_param(params, 'chain_letter', chain_letter)
            # for each best site, see if any of the true sites are within that cutoff
            true_sites = self.get_var_or_file(anW, params, False, True, True, False)
            some_dists = self.get_var_or_file(amW, params, False, True, True, False)
            total_num_sites = total_num_sites + len(true_sites)

#           pdb.set_trace()
            num_true_sites = len(true_sites)
            for j in range(num_true_sites):
                true_site = true_sites[j]
                for l in range(len(cutoffs)):
                    for k in range(num_true_sites):
                        ranked_site = sorted_temp[k][1]
                        if some_dists[j][ranked_site] < cutoffs[l]:
                            these_results[l] = these_results[l] + 1
                            break
 #           pdb.set_trace()
        these_results = [float(these_results[i]) / float(total_num_sites) for i in range(len(cutoffs))]
        #pdb.set_trace()
        to_write = [self.get_param(params, 'iter', record=False), total_num_sites] + these_results
#        pdb.set_trace()
        if old_obj == None:
            return [to_write]
        else:
 #           pdb.set_trace()
            old_obj.append(to_write)
            return old_obj
            
                        
class aiW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        # assume params contains 2 distributions: d1 and d2
        d1 = self.get_param(params, 'd1')
        d2 = self.get_param(params, 'd2')
        # keep dictionary of counts for each distribution
        d1_dict = {}
        for i in range(len(d1)):
            if d1[i] in d1_dict.keys():
                d1_dict[i] = d1_dict[i] + 1
            else:
                d1_dict[i] = 1
        for k in d1_dict.keys():
            d1_dict[k] = float(d1_dict[k]) / float(len(d1))

        d2_dict = {}
        for i in range(len(d2)):
            if d2[i] in d2_dict.keys():
                d2_dict[i] = d2_dict[i] + 1
            else:
                d2_dict[i] = 1
        for k in d2_dict.keys():
            d2_dict[k] = float(d2_dict[k]) / float(len(d2))

        ans = 0
        for k in d1_dict.keys():
            if d1_dict[k] > 0:
                ans += d1_dict[k] * math.log(d1_dict[k] / d2_dict[k])
                
        return ans
        
# returns sorted dists, sorted pos for all position in a chain
class aoW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        dists = self.get_var_or_file(gW, params, recalculate, True, True, False)

        all_sorted_dists = []
        all_sorted_pos = []

        for j in range(len(dists)):
        
            dist = dists[j]
            temp = [ [dist[i],i] for i in range(len(dist)) ]
            sorted_temp = sorted(temp, key = lambda x: x[0])
            all_sorted_dists.append([sorted_temp[i][0] for i in range(len(sorted_temp))])
            all_sorted_pos.append([sorted_temp[i][1] for i in range(len(sorted_temp))])

        return [all_sorted_dists, all_sorted_pos]

class alW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        all_sorted_dists, all_sorted_pos = self.get_var_or_file(aoW, params, recalculate, True, True, False)
        aa = self.get_param(params, 'aa')
        return [all_sorted_dists[aa], all_sorted_pos[aa]]

# returns list of the true states in a chain
class anW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        true_states = self.get_var_or_file(oW, params, recalculate, True, True, False)
        true_list = []
        for i in range(len(true_states)):
            if true_states[i] == 1:
                true_list.append(i)
        return true_list

# returns distances rows of only the true states
class amW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        all_dists = self.get_var_or_file(gW, params, recalculate, True, True, False)
        true_list = self.get_var_or_file(anW, params, recalculate, True, True, False)
        ans = []
        for pos in true_list:
            ans.append(all_dists[pos])
        return ans

# returns sorted distances for every column, but truncated
class apW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = False, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        all_sorted_dists, all_sorted_pos = self.get_var_or_file(aoW, params, recalculate, True, False, False)
        dists = []
        pos = []
        trun = self.get_param(params, "trun", record = False)
        trun = min(trun, len(all_sorted_dists))
        for i in range(len(all_sorted_dists)):
            dists.append(all_sorted_dists[i][0:trun])
            pos.append(all_sorted_pos[i][0:trun])
        #pdb.set_trace()
        return [dists, pos]


# for each site in chain, returns distance of the closest active site
class aqW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        all_dists, all_pos = self.get_var_or_file(aoW, params, recalculate, True, True, False)
        true_states = self.get_var_or_file(oW, params, recalculate, True, True, False)
        num_sites = len(all_dists)
        closest_dists = [0 for i in range(num_sites)]
        for i in range(num_sites):
            closest = -1
            for j in range(num_sites):
                if true_states[all_pos[i][j]] == 1:
                    closest = all_dists[i][j]
                    break
            assert closest != -1
            closest_dists[i] = closest
#        pdb.set_trace()
        return closest_dists

# given a data_list_file, returns a huge list of features/true classes for every site in every chain in that list
class arW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        data_list = self.get_var_or_file(mW, params, recalculate, False, True, False)
        all_features = []
        all_classes = []
        for i in range(len(data_list)):
            pdb_name = data_list[i][0]
            chain_letter = data_list[i][1]
            this_features = self.get_var_or_file(jW, params, global_stuff.recalculate, True, True, False)
            all_features = all_features + this_features
            this_classes = self.get_var_or_file(oW, params, recalculate, True, True, False)
            all_classes = all_classes + this_classes
        return all_features, all_classes


# calls dssp
class avW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        input_f = self.get_var_or_file(fW, params, recalculate, False, False, False)
        subprocess.call([global_stuff.DSSP_PATH, input_f.name, self.get_holding_location()])
#        pdb.set_trace()
        return open(self.get_holding_location())

# parses dssp output into a dictionary
class awW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        f = self.get_var_or_file(avW, params, recalculate, False, False, False)
        lines = f.readlines()
        idx = 0
        while lines[idx][2] != '#':
            idx = idx + 1
        idx = idx + 1
        ans = {}

        for i in range(idx,len(lines)):
            s = string.split(lines[i])
            try:
                ans[(int(s[1]),s[2])] = lines[i][16]
            except:
                print 'error in reading dssp for ', self.get_param(params, 'pdb_name')

        return ans
                       
# creates raw naccess file output
class azW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        # first copy pdb file to have more manageable file name

        f = self.get_var_or_file(fW, params, recalculate, False, False, False)
        temp_location = self.get_holding_folder() + self.get_param(params, 'pdb_name') + '.pdb'
        subprocess.call(['cp', f.name, temp_location])
        # run Naccess which then puts the stuff in current working directory
        subprocess.call([global_stuff.NACCESS_PATH, temp_location])
        # delete the temp pdb file
        subprocess.call(['rm', temp_location])
        # move the newly created file to holding location
        rsa_location = self.get_param(params, 'pdb_name') + '.rsa'
        subprocess.call(['mv', rsa_location, self.get_holding_location()])
        # delete the other 2 files
        subprocess.call(['rm', self.get_param(params, 'pdb_name') + '.log'])
        subprocess.call(['rm', self.get_param(params, 'pdb_name') + '.asa'])
        return open(self.get_holding_location())

# dictionary for processed NACCESS values.  will also store the residue for error checking
class baW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        # a line is relevant only if first 3 characters are RES

        f = self.get_var_or_file(azW, params, recalculate, False, False, False)

        lines = f.readlines()
        ans_dict = {}
        for i in range(len(lines)):
            line = lines[i]
            if line[0:3] == 'RES':
                s = string.split(line)
                key = (int(s[3]), s[2])
                vals = [s[1]] + [float(x) for x in s[4:]]
                ans_dict[key] = vals
        return ans_dict


# gets the raw ligsite output
class bcW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        f = self.get_var_or_file(fW, params, recalculate, False, False, False)
        #pdb.set_trace()
        subprocess.call(['lcs','-i', f.name, '-o', self.get_holding_location()])
        #pdb.set_trace()
        return open(self.get_holding_location())
    
# gets processed ligsite output
class bdW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        
        f = self.get_var_or_file(bcW, params, recalculate, False, False, False)
        #pdb.set_trace()
        lines = f.readlines()
        ans = []
        for i in range(len(lines)):
            line = lines[i]
            s = string.split(line)
            ll = len(s)
            ans.append([float(s[ll-3]), float(s[ll-2]), float(s[ll-1])])
        #pdb.set_trace()
        return ans

# computes b-factor for a chain
class bfW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        chain = self.get_var_or_file(cW, params, recalculate, True, False, always_recalculate)
        ans = []
        for i in range(len(chain)):
            ans.append(global_stuff.get_representative_atom(chain[i]).get_bfactor())
        return ans


# list of the coordinates corresponding to each atom
class bgW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        chain = self.get_var_or_file(cW, params, recalculate, True, False, False)
        ans = []
        for i in range(len(chain)):
            ans.append(global_stuff.get_representative_atom(chain[i]).coord)
        return ans


# returns for each site a number between 0 and 1 representing how 'active' it is.
class bhW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        # get, for each site, the distance to closest active site
#        pdb.set_trace()
        closest_dists = self.get_var_or_file(aqW, params, recalculate, True, True, False)
        which_f = self.get_param(params, 'wtpr')
        if which_f == 0:
            def f(x):
                taper_c = self.get_param(params, 'nwc')
                return math.exp(taper_c * x)
        elif which_f == 1:
            def f(x):
                if abs(x) > .001:
                    return 0
                else:
                    return 1
        ans = [f(x) for x in closest_dists]

        return ans


# calculates pairwise KL for all residues within a certain cutoff, then averages them
class biW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        trun_sorted_dist, trun_sorted_pos = self.get_var_or_file(apW, params, recalculate, True, False, False)
        cutoff = self.get_param(params, 'micut')
        #find to which index to calculate stuff
        #pdb.set_trace()
        msa = self.get_var_or_file(agW, params, recalculate, True, True, False)
        ans = []
        for k in range(len(msa[0,:])):
            count = 0
            total = 0
            stop = 0
            for i in range(len(trun_sorted_dist[k])):
                stop = i
                if trun_sorted_dist[k][i] > cutoff:
                    break
            if stop <= 4:
                stop = 5
            for i in range(stop):
                for j in range(stop):
                    if i != j:
                        aa1 = trun_sorted_pos[k][i]
                        aa2 = trun_sorted_pos[k][j]
                        col1 = msa[:,aa1]
                        col2 = msa[:,aa2]
                        d1 = re.sub(r'-','',col1)
                        d2 = re.sub(r'-','',col2)
                        total = total + global_stuff.get_KL(d1,d2)
                        count = count + 1
            ans.append(total / float(count))
        return ans

# returns big matrix of features + labels + distance to closest site
class bkW(wrapper.mat_obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        data_list = self.get_var_or_file(mW, params, recalculate, False, True, False)
        ans = []
        j = 0
        for line in data_list:
            print line, j
            self.set_param(params, 'pdb_name', line[0])
            self.set_param(params, 'chain_letter', line[1])
            node_features = self.get_var_or_file(jW, params, recalculate, True, True, False)
            true_classes = self.get_var_or_file(oW, params, recalculate, True, True, False)
            dist = self.get_var_or_file(aqW, params, recalculate, True, True, False)
            num_sites = len(node_features)
            temp = [node_features[i] + [true_classes[i]] + [dist[i]] for i in range(num_sites)]
            ans = ans + temp
            j = j + 1
        return ans


# file written will contain param as a string, auroc, precrec area, iterations
class blW(wrapper.file_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        # first get instance of pW
        from wc import wc
        # figure out which roc file input wrapper to use
        which_classifier_roc_input = self.get_param(params, 'wclf')
        # if using svm, treat 'iter' as the fold number just so i don't have to change roc curve code yet
        #import try_svm
        # if which_classifier_roc_input == try_svm.atW:
        #    self.set_param(params, 'iter', self.get_param(params, 'wfld'))
        self.set_param(params, 'which_wrapper_class', which_classifier_roc_input)
        pW_instance = self.old_get_var_or_file(wc, params, True, False, False)
        if always_recalculate:
            f  = self.old_get_var_or_file(pW_instance.cache.file_dumper_wrapper, params, recalculate, False, False, 2)
        else:
            f  = self.old_get_var_or_file(pW_instance.cache.file_dumper_wrapper, params, recalculate, False, False, False)
        # call R roc area calculating script, specifying input location and where to write results
        source = f.name
        destination = self.get_holding_location()
        assert os.path.isfile(source)
        iteration = self.get_param(params, 'iter', record=False)
        obj_val = self.get_param(params, 'obj_val', record=False)
        subprocess.call(['Rscript', constants.ROC_INFO_SCRIPT, source, destination, str(iteration), str(obj_val)])
        return open(destination)


def print_stuff(x):
    #pdb.set_trace()
    print 'printing!!!!!!!!'
    print x


