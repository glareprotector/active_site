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
import helper
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


import cross_validation_pseudo as cv


# doesn't matter if pdb_name and c are capitalized
class fW(wrapper.file_wrapper):

    def get_folder(self, object_key):
        return global_stuff.BIN_FOLDER + 'fWs' + '/'

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        pdb_file_name = self.get_param(params, 'p')
        pdbl = Bio.PDB.PDBList()
        pdbl.retrieve_pdb_file(pdb_file_name, pdir=self.get_holding_folder())
        subprocess.call(['mv', self.get_holding_folder() + string.lower('pdb'+pdb_file_name+'.ent'), self.get_holding_location()])

        return open(self.get_holding_location(), 'r')


        
class cW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        f = self.get_var_or_file(fW, params, recalculate, to_pickle)
        structure = Bio.PDB.PDBParser().get_structure(self.get_param(params, 'p'), f)
        chain = Bio.PDB.PPBuilder().build_peptides(structure[0][self.get_param(params, 'c')])
        to_return = []
        for chain_frag in chain:
            to_return = to_return + chain_frag
        # should raise exception here if error
        # now, only keep the parts of the chain whose positions are in the right range
        start_pos = self.get_param(params, 'st')
        end_pos = self.get_param(params, 'en')
        # search for corresponding positions
        start_aa = -1
        end_aa = -1

        for i in range(len(to_return)):
            if to_return[i].get_id()[1] == start_pos:
                start_aa = i
        for i in reversed(range(len(to_return))):
            if to_return[i].get_id()[1] == end_pos:
                end_aa = i


        assert start_aa != -1
        assert end_aa != -1

        #


        node_features_mode = self.get_param(params, 'nfm')

        if node_features_mode == 1:
            real_to_return = to_return[start_aa:end_aa+1]
        elif node_features_mode == 2:
            raw_text = self.get_var_or_file(chW, params, recalculate, True, False, False)
            # get a set of all of the positions for which there are features
            positions = set([int(row[0]) for row in raw_text])
            real_to_return = []
            for i in range(len(to_return)):
                if to_return[i].get_id()[1] in positions:
                    real_to_return.append(to_return[i])

        #pdb.set_trace()

        # at this point, there might still be the same position represented in 2 positions in the chain
        seen_pos = set()
        real_real_to_return = []
        for i in range(len(real_to_return)):
            if real_to_return[i].get_id()[1] not in seen_pos:
                real_real_to_return.append(real_to_return[i])
                seen_pos.add(real_to_return[i].get_id()[1])
            else:
                #pdb.set_trace()
                print 'double parking'


        # also filter out the positions that are not standard residues

        #pdb.set_trace()
        
        real_real_real_to_return = []
        for i in range(len(real_real_to_return)):
            try:
                Polypeptide.three_to_one(real_real_to_return[i].resname)
            except Exception, err:
                print err
                print real_real_to_return[i], i
                #pdb.set_trace()
                print 'weird AA'
            else:
                real_real_real_to_return.append(real_real_to_return[i])

        return real_real_real_to_return

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

        return pos_to_aa_dict

class gW(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        residues = self.get_var_or_file(cW, params, recalculate, True)
        rep_atoms = [helper.get_representative_atom(res) for res in residues]
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
                elif res_dists[i][j] < self.get_param(params, 'co') and math.fabs(i-j) > 5:
                    edges.append([int(i),int(j)])

        return edges

# returns the raw file associated with a catres-fischer data_set
class chW(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        pdb_name = self.get_param(params, 'p')
        c = self.get_param(params, 'c')
        start = self.get_param(params, 'st')
        end = self.get_param(params, 'en')
        aux_folder = helper.get_aux_folder(pdb_name, c, start, end)
        node_feature_file = aux_folder + pdb_name + '.features'
        f = open(node_feature_file, 'r')
        to_return = []
        for line in f:
            s = line.strip().split('\t')
            pos = s[0]
            feature_string = s[1]
            this_features = feature_string.split(',')
            to_return.append([pos] + this_features + [s[-1]])
        return to_return
            

# returns data list but with samples lacking any true sites omitted
class ciW(wrapper.obj_wrapper):
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        data_list_file = self.get_var_or_file(mW, params, recalculate, True, False)
        to_return = []
        pdb.set_trace()
        num_good = 0
        num_bad = 0
        asdf = 0
        for sample in data_list_file:
            print sample
            print asdf, len(data_list_file)
            asdf += 1

            self.set_param(params, 'p', sample.pdb_name)
            self.set_param(params, 'c', sample.chain_letter)
            self.set_param(params, 'st', sample.start)
            self.set_param(params, 'en', sample.end)
            #true_states = self.get_var_or_file(oW, params, recalculate, False, False)
            #num_true = sum(true_states)
            feat = self.get_var_or_file(chW, params, recalculate, False, False)

            # get the number of true sites
            states = [int(row[-1]) for row in feat]
            pos_pos = []
            for i in range(len(states)):
                if states[i] == -1:
                    states[i] = 0
                else:
                    pos_pos.append(int(feat[i][0]))
            num_true = sum(states)
            
            # see if any true sites are heteroatoms.  do so by looking at
            f = self.get_var_or_file(fW, params, False, to_pickle)
            structure = Bio.PDB.PDBParser().get_structure(self.get_param(params, 'p', False), f)


            passes_hetero = True
            for pos in pos_pos:
                try:
                    Polypeptide.three_to_one(structure[0][sample.chain_letter][pos].resname)
                except:
                    print 'true site was hetero', pos

                    passes_hetero = False
                

            if num_true > 0 and passes_hetero:
                to_return.append(sample)
                num_good += 1
            else:
                for i in range(len(data_list_file)):
                    if data_list_file[i].pdb_name == sample.pdb_name:
                        print data_list_file[i], 'OTHER'
                print sample, 'bAD'
                num_bad += 1

        print num_good, num_bad
        import time
        time.sleep(10)

        return to_return


# there are now 2 modes for this.  read from file, or calculate, depending on nmd: node feature mode
class jW(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    # params will contain node_feature_list
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = True, always_recalculate = False, old_obj = None):

        node_feature_mode = self.get_param(params, 'nfm')

        if node_feature_mode == 1:

            feature_list = self.get_param(params, 'n')
            aa_to_pos = self.get_var_or_file(aW, params, recalculate, True)
            node_features = []
            for i in range(len(aa_to_pos)):
                pos = aa_to_pos[i]
                self.set_param(params, 'pos', pos)
                temp = []
                for fxn_w in feature_list:
                    temp = temp + self.get_var_or_file(fxn_w, params, recalculate, False)
                    node_features.append(temp)
            return node_features

        elif node_feature_mode == 2:

            #reading from file.
            pdb_name = self.get_param(params, 'p')
            c = self.get_param(params, 'c')
            start = self.get_param(params, 'st')
            end = self.get_param(params, 'en')
            aux_folder = helper.get_aux_folder(pdb_name, c, start, end)
            node_feature_file = aux_folder + pdb_name + '.features'
            f = open(node_feature_file, 'r')
            # afraid that not all positions got features assigned
            pos_to_aa = self.get_var_or_file(eW, params, recalculate, False, False, False)
            set_pos = set()
            num_nodes = len(pos_to_aa.keys())
            node_features = [-1 for i in range(num_nodes)]

            for line in f:

                s = line.strip().split('\t')
                pos = int(s[0])
                feature_string = s[1]
                this_features = [1.0] + [float(x) for x in feature_string.split(',')]
                # some positions in file will not be in the chain (due to say weird AA)
                try:
                    node_features[pos_to_aa[pos]] = this_features
                except Exception, err:

                    print pos

                # keep track of which positions i have features for.  at the end, this set should include all of the keys in pos_to_aa.
                # i got pos_to_aa by taking all positions in the specified range that worked
                set_pos.add(pos)

            # assert that all positions have been accounted for

            assert(len(set(pos_to_aa.keys()) - set_pos) == 0)

            return node_features
                
            

class kW(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    # params will contain edge_feature_list
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        feature_list = self.get_param(params, 'e')
        aa_to_pos = self.get_var_or_file(aW, params, recalculate, True)
        edge_list = self.get_var_or_file(iW, params, recalculate, True)
        edge_features = []
        for i in range(len(edge_list)):
            try:
                edge = edge_list[i]
                pos1 = aa_to_pos[edge[0]]
                pos2 = aa_to_pos[edge[1]]
                self.set_param(params, 'pos1', pos1)
                self.set_param(params, 'pos2', pos2)
                temp = []
                for fxn_w in feature_list:
                    temp = temp + self.get_var_or_file(fxn_w, params, recalculate, False)
                edge_features.append(temp)
            except Exception, err:
                print 'error in getting edge features'
                raise Exception

        return edge_features

class lW(wrapper.obj_wrapper, wrapper.shorten_name_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        # this simply extracts every single parameter, puts them into a dict, and prints the dict.  it also sets params if they weren't gotten here
        print params
        #pdb.set_trace()
        the_dict = {}
        the_dict["dist_cut_off"] = self.get_param(params, "co");
        # data_list will be list of tuples of (pdb_name, c)
        the_dict["d"] = self.get_param(params, "d")
        the_dict['n'] = self.get_param(params, "n")
        the_dict['e'] = self.get_param(params, "e")
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

class mW(wrapper.obj_wrapper):

    # params should include file name want to read - 'data_list_file'
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(self.get_param(params, 'd'), 'r')
        ans = []
        for line in f:
            line = string.split(string.strip(line), sep=',')
            ans.append(cv.pdb_name_struct(line[0], line[1], int(line[2]), int(line[3])))
            # pdb.set_trace()
        print ans
        return ans

# this takes merged regular results folder and puts it in the form needed by roc input and alternative measure
# so it has to get the results object somehow.  there will be 2 ways.  there will be 2 descendants.  1 will have results struct in params, and simply get it and return it.  another will call c++ function to get the wrapper.
# there will be a parameter specifying which descendant to call
class nW(wrapper.mat_obj_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

    # params will be [(scores, sizes, pdb_name, c),   ]
    # params which which it is stored does NOT include these things...only data_list, params for getting features
    # actually, now call the interface library to get a pdb_results_struct
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):


        mode = self.get_param(params, 'md')



        if mode == 0:
            # in perfect world, this would have been generated model file, which then generates results.  so simulate this by getting parameters from experiment info


            #experiment_info_pretending_to_be_model_params = self.get_var_or_file(lW, params, recalculate, False, True)
            results = self.get_param(params, 'results', False)
            

            """
            scores = self.get_param(params, 'scores', False)
            true_states = self.get_param(params, 'true_states', False)
            sizes = self.get_param(params, 'sizes', False)
            pdb_names = self.get_param(params, 'pdb_names', False)
            cs = self.get_param(params, 'cs', False)
            """

        elif mode == 1:
            # trying to get outer cv results
            the_overall_results = self.get_var_or_file(btW, params, recalculate, False, False, False)
            results = the_overall_results.get_raw_results()

        elif mode == 2:
            # trying to get one train_test results

            the_train_test_results = self.get_var_or_file(bwW, params, recalculate, True, True, False)
            results = the_train_test_results.get_raw_results()

        elif mode == 3:
            # trying to get inner cv results
            # need to have ifold specified
            the_cv_results = self.get_var_or_file(cdW, params, recalculate, False, False, False)

            results = the_cv_results.get_raw_results()



        scores = results.scores
        true_states = results.true_classes
        sizes = results.sample_lengths
        pdb_struct_names = results.pdb_structs
        pdb_names = [x.pdb_name for x in pdb_struct_names]
        cs = [x.chain_letter for x in pdb_struct_names]
        sample_starts = [x.start for x in pdb_struct_names]
        sample_ends = [x.end for x in pdb_struct_names]
        
        num_samples = len(pdb_names)
        pos = 0
        mat = []
        # alternate between pdb_names, c and scores
        for i in range(num_samples):
            mat.append([pdb_names[i], cs[i], sizes[i], sample_starts[i], sample_ends[i]])
            print [pdb_names[i], cs[i], sizes[i], sample_starts[i], sample_ends[i]]
            mat.append(scores[pos:pos+sizes[i]])
            mat.append(true_states[pos:pos+sizes[i]])
            pos = pos + sizes[i]



        return mat

# for features i can't get automatically, assume that there is folder for each chain containing the features.  wrappers will then read from those folders
# true_states file depends on which dataset i am working with.  store this as parameter that i don't keep track of
class oW(wrapper.vect_obj_wrapper, wrapper.by_pdb_folder_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        pos_to_aa = self.get_var_or_file(eW, params, recalculate, True, False)

        pdb_name = self.get_param(params, 'p')
        c = self.get_param(params, 'c')
        start = self.get_param(params, 'st')
        end = self.get_param(params, 'en')
        aux_folder = helper.get_aux_folder(pdb_name, c, start, end)
        true_states = [0 for i in range(len(pos_to_aa.keys()))]
        chain_seq_in_one = self.get_var_or_file(dW, params, recalculate, True, False)


        true_states_mode = self.get_param(params, 'tsmd', False)

        if true_states_mode == 1:

            states_file = aux_folder + string.lower(pdb_name) + '.catres'
            f = open(states_file, 'r')

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

        elif true_states_mode == 2:
            # will just read from the node features file
            features_file = aux_folder + pdb_name + '.features'
            f = open(features_file, 'r')
            for line in f:
                s = line.strip().split('\t')
                pos = int(s[0])
                state = int(s[2])
                try:
                    if state == 1:
                        true_states[pos_to_aa[pos]] = 1
                except Exception, err:
                    pdb.set_trace()
                    print err
            return true_states


class pW(wrapper.mat_obj_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        #self.get_param(params, 'md')

        results = self.get_var_or_file(nW, params, recalculate, to_pickle, False, always_recalculate)

        assert len(results) % 3 == 0
        num_samples = len(results) / 3
        roc_classes = []
        roc_scores = []
        # if using regular results struct
        for i in range(num_samples):
            pdb_name = results[3 * i][0]
            c = results[3 * i][1]
            start = results[3 * i][3]
            end = results[3 * i][4]
            scores = results[(3 * i) + 1]
            true_states_from_cpp = results[(3 * i) + 2]
            self.set_param(params, 'p', pdb_name)
            self.set_param(params, "c", c)
            self.set_param(params, 'st', start)
            self.set_param(params, 'en', end)

            true_states = self.get_var_or_file(oW, params, global_stuff.recalculate, True, True)
            roc_classes = roc_classes + true_states
            roc_scores = roc_scores + scores
            assert(len(roc_classes) == len(roc_scores))

        print [(roc_classes[i],roc_scores[i]) for i in range(len(roc_scores)) if roc_classes[i] == 1 ]
            
        ans = helper.get_transpose([roc_classes, roc_scores])
        return ans


class qW(wrapper.file_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

    # this part does not work
    # params will be those required by experiment_results, which are those required by experiment_info.
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = True, old_obj = None):


        #self.get_param(params, 'md')


        # first get instance of pW
        from wc import wc
        # figure out which roc file input wrapper to use
        #which_classifier_roc_input = self.get_param(params, 'wclf')
        which_classifier_roc_input = pW
        # if using svm, treat 'iter' as the fold number just so i don't have to change roc curve code yet
        #import try_svm
        #if which_classifier_roc_input == try_svm.atW:
        #    self.set_param(params, 'iter', self.get_param(params, 'wfld'))
        self.set_param(params, 'which_wrapper_class', which_classifier_roc_input)
        pW_instance = self.old_get_var_or_file(wc, params, True, False, False)

        print "OOOOOHHHHHH"
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
        #iteration = self.get_param(params, 'mx', record=False)
        #obj_val = self.get_param(params, 'obj_val', record=False)
        iteration = '.'
        obj_val = '.'
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
        pdb_name = self.get_param(params, 'p')
        c = self.get_param(params, 'c')
        aux_folder = helper.get_aux_folder(pdb_name, c)
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

            while 1:
                try:
                    pdb_name = data_list[i].pdb_name
                    #c = data_list[i].c
                    print pdb_name
                    self.set_param(params, 'p', pdb_name)
                    self.get_var_or_file(which_wrapper, params, recalculate, False, False)
                except:
                    print 'error while downloading ', data_list[i]

                else:
                    break
    
# blast results file wrapper(xml format)
class adW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):



    def whether_to_override(self, location):
        
        #if the file size is too small, we know there was something wrong
        import os
        if os.path.getsize(location) < 1:
            return True

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
                if hsp.expect < self.get_param(params, 'ev') and not hsp.sbjct in seen:
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
class ahW(wrapper.mat_obj_wrapper, wrapper.experiment_type_wrapper, wrapper.shorten_name_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = True, old_obj = None):

#        pdb.set_trace()
        #self.get_param(params, 'md')

        results = self.get_var_or_file(nW, params, recalculate, False, False, always_recalculate)
        num_samples = len(results) / 3
        cutoffs = global_stuff.metric_cutoffs
        these_results = [0 for i in range(len(cutoffs))]
        total_num_sites = 0
#        pdb.set_trace()
        for i in range(num_samples):
            pdb_name = results[3*i][0]

            c = results[3*i][1]
            start = results[3*i][3]
            end = results[(3*i)][4]
            scores = results[(3*i)+1]
            true_states = results[(3*i)+2]
            # sort scores by marginals
            temp = [ [scores[i], i] for i in range(len(scores))]


            sorted_temp = sorted(temp, key = lambda x: x[0], reverse = True)
            self.set_param(params, 'p', pdb_name)
            self.set_param(params, 'c', c)
            self.set_param(params, 'st', start)
            self.set_param(params, 'en', end)

            # for each best site, see if any of the true sites are within that cutoff
            true_sites = self.get_var_or_file(anW, params, recalculate, True, True, False)
            some_dists = self.get_var_or_file(amW, params, recalculate, True, True, False)
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

        to_write = [self.get_param(params, 'mx', record=False), total_num_sites] + these_results

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
            pdb_name = data_list[i].pdb_name
            c = data_list[i].c
            this_features = self.get_var_or_file(jW, params, global_stuff.recalculate, True, True, False)
            all_features = all_features + this_features
            this_classes = self.get_var_or_file(oW, params, recalculate, True, True, False)
            all_classes = all_classes + this_classes
        return all_features, all_classes


# calls dssp
class avW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):

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
                print 'error in reading dssp for ', self.get_param(params, 'p')

        return ans
                       
# creates raw naccess file output
class azW(wrapper.file_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        # first copy pdb file to have more manageable file name
        #pdb.set_trace()
        f = self.get_var_or_file(fW, params, recalculate, False, False, False)
        temp_location = self.get_holding_folder() + self.get_param(params, 'p') + '.pdb'
        temp_location = global_stuff.NACCESS_FOLDER + self.get_param(params, 'p') + '.pdb'
        subprocess.call(['cp', f.name, temp_location])
        # run Naccess which then puts the stuff in current working directory
        subprocess.call([global_stuff.NACCESS_PATH, temp_location])
        # delete the temp pdb file
        subprocess.call(['rm', temp_location])
        # move the newly created file to holding location
        rsa_location = self.get_param(params, 'p') + '.rsa'
        subprocess.call(['mv', rsa_location, self.get_holding_location()])
        # delete the other 2 files
        subprocess.call(['rm', self.get_param(params, 'p') + '.log'])
        subprocess.call(['rm', self.get_param(params, 'p') + '.asa'])
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

        lines = f.readlines()
        ans = []
        for i in range(len(lines)):
            line = lines[i]
            s = string.split(line)
            ll = len(s)
            ans.append([float(s[ll-3]), float(s[ll-2]), float(s[ll-1])])

        return ans

# computes b-factor for a chain
class bfW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        #pdb.set_trace()
        chain = self.get_var_or_file(cW, params, recalculate, True, False, always_recalculate)
        ans = []
        for i in range(len(chain)):
            ans.append(helper.get_representative_atom(chain[i]).get_bfactor())
        return ans


# list of the coordinates corresponding to each atom
class bgW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        chain = self.get_var_or_file(cW, params, recalculate, True, False, False)
        ans = []
        for i in range(len(chain)):
            ans.append(helper.get_representative_atom(chain[i]).coord)
        return ans


# returns for each site a number between 0 and 1 representing how 'active' it is.
class bhW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        closest_dists = self.get_var_or_file(aqW, params, recalculate, True, True, False)
        which_f = self.get_param(params, 'wnv')
        if which_f == 0:
            c = self.get_param(params, 'nvec')
            def f(x):
                return math.exp(taper_c * x)
        elif which_f == 1:
            cut_off = self.get_param(params, 'nvjd')
            def f(x):
                if abs(x) > cut_off:
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
                        total = total + helper.get_KL(d1,d2)
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
            self.set_param(params, p, line.pdb_name)
            self.set_param(params, 'c', line.c)
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
        pdb.set_trace()
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



class blW2(wrapper.file_wrapper, wrapper.experiment_results_wrapper, wrapper.shorten_name_wrapper):

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
        #iteration = self.get_param(params, 'mx', record=False)
        #obj_val = self.get_param(params, 'obj_val', record=False)
        
        subprocess.call(['Rscript', constants.PREC_POINTS_SCRIPT, source, destination, '.', '.'])
        return open(destination)




# returns node_features for a sample, but normalized within the sample
class bmW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        
        node_features = self.get_var_or_file(jW, params, recalculate, True, True, False)
        return helper.normalize_mat(node_features)
                
# returns edge_features for a sample, but normalized within the sample
class bnW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        edge_features = self.get_var_or_file(kW, params, recalculate, True, True, False)
        return helper.normalize_mat(edge_features)



# returns a data.  data_list_file should be specified
class brW(wrapper.obj_wrapper, wrapper.shorten_name_wrapper):

    def whether_to_override(self):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):

        data_list_file = self.get_param(params, 'd')
        return cv.data(self, params, recalculate, data_list_file)

# returns a fold.  source instance and num folds and which fold should be specified in params
class buW(wrapper.obj_wrapper, wrapper.shorten_name_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        source = self.get_param(params, 's')
        num_folds = self.get_param(params, 'm')
        which_fold = self.get_param(params, 'k')
        return cv.fold(self, params, recalculate, source, which_fold, num_folds)


# returns sorted distances(previous thing returned a tuple of distances and positions)
class bxW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        trun_sorted_dist, trun_sorted_pos = self.get_var_or_file(apW, params, recalculate, True, False, False)
        return trun_sorted_dist


# returns positions corresponding to closest positions
class byW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        trun_sorted_dist, trun_sorted_pos = self.get_var_or_file(apW, params, recalculate, True, False, False)
        return trun_sorted_pos


# for a specified fold, returns the result of training on fold's training, testing on fold's testing.  actually, returns (fold, results, object_key) list
# hyperparams should also be passed as separate parameter
class bwW(wrapper.obj_wrapper, wrapper.shorten_name_wrapper, wrapper.experiment_type_wrapper):

    #def get_folder(self, object_key):
    #    folder = str(self.get_param(object_key, 'wif', False)) + '_' + str(self.get_param(object_key, 'wob', False)) + '/'
    #    return global_stuff.RESULTS_FOLDER + folder

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):

        print self.get_param(params, 'hp', False)
        import time

        return cv.train_test_result(self, params, recalculate)




# returns the results of outermost cv
class btW(wrapper.obj_wrapper, wrapper.shorten_name_wrapper, wrapper.experiment_results_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
#        pdb.set_trace()
        return cv.overall_results(self, params, recalculate)














class test(wrapper.obj_wrapper):
    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        a_data = self.get_var_or_file(brW, params, True, False, False)
        self.set_param(params, 'f')
        folds = a_data.get_folds(self, params, recalculate, 4)
        return folds, folds[0].get_folds(self, params, recalculate, 4)


# ideally, param contains an object that implements get grid points method
class ccW(wrapper.obj_wrapper):
    
    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):

        hp_values = helper.read_hp_values(self.get_param(params, 'hpvf'))
        return hp_values
        """
        the_file = constants.HP_VALUES_FOLDER + self.get_param(params, 'hpvf')
        # format is param name followed by the values
        f = open(the_file, 'r')
        the_dict = {}
        for line in f:
            #s = line.split(line, ',')
            s = line.strip().split(',')
            param_name = s[0]
            param_values = []
            for i in range(1, len(s)):
                param_values.append(float(s[i]))
            the_dict[param_name] = param_values
        return param.param(the_dict)
        """

class caW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):



        # instead, will read hp_values from another wrapper that just reads files
        hp_values = self.get_var_or_file(ccW, params, False, False, False, False)
        #hp_values = self.get_param(params, 'hpv')
        total_jobs = self.get_param(params, 'tj')
        which_job = self.get_param(params, 'wj')

        param_combos = cv.cross_product(hp_values.sorted_listify())
        sorted_keys = hp_values.get_sorted_keys()
        hp_stash = []
        for i in range(len(param_combos)):
            if i % total_jobs == which_job:
                combo = param_combos[i]
                temp = {}
                for i in range(len(sorted_keys)):
                    temp[sorted_keys[i]] = combo[i]
                hp_stash.append(param.param(temp))
        return hp_stash
        

# hp_values, total_jobs, which_job should be specified in params, as well as which fold
class cbW(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):

        a_hp_searcher = cv.hp_searcher(self, params, recalculate)
        self.set_param(params, 'hps', a_hp_searcher)

        return cv.hp_search_results(self, params, recalculate)
        

# make wrapper that returns cv search results
class cdW(wrapper.obj_wrapper, wrapper.experiment_results_wrapper):
    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):

        return cv.cv_results(self, params, recalculate)

class ceW(wrapper.mat_obj_wrapper, wrapper.experiment_results_wrapper):
    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):

        a_data = self.get_var_or_file(brW, params, True, False, False)
        self.set_param(params, 'f', a_data)
        self.set_param(params, 'md', 1)
        the_overall_results = self.get_var_or_file(btW, params, recalculate, False, False, False)

        self.get_var_or_file(qW, params, False, False, True)

        metric = self.get_var_or_file(ahW, params, False, True, True)
        results = metric , the_overall_results.best_hps

        return results


# returns a positive weight for each position
class cfW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):

        closest_dists = self.get_var_or_file(aqW, params, recalculate, True, True, False)

        which = self.get_param(params, 'ww')
        if which == 0:
            c = self.get_param(params, 'wec')
            def f(x):
                return math.exp(c*x)

        elif which == 1:
            cut_off = self.get_param(params, 'wjd')
            posw = self.get_param(params, 'wpw')
            def f(x):
                if x <= cut_off:
                    return posw
                else:
                    return 1.0

        return [f(x) for x in closest_dists]

# for a chain, returns the min and max positions
class cgW(wrapper.obj_wrapper):

    def get_folder(self, object_key):

        return constants.BIN_FOLDER + 'cgWs' + '/'

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, old_obj = None):
        f = self.get_var_or_file(fW, params, recalculate, to_pickle)
        structure = Bio.PDB.PDBParser().get_structure(self.get_param(params, 'p'), f)
        chain = Bio.PDB.PPBuilder().build_peptides(structure[0][self.get_param(params, 'c')])
        to_return = []
        for chain_frag in chain:
            to_return = to_return + chain_frag
        pdb.set_trace()
        min_pos = to_return[0].get_id()[1]
        max_pos = to_return[-1].get_id()[1]
        return [min_pos, max_pos]

# for chain segment, looks at the supplied feature file and returns whether 


def print_stuff(x):
    #pdb.set_trace()
    print 'printing!!!!!!!!'
    print x




