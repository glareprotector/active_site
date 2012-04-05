'''
Created on Mar 10, 2012

@author: glareprotector
'''


import Bio.PDB
import constants
import global_stuff

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.PDB import Polypeptide
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline

import subprocess
import string
import os

class param(object):
    
    def __init__(self, param_dict):
        self.param_dict = param_dict
        
    def get_param(self, which):
        return self.param_dict[which]
    
    def get_keys(self):
        return self.param_dict.keys()
    
    def set_param(self,key,val):
        self.param_dict[key] = val
        
    def has_param(self,key):
        return key in self.param_dict.keys()
        
    def __add__(self, other):
        return param(dict(self.param_dict.items() + other.param_dict.items()))
    
    def __str__(self):
        return str(self.param_dict)
    
    def get_copy(self):
        to_return = param({})
        for key in self.param_dict.keys():
            to_return.set_param(key, self.get_param(key))
        return to_return
    
    # returns A U B, with values in A taking precedence
    @staticmethod
    def merge(A, B):
        toReturn = A.get_copy()
        for key in B.get_keys():
            if key not in A.get_keys():
                toReturn.set_param(key, B.get_param(key))
        return toReturn
        


class param_tree(object):
    
    def __init__(self, node_params, child_params_dict):
        self.node_params = node_params
        self.child_param_list = child_params_dict
        

'''
3 member variables: self.params, which will store all parameters, self.child_wrappers, and param_tree
'''
class wrapper(object):
    
    def _get_child_wrappers(self, self_params, leaf_params):
        raise NotImplementedError
                
    def get_name(self):
        return self.get_param('name')
    
    def get_hard_coded_params(self):
        raise NotImplementedError
    
    def get_param(self,key):
        return self.params.get_param(key)
    
    def set_param(self,key,val):
        self.params.set_param(key, val)
        
    def has_param(self, key):
        return self.params.has_param(key)
    # stores which params are relevant to the function.  
    def get_self_param_keys(self):
        raise NotImplementedError
    
    def get_self_params(self, inherited_params):
        hard_coded_params = self.get_hard_coded_params()
        return param.merge(inherited_params,hard_coded_params)
                
    def set_name(self):
        bits = [str(param) + '-' + str(self.get_param(param)) for param in self.get_self_param_keys()]
        bits.append(self.__class__.__name__)
        val = string.join(bits, sep='@')
        self.set_param('name', val)
    
    def __init__(self, inherited_params=param({}), recalculate = False):
        self.params = self.get_self_params(inherited_params)
        self.set_name() 
        self.recalculate = recalculate
        print self.params, self.get_name()

    def constructor(self, recalculate):
        raise NotImplementedError


class file_wrapper(wrapper):
        
    def get_file_location(self):
        return self.get_param('location') + self.get_name()
    

    
class obj_wrapper(wrapper):
        
    def get_pickle_location(self):
        return self.get_param('location') + self.get_name() + '.pickle'
        
'''
this is a leaf object.  it should always be obvious when something is a leaf, since it will always be a file
and it can either be assumed to be already present(file manually downloaded in which case constructor doesn't 
have to do anything, or can be downloaded from internet by some program
'''    
class pdb_file_wrapper(file_wrapper):
    
    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})
    
    def get_self_param_keys(self):
        return ['pdb_name']
    
    def constructor(self, recalculate):
        pdb_file_name = self.get_param('pdb_name')
        pdbl = Bio.PDB.PDBList()
        print pdb_file_name, self.get_param('location')
        pdbl.retrieve_pdb_file(pdb_file_name, pdir=self.get_param('location'))
        subprocess.call(['mv', self.get_param('location') + string.lower('pdb'+pdb_file_name+'.ent'), self.get_file_location()])
        
class pdb_chain_wrapper(obj_wrapper):
    
    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})
    
    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter']
    
    def constructor(self, recalculate):
        f = global_stuff.the_file_manager.get_file_handle(pdb_file_wrapper(self.params, self.recalculate))
        structure = Bio.PDB.PDBParser().get_structure(self.get_param('name'), f)
        chain = Bio.PDB.PPBuilder().build_peptides(structure[0][self.get_param('chain_letter')])
        to_return = []
        for chain_frag in chain:
            to_return = to_return + chain_frag
        # should raise exception here if error
        return to_return

class pdb_chain_seq_file_wrapper(file_wrapper):
    
    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})
    
    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter']
        
    def constructor(self, recalculate):
        chain_obj = global_stuff.the_obj_manager.get_variable(pdb_chain_wrapper(self.params, self.recalculate))
        # write the seq file at location + name
        raw_seq_string = ''.join([Polypeptide.three_to_one(res.resname) for res in chain_obj])
        seq = Bio.Seq.Seq(raw_seq_string)
        seq_record = Bio.SeqRecord.SeqRecord(seq)
        SeqIO.write(seq_record, self.get_file_location(), 'fasta')
        return open(self.get_file_location(),'r')

class pdb_chain_seq_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})
    
    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter']

    def constructor(self, recalculate):
        seq_file_handle = global_stuff.the_file_manager.get_file_handle(pdb_chain_seq_file_wrapper(self.params, self.recalculate))
        return SeqIO.read(seq_file_handle, 'fasta')
    
class pdb_chain_blast_results_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER,'evalue':10})

    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter', 'evalue']

    def constructor(self, recalculate):
        seq_records = []
        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_seq_file_wrapper(self.params, self.recalculate))
        query = SeqIO.parse(f, 'fasta')
        seq_records.append(query)
        print 'RUNNING BLAST!!!!!!!'
        psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLAST_PATH, outfmt = 5, query = f.name, db = 'nr', out = self.get_file_location())
        subprocess.Popen(str(psi_blast_cline), shell=True, executable='/bin/bash')

    
class pdb_chain_aa_to_pos_obj_wrapper(obj_wrapper):
    
    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})
        
    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter']
    
    def constructor(self, recalculate):
        chain_obj = global_stuff.the_obj_manager.get_variable(pdb_chain_wrapper(self.params, self.recalculate))
        chain_positions = [chain_obj[j].get_id()[1] for j in range(len(chain_obj))]
        return chain_positions 

class pdb_chain_msa_input_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER, 'evalue':1e-10, 'msa_input_max_num':100})

    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter', 'evalue', 'msa_input_max_num']

    def constructor(self, recalculate):
        # parse blast xml file, then do processing
        blast_xml_handle = global_stuff.the_file_manager.get_file_handle(pdb_chain_blast_results_file_wrapper(self.params, recalculate))
        record = NCBIXML.read(blast_xml_handle)
        seen = set()
        seq_records = []
        # add the query sequence, and have a set so that only add each sequence once
        query = global_stuff.the_obj_manager.get_variable(pdb_chain_seq_obj_wrapper(self.params, self.recalculate))
        seq_records.append(query)
        seen.add(query.seq.tostring())
        # add high scoring pairs in alignments with sufficiently low evalue that have not been seen
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < self.get_param('evalue') and not hsp.sbjct in seen:
                    seq_records.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(hsp.sbjct), id = alignment.hit_id))
        # write hits to fasta file
        output_handle = open(self.get_file_location(), 'w')
        SeqIO.write(seq_records, output_handle, 'fasta')
        print 'WROTE ', self.get_file_location()

class pdb_chain_pairwise_distance_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter']
    
    def get_representative_atom(self, res):
        if 'CA' in res.child_dict.keys():
            return res['CA']
        elif 'CB' in res.child_dict.keys():
            return res['CB']
        else:
            print 'no CA or CB atom in residue'
            raise RuntimeError

    def constructor(self, recalculate):
        residues = global_stuff.the_obj_manager.get_variable(pdb_chain_wrapper(self.params, self.recalculate))
        rep_atoms = [self.get_representative_atom(res) for res in residues]
        num_res = len(residues)
        dists = [[-1 for i in range(num_res)] for j in range(num_res)]
        for i in range(num_res):
            for j in range(num_res):
                dists[i][j] = rep_atoms[i] - rep_atoms[j]
        return dists

class pdb_chain_inverse_average_distances_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter']

    def constructor(self, recalculate):
        dists = global_stuff.the_obj_manager.get_variable(pdb_chain_pairwise_distance_obj_wrapper(self.params, self.recalculate))
        inv_avg_dists = [-1 for i in range(len(dists))]
        for i in range(len(dists)):
            val = 0;
            for j in range(len(dists)):
                val = val + dists[i][j]
            inv_avg_dists[i] = 1.0 / (val / len(dists))
        return inv_avg_dists

class pdb_chain_msa_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER, 'evalue':1e-10, 'msa_maxiter':4, 'msa_input_max_num':100})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter', 'evalue', 'msa_maxiter', 'msa_input_max_num']

    def constructor(self, recalculate):
        msa_input_handle = global_stuff.the_file_manager.get_file_handle(pdb_chain_msa_input_file_wrapper(self.params, recalculate))
        cline = MuscleCommandline(cmd = global_stuff.MUSCLE_PATH, input = msa_input_handle.name, out = self.get_file_location(), clw = False, maxiters = 2)
        cline()

class pdb_chain_conservation_score_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER, 'evalue':1e-10, 'msa_maxiter':4, 'msa_input_max_num':100})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter', 'evalue', 'msa_maxiter', 'msa_input_max_num']

    def constructor(self, recalculate):
        msa_file_handle = global_stuff.the_file_manager.get_file_handle(pdb_chain_msa_file_wrapper(self.params, recalculate))
        args = ['python', global_stuff.CONSERVATION_FOLDER+'score_conservation.py', '-m', global_stuff.CONSERVATION_FOLDER+'matrix/'+'blosum50.bla', '-o', self.get_file_location(), msa_file_handle.name]
        subprocess.Popen(args)


class crf_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter', 'dist_cut_off']

    def constructor(self, recalculate):
        return crf(self.params, recalculate)
