from wrapper import *
    
class pdb_file_wrapper(file_wrapper):
    
    def constructor(self, params, recalculate, to_pickle):
        pdb_file_name = self.get_param(params, 'pdb_name')
        pdbl = Bio.PDB.PDBList()
        pdbl.retrieve_pdb_file(pdb_file_name, pdir=self.get_holding_folder(params))
        subprocess.call(['mv', self.get_holding_folder(params) + string.lower('pdb'+pdb_file_name+'.ent'), self.get_holding_location(params)])

the_pdb_file_wrapper = pdb_file_wrapper()

f = the_pdb_file_wrapper.constructor(param({'pdb_name':'1asy'}), True, False)

# note that the file version is constrained to be in the same folder as the object version
class file_cache_decorator(object):

    def __init__(self, f):
        self.f = f

    def __get__(self, inst, cls):
        self.obj = inst

    def write_result_to_file(self, result, file_name):
        pass

    def __call__(self, *args):
        file_name = self.obj.get_location() + self.obj.get_name() + ".file_version"
        result = self.f(*args)
        self.write_result_to_file(result, file_name)
        return result

class vect_file_cache_decorator(file_cache_decorator):

    def write_result_to_file(self, result, file_name):
        global_stuff.write_vect(result, file_name)

class mat_file_cache_decorator(file_cache_decorator):

    def write_result_to_file(self, result, file_name):
        global_stuff.write_mat(result, file_name)
        

#this is wrapper that knows how to write a mat.  constructor is not defined
class vect_file_wrapper(file_wrapper):

    # vect to write is stored as "vect" param
    def write_vect(self, vect):
        global_stuff.write_vect(vect, self.get_file_location())

class mat_file_wrapper(file_wrapper):

    # the matrix to write is stored in params as 'mat'
    def write_mat(self, mat):
        global_stuff.write_mat(mat, self.get_file_location())

class pdb_chain_score_file_wrapper(vect_file_wrapper, pdb_chain_specific_wrapper, in_folder_wrapper):

    # scores to write are in the "scores" param
    def constructor(self, recalculate):
        self.write_vect(self.get_param('score'))

class crf_obj_wrapper(obj_wrapper, pdb_chain_specific_wrapper):

    @classmethod
    def get_self_attr_keys(cls):
        return ['dist_cut_off']

    # make the objects I make have constructors that support the params input for __init__ so i don't have to rename stuff
    def constructor(self, recalculate):
        return crf(self.params, recalculate)

# get node features file.  current convention is to force creating file version of any object before getting actual object
# this forces there to be a text backup of object
class pdb_chain_node_features_file_wrapper(pdb_chain_specific_wrapper, in_folder_wrapper, file_wrapper, mat_file_wrapper):

    @classmethod
    def get_dependent_wrappers(cls):
        return [crf_obj_wrapper]

    def constructor(self, recalculate):
        the_crf = global_stuff.the_obj_manager.get_variable(crf_obj_wrapper(self.params), recalculate)
        Xnode = the_crf.Xnode
        self.write_mat(Xnode)

class pdb_chain_edge_features_file_wrapper(pdb_chain_specific_wrapper, in_folder_wrapper, file_wrapper, mat_file_wrapper):

    @classmethod
    def get_dependent_wrappers(cls):
        return [crf_obj_wrapper]

    def constructor(self, recalculate):
        the_crf = global_stuff.the_obj_manager.get_variable(crf_obj_wrapper(self.params), recalculate)
        Xedge = the_crf.Xedge
        self.write_mat(Xedge)

class pdb_chain_edge_list_file_wrapper(pdb_chain_specific_wrapper, in_folder_wrapper, file_wrapper, mat_file_wrapper):

    @classmethod
    def get_dependent_wrappers(cls):
        return [crf_obj_wrapper]

    def constructor(self, recalculate):
        the_crf = global_stuff.the_obj_manager.get_variable(crf_obj_wrapper(self.params), recalculate)
        edges = the_crf.edges
        self.write_mat(edges)

class pdb_chain_true_states_file_wrapper(pdb_chain_specific_wrapper, in_folder_wrapper, file_wrapper, mat_file_wrapper):

    @classmethod
    def get_dependent_wrappers(cls):
        return [crf_obj_wrapper]

    def constructor(self, recalculate):
        the_crf = global_stuff.the_obj_manager.get_variable(crf_obj_wrapper(self.params), recalculate)
        true_states = the_crf.true_y
        self.write_mat(true_y)

class pdb_chain_info_file_wrapper(pdb_chain_specific_wrapper, in_folder_wrapper, file_wrapper, vect_file_wrapper):

    @classmethod
    def get_dependent_wrappers(cls):
        return [crf_obj_wrapper]

    def constructor(self, recalculate):
        the_crf = global_stuff.the_obj_manager.get_variable(crf_obj_wrapper(self.params), recalculate)
        to_write = [the_crf.numNodes, the_crf.numEdges, the_crf.numStates, the_crf.numFeatures, the_crf.numParams]
        self.write_vect(to_write)

# all of the info needed to make a crf returned as one long tuple of lots of stuff
class pdb_chain_data_obj_wrapper(pdb_chain_specific_wrapper, obj_wrapper):

    @classmethod
    def get_dependent_wrappers(cls):
        return [pdb_chain_true_states_file_wrapper, pdb_chain_node_features_file_wrapper, pdb_chain_edge_features_file_wrapper, pdb_chain_edge_list_file_wrapper]

    def constructor(self, recalculate):
        # read in all of the stuff
        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_true_states_file_wrapper(self.params), recalculate)
        true_states = global_stuff.read_vect_to_int(f)
        
        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_node_features_file_wrapper(self.params), recalculate)
        node_features = global_stuff.read_vect_to_float(f)

        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_edge_features_file_wrapper(self.params), recalculate)
        edge_features = global_stuff.read_vect_to_float(f)

        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_edge_list_file_wrapper(self.params), recalculate)
        edge_list = global_stuff.read_vect_to_float(f)

        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_info_file_wrapper(self.params), recalculate)
        info = global_stuff.read_vect_to_float(f)

        return [node_features, edge_features, edge_list, edge_list, info[0], info[1]]

class pdb_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})
    
    def get_self_param_keys(self):
        return ['pdb_name']

    def constructor(self, recalculate):
        f = global_stuff.the_file_manager.get_file_handle(pdb_file_wrapper(self.params), recalculate)    
        structure = Bio.PDB.PDBParser().get_structure(self.get_param('name'), f)
        return structure
        
class pdb_chain_wrapper(obj_wrapper):
    
    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})
    
    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter']
    
    def constructor(self, recalculate):
        f = global_stuff.the_file_manager.get_file_handle(pdb_file_wrapper(self.params), recalculate)
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
        chain_obj = global_stuff.the_obj_manager.get_variable(pdb_chain_wrapper(self.params), recalculate)
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
        seq_file_handle = global_stuff.the_file_manager.get_file_handle(pdb_chain_seq_file_wrapper(self.params), recalculate)
        #pdb.set_trace()
        asdf = SeqIO.read(seq_file_handle, 'fasta')
        return asdf
    
class pdb_chain_blast_results_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER,'evalue':10})

    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter', 'evalue']

    def constructor(self, recalculate):
        seq_records = []
        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_seq_file_wrapper(self.params), recalculate)
        query = SeqIO.parse(f, 'fasta')
        seq_records.append(query)
        print 'RUNNING BLAST!!!!!!!'
        psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLAST_PATH, outfmt = 5, query = f.name, db = 'nr', out = self.get_file_location())
        subprocess.call(str(psi_blast_cline), shell=True, executable='/bin/bash')
        #pdb.set_trace()

    
class pdb_chain_aa_to_pos_obj_wrapper(obj_wrapper):
    
    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})
        
    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter']
    
    def constructor(self, recalculate):
        chain_obj = global_stuff.the_obj_manager.get_variable(pdb_chain_wrapper(self.params), recalculate)
        chain_positions = [chain_obj[j].get_id()[1] for j in range(len(chain_obj))]
        return chain_positions


class pdb_chain_pos_to_aa_dict_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter']

    def constructor(self, recalculate):
        aa_to_pos = global_stuff.the_obj_manager.get_variable(pdb_chain_aa_to_pos_obj_wrapper(self.params), recalculate)
        pos_to_aa_dict = {}
        for i in range(len(aa_to_pos)):
            pos_to_aa_dict[aa_to_pos[i]] = i
        return pos_to_aa_dict

class pdb_chain_msa_input_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER, 'evalue':1e-10, 'msa_input_max_num':100})

    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter', 'evalue', 'msa_input_max_num']

    
    def constructor(self, recalculate):
        # parse blast xml file, then do processing
        blast_xml_handle = global_stuff.the_file_manager.get_file_handle(pdb_chain_blast_results_file_wrapper(self.params), recalculate)
        record = NCBIXML.read(blast_xml_handle)
        seen = set()
        seq_records = []
        # add the query sequence, and have a set so that only add each sequence once
        query = global_stuff.the_obj_manager.get_variable(pdb_chain_seq_obj_wrapper(self.params), recalculate)
        query.id = 'QUERY'
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
        return ['pdb_name', 'chain_letter']
    
    def constructor(self, recalculate):
        residues = global_stuff.the_obj_manager.get_variable(pdb_chain_wrapper(self.params), recalculate)
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
        #print dists
        #pdb.set_trace()
        return dists

# for a site, tuple consisting of the sites indicies/distances sorted by distance from the site
class pdb_chain_site_sorted_distances_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter', 'aa']

    def constructor(self, recalculate):
        dists = global_stuff.the_obj_manager.get_variable(pdb_chain_pairwise_distance_obj_wrapper(self.params), recalculate)
        dist = dists[self.get_param('aa')]
        temp = [ [dist[i],i] for i in range(len(dist)) ]
        sorted_temp = sorted(temp, key = lambda x: x[0])
        return [sorted_temp[i][0] for i in range(len(sorted_temp))], [sorted_temp[i][1] for i in range(len(sorted_temp))]
        

class pdb_chain_inverse_average_distances_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return ['pdb_name', 'chain_letter']

    def constructor(self, recalculate):
        dists = global_stuff.the_obj_manager.get_variable(pdb_chain_pairwise_distance_obj_wrapper(self.params), recalculate)
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

class pdb_chain_msa_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER, 'evalue':1e-10, 'msa_maxiter':4, 'msa_input_max_num':100})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter', 'evalue', 'msa_maxiter', 'msa_input_max_num']

    def constructor(self, recalculate):
        msa_input_handle = global_stuff.the_file_manager.get_file_handle(pdb_chain_msa_input_file_wrapper(self.params), recalculate)
        cline = MuscleCommandline(cmd = global_stuff.MUSCLE_PATH, input = msa_input_handle.name, out = self.get_file_location(), clw = False, maxiters = 2)
        cline()

class pdb_chain_processed_msa_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER, 'evalue':1e-10, 'msa_maxiter':4, 'msa_input_max_num':100})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter', 'evalue', 'msa_maxiter', 'msa_input_max_num']

    def constructor(self, recalculate):
        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_msa_file_wrapper(self.params), recalculate)
        msa = AlignIO.read(f.name, 'fasta')
        # search for the query sequence
        idx = -1
        for i in range(len(msa)):
            if msa[i].id == 'QUERY':
                idx = i
                break
        #print 'AAAAAAAAAAAAAAAAAAAAAAAAAAAA', idx
        #pdb.set_trace()
        # find the first non-insertion column
        i = 0
        while msa[idx,i] == '-':
            #print msa[idx,i]
            i = i + 1
            #print idx, i
        to_return = msa[:,i:(i+1)]
        print 'EEEEEEEEEEEEEEE'
        # add in all the other columns
        for k in range(i+1, msa.get_alignment_length()):
            if msa[idx,k] != '-':
                #print k
                to_return = to_return + msa[:,k:(k+1)]
        AlignIO.write(to_return, open(self.get_file_location(),'w'), 'fasta')

class pdb_chain_conservation_score_file_wrapper(file_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER, 'evalue':1e-10, 'msa_maxiter':4, 'msa_input_max_num':100, 'conservation_method':'shannon_entropy', 'use_gap_penalty':False})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter', 'evalue', 'msa_maxiter', 'msa_input_max_num', 'conservation_method', 'use_gap_penalty']

    def constructor(self, recalculate):
        msa_file_handle = global_stuff.the_file_manager.get_file_handle(pdb_chain_processed_msa_file_wrapper(self.params), recalculate)
        args = ['python', global_stuff.CONSERVATION_FOLDER+'score_conservation.py', '-m', global_stuff.CONSERVATION_FOLDER+'matrix/'+'blosum50.bla', '-o', self.get_file_location(), '-s', self.get_param('conservation_method'), '-p', str(self.get_param('use_gap_penalty')), msa_file_handle.name]
        subprocess.call(args)

class pdb_chain_conservation_score_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER, 'evalue':1e-10, 'msa_maxiter':4, 'msa_input_max_num':100, 'conservation_method':'shannon_entropy', 'use_gap_penalty':False})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter', 'evalue', 'msa_maxiter', 'msa_input_max_num', 'conservation_method', 'use_gap_penalty']

    def constructor(self, recalculate):
        #pdb.set_trace()
        f = global_stuff.the_file_manager.get_file_handle(pdb_chain_conservation_score_file_wrapper(self.params), recalculate)
        f.readline()
        f.readline()
        to_return = []
        for line in f:
            s = string.split(line, sep='\t')
            to_return.append(float(s[1]))
        return to_return
