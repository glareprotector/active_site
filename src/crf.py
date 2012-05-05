"""
created 3-26-2012
CRF will consist of potentials(node_id's, list of feature functions), values of feature/class associated weights, values of node_id associated variables.  potentials are really just a functions(the node_id's associated with function determine its C, and the values of node_id associated variables determine the values of y_C).  potentials can evaluate themselves.  also need a way to evaluate.  need to create the CRF by taking in a pdb chain
"""

import constants
import global_stuff
from objects import *
from cache import cache
import operator
import random
from features import *
import pdb
import math
import random

class node(object):

    # params should have pdb_name, chain_letter, pos
    def __init__(self, params):
        self.params = params

class edge(object):

    def __init__(self, res1, res2):
        self.res1 = res1
        self.res2 = res2

# needs to be able to generate for a single example: edgeStruct, nodeMap, edgeMap, Xnode, Xedge
class crf(object):

    def get_feature_function_wrappers(self):
        feature_wrappers = [ones_function_wrapper(), inverse_average_distance_feature_function_wrapper(), conservation_feature_function_wrapper(), b_factor_feature_function_wrapper()]
        #feature_wrappers = [b_factor_feature_function_wrapper()]
        # add categorical variable for the nucleotide
        temp_params = param({'g_wrapper':get_residue_function_wrapper(), 'indicator_values_list':[ [x] for x in constants.AMINO_ACID_LIST]})
        res_cat = global_stuff.the_obj_manager.get_variable(categorical_function_wrapper(temp_params), self.recalculate)
        feature_wrappers.extend(res_cat)
        # add feature for which of charged, polar, hydrophobic the amino acid belongs to
        temp_params = param({'g_wrapper':get_residue_function_wrapper(), 'indicator_values_list':constants.AMINO_ACID_CATEGORIES})
        res_type_cat = global_stuff.the_obj_manager.get_variable(categorical_function_wrapper(temp_params), self.recalculate)
        feature_wrappers.extend(res_type_cat)
        return feature_wrappers
        

    def get_res_features(self, res):
        features = []
        for f_wrapper in self.get_feature_function_wrappers():
            f = global_stuff.the_obj_manager.get_variable(f_wrapper, self.recalculate, to_pickle=False, use_pickle=False)
            features.append(f(res.params))
        return features

    # returns the node_map for a single node, edge_map for a single edge
    def get_maps(self):
        k = 1 # matlab indices are 1-based
        num_states = len(self.y_range)
        num_features = len(self.get_feature_function_wrappers())
        #pdb.set_trace()
        node_map = [ [-1 for i in range(num_features)] for j in range(num_states)]
        for s in range(num_states):
            for f in range(num_features):
                node_map[s][f] = k
                k = k + 1

        # for now, only 1 feature per edge, so edge_map is 2d
        edge_map = [[-1 for i in range(num_states)] for j in range(num_states)]
        for s1 in range(num_states):
            for s2 in range(num_states):
                #print s1
                #print s2
                edge_map[s1][s2] = k
                k = k + 1
        
        return node_map, edge_map, k-1

    # returns a num_feature by num_nodes matrix of the node features
    # this will be replicated for each example in the matlab code
    def get_Xnode(self):
        Xnode = []
        num_nodes = len(self.residues)
        for i in range(num_nodes):
            Xnode.append(self.get_res_features(self.residues[i]))
        return global_stuff.get_transpose(Xnode)

    # returns a (for now) num_feature by num_nodes matrix of edge features
    # for now this 1 by num_edges matrix is all 1's. 
    def get_Xedge(self):
        num_edges = len(self.edges)
        Xedge = [[1 for i in range(num_edges)]]
        return Xedge

    # reads data file to get the true label for each position
    def get_true_y(self):
        true_y = [1 for i in range(len(self.residues))]


        aa_to_pos = global_stuff.the_obj_manager.get_variable(pdb_chain_aa_to_pos_obj_wrapper(self.params), self.recalculate)
        g = global_stuff.the_obj_manager.get_variable(b_factor_feature_function_wrapper(self.params), self.recalculate)

        #for i in range(len(true_y)):
        #    print i
        #    temp_params = self.params.get_copy()
        #    temp_params.set_param('pos', aa_to_pos[i])
        #    bf = g(temp_params)
        #    if random.uniform(0,100) > bf:
        #        true_y[i] = 2
        #return true_y




        m = global_stuff.the_obj_manager.get_variable(pdb_chain_pos_to_aa_dict_obj_wrapper(self.params), self.recalculate)
        f = open(global_stuff.CSA_FILE,'r')
        f.readline()
        #print m.keys()
        for line in f:
            s = string.split(line, sep=',')
            pos = int(s[2])
            if string.upper(s[0]) == string.upper(self.params.get_param('pdb_name')) and string.upper(s[1]) == string.upper(self.params.get_param('chain_letter')):
                #try:
                #    if pos not in m.keys():
                #    pdb.set_trace()

                true_y[m[pos]] = 2
                #except Exception as e:
                #    print 'ERROR: map key out of bounds', pos, m.keys()



        return true_y

    # writes to file all the info needed for training
    def write_info(self, folder):
        global_stuff.write_mat(self.node_map, folder + 'node_map.csv', ',')
        global_stuff.write_mat(self.edge_map, folder + 'edge_map.csv', ',')
        global_stuff.write_mat([self.true_y], folder + 'true_y.csv', ',')
        global_stuff.write_mat(self.Xnode, folder + 'Xnode.csv', ',')
        global_stuff.write_mat(self.Xedge, folder + 'Xedge.csv', ',')
        global_stuff.write_mat(self.adj_mat, folder + 'adj_mat.csv', ',')
        info = [self.numNodes, self.numEdges, self.numStates, self.numFeatures, self.numParams]
        global_stuff.write_mat([info], folder + 'info.txt', ' ')

    # params contains 'pdb_name' and 'chain_letter' and 'dist_cut_off'
    def __init__(self, params, recalculate):
        
        self.params = params
        self.recalculate = recalculate

        # if chain letter is not specified in params, figure out what it is
        if self.params.get_param('chain_letter') == '-1':
            structure = global_stuff.the_obj_manager.get_variable(pdb_obj_wrapper(self.params))
            self.params.set_param('chain_letter',structure[0].child_dict.keys()[0])



        
        res_dists = global_stuff.the_obj_manager.get_variable(pdb_chain_pairwise_distance_obj_wrapper(self.params), self.recalculate)
        aa_to_pos = global_stuff.the_obj_manager.get_variable(pdb_chain_aa_to_pos_obj_wrapper(self.params), self.recalculate)
        chain_seq = global_stuff.the_obj_manager.get_variable(pdb_chain_seq_obj_wrapper(self.params), self.recalculate)

        self.residues = []
        self.edges = []
        self.adj_mat = [ [0 for i in range(len(chain_seq))] for j in range(len(chain_seq))]
        self.y_range = [0,1]


        # add all residues to residues
        for i in range(len(chain_seq)):
            res_param = self.params.get_copy()
            res_param.set_param('pos', aa_to_pos[i])
            self.residues.append(node(res_param))
        
        # add residues within dist_cut_off to edges.  also create adj_mat matlab will use for edge_struct
        for i in range(len(chain_seq)):
            for j in range(i):
                if math.fabs(i-j) == 1:
                    self.edges.append(edge(self.residues[i], self.residues[j]))
                    self.adj_mat[i][j] = 1
                    self.adj_mat[j][i] = 1

                elif res_dists[i][j] < self.params.get_param('dist_cut_off') and math.fabs(i-j) > 5:
                    self.edges.append(edge(self.residues[i], self.residues[j]))
                    self.adj_mat[i][j] = 1
                    self.adj_mat[j][i] = 1
        


        # need map from position in residues to true class.  can do this by reading in file, but for now assign randomly
        self.node_map, self.edge_map, self.numParams = self.get_maps()
        self.Xnode = self.get_Xnode()
        self.Xedge = self.get_Xedge()
        self.true_y = self.get_true_y()

        self.numNodes = len(self.residues)
        self.numEdges = len(self.edges)
        self.numStates = len(self.y_range)
        self.numFeatures = len(self.get_feature_function_wrappers())

class crf_obj_wrapper(obj_wrapper):

    def get_hard_coded_params(self):
        return param({'location':constants.BIN_FOLDER})

    def get_self_param_keys(self):
        return['pdb_name', 'chain_letter', 'dist_cut_off']

    def constructor(self, recalculate):
        return crf(self.params, recalculate)
