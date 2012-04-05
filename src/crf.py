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
        return [conservation_feature_function_wrapper(), inverse_average_distances_feature_function_wrapper()]

    def get_res_features(self, res):
        features = []
        for f_wrapper in self.get_feature_function_wrappers():
            f = global_stuff.the_obj_manager.get_variable(f_wrapper)
            features.append(f(res.params))

    # returns the node_map for a single node, edge_map for a single edge
    def get_maps(self):
        k = 1 # matlab indices are 1-based
        num_states = len(self.y_range)
        num_features = len(self.get_feature_function_wrapper())
        node_map = [ [-1 for i in range(num_features)] for j in range(num_states)]
        for s in range(num_states):
            for f in range(num_features):
                node_map[s][f] = k
                k = k + 1

        # for now, only 1 feature per edge, so edge_map is 2d
        edge_map = [[-1 for i in range(num_states)] for j in range(num_states)]
        for s1 in range(num_states):
            for s2 in range(num_states):
                edge_map[s1, s2] = k
                k = k + 1
        
        return node_map, edge_map

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
        Xedge = [1 for i in range(num_edges)]
        return Xedge

    # reads data file to get the true label for each position
    def get_true_y(self):
        true_y = [0 for i in range(len(self.residues))]
        for i in range(len(self.residues)):
            true_y = random.randint(0,1)

    # writes to file all the info needed for training
    def write_info(self, folder):
        global_stuff.write_mat(self.node_map, folder + 'node_map.csv', ',')
        global_stuff.write_mat(self.edge_map, folder + 'edge_map.csv', ',')
        global_stuff.write_mat([self.true_y], folder + 'true_y.csv', ',')
        global_stuff.write_mat(self.Xnode, folder + 'Xnode.csv', ',')
        global_stuff.write_mat(self.Xedge, folder + 'Xedge.csv', ',')
        global_stuff.write_mat(self.adj_mat, folder + 'adj_mat.csv', ',')

    # params contains 'pdb_name' and 'chain_letter' and 'dist_cut_off'
    def __init__(self, params, recalculate):
        
        self.residues = []
        self.edges = []
        self.adj_mat = [ [0 for i in range(len(chain_seq))] for j in range(len(chain_seq))]
        self.params = params
        self.y_range = [0,1]
        
        res_dists = global_stuff.the_obj_manager.get_variable(pdb_chain_pairwise_distance_obj_wrapper(self.params, recalculate))
        aa_to_pos = global_stuff.the_obj_manager.get_variable(pdb_chain_aa_to_pos_obj_wrapper(self.params, recalculate))
        chain_seq = global_stuff.the_obj_manager.get_variable(pdb_chain_seq_obj_wrapper(self,params, recalculate))

        # add all residues to resides
        for i in range(len(chain_seq)):
            res_param = global_stuff.deep_copy(self.params)
            res_param['pos'] = aa_to_pos[i]
            self.residues.append(node(res_param))
        
        # add residues within dist_cut_off to edges.  also create adj_mat matlab will use for edge_struct
        for i in range(len(chain_seq)):
            for j in range(i):
                if res_dists[i][j] < self.params['dist_cut_off']:
                    self.edges.append(edge(self.residues[i], self.residues[j]))
                    adj_mat[i][j] = 1
                    adj_mat[j][i] = 1

        # need map from position in residues to true class.  can do this by reading in file, but for now assign randomly
        self.true_y = self.get_true_y()
        self.node_map, self.edge_map = self.get_maps()
        self.Xnode = self.get_Xnode()
        self.Xedge = self.get_Xedge()

# reads in list of pdb chains, and makes crf, writes down info for them
