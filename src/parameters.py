import f as features
import param
the_dict = {'data_list_file':'catres.pdb_list', 'edge_feature_list':[features.xW], 'node_feature_list':[features.vW, features.uW, features.wW, features.zW], 'dist_cut_off':5, 'pdb_name':'1afw', 'chain_letter':'A', 'reg':10, 'mfmi':100, 'wif':1, 'wfld':0, 'nfld':2, 'ns':2, 'wob':2}

#the_dict = {'data_list_file':'catres_six.pdb_list', 'edge_feature_list':[features.xW], 'node_feature_list':[features.aaW], 'dist_cut_off':5, 'pdb_name':'12as', 'chain_letter':'A', 'reg':2, 'mfmi':15, 'obj':0, 'wfld':0, 'nfld':0, 'ns':2}

the_params = param.param(the_dict)
