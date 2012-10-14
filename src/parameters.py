import f as features
import new_new_objects as objects
import param
#import try_svm

# [features.xW, features.akW]
 
#the_dict = {'data_list_file':'q.pl', 'edge_feature_list':[features.xW, features.akW], 'node_feature_list':[features.ayW, features.xW, features.vW, features.uW, features.wW, features.zW], 'dist_cut_off':5, 'pdb_name':'1b6b', 'chain_letter':'B', 'reg':10, 'mfmi':100, 'wif':0, 'wfld':0, 'nfld':2, 'ns':2, 'wob':0, 'evalue':1e-10, 'pos1':1, 'pos2':2, 'wreg':1, 'trun':99, 'wclf':objects.pW, 'svmC':10}

the_dict = {'hp':param.param(), 'd':'cw.pl', 'e':[features.xW, features.akW], 'n':[features.beW, features.bbW, features.ayW, features.xW, features.vW, features.uW, features.wW, features.zW], 'co':5, 'pdb_name':'12as', 'chain_letter':'A', 'reg':100, 'mfmi':20, 'wif':0, 'wfld':0, 'nfld':2, 'ns':2, "wob":0, 'wob2':2, 'evalue':1e-10, 'pos1':1, 'pos2':2, 'wreg':1, 'trun':99, 'wclf':objects.pW, 'svmC':10, 'lgn':9, 'lgc':5, 'pos':45, 'nwc':-1.0, 'micut':5, 'wtpr':0, 'posw':150, 'sfc':1, 'self':False, 'mx':50, 'ok':3, 'ik':2, 'md':1, 'tj':1, 'wj':0}


#the_dict = {'data_list_file':'catres.pdb_list', 'edge_feature_list':[], 'node_feature_list':[features.xW, features.vW, features.uW, features.wW, features.zW], 'dist_cut_off':5, 'pdb_name':'2jcw', 'chain_letter':'A', 'reg':0, 'mfmi':100, 'wif':2, 'wfld':0, 'nfld':2, 'ns':2, 'wob':0, 'evalue':1e-10, 'pos1':1, 'pos2':2, 'wreg':1, 'trun':99, 'wclf':objects.pW, 'svmC':10}

#asdf = features.xW


#the_dict = {'data_list_file':'catres_one.pdb_list', 'edge_feature_list':[features.xW], 'node_feature_list':[features.xW], 'dist_cut_off':5, 'pdb_name':'2jcw', 'chain_letter':'A', 'reg':10, 'mfmi':100, 'wif':1, 'wfld':0, 'nfld':2, 'ns':2, 'wob':0, 'evalue':1e-10, 'pos1':1, 'pos2':2, 'wreg':1, 'trun':20, 'wclf':try_svm.atW, 'svmC':10}

#the_dict = {'data_list_file':'catres_six.pdb_list', 'edge_feature_list':[], 'node_feature_list':[features.xW, features.vW, features.uW, features.wW, features.zW], 'dist_cut_off':5, 'pdb_name':'2jcw', 'chain_letter':'A', 'reg':10, 'mfmi':100, 'wif':2, 'wfld':0, 'nfld':2, 'ns':2, 'wob':2, 'evalue':1e-10, 'pos1':1, 'pos2':2, 'wreg':0}

#the_dict = {'data_list_file':'catres_six.pdb_list', 'edge_feature_list':[features.xW], 'node_feature_list':[features.xW], 'dist_cut_off':5, 'pdb_name':'1afw', 'chain_letter':'A', 'reg':10, 'mfmi':100, 'wif':1, 'wfld':0, 'nfld':2, 'ns':2, 'wob':2, 'wreg':0}

#the_dict = {'data_list_file':'catres_six.pdb_list', 'edge_feature_list':[features.xW], 'node_feature_list':[features.aaW], 'dist_cut_off':5, 'pdb_name':'12as', 'chain_letter':'A', 'reg':2, 'mfmi':15, 'obj':0, 'wfld':0, 'nfld':0, 'ns':2, 'wreg':0}

the_params = param.param(the_dict)


the_params.set_param('hpv',param.param({'reg':[0,10,50,100,400,900]}))
