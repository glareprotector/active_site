'''
Created on Mar 18, 2012

@author: glareprotector
'''
INFO_FOLDER = '../info_folder/'
ROC_CURVE_SCRIPT = './get_roc_from_python.r'
ROC_INFO_SCRIPT = './get_roc_text.r'
BIN_FOLDER = '../test8/'
HOLDING_FOLDER = '../holding_folder/'
AUX_FOLDER = '../data/catres-fam/'
EVALUE = 1e-10
AMINO_ACID_LIST = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','Y']
AMINO_ACID_CATEGORIES = [ ['D','E','H','K','R'], ['Q','T','S','N'], ['A','F','G','I','L','M','P','V','W']]
DSSP_CATEGORIES = ['H','B','E','G','I','T','S']
HP_VALUES_FOLDER = INFO_FOLDER

CACHE_MAX_SIZE = 100




def get_master_node_feature_list():
    import f
    return [f.beW, f.bbW, f.ayW, f.xW, f.vW, f.uW, f.wW, f.zW]

def get_master_edge_feature_list():
    import f
    return [f.xW, f.akW]
