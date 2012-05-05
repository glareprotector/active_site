# will, for each protein chain, generate the input needed for training

from objects import *
from crf import *
import global_stuff
import constants
import pdb
import string

#pdb_chains = [['112l','A'],['119l','A']]
pdb_chains = []

f = open(global_stuff.ORIG_CHAINS, 'r')
for line in f:
    pdb_chains.append(string.split(string.strip(line), '_'))

print pdb_chains

success_f = open(global_stuff.success_file, 'w')
fail_f = open(global_stuff.fail_file, 'w')

all_node_features_file = 'all_node_features.csv'
true_class_file = 'true_class'

for pdb_chain in pdb_chains:
    temp_params = param({'pdb_name':pdb_chain[0], 'chain_letter':pdb_chain[1], 'dist_cut_off':8, 'evalue':1e-4})
    chain_str = pdb_chain[0] + '_' + pdb_chain[1]
    print 'CURRENT CHAIN: ', chain_str
    try:
    
        temp_crf = global_stuff.the_obj_manager.get_variable(crf_obj_wrapper(temp_params), recalculate=False, use_pickle = False)
        success_f.write(chain_str + '\n')
        temp_folder = constants.BIN_FOLDER + chain_str + '/'
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)
        temp_crf.write_info(temp_folder)

        #crf_node_features = global_stuff.get_transpose(temp_crf.Xnode)
        #global_stuff.write_mat(crf_node_features, all_node_features_file, ',', 'a')
        #global_stuff.write_mat([[temp_crf.true_y[i]] for i in range(len(temp_crf.true_y))], true_class_file,',', 'a')

        
    except Exception as e:
        print e
        fail_f.write(chain_str + '\n')
        fail_f.write(str(e))

    success_f.flush()
    fail_f.flush()


    
    
