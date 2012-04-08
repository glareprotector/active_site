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

for pdb_chain in pdb_chains:
    temp_params = param({'pdb_name':pdb_chain[0], 'chain_letter':pdb_chain[1], 'dist_cut_off':8, 'evalue':1e-4})
    chain_str = pdb_chain[0] + '_' + pdb_chain[1]
    try:
    
        temp_crf = global_stuff.the_obj_manager.get_variable(crf_obj_wrapper(temp_params), recalculate=False, use_pickle = False)
        success_f.write(chain_str + '\n')
        temp_folder = constants.BIN_FOLDER + chain_str + '/'
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)
        temp_crf.write_info(temp_folder)
    
    except Exception as e:
        print e
        fail_f.write(chain_str + '\n')
        fail_f.write(str(e))

    success_f.flush()
    fail_f.flush()


    
    
