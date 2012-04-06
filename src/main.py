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

for pdb_chain in pdb_chains:
    temp_params = param({'pdb_name':pdb_chain[0], 'chain_letter':pdb_chain[1], 'dist_cut_off':8})
    temp_crf = global_stuff.the_obj_manager.get_variable(crf_obj_wrapper(temp_params), recalculate=False, use_pickle = False)
    temp_folder = constants.BIN_FOLDER + pdb_chain[0] + '_' + pdb_chain[1] + '/'
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    #pdb.set_trace()
    temp_crf.write_info(temp_folder)
    
