import pdb
import wc
import new_new_objects as objects
import param
import global_stuff
import helper
import wrapper
import os
from wrapper_decorator import dec
import f as features
import constants


class generate_old_input_files(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate = False, to_pickle = False, to_filelize = False, always_recalculate = False):
        pdb_names = self.get_param(params, 'pdb_names')
        chain_letters = self.get_param(params, 'chain_letters')
        for i in range(len(pdb_names)):

            pdb_name = pdb_names[i]
            chain_letter = chain_letters[i]
            self.set_param(params, 'pdb_name', pdb_name)
            self.set_param(params, 'chain_letter', chain_letter)

            the_folder = constants.BIN_FOLDER + pdb_name + '/'

            if not os.path.exists(the_folder):
                os.makedirs(the_folder)
            pdb.set_trace()
            node_features = self.get_var_or_file(objects.jW, params, False, False, False)
            transposed_node_features = helper.get_transpose(node_features)
            pdb.set_trace()
            helper.write_mat(transposed_node_features, the_folder + 'Xnode.csv')
            pdb.set_trace()
            edge_features = self.get_var_or_file(objects.kW, params, False, False, False)
            edge_features_transposed = helper.get_transpose(edge_features)
            helper.write_mat(edge_features_transposed, the_folder + 'Xedge.csv')

            true_states = self.get_var_or_file(objects.oW, params, False, False, False)
            helper.write_vect(true_states, the_folder + 'true_y.csv')

            edge_list = self.get_var_or_file(objects.iW, params, False, False, False)
            helper.write_mat(edge_list, the_folder + 'edge_list.csv')

            info = [str(len(node_features)), str(len(edge_features))]
            helper.write_vect(info, the_folder + 'info.txt', the_sep = ' ')
            
            
        return None

#pdb_names = ['12as','2jcw','13pk','1a4i','1a4s','1ab8']
pdb_names = ['2jcw']
#chain_letters = ['A','A','A','A','A','A']
chain_letters = ['A']

from parameters import the_params

the_params.set_param('pdb_names', pdb_names)
the_params.set_param('chain_letters', chain_letters)

wc.get_stuff(generate_old_input_files,the_params, False, False, False, False)
