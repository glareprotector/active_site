import pdb
import f as features
import new_new_objects as objects
import param
import wc
import global_stuff

# hardcode parameters for the experiment here for now.

#the_dict = {'pdb_name':'1asy', 'chain_letter':'A', 'edge_feature_list':[features.xW], 'node_feature_list':[features.vW, features.uW, features.wW], 'dist_cut_off':5}

#the_params = param.param(the_dict)


import helper

file_location = 'mf_nodewise_0'

folder_name, the_params = helper.read_param(file_location)

the_params.set_param('p', '1p3d')
the_params.set_param('c', 'A')
the_params.set_param('st', 322)
the_params.set_param('en', 473)


ans = wc.get_stuff(objects.ciW, the_params, False, False, False)

print ans
