import new_features as features
import new_new_objects as objects
import param
import wc
import global_stuff

the_dict = {'pdb_name':'1asy', 'chain_letter':'A', 'edge_feature_list':[features.ones_fxn_w], 'node_feature_list':[features.residue_class_categorical_fxn_w, features.residue_categorical_fxn_w, features.inverse_avg_dist_fxn_w], 'dist_cut_off':5}

the_params = param.param(the_dict)

wc.get_stuff(objects.node_features_obj_w, the_params, global_stuff.recalculate, True, True)
