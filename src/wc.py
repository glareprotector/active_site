import wrapper

import new_features as features
import param
import pdb
import new_new_objects as objects

pdb.set_trace()
useless = wrapper.famished_wrapper()
wc = wrapper.wrapper_catalog(useless, param.param({}))

def get_stuff(wrapper_class, params, recalculate, to_pickle, to_filelize):
    #pdb.set_trace()
    params.set_param('which_wrapper_class', wrapper_class)
    wc_used_keys, wc_all_keys, wrapper_instance = wc.constructor(params, True, False, False)
    stuff_used_keys, stuff_all_keys, stuff = wrapper_instance.constructor(params, recalculate, to_pickle, to_filelize)
    return stuff

import tester
