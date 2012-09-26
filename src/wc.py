import wrapper

import f as features
import param
import pdb
import new_new_objects as objects
import global_stuff

#pdb.set_trace()
useless = wrapper.famished_wrapper()
wc = wrapper.wrapper_catalog(useless, param.param({}))

def get_stuff(wrapper_class, params, recalculate, to_pickle, to_filelize, always_recalculate = False):
    params.set_param('which_wrapper_class', wrapper_class)
    pdb.set_trace()
    wc_used_keys, wc_all_keys, wrapper_instance = wc.constructor(params, True, False, False)

    try:
        stuff_used_keys, stuff_all_keys, stuff = wrapper_instance.constructor(params, recalculate, to_pickle, to_filelize, always_recalculate = always_recalculate)
    except Exception, err:
        print 'ERROR when calling get_stuff with this error', err, params
    return stuff


