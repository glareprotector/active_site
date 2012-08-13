import constants
import pdb
from param import *

def dec(f):

    def affects_which_wrappers_are_called(key, val):
        if val.__class__.__name__ in ['int', 'float', 'str']:
            return False
        else:
            return True

    # a key/value should go into the key only if the value determines which wrappers will be called
    # assumptions: if a param does not affect which wrappers are called at this node, then it does not affect which wrappers are called
    # at any descendant node.  as long as the param is used in the same way in every node, this will be the case
    #@print_stuff_dec
    def get_all_keys_key(self, params, used_keys):
        temp = {}
        for key in used_keys:
            if affects_which_wrappers_are_called(key, params.get_param(key)):
                temp[key] = params.get_param(key)
        return param(temp)

    def get_all_keys(self, params, recalculate):
        return self.used_keys_cache.get(recalculate).union(self.temp_dependents_keys) - self.temp_new_param_keys

    def get_objects_key(self, params, all_keys):
        return param.restriction(params, all_keys)

    def cache_everything_f_poster(self, params, recalculate, to_pickle, to_filelize):
        object = f(self, params, recalculate, to_pickle, to_filelize)
        # always pickle used_keys and all_keys
        self.used_keys_cache.set(self.temp_used_keys, True)
        all_keys_key = get_all_keys_key(self, params, self.used_keys_cache.get(recalculate))
        all_keys = get_all_keys(self, params, recalculate)
        self.temp_dependents_keys = set()
        self.all_keys_cache.set(all_keys_key, all_keys, True)
        object_key = get_objects_key(self, params, all_keys)
        object = self.cache.set(object_key, object, to_pickle, params, to_filelize)
        return self.used_keys_cache.get(recalculate), all_keys, object
    
    def h(self, params, recalculate, to_pickle, to_filelize = False):
#        pdb.set_trace()
        if self.used_keys_cache.has(recalculate):
            used_keys = self.used_keys_cache.get(recalculate)
            all_keys_key = get_all_keys_key(self, params, used_keys)
            if self.all_keys_cache.has(all_keys_key, recalculate):
                all_keys = self.all_keys_cache.get(all_keys_key)
                object_key = get_objects_key(self, params, all_keys)
                if self.cache.has(object_key, recalculate):
                    object = self.cache.get(object_key, recalculate)
                    return used_keys, all_keys, object
        return cache_everything_f_poster(self, params, recalculate, to_pickle, to_filelize)

    return h
