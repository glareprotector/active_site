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
        return self.temp_used_keys[-1].union(self.temp_dependents_keys[-1]) - self.temp_new_param_keys[-1]

    def get_objects_key(self, params, all_keys):
        return param.restriction(params, all_keys)

    def cache_everything_f_poster(self, params, recalculate, to_pickle, to_filelize):
        pdb.set_trace()
        # add empty sets to the used_keys and dependents keys and set params
        self.temp_used_keys.append(set())
        self.temp_dependents_keys.append(set())
        self.temp_new_param_keys.append(set())
        object = f(self, params, recalculate, to_pickle, to_filelize)
        # always pickle used_keys and all_keys
        self.used_keys_cache.set(self.temp_used_keys, True)
        all_keys_key = get_all_keys_key(self, params, self.used_keys_cache.get(recalculate))
        all_keys = get_all_keys(self, params, recalculate)
        self.temp_dependents_keys = set()
        self.all_keys_cache.set(all_keys_key, all_keys, True)
        object_key = get_objects_key(self, params, all_keys)
        # wrapper decides whether to make index.  index provided to object constructor in params.  object can decide whether to do anything with it.
        if self.makes_index():
            # this index may already be in the index, or may not
            # NEED TO INCLUDE TO_REINDEX
            index = self.object_key_to_index.get_and_set_index(object_key, to_reindex)
            # self.object_key_to_index.set(object_key, index, to_reindex)
            if hasattr(object, 'process_index'):
                object.process_index(index)
        object = self.cache.set(object_key, object, to_pickle, params, to_filelize)
        return self.used_keys_cache.get(recalculate), all_keys, object
    
    def h(self, params, recalculate, to_pickle, to_filelize = False, to_reindex = False):
#        pdb.set_trace()
        if self.used_keys_cache.has(recalculate):
            used_keys = self.used_keys_cache.get(recalculate)
            all_keys_key = get_all_keys_key(self, params, used_keys)
            if self.all_keys_cache.has(all_keys_key, recalculate):
                all_keys = self.all_keys_cache.get(all_keys_key)
                object_key = get_objects_key(self, params, all_keys)
                # if are making index, object is stored under index(object_key) if it exists
                if self.makes_index():
                    if self.object_key_to_index.has(object_key, to_reindex):
                        object_index = self.object_key_to_index.get(object_key)
                        if self.cache.has(object_index, recalculate):
                            object = self.cache.get(object_index, recalculate)
                            return used_keys, all_keys, object
                        
                if self.cache.has(object_key, recalculate):
                    object = self.cache.get(object_key, recalculate)
                    return used_keys, all_keys, object
        return cache_everything_f_poster(self, params, recalculate, to_pickle, to_filelize)

    return h
