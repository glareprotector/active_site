import constants
import pdb
from param import *
import global_stuff

def dec(f):

    def affects_which_wrappers_are_called(self, key, val):
        if key == 'md':
            
            return True
        if val.__class__.__name__ in ['int', 'float', 'str']:
            return False
        import cross_validation_pseudo as cv
        if type(val) == type(cv.data) or type(val) == cv.fold:
            return False
        elif type(val) == param:
            return False
        else:
            return True

    # a key/value should go into the key only if the value determines which wrappers will be called
    # assumptions: if a param does not affect which wrappers are called at this node, then it does not affect which wrappers are called
    # at any descendant node.  as long as the param is used in the same way in every node, this will be the case
    #@print_stuff_dec
    """
    def get_all_keys_key(self, params, used_keys, set_keys):
        temp = {}
        for key in used_keys:
            if key not in set_keys:
                if affects_which_wrappers_are_called(self, key, params.get_param(key)):
                    temp[key] = params.get_param(key)
        return param(temp)
    """

    def get_all_keys_key_key_set(self, params, used_keys, set_keys):
        ans = set()
        for key in used_keys.union(self.temp_dependents_all_keys_key_key_sets[-1]):
            if key not in set_keys:
                if affects_which_wrappers_are_called(self, key, params.get_param(key)):
                    ans.add(key)
        #return ans.union(self.temp_dependents_all_keys_key_key_sets[-1])
        return ans


    def get_all_keys(self, params, recalculate):
        A = self.temp_used_keys[-1]
        B = self.temp_dependents_keys[-1]
        C = self.temp_new_param_keys[-1]
        return ((A | B) - C) | (A & B & C)
        #return (self.temp_used_keys[-1].union(self.temp_dependents_keys[-1]) - self.temp_new_param_keys[-1]) + 


    def get_all_keys_key(self, params, all_keys_key_key_set):
        return param.restriction(params, all_keys_key_key_set)

    def get_objects_key(self, params, all_keys):
        return param.restriction(params, all_keys)

    def cache_everything_f_poster(self, params, recalculate, to_pickle, to_filelize, always_recalculate, old_obj = None):
#        pdb.set_trace()
        # add empty sets to the used_keys and dependents keys and set params

        self.temp_used_keys.append(set())
        self.temp_dependents_keys.append(set())
        self.temp_new_param_keys.append(set())
        self.temp_dependents_all_keys_key_key_sets.append(set())
        # added maker (pointer to self) to params, in case it is needed when making my own classes
 #       print '           starting to calculate NEW stuff in wrapper ', self, recalculate, to_pickle, to_filelize, always_recalculate
        try:
            if old_obj == None:
                obj = f(self, params, recalculate, to_pickle, to_filelize, always_recalculate)
            else:
                obj = f(self, params, recalculate, to_pickle, to_filelize, always_recalculate, old_obj)
        except Exception, err:
            print 'ERROR DOWNSTREAM OF CALLING ', self, err
            print params
            self.temp_used_keys.pop()
            self.temp_dependents_keys.pop()
            self.temp_new_param_keys.pop()
            self.temp_dependents_all_keys_key_key_sets.pop()

            raise
#        print '        finished to calculate NEW stuff in wrapper ', self, recalculate, to_pickle, to_filelize, always_recalculate
#        print '                                                             ', len(self.temp_used_keys)
        # always pickle used_keys and all_keys
        self.used_keys_cache.set(self.temp_used_keys[-1], True)
        self.set_keys_cache.set(self.temp_new_param_keys[-1], True)
        try:
            asdf = 1
            #assert len(self.temp_used_keys) > 0
            #assert len(self.temp_used_keys) <= 1
        except:
            print 'TOO MANY THINGS?FEW ', self, len(self.temp_used_keys)

        #pdb.set_trace()
        # split below step into pieces - figuring out what keys go into the key, and then getting actual value of key
        all_keys_key_key_set = get_all_keys_key_key_set(self, params, self.used_keys_cache.get(recalculate), self.set_keys_cache.get(recalculate))
        self.all_keys_key_key_set_cache.set(all_keys_key_key_set, True)
        all_keys_key = get_all_keys_key(self, params, all_keys_key_key_set)
        #all_keys_key = get_all_keys_key(self, params, self.used_keys_cache.get(recalculate), self.set_keys_cache.get(recalculate))
        all_keys = get_all_keys(self, params, recalculate)
        # pop the used_keys and dependents keys and set params
        object_key = get_objects_key(self, params, all_keys)
        # wrapper decides whether to make index.  index provided to object constructor in params.  object can decide whether to do anything with it.
        #pdb.set_trace()
        if self.makes_index():
            # this index may already be in the index, or may not
            # NEED TO INCLUDE TO_REINDEX
            #pdb.set_trace()
            index = self.object_key_to_index.get_and_set_index(object_key, global_stuff.to_reindex)
            # self.object_key_to_index.set(object_key, index, to_reindex)
            if hasattr(obj, 'process_index'):
                obj.process_index(index)
            object_key = index
        self.all_keys_cache.set(all_keys_key, all_keys, True, recalculate)
        # set always recalculate to false here? if call dumper wrapper directly can set it to true
        obj = self.cache.set(object_key, obj, to_pickle, params, to_filelize, always_recalculate)

        self.temp_used_keys.pop()
        self.temp_dependents_keys.pop()
        self.temp_new_param_keys.pop()
        self.temp_dependents_all_keys_key_key_sets.pop()
        

                    
        return self.used_keys_cache.get(recalculate), all_keys, obj, all_keys_key_key_set
    
    def h(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        # make a copy of params
        
        params = params.get_copy()
#        print self, recalculate, params
        if self.used_keys_cache.has(recalculate) and self.set_keys_cache.has(recalculate) and self.all_keys_key_key_set_cache.has(recalculate):
            used_keys = self.used_keys_cache.get(recalculate)
            set_keys = self.set_keys_cache.get(recalculate)
            all_keys_key_key_set = self.all_keys_key_key_set_cache.get(recalculate)
            all_keys_key = get_all_keys_key(self, params, all_keys_key_key_set)
            # so i don't have to do any processing that uses the used_keys and set_keys to get all_keys_key_keyset.  used_keys is still useful for getting the keys if i don't have all_keys_key_key_set
            #all_keys_key = get_all_keys_key(self, params, used_keys, set_keys)
            if self.all_keys_cache.has(all_keys_key, recalculate):
                all_keys = self.all_keys_cache.get(all_keys_key)
                object_key = get_objects_key(self, params, all_keys)
                # if are making index, object is stored under index(object_key) if it exists
                if self.makes_index():
                    if self.object_key_to_index.has(object_key, global_stuff.to_reindex):
                        object_index = self.object_key_to_index.get(object_key)
                        if self.cache.has(object_index, recalculate):
                            #print '           getting cached value in wrapper ', self
                            obj = self.cache.get(object_index, recalculate)
                            if always_recalculate:
                                return cache_everything_f_poster(self, params, recalculate, to_pickle, to_filelize, always_recalculate, obj)
                            return used_keys, all_keys, obj, all_keys_key_key_set
                        
                if self.cache.has(object_key, recalculate):
                    #print '           getting cached value in wrapper ', self
                    obj = self.cache.get(object_key, recalculate)
                    if always_recalculate:
                        return cache_everything_f_poster(self, params, recalculate, to_pickle, to_filelize, always_recalculate, obj)
                    return used_keys, all_keys, obj, all_keys_key_key_set
        return cache_everything_f_poster(self, params, recalculate, to_pickle, to_filelize, always_recalculate)

    return h
