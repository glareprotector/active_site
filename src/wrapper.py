from param import param

import pdb, os, subprocess, constants, pickle
from global_stuff import print_stuff_dec, write_mat, write_vect


import caches
from wrapper_decorator import *
#from wc import wc

# since wrappers and caches are objects that can be cached, should record the wrapper that is constructing them.

# note for all caches: recalculate specifies whether you want to use the files/pickles that are present when cache is created

# act of caching: moving file from temporary holding spot to proper folder.  cache structure: file system

# this is the function wrapper

class wrapper(object):

    def is_indexed(self):
        return True

    def makes_index(self):
        return False

    def __repr__(self):
        return self.__class__.__name__

    def get_holding_location(self):
        return constants.BIN_FOLDER + str(id(self))

    def get_folder(self, object_key):
        return constants.BIN_FOLDER


    # if the wrapper index stuff, its objects can obtain indicies thru maker->object_key_to_index
    # in general, the methods i make should have a maker pointer so that params can be accessed thru the maker
    def get_name(self, object_key, to_reindex = False):
        if not to_reindex:
            return self.__repr__() + str(object_key)
        else:
            assert self.maker.object_key_to_index.has(object_key, to_reindex) == True
            return self.__repr__() + str(self.maker.object_key_to_index.get(object_key))

    # can't get and set same parameter in the same run of constructor.  at the beginning of constructor, set maker.  in init of whatever object, would fetch maker.  can NOT get, then set.  however, can set and THEN get.
    def basic_init(self, maker, params):
        self.temp_used_keys = []
        self.temp_dependents_keys = []
        self.temp_new_param_keys = []
        maker.set_param(params, "source_instance", self, self)
        self.used_keys_cache = caches.ukcO(maker, params)
        self.set_keys_cache = caches.skcO(maker, params)
        self.all_keys_cache = caches.akcO(maker, params)
        if self.makes_index():
            self.object_key_to_index = caches.index_cache(maker, params)
        self.maker = maker

    def before_init(self, maker, params):
        maker.set_param(params, "dumper_source_instance", self, self)
        return maker, params
    
# if set param, means its value was based on other params, so only other params matter.
    def __init__(self, maker, params = param({})):
#        pdb.set_trace()
        maker, params = self.before_init(maker, params)
        self.other_init(maker, params)
        self.basic_init(maker, params)

    def other_init(self, maker, params):
        pass

    def set_param(self, params, key, val, old_self = None):
        if old_self != self:
            self.temp_new_param_keys[-1].add(key)
        params.set_param(key, val)
        return params

    def get_param(self, params, key, record = True, old_self = None):
        if record and old_self != self:
            self.temp_used_keys[-1].add(key)
        return params.get_param(key)

    # wrapper refers to which class, not a specific wrapper instance
    @print_stuff_dec
    def get_var_or_file(self, wrapper, params, recalculate, to_pickle, to_filelize = False, old_self = None):
        print '                      called  ', self, wrapper, recalculate, to_pickle, to_filelize
        #pdb.set_trace()
        self.set_param(params, "which_wrapper_class", wrapper)
        from wc import wc
        the_wrapper = self.old_get_var_or_file(wc, params, True, False, False, old_self)
        used_keys, all_keys, x = the_wrapper.constructor(params, recalculate, to_pickle, to_filelize)
        self.temp_dependents_keys[-1]  = self.temp_dependents_keys[-1].union(all_keys)
        print '                      returned ', self, wrapper, recalculate, to_pickle, to_filelize
        return x

    def old_get_var_or_file(self, the_wrapper, params, recalculate, to_pickle, to_filelize = False, old_self = None):

        used_keys, all_keys, x = the_wrapper.constructor(params, recalculate, to_pickle, to_filelize)
        if old_self != self:
            self.temp_dependents_keys[-1]  = self.temp_dependents_keys[-1].union(all_keys)
        return x




    # returns only the object, but after decorating, will return used_keys, all_keys, object
    @dec
    @print_stuff_dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        pass

    def get_cache(self):
        pass

    def has(self, object_key, recalculate):
        return self.cache.has(object_key, recalculate)

    def get(self, object_key, recalculate):
        return self.cache.get(object_key, recalculate)

    def set(self, object_key, object, to_pickle, params):
        return self.cache.set(object_key, object, to_pickle, params)

class by_pdb_folder_wrapper(wrapper):

    def get_folder(self, object_key):
        return constants.BIN_FOLDER + '/' + object_key.get_param('pdb_name') + '/'

class indexing_wrapper(wrapper):

    def makes_index(self):
        return True

class indexed_wrapper(wrapper):

    def process_index(self, index):
        self.name = str(index)

class obj_wrapper(wrapper):

    def get_file_dumper(self, maker, params):
        return None

    def get_file_location(self, object_key):
        return self.get_folder(object_key) + self.get_name(object_key) + '.pk'

    def other_init(self, maker, params):
        maker.set_param(params, "source_instance", self)
        self.cache = caches.object_cache_for_wrapper(maker, params)

class mat_obj_wrapper(obj_wrapper):
    
    def get_file_dumper(self, maker, params):
        return mfdW(maker, params)
        maker.set_param(params, "which_wrapper_class", wrapper.mfdW)
        return maker.get_var_or_file(wrapper_catalog, params, True, False, False)

class vect_obj_wrapper(obj_wrapper):

    def get_file_dumper(self, maker, params):
        return generic_vect_file_dumper_wrapper(maker, params)
        maker.set_param(params, "which_wrapper_class", wrapper.generic_vect_file_dumper_wrapper)
        return maker.get_var_or_file(wrapper_catalog, params, True, False, False)

class file_wrapper(wrapper):

    def get_file_location(self, object_key):
        return self.get_folder(object_key) + self.get_name(object_key)

    def other_init(self, maker, params):
        maker.set_param(params, "source_instance", self)
        self.cache = caches.file_cache_for_wrapper(maker, params)
        
    @print_stuff_dec
    def get_holding_folder(self):
        return constants.BIN_FOLDER

    
# for now, assume source_wrapper instance is available.  if instance is available, that means it wasn't created yet, so wouldn't be in wrapper registry, so don't need to fetch it from there.  in the future, might separately create source instance first, in which case might instead pass in the parameters needed to fetch the param.  or could fetch the instance first externally, then pass to this __init__
class generic_dumper_wrapper(file_wrapper):

    def __repr__(self):
        return self.__class__.__name__  + '-' + self.source_wrapper.__repr__()

    def other_init(self, maker, params):
        self.source_wrapper = maker.get_param(params, "dumper_source_instance")
        #print self, self.source_wrapper
        #pdb.set_trace()
        self.cache = caches.file_cache_for_wrapper(maker, params)

    def before_init(self, maker, params):
        return maker, params

    def dump_object(self, object):
        pass

    def get_folder(self, object_key):
        return self.source_wrapper.get_folder(object_key)

    # never pickle, since object returned is a file handle
    @dec
    @print_stuff_dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        #pdb.set_trace()
        assert to_pickle == False
        assert recalculate == True
        object = self.old_get_var_or_file(self.source_wrapper, params, global_stuff.recalculate, to_pickle)
        self.dump_object(object)
        return open(self.get_holding_location(), 'rb')

class pkdW(generic_dumper_wrapper):
    
    def dump_object(self, object):
        pickle.dump(object, open(self.get_holding_location(), 'wb'))

class mfdW(generic_dumper_wrapper):

    # do i have to create the folder first?
    def dump_object(self, object):
        write_mat(object, self.get_holding_location())

class generic_vect_file_dumper_wrapper(generic_dumper_wrapper):

    # do i have to create the folder first?
    def dump_object(self, object):
        write_mat(object, self.get_holding_location())

class wrapper_catalog(obj_wrapper, indexing_wrapper):

    # since there is no maker, hackishly set it to self
    def other_init(self, maker, params):
        maker.set_param(params, "source_instance", self)
        self.cache = caches.object_cache_for_wrapper(maker, params)

    def is_indexed(self):
        return False

    # params contains which_wrapper, and if which_wrapper is a generic_dumper_wrapper, contains 
    @dec
    def constructor(self, params, recalculate = False, to_pickle = False, to_filelize = False):
        wrapper_instance = self.get_param(params, "which_wrapper_class")(self, params)
        return wrapper_instance

# only purpose of this class is to give wrapper_catalog instance a "maker"
class famished_wrapper(object):

    def set_param(self, params, key, val, old_self = None):
        params.set_param(key, val)
        return params

    def get_param(self, params, key, record = True, old_self = None):
        return params.get_param(key)

    def get_var_or_file(self, wrapper, params, recalculate, to_pickle, to_filelize = False, old_self = None):
        #print '                       ', self, wrapper

        self.set_param(params, "which_wrapper_class", wrapper)
        from wc import wc
        wrapper_used_keys, wrapper_all_keys, the_wrapper = old_get_var_or_file(wc, params, True, False, False, old_self)
        used_keys, all_keys, x = the_wrapper.constructor(params, recalculate, to_pickle, to_filelize)
        #self.temp_dependents_keys[-1]  = self.temp_dependents_keys[-1].union(all_keys)
        return x

    def old_get_var_or_file(self, the_wrapper, params, recalculate, to_pickle, to_filelize = False, old_self = None):
        #print '                       ', self, wrapper
        used_keys, all_keys, x = the_wrapper.constructor(params, recalculate, to_pickle, to_filelize)
        return x

