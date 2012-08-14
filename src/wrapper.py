from param import param

import pdb, os, subprocess, constants, pickle
from global_stuff import print_stuff_dec, write_mat, write_vect

import caches
from wrapper_decorator import *

# note for all caches: recalculate specifies whether you want to use the files/pickles that are present when cache is created

# act of caching: moving file from temporary holding spot to proper folder.  cache structure: file system

# this is the function wrapper

class wrapper(object):

    def is_indexed(self):
        return True

    def makes_index(self):
        return False

    def __repr__(self):
        return self.name

    def get_holding_location(self):
        return constants.BIN_FOLDER + str(id(self))

    def get_folder(self, object_key):
        return constants.BIN_FOLDER

    def get_wrapper_name(self):
        return self.__repr__()

    def get_name(self, object_key):
        return self.get_wrapper_name() + str(object_key)
    
    def basic_init(self, maker, params):
        maker.set_param(params, "source_instance", self)
        self.used_keys_cache = caches.used_keys_obj_cache(self, maker, params)
        self.all_keys_cache = caches.all_keys_obj_cache(self, maker, params)
        self.temp_used_keys = set()
        self.temp_dependents_keys = set()
        self.temp_new_param_keys = set()
        if self.makes_index():
            self.object_key_to_index = index_cache(self, maker, params)
        self.name = self.__class__.__name__
        self.maker = maker
        self.construction_params = params

    def __init__(self, maker, params = param({})):
        self.basic_init(maker, params)

    def set_param(self, params, key, val):
        self.temp_new_param_keys.add(key)
        params.set_param(key, val)
        return params

    def get_param(self, params, key, record = True):
        if record:
            self.temp_used_keys.add(key)
        return params.get_param(key)

    @print_stuff_dec
    def get_var_or_file(self, wrapper, params, recalculate, to_pickle, to_filelize = False):
        print '                       ', self, wrapper
        used_keys, all_keys, x = wrapper.constructor(params, recalculate, to_pickle, to_filelize)
        self.temp_dependents_keys  = self.temp_dependents_keys.union(all_keys)
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

class indexing_wrapper(wrapper):

    def makes_index(self):
        return True

class indexed_wrapper(wrapper):

    def process_index(self, index):
        self.name = str(index)

class obj_wrapper(wrapper):

    def get_file_dumper(self):
        return None

    def get_file_location(self, object_key):
        return self.get_folder(object_key) + self.get_name(object_key) + '.pickle'

    def __init__(self):
        self.basic_init()
        self.cache = caches.object_cache_for_wrapper(self)

class mat_obj_wrapper(obj_wrapper):
    
    def get_file_dumper(self):
        return generic_mat_file_dumper_wrapper(self)

class vect_obj_wrapper(obj_wrapper):

    def get_file_dumper(self):
        return generic_vect_file_dumper_wrapper(self)

class file_wrapper(wrapper):

    def get_file_location(self, object_key):
        return self.get_folder(object_key) + self.get_name(object_key)

    def __init__(self):
        self.basic_init()
        self.cache = caches.file_cache_for_wrapper(self)
        
    @print_stuff_dec
    def get_holding_folder(self):
        return constants.BIN_FOLDER

    
# normal file_wrappers don't have a source, but a dumper does.  dumper's folder should be same as its source
# what if you have 2 generic dumpers with different sources.  then they're different.  what if same wrapper with 2 different sources? still different.
class generic_dumper_wrapper(file_wrapper):

    # same as regular file wrapper, except that it is initialized with a object_wrapper - the wrapper that created it
    def __init__(self, maker, params):
        self.basic_init()
        #self.source_wrapper = source_wrapper # source makes the object
        self.source_wrapper = self.maker.get_param("source_instance")
        self.cache = caches.file_cache_for_wrapper(self)

    def get_wrapper_name(self):
        return self.__class__.__name__ + '-' + self.source_wrapper.get_wrapper_name()

    def dump_object(self, object):
        pass

    def get_folder(self, object_key):
        return self.source_wrapper.get_folder(object_key)

    # never pickle, since object returned is a file handle
    @dec
    @print_stuff_dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False):
        object = self.get_var_or_file(self.source_wrapper, params, recalculate, to_pickle)
        self.dump_object(object)
        return open(self.get_holding_location(), 'rb')

class generic_pickle_dumper_wrapper(generic_dumper_wrapper):
    
    def dump_object(self, object):
        pickle.dump(object, open(self.get_holding_location(), 'wb'))

class generic_mat_file_dumper_wrapper(generic_dumper_wrapper):

    # do i have to create the folder first?
    def dump_object(self, object):
        write_mat(object, self.get_holding_location())

class generic_vect_file_dumper_wrapper(generic_dumper_wrapper):

    # do i have to create the folder first?
    def dump_object(self, object):
        write_mat(object, self.get_holding_location())

    

