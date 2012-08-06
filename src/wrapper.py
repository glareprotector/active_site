from param import param

import pdb, os, subprocess, constants, pickle
from global_stuff import print_stuff_dec

# note for all caches: recalculate specifies whether you want to use the files/pickles that are present when cache is created

# act of caching: moving file from temporary holding spot to proper folder.  cache structure: file system
class file_cache_for_wrapper(object):

    # files_created is the equivalent of the dump in the objects_cache
    def __init__(self, the_wrapper):
        self.the_wrapper = the_wrapper
        self.files_created = set()


    # if we are not trusting the files previously in file system, then we only say a file is there if it was created this round
    @print_stuff_dec
    def has(self, object_key, recalculate):
        #pdb.set_trace()
        if object_key in self.files_created:
            return True
        elif not recalculate:
            if os.path.isfile(self.the_wrapper.get_file_location(object_key)):
                return True

        return False

    # get is only called if has was just called and returned true
    # actually, can also be used for object creation
    def get(self, object_key, mode = 'r'):
        self.allocate(object_key)
        return open(self.the_wrapper.get_file_location(object_key), mode)

    # moves object from temporary holding place to where it should be stored permenantly.  makes destination folder if needed. keep track of which files were get'ed in this execution
    # set is only called right after get is called
    @print_stuff_dec
    def set(self, object_key, object, to_pickle):
        #pdb.set_trace()
        self.allocate(object_key)
        subprocess.call(['mv', object.name, self.the_wrapper.get_file_location(object_key)])
        self.files_created.add(object_key)
        return open(self.the_wrapper.get_file_location(object_key), 'r')

    # creates folder for object whose type is specified by the wrapper
    # this step does the prep work allowing the object to be put in the cache, which in this case is the folder in the file system
    @print_stuff_dec
    def allocate(self, object_key):
        if not os.path.exists(self.the_wrapper.get_folder(object_key)):
            os.makedirs(self.the_wrapper.get_folder(object_key))

# act of caching: putting object into cache's dictionary
class object_cache_for_wrapper(object):

    def __init__(self, the_wrapper):
        self.dump = {}
        self.the_wrapper = the_wrapper
        self.a_file_cache = file_cache_for_wrapper(the_wrapper)

        self.pickles_created = set()

        self.finish_init()

    def has(self, object_key, recalculate):
        if object_key in self.dump:
            return True
        elif not recalculate:
            # if self.the_wrapper.pickle_wrapper.has(object_key, recalculate):
            if self.a_file_cache.has(object_key, recalculate):
                return True
    
    # if we don't trust pickles, only trust it if it was created this round.  stuff created this round is in self.dump
    def get(self, object_key, recalculate):
        pdb.set_trace()
        if object_key in self.dump:
            return dump[object_key]
        elif self.a_file_cache.has(object_key, recalculate):
            # if self.the_wrapper.pickle_wrapper.has(object_key, recalculate):
            object = pickle.load(self.a_file_cache.get(object_key, 'rb'))
            self.dump[object_key] = object
            return object
        pdb.set_trace()
        raise KeyError
    
    # is it possible that set is called when the pickled object is accurate.  yes, only if you created the object, pickled, then deleted the used_keys and all_keys_cache pickles
    # self.pickle_wrapper where constructor makes the pickle.  params would include object_key and object
    @print_stuff_dec
    def set(self, object_key, object, to_pickle, params):
        self.dump[object_key] = object
        if to_pickle and not self.a_file_cache.has(object_key, True):
        #if to_pickle and object_key not in self.pickles_created:
            #f = self.a_file_cache.get(object_key, 'wb')
            ##f = open(self.the_wrapper.get_holding_location(), 'rb')
            ##pickle.dump(object, f)
            ##self.a_file_cache.set(object_key, f)
            #self.pickles_created.add(object_key)
        if to_pickle: 
            temp_f = self.the_wrapper.get_var_or_file(pickle_wrapper, params, True, None)
            
        return object

    # don't have to do anything to put the object in the cache
    def allocate(self, object_key):
        pass

# is like a regular objects cache, but pickles entire dictionary instead of individual items
class all_keys_obj_cache(object):

    def __init__(self, the_wrapper):
        self.dump = {}
        self.the_wrapper = the_wrapper
        self.pickle_location = constants.BIN_FOLDER + self.__class__.__name__ + '-' + the_wrapper.__class__.__name__ + '.pickle'
        if os.path.isfile(self.pickle_location):
            self.existing_dump = pickle.load(open(self.pickle_location, 'rb'))
        else:
            self.existing_dump = {}
        
        self.pickles_created = set()

    def has(self, key, recalculate):
        if key in self.dump:
            return True
        elif not recalculate:
            if key in self.existing_dump:
                return True
        return False
    
    # only called if has returned true.  if has was called with recalculate = True and has was true, means key was in self.dump.
    # if has was called with recalculate = False and was true, key could be in both self.dump or self.existing dump.
    def get(self, key):
        if key in self.dump:
            return self.dump[key]
        else:
            assert key in self.existing_dump
            self.dump[key] = self.existing_dump[key]
            return self.existing_dump[key]
        assert False

    # if item was pickled this round, don't need to repickle it.  so, keep track of the stuff that has been pickled
    def set(self, key, val, to_pickle):
        self.dump[key] = val
        if to_pickle and key not in self.pickles_created:
            pickle.dump(self.dump, open(self.pickle_location, 'wb'))
            self.pickles_created.add(key)

class used_keys_obj_cache(object):

    def __init__(self, the_wrapper):
        self.dump = None
        self.the_wrapper = the_wrapper
        self.pickle_location = constants.BIN_FOLDER + self.__class__.__name__ + '-' + the_wrapper.__class__.__name__ + '.pickle'
        
        self.pickle_created = False

    # recalculate specifies whether to look at possible pickle file and use it
    def has(self, recalculate):
        if self.dump != None:
            return True
        elif not recalculate:
            if os.path.isfile(self.pickle_location):
                return True
        return False

    # recalculate determines whether to look at pickle
    def get(self):
        if self.dump != None:
            return self.dump
        elif os.path.isfile(self.pickle_location):
            val = pickle.load(open(self.pickle_location, 'rb'))
            self.dump = val
            return val
        raise KeyError

    def set(self, val, to_pickle):
        self.dump = val
        if to_pickle and not self.pickle_created:
            pickle.dump(val, open(self.pickle_location, 'wb'))
            self.pickle_created = True



# this is the decorator
#@print_stuff_dec 
def dec(f):

    def affects_which_wrappers_are_called(key, val):
        if val.__class__.__name__ in ['int', 'float', 'str']:
            return False
        else:
            return True

    # a key/value should go into the key only if the value determines which wrappers will be called
    # assumptions: if a param does not affect which wrappers are called at this node, then it does not affect which wrappers are called
    # at any descendant node.  as long as the param is used in the same way in every node, this will be the case
    @print_stuff_dec
    def get_all_keys_key(self, params, used_keys):
        temp = {}
        for key in used_keys:
            if affects_which_wrappers_are_called(key, params.get_param(key)):
                temp[key] = params.get_val(key)
        return param(temp)

    def get_all_keys(self, params):
        return self.used_keys_cache.get().union(self.temp_dependents_keys)

    def get_objects_key(self, params, all_keys):
        return param.restriction(params, all_keys)

    @print_stuff_dec
    def cache_everything_f_poster(self, params, recalculate, to_pickle):
        object = f(self, params, recalculate, to_pickle)
        self.used_keys_cache.set(self.temp_used_keys, to_pickle)
        all_keys_key = get_all_keys_key(self, params, self.used_keys_cache.get())
        all_keys = get_all_keys(self, params)
        self.temp_dependents_keys = set()
        self.all_keys_cache.set(all_keys_key, all_keys, to_pickle)
        object_key = get_objects_key(self, params, all_keys)
        object = self.cache.set(object_key, object, to_pickle, params)
        return self.used_keys_cache.get(), all_keys, object
    
    @print_stuff_dec    
    def h(self, params, recalculate, to_pickle):
        #pdb.set_trace()
        if self.used_keys_cache.has(recalculate):
            used_keys = self.used_keys_cache.get()
            all_keys_key = get_all_keys_key(self, params, used_keys)
            if self.all_keys_cache.has(all_keys_key, recalculate):
                all_keys = self.all_keys_cache.get(all_keys_key)
                object_key = get_objects_key(self, params, all_keys)
                if self.cache.has(object_key, recalculate):
                    object = self.cache.get(object_key)
                    return used_keys, all_keys, object
        return cache_everything_f_poster(self, params, recalculate, to_pickle)

    return h


# this is the function wrapper

class wrapper(object):

    def get_holding_location(self):
        return constants.BIN_FOLDER + str(id(self))

    def get_folder(self, object_key):
        return constants.BIN_FOLDER

    def get_wrapper_name(self):
        return self.__class__.__name__

    def get_name(self, object_key):
        return self.get_wrapper_name + str(object_key)
    
    def __init__(self):

        self.used_keys_cache = used_keys_obj_cache(self)
        self.all_keys_cache = all_keys_obj_cache(self)

        self.temp_used_keys = set()
        self.temp_dependents_keys = set()

        self.finish_init()

    def get_param(self, params, key):
        self.temp_used_keys.add(key)
        return params.get_param(key)
    @print_stuff_dec
    def get_var_or_file(self, wrapper, params, recalculate, to_pickle):
        used_keys, all_keys, x = wrapper.constructor(params, recalculate, to_pickle)
        self.temp_dependents_keys  = self.temp_dependents_keys.union(all_keys)
        return x

    # returns only the object, but after decorating, will return used_keys, all_keys, object
    @dec
    def constructor(self, params, recalculate, to_pickle):
        pass

    def get_cache(self):
        pass

class obj_wrapper(wrapper):


    def get_file_location(self, object_key):
        return self.get_folder(object_key) + self.get_name(object_key) + '.pickle'

    def finish_init(self):
        self.cache = object_cache_for_wrapper(self)
        self.pickle_wrapper = generic_pickle_wrapper(self)

class file_wrapper(wrapper):
    @print_stuff_dec
    def get_file_location(self, object_key):
        return self.get_folder(object_key) + self.get_name(object_key)

    def finish_init(self)
        self.cache = file_cache_for_wrapper(self)
        
    #@print_stuff_dec
    #def get_holding_location(self, params):
    #    return constants.BIN_FOLDER + self.__class__.__name__ + '-' + str(params)
    @print_stuff_dec
    def get_holding_folder(self):
        return constants.BIN_FOLDER

class base_wrapper(obj_wrapper):

    @dec
    @print_stuff_dec
    def constructor(self, params, recalculate, to_pickle):
        x = wrapper.get_param(self, params, 'x')
        return x

the_base_wrapper = base_wrapper()

class squared_wrapper(obj_wrapper):

    @dec
    @print_stuff_dec
    def constructor(self, params, recalculate, to_pickle):
        return self.get_var_or_file(the_base_wrapper, params, recalculate, to_pickle) + 1
        
#the_squared_wrapper = squared_wrapper()


#print 'aff'
#params = param({'x':5})
#pdb.set_trace()
#a = the_squared_wrapper.constructor(params, False, True)
