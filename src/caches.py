import constants
import os
import cPickle as pickle
import pdb
#import wrapper
import subprocess
import wrapper
import global_stuff

#from import *

# act of caching: should be putting key, location of file in dictionary, but i don't have that.
class file_cache_for_wrapper(object):

    # files_created is the equivalent of the dump in the objects_cache.
    def __init__(self, maker, params):
        self.maker = maker
        self.the_wrapper = self.maker.get_param(params, "source_instance")
        self.files_created = set()


    # if we are not trusting the files previously in file system, then we only say a file is there if it was created this round
    def has(self, object_key, recalculate):
#        pdb.set_trace()
        if object_key in self.files_created:
            return True
        elif not recalculate:

            location = self.the_wrapper.get_file_location(object_key)
            # check if we want to override
            

            if os.path.isfile(location):
                if self.the_wrapper.whether_to_override(location):
                    return False
                return True
        return False

    # get is only called if has was just called and returned true
    def get(self, object_key, recalculate, mode = 'r'):
        self.allocate(object_key)
        return open(self.the_wrapper.get_file_location(object_key), mode)

    # moves object from temporary to permanent location.  makes destination folder if needed. keep track of which files were get'ed in this execution.  set is only called right after get is called
    # if this is not chain specific, only do this
    def set(self, object_key, object, to_pickle, params, to_filelize = None, always_recalculate = False):
        self.allocate(object_key)
        subprocess.call(['mv', object.name, self.the_wrapper.get_file_location(object_key)])
        self.files_created.add(object_key)
        #pdb.set_trace()
        return open(self.the_wrapper.get_file_location(object_key), 'r')

    # creates folder for object whose type is specified by the wrapper. this step does the prep work allowing the object to be put in the cache, which in this case is the folder in the file system
    def allocate(self, object_key):
        #if not os.path.exists(self.the_wrapper.get_folder(object_key)):
        try:
            os.makedirs(self.the_wrapper.get_folder(object_key))
        except:
            pass

# act of caching: putting object into cache's dictionary
class object_cache_for_wrapper(object):

    # wrapper is instance.  set instance if want dumper to use it.
    def __init__(self, maker, params):
        self.maker = maker
        self.dump = {}
        self.the_wrapper = self.maker.get_param(params, "source_instance")
        # it's not possible that the dumper exists already bc that would have required an instance of the_wrapper, which can't exist since we are creating it.  give dumper_wrapper the_wrapper instance.
        # POSSIBLY PROBLEM AREA
        #self.maker.set_param(params, "which_wrapper_class", wrapper.pkdW)
        #pdb.set_trace()
        self.maker.set_param(params, "dumper_source_instance", self.the_wrapper)
        self.pickle_dumper_wrapper = wrapper.pkdW(maker, params)
        self.maker.set_param(params, "dumper_source_instance", self.the_wrapper)
        #pdb.set_trace()
        self.file_dumper_wrapper = self.the_wrapper.get_file_dumper(maker, params)


    def has(self, object_key, recalculate):
        #print self, self.the_wrapper
        if object_key in self.dump:
            
            return True
        else:

            if self.the_wrapper.whether_to_override():

                return False

            
            return self.pickle_dumper_wrapper.has(object_key, recalculate)
    
    # if we don't trust pickles, only trust it if it was created this round.  stuff created this round is in self.dump
    def get(self, object_key, recalculate):
        if object_key in self.dump:
            return self.dump[object_key]
        elif self.pickle_dumper_wrapper.has(object_key, recalculate):
            assert recalculate == False
            f = self.pickle_dumper_wrapper.get(object_key, recalculate)

            obj = pickle.load(open(f.name, 'rb'))
            return obj
        raise KeyError
    
    # is it possible that set is called when the pickled object is accurate.  yes, only if you created the object, pickled, then deleted the used_keys and all_keys_cache pickles
    # self.pickle_wrapper where constructor makes the pickle.  params would include object_key and object
    #@print_stuff_dec
    def set(self, object_key, object, to_pickle, params, to_filelize, always_recalculate = False):
        #print self.the_wrapper
        #pdb.set_trace()


        if len(self.dump) > constants.CACHE_MAX_SIZE:
            self.dump.clear()

        
        self.dump[object_key] = object
        if to_pickle: 
            #pdb.set_trace()
            # here, already stored reference to a wrapper instance.  so just call it directly.
            temp_f = self.the_wrapper.old_get_var_or_file(self.pickle_dumper_wrapper, params, True, False, False, always_recalculate)
        if to_filelize:
            #pdb.set_trace()
            assert self.file_dumper_wrapper != None
            #pdb.set_trace()
            temp_f = self.the_wrapper.old_get_var_or_file(self.file_dumper_wrapper, params, True, False, False, always_recalculate)

        try:
            object.set_object_key(str(object_key))
        except:
            pass

            
        return object

    # don't have to do anything to put the object in the cache
    def allocate(self, object_key):
        pass

# is like a regular objects cache, but pickles entire dictionary instead of individual items
class akcO(object):

    def __repr__(self):
        return self.__class__.__name__ + '-' + self.the_wrapper.__repr__()

    def __init__(self, maker, params):
        self.dump = {}
        self.the_wrapper = maker.get_param(params, "source_instance")
        self.pickle_location = constants.BIN_FOLDER + self.__repr__() + '.pk'
        self.existing_dump = {}
        if global_stuff.recalculate == False:
            try:
                self.existing_dump = pickle.load(open(self.pickle_location, 'rb'))
            except:
                self.existing_dump = {}
        
        self.pickles_created = set()
        self.maker = maker

    def get_keys(self, recalculate):
        if not recalculate:
            return set(self.dump.keys()).union(set(self.existing_dump.keys()))
        else:
            return self.dump.keys()
            

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
    # if not recalculating, should i pickle existing dump?
    def set(self, key, val, to_pickle, recalculate):
        self.dump[key] = val
        if to_pickle and key not in self.pickles_created:
            if recalculate:
                pickle.dump(self.dump, open(self.pickle_location, 'wb'))
            else:
                self.existing_dump.update(self.dump)
                pickle.dump(self.existing_dump, open(self.pickle_location, 'wb'))
            self.pickles_created.add(key)

class index_cache(akcO):

    # for now, decide that if i'm going to index, then i'm also going to pickle.
    # it's possible that you are caching an object(it wasn't in object cache), but its index is already here
    def get_and_set_index(self, key, to_reindex = global_stuff.to_reindex):
        #print self.dump
        if self.has(key, to_reindex):
            return self.get(key)
        else:
            cur_keys = self.get_keys(to_reindex)
            if len(cur_keys) == 0:
                
                new_val = 0
            else:
                new_val = len(cur_keys) + 1
            self.set(key, new_val, True, to_reindex)
            return new_val
            
        
        

class ukcO(object):

    def __repr__(self):
        return self.__class__.__name__ + '-' + self.the_wrapper.__repr__()

    def __init__(self, maker, params):

        self.dump = None
        self.val = None
        self.the_wrapper = maker.get_param(params, "source_instance")
        self.pickle_location = constants.BIN_FOLDER + self.__repr__() + '.pk'
        self.pickle_created = False
        if global_stuff.recalculate == False:
            try:
                self.dump = pickle.load(open(self.pickle_location, 'rb'))

            except:
                self.dump = None
        self.maker = maker

    # recalculate specifies whether to look at possible pickle file and use it
    def has(self, recalculate):
        if self.val != None:
            return True
        elif not recalculate:
            if self.dump != None:
                return True
        elif self.pickle_created == True:
            assert False
        return False

    # recalculate determines whether to look at pickle
    def get(self, recalculate):
        if self.val != None:
            return self.val
        elif not recalculate:
            if self.dump != None:
                return self.dump
        elif self.pickle_created == True:
            assert False
            return self.dump
        raise KeyError

    def set(self, val, to_pickle):
        self.val = val
        if to_pickle and not self.pickle_created:
            pickle.dump(val, open(self.pickle_location, 'wb'))
            self.pickle_created = True

class skcO(ukcO):
    pass

class akkkscO(ukcO):
    pass

import wrapper
