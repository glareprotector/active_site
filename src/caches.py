import constants
import os
import pickle
import pdb
import wrapper
import subprocess

from new_features import *
#import new_features as features

class file_cache_for_wrapper(object):

    # files_created is the equivalent of the dump in the objects_cache
    def __init__(self, the_wrapper):
        self.the_wrapper = the_wrapper
        self.files_created = set()


    # if we are not trusting the files previously in file system, then we only say a file is there if it was created this round
    #@print_stuff_dec
    def has(self, object_key, recalculate):
        #print self,self.the_wrapper
        #pdb.set_trace()
        if object_key in self.files_created:
            return True
        elif not recalculate:
            if os.path.isfile(self.the_wrapper.get_file_location(object_key)):
                return True

        return False

    # get is only called if has was just called and returned true
    # actually, can also be used for object creation
    def get(self, object_key, recalculate, mode = 'r'):
        self.allocate(object_key)
        return open(self.the_wrapper.get_file_location(object_key), mode)

    # moves object from temporary holding place to where it should be stored permenantly.  makes destination folder if needed. keep track of which files were get'ed in this execution
    # set is only called right after get is called
    #@print_stuff_dec
    def set(self, object_key, object, to_pickle, params, to_filelize = None):
        #pdb.set_trace()
        self.allocate(object_key)
        subprocess.call(['mv', object.name, self.the_wrapper.get_file_location(object_key)])
        self.files_created.add(object_key)
        return open(self.the_wrapper.get_file_location(object_key), 'r')

    # creates folder for object whose type is specified by the wrapper
    # this step does the prep work allowing the object to be put in the cache, which in this case is the folder in the file system
    #@print_stuff_dec
    def allocate(self, object_key):
        if not os.path.exists(self.the_wrapper.get_folder(object_key)):
            os.makedirs(self.the_wrapper.get_folder(object_key))

# act of caching: putting object into cache's dictionary
class object_cache_for_wrapper(object):

    def __init__(self, the_wrapper):
        self.dump = {}
        self.the_wrapper = the_wrapper
        self.pickle_dumper_wrapper = wrapper.generic_pickle_dumper_wrapper(self.the_wrapper)
        self.file_dumper_wrapper = self.the_wrapper.get_file_dumper()


    def has(self, object_key, recalculate):
        #print self, self.the_wrapper
        #pdb.set_trace()
        if object_key in self.dump:
            return True
        else:
            return self.pickle_dumper_wrapper.has(object_key, recalculate)
    
    # if we don't trust pickles, only trust it if it was created this round.  stuff created this round is in self.dump
    def get(self, object_key, recalculate):
        if object_key in self.dump:
            return self.dump[object_key]
        elif self.pickle_dumper_wrapper.has(object_key, recalculate):
            assert recalculate == False
            f = self.pickle_dumper_wrapper.get(object_key, recalculate)
            object = pickle.load(open(f.name, 'rb'))
            return object
        raise KeyError
    
    # is it possible that set is called when the pickled object is accurate.  yes, only if you created the object, pickled, then deleted the used_keys and all_keys_cache pickles
    # self.pickle_wrapper where constructor makes the pickle.  params would include object_key and object
    #@print_stuff_dec
    def set(self, object_key, object, to_pickle, params, to_filelize = False):
        print self.the_wrapper
        #pdb.set_trace()
        self.dump[object_key] = object
        if to_pickle: 
            #pdb.set_trace()
            temp_f = self.the_wrapper.get_var_or_file(self.pickle_dumper_wrapper, params, True, None, False)
        if to_filelize:
            #pdb.set_trace()
            assert self.file_dumper_wrapper != None
            temp_f = self.the_wrapper.get_var_or_file(self.file_dumper_wrapper, params, True, None, False)
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
            print self.pickle_location, self, the_wrapper
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
