'''
Created on Mar 9, 2012

@author: glareprotector
'''

from constants import *
import global_stuff
import cPickle as pickle
import cache
import os
import string
import pdb

class obj_manager(object):
    '''
    responsible for ensuring presence of objects in cache
    '''


    def __init__(self):
        self.the_cache = cache.cache()
    
    '''
    returns the variable with the desired name.  pickles and caches if desired
    '''
    def get_variable(self, the_obj_wrapper, recalculate = False, to_pickle=True, use_pickle=True):
        #to_pickle = False
        #use_pickle = False
        pickle_file = the_obj_wrapper.get_pickle_location()
        #print 'AAAAAAAAAAAAAAAAA', the_obj_wrapper.get_pickle_location(), recalculate
        if self.has_in_cache(the_obj_wrapper.get_name()):
            #pdb.set_trace()
            #print the_obj_wrapper.get_name(), 'has in cache'
            return self.get_from_cache(the_obj_wrapper.get_name())
        elif not recalculate and use_pickle and global_stuff.the_file_manager.has_file(the_obj_wrapper.get_pickle_location()):
            #pdb.set_trace()
            #print the_obj_wrapper.get_name(), 'has pickle'
            pickle_handle = global_stuff.the_file_manager.get_pickle_handle(the_obj_wrapper.get_pickle_location())
            obj = pickle.load(pickle_handle)
            self.cache_variable(obj, the_obj_wrapper.get_name())
        else:
            #print the_obj_wrapper.get_name(), 'constructing from scratch'
            obj = the_obj_wrapper.constructor(recalculate)
            self.cache_variable(obj, the_obj_wrapper.get_name())
            if to_pickle:
                global_stuff.the_file_manager.pickle_obj(obj, pickle_file)
                
        return obj
    
    '''
    caches the specified variable
    '''            
    def cache_variable(self, obj, cache_name):
        self.the_cache.put(obj, cache_name)
        
    
    '''
    indicates whether variable of specified name is in cache
    '''
    def has_in_cache(self, cache_name):
        return self.the_cache.has(cache_name)
    
    def get_from_cache(self, cache_name):
        return self.the_cache.get(cache_name)
    
    
    
class file_manager(object):
    '''
    responsible for ensuring the presence of files on disk
    '''

    
    def __init__(self):
        pass
    
    def has_file(self,f):
        return os.path.isfile(f)
    
    '''
    returns a handle to thec desired (non-pickle) file.  if it isn't there, call the file_wrapper's constructor
    '''
    def get_file_handle(self, the_file_wrapper, recalculate = False, option = 'r'):
        full_path = the_file_wrapper.get_file_location()
        #pdb.set_trace()
        if not recalculate and os.path.isfile(full_path):
            #print the_file_wrapper.get_name(), 'has file'
            #pdb.set_trace()
            return open(full_path)
        else:
            #print the_file_wrapper.get_name(), 'NO file'
            file_folder = string.join(full_path.split('/')[:-1],'/')
            if not os.path.exists(file_folder):
                os.makedirs(file_folder)
            the_file_wrapper.constructor(recalculate)
            #pdb.set_trace()
            return open(full_path, option)
    
    '''
    returns a handle to a pickle file.  if the pickle isn't there, don't try to create it
    '''
    def get_pickle_handle(self, the_pickle_file):
        return open(the_pickle_file, 'rb')
    
    def pickle_obj(self, obj, pickle_file):
        handle = open(pickle_file,'wb')
        pickle.dump(obj, handle)    
