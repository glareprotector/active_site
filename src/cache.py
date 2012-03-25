'''
Created on Mar 9, 2012

@author: glareprotector
'''


class cache(object):

    def __init__(self):
        self.the_dict = {}
        
    def has(self, cache_name):
        return cache_name in self.the_dict.keys()
    
    def put(self, obj, cache_name):
        self.the_dict[cache_name] = obj
    
    def get(self, cache_name):
        return self.the_dict[cache_name]