from wrapper_decorator import dec
import new_new_objects as objects
import new_features as features
import caches
import param

# will contain registry of wrappers.  if it is a wrapper, in constructor will obtain current constructor's number.  the constructor returns wrappers.  make it a wrapper even though not taking advantage of pickling, still need caching (wrappers string to idx doesn't store actual wrappers).  if you want to create a wrapper with parameters, would pass that in as a parameter in params  

class wrapper_catalog(wrapper, indexing_wrapper):

    def __init__(self):
        self.basic_init()


    # params contains which_wrapper, and if which_wrapper is a generic_dumper_wrapper, contains 
    @dec
    def constructor(self, params, recalculate = False, to_pickle, to_filelize = False):
        wrapper_instance = self.get_param(params, "which_wrapper")(self, params)
        
    

    # wrapper for constructor.  recalculate should always be the same value or else get with recalculate = false or true could return 2 different instances of wrapper.  saying recalculate = false means you want to keep the current indexing.  the_wrapper is a class.  this returns an instance
    def get_wrapper(self, the_wrapper, params = param.param({}), recalculate = False):
        params.set_param("which_wrapper", the_wrapper)
        used_keys, all_keys, wrapper_instance = constructor(params, recalculate, False, False)
        return wrapper_instance
