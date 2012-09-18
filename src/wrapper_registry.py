from wrapper_decorator import dec
import wrapper
import param
import pdb

# will contain registry of wrappers.  if it is a wrapper, in constructor will obtain current constructor's number.  the constructor returns wrappers.  make it a wrapper even though not taking advantage of pickling, still need caching (wrappers string to idx doesn't store actual wrappers).  if you want to create a wrapper with parameters, would pass that in as a parameter in params  

class wrapper_catalog(wrapper.obj_wrapper, wrapper.indexing_wrapper):

    # since there is no maker, hackishly set it to self
    def __init__(self, maker, params):
        #pdb.set_trace()
        maker = self
        self.basic_init(maker, params)
        self.maker.set_param(params, "source_instance", self)
        self.cache = caches.object_cache_for_wrapper(maker, params)


    def is_indexed(self):
        return False

    # params contains which_wrapper, and if which_wrapper is a generic_dumper_wrapper, contains 
    @dec
    def constructor(self, params, recalculate = False, to_pickle = False, to_filelize = False):
        wrapper_instance = self.get_param(params, "which_wrapper_class")(self, params)
        return wrapper_instance
    
#pdb.set_trace()
wc = wrapper_catalog(None, param.param({}))
