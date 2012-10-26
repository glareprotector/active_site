import pdb
import global_stuff
import helper
import string

class param(object):
    
    def __init__(self, param_dict={}):
        self.param_dict = param_dict
        
    def get_param(self, which):
        return self.param_dict[which]
    
    def get_keys(self):
        return self.param_dict.keys()

    def get_sorted_keys(self):
        return sorted(self.param_dict.keys())
    
    def set_param(self,key,val):
        self.param_dict[key] = val
        
    def has_param(self,key):
        return key in self.param_dict.keys()
        
    def __add__(self, other):
        return param(dict(self.param_dict.items() + other.param_dict.items()))
    
    def __str__(self):
        to_join = []
        for key in sorted(self.param_dict.keys()):
            temp = str(key) + str(self.param_dict[key])
            to_join.append(temp)
            
        return '(' + string.join(to_join, '_') + ')'
        
        return helper.shorten(str(sorted(self.param_dict.iteritems())))

    # returns a list of the keys in sorted key order
    def sorted_listify(self):
        res = []
        for key in self.get_sorted_keys():
            res.append(self.get_param(key))
        return res
    
    def get_copy(self):
        to_return = param({})
        for key in self.param_dict.keys():
            to_return.set_param(key, self.get_param(key))
        return to_return

    def remove_param(self, key):
        del self.param_dict[key]

    def rename(self, old_key, new_key):
        temp = self.get_param(old_key)
        self.remove_param(old_key)
        self.set_param(new_key, temp)
    
    @classmethod
    def restriction(cls, params, A):
        #pdb.set_trace()
        temp = {}
        for key in A:
            temp[key] = cls.get_param(params, key)
        return param(temp)

    def __hash__(self):
        return hash(self.__str__())

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    # returns A U B, with values in A taking precedence

    @classmethod
    def merge(cls, A, B):
        import pdb
        pdb.set_trace()
        toReturn = cls.get_copy(A)
        for key in cls.get_keys(B):
            toReturn.set_param(key, B.get_param(key))
        return toReturn


    def merge_non_class(self, B):
        toReturn = self.get_copy()
        for key in B.get_keys():
            toReturn.set_param(key, B.get_param(key))
        return toReturn
        
    def __repr__(self):
        return self.__str__()

    def flatten_hp(self, maker):

        hp = maker.get_param(self,'hp')

        for key in hp.get_keys():
            maker.set_param(self, key, hp.get_param(key))
