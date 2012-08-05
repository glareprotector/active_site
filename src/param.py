import pdb

class param(object):
    
    def __init__(self, param_dict):
        self.param_dict = param_dict
        
    def get_param(self, which):
        return self.param_dict[which]
    
    def get_keys(self):
        return self.param_dict.keys()
    
    def set_param(self,key,val):
        self.param_dict[key] = val
        
    def has_param(self,key):
        return key in self.param_dict.keys()
        
    def __add__(self, other):
        return param(dict(self.param_dict.items() + other.param_dict.items()))
    
    def __str__(self):
        return str(self.param_dict)
    
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



    # returns A U B, with values in A taking precedence

    @classmethod
    def merge(cls, A, B):
        toReturn = cls.get_copy(A)
        for key in cls.get_keys(B):
            toReturn.set_param(key, B.get_param(key))
        return toReturn
        
    
