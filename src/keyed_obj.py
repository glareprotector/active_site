

class keyed_obj(object):

    def set_object_key(self, object_key):
        self.object_key = object_key

    def __repr__(self):
        return self.object_key
