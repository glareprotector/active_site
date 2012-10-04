import new_new_objects as objects



class pdb_name_struct(object):

    def __init__(self, pdb_name, chain_letter):
        self.pdb_name = pdb_name
        self.chain_letter = chain_letter

class pdb_results_struct(object):

    def __init__(self, scores=[], true_classes=[], pdb_structs=[], sample_lengths=[]):
        self.scores = scores
        self.true_classes = true_classes
        self.pdb_structs = pdb_structs
        self.sample_lengths = sample_lengths;


    


    def __add__(self, other):
        new_scores = self.scores + other.scores
        new_true_classes = self.true_classes + other.true_classes
        new_pdb_structs = self.pdb_structs + other.pdb_structs
        new_sample_lengths = self.sample_lengths + other.sample_lengths
        return pdb_results_struct(new_scores, new_true_classes, new_pdb_structs, new_sample_lengths)




class data(object):

    def __init__(self, maker, params, recalculate):
        self.data_list_file = maker.get_param(params, 'data_list_file')
        # ERROR objects.mW should return list of pdb_name_structs
        self.names = maker.get_var_or_file(params, objects.mW, recalculate, True, True)

    def get_training_names(self):
        return self.names

    def get_testing_names(self):
        return []

    # what kind of objects does iterator return? returns fold objects, which will have a source of self
    def get_folds(self, maker, params, recalculate):
        maker.set_param(params, 'fsc', self)
        maker.set_param(params, 'fldk', self.k)
        maker.set_param(params, 'wfldk', self.cur_k)
        return helper.get_var_or_file(params, objects.brW, recalculate)


class fold(data):

    def __init__(self, maker, params, recalculate):
        self.source = maker.get_param('fsc')
        self.k = maker.get_param(params, 'fldk')
        self.which_fold = maker.get_param(params, 'wfldk')
        self.params = params.get_copy()

        self.cur_k = 0

    def get_training_names(self):
        source_training_names = self.source.get_training_names()
        return [source_training_names[i] for i in range(len(source_training_names)) if i % self.k == self.which_fold]

    def get_testing_names(self):
        source_training_names = self.source.get_training_names()
        return [source_training_names[i] for i in range(len(source_training_names)) if i % self.k != self.which_fold]

class cv_results_struct(object):

    def __init__(self, helper, params):
        self.data = helper.get_param(params, 'data')
        self.H = helper.get_param(params, 'H')
        self.k = helpers.get_param(params, 'fldk')
        self.cache = []
    
    def add(self, fold, theta, results):
        self.cache.append([fold, theta, results])

    def get_regular_results(self):
        merged_results = pdb_results_struct()
        if len(self.cache == 0):
            return pdb_results_struct()
        else:
            merged_results = self.cache[0][2]
            for i in range(2,len(self.cache)):
                merged_results = merged_results + self.cache[i][2]

    def get_ROC(self):
        merged_results = get_regular_results(self)



class H_searcher(object):

    def __init__(self, helper, params):
        self.data = helper.get_param(params, 'data')
        if helper.get_param(params, 'wobj') == 0:
            self.param_names = ['asdf', 'asdf']
        self.count = 0
        self.cache = []

    def done(self):
        return self.count > global_stuff.max_search_num

    def get_next_H(self):
        import random
        return [random.random() for i in len(self.param_names)]

    def add(self, H, cv_results):
        self.cache.add(H, cv_results)


