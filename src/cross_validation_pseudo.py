import new_new_objects as objects

import pdb

from param import param

class keyed_object(object):

    def __repr__(self):
        try:
            return self.object_key
        except:
            return self.__class__.__name__

    def set_object_key(self, object_key):
        self.object_key = object_key

    def get_object_key(self):
        return self.object_key


class pdb_name_struct(object):

    def __init__(self, pdb_name, chain_letter, start, end):
        self.pdb_name = pdb_name
        self.chain_letter = chain_letter
        self.start = start
        self.end = end

    def __repr__(self):
        return self.pdb_name + '_' + self.chain_letter

class pdb_results_struct(keyed_object):

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




class data(keyed_object):

    def __init__(self, maker, params, recalculate, data_list_file):
        self.data_list_file = data_list_file
        # ERROR objects.mW should return list of pdb_name_structs

        maker.set_param(params, 'data_list_file', data_list_file)
        self.names = maker.get_var_or_file(objects.ciW, params, recalculate, False, False)

    def get_training_names(self):
        return self.names

    def get_testing_names(self):
        return []

    # what kind of objects does iterator return? returns fold objects, which will have a source of self
    def get_fold(self, maker, params, recalculate, which_fold, num_folds):
        maker.set_param(params, 's', self)
        maker.set_param(params, 'm', num_folds)
        maker.set_param(params, 'k', which_fold)
        return maker.get_var_or_file(objects.buW, params, recalculate, False)

    def get_folds(self, maker, params, recalculate, num_folds):
        folds = []
        for i in range(num_folds):
            folds.append(self.get_fold(maker, params, recalculate, i, num_folds))
        return folds


class fold(data):

    def __init__(self, maker, params, recalculate, source, which_fold, num_folds):
        self.source = source
        self.k = num_folds
        self.which_fold = which_fold

        self.cur_k = 0

    def get_training_names(self):
        source_training_names = self.source.get_training_names()
        return [source_training_names[i] for i in range(len(source_training_names)) if i % self.k != self.which_fold]

    def get_testing_names(self):
        source_training_names = self.source.get_training_names()
        return [source_training_names[i] for i in range(len(source_training_names)) if i % self.k == self.which_fold]


class hp_searcher(keyed_object):

    # assume hp_values comes in form of param where each value is the points in the grid for that parameter
    def __init__(self, maker, params, recalculate):

        self.the_hp_search_results = None

        self.hp_stash = maker.get_var_or_file(objects.caW, params, recalculate, False, False)
        
        self.pos = -1

    def turn_on_search(self):
        self.pos = 0

    def done_searching(self):
        if self.pos == -1:
            assert False
        else:
            return self.pos >= len(self.hp_stash)

    def get_next_hp(self):
        ans = self.hp_stash[self.pos]
        self.pos = self.pos + 1
        return ans

    def set_hp_search_results(self, the_hp_search_results):
        self.the_hp_search_results = the_hp_search_results

class train_test_result(keyed_object):

    def __init__(self, maker, params, recalculate):



        the_fold = maker.get_param(params, 'f')
        hp = maker.get_param(params, 'hp')
        num_states = maker.get_param(params, 'ns')
        max_iter = maker.get_param(params, 'mx')
        which_obj = maker.get_param(params, 'wob')
        which_infer = maker.get_param(params, 'wif')
        which_reg = maker.get_param(params, 'wreg')


        self.hp = hp
        
        self.the_fold = the_fold

        training_list = self.the_fold.get_training_names()
        testing_list = self.the_fold.get_testing_names()
        import _test

        try:
            self.theta = _test.get_theta_given_training_data_and_hypers(maker, params, recalculate, num_states, max_iter, which_obj, which_reg, which_infer, training_list)
            if len(self.theta) == 0:
                raise Exception("couldn't train")
            self.results = _test.get_results_given_testing_data_and_theta(maker, params, recalculate, self.theta, which_infer, testing_list, num_states)
            if len(self.results.scores) == 0:
                raise Exception("couldn't test")
        except Exception as e:
            print e
            self.set_faulty()

        else:
            self.set_not_faulty()
        #pdb.set_trace()

    def is_faulty(self):
        return self.is_broken

    def set_faulty(self):
        self.is_broken = True

    def set_not_faulty(self):
        self.is_broken = False

    def get_hp(self):
        return self.hp

    def get_fold(self):
        return self.the_fold

    def get_results(self):
        if self.is_faulty():
            raise Exception("faulty train_test results with hyperparameters: " + str(self.get_hp()))
        return self.results

    def get_theta(self):
        return self.theta

# have to set ifold when calling it
# had ifold
class cv_results(keyed_object):

    def __init__(self, maker, params, recalculate):


        self.the_data = maker.get_param(params,'f')
        self.k = maker.get_param(params, 'ik')

        self.train_test_result_list = []
        the_folds = self.the_data.get_folds(maker, params, recalculate, self.k)
        for a_fold in the_folds:

            print a_fold, maker.get_param(params, 'hp')
            maker.set_param(params, 'f', a_fold)
#            pdb.set_trace()
            result = maker.get_var_or_file(objects.bwW, params, recalculate, True, False, False)



            
            self.train_test_result_list.append(result)
        #pdb.set_trace()

    def is_faulty(self):
        for result in self.train_test_result_list:
            if result.is_faulty():
                return True
        return False

    def get_hp(self):
        return self.train_test_result_list[0].get_hp()

    def get_k(self):
        return self.k

    def get_data(self):
        return self.data

    def get_train_test_result_list(self):
        return self.get_train_test_result_list

    def get_raw_results(self):

        raw = pdb_results_struct()
        for result in self.train_test_result_list:
            raw = raw + result.get_results()
        return raw


    def get_metric(self, maker, params, recalculate):
        #raw = self.get_raw_results()
        maker.set_param(params, 'md', 3)
        #maker.set_param(params, 'results', raw)
        
        maker.set_param(params, 'hp', self.get_hp())
        metric_data = maker.get_var_or_file(objects.ahW, params, False, False, False, False)
        #assert(len(metric_data) == 1)


        which_distance = maker.get_param(params, 'wdm')

        metric = metric_data[0][which_distance]

        #pdb.set_trace()
        return metric
        

def cross_product(stuff):
    if len(stuff) == 1:
        return [ [x] for x in stuff[0]]
    else:
        first = stuff[0]
        second = cross_product(stuff[1:])
        results = []
        for it in second:
            for a in first:
                results.append([a] + it)
        return results
    

class hp_search_results(keyed_object):

    def add_cv_result(self, hp, a_cv_result):
        self.result_keys.append(hp)
        self.cv_result_list.append(a_cv_result)

    def __init__(self, maker, params, recalculate):
        #pdb.set_trace()
        self.result_keys = []
        self.cv_result_list = []
        self.the_hp_searcher = maker.get_param(params, 'hps')
        
        self.the_hp_searcher.set_hp_search_results(self)
        self.the_hp_searcher.turn_on_search()
        i = 0
        while not self.the_hp_searcher.done_searching():
            #pdb.set_trace()
            next_hp = self.the_hp_searcher.get_next_hp()
            maker.set_param(params, 'hp', next_hp)
            print i, next_hp
            i += 1
            #a_cv_result = cv_results(maker, params, recalculate)
            try:
                a_cv_result = maker.get_var_or_file(objects.cdW, params, recalculate, False, False, False)
                self.add_cv_result(next_hp, a_cv_result)
                # also get ahW corresponding to this cv_result so it is easy to get the best one later on
                maker.set_param(params, 'md', 3)
                temp = maker.get_var_or_file(objects.ahW, params, recalculate, True, True, False)
            except:
                print 'couldnt gets results', next_hp
        #pdb.set_trace()

    def get_cv_result(self, hp):
        pass

    def get_best_hp(self, maker, params, recalculate):
        best = -1
        best_idx = -1

        for i in range(len(self.result_keys)):
            print i, 'getting best'

            try:
                temp = self.cv_result_list[i].get_metric(maker, params, recalculate)
                if temp > best:
                    best_idx = i
                    best = temp
            except:
                print "failed at ", i
        assert best != -1

        return self.result_keys[best_idx]


# probably set outer_k to 3 and inner_k to 2
class overall_results(keyed_object):

    def __init__(self, maker, params, recalculate):



#        pdb.set_trace()


        
        self.the_data = maker.get_var_or_file(objects.brW, params, recalculate, True, False, False)
        self.train_test_result_list = []
        self.best_hps = []
        outer_k = maker.get_param(params, 'ok')
        folds = self.the_data.get_folds(maker, params, recalculate, outer_k)
        for a_fold in folds:

            maker.set_param(params, 'f', a_fold)
            maker.set_param(params, 'wj', 0)
            maker.set_param(params, 'nj', 1)
            the_searcher = hp_searcher(maker, params, recalculate)
            maker.set_param(params, 'hps', the_searcher)

            a_hp_search_result = hp_search_results(maker, params, recalculate)

            best_hp = a_hp_search_result.get_best_hp(maker, params, recalculate)

            self.best_hps.append(best_hp)

            maker.set_param(params, 'hp', best_hp)
            result = maker.get_var_or_file(objects.bwW, params, recalculate, True, False, False)
            self.train_test_result_list.append(result)
        #pdb.set_trace()

    def get_raw_results(self):

        raw = pdb_results_struct()
        for result in self.train_test_result_list:
            raw = raw + result.get_results()
        return raw

# make sure that every time i set a fold, the wrapper that created it doesn't set any other folds

