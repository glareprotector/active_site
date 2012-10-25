import new_new_objects as objects


class data(object):

    def __init__(self, maker, params, recalculate):
        self.data_list_file = maker.get_param(params, 'data_list_file')
        self.names = maker.get_var_or_file(params, objects.mW, recalculate, True, True)

    def get_training_names(self):
        return self.names


class fold(data):

    def __init__(self, maker, params, recalculate):
        self.source = maker.get_param('fsc')
        self.k = maker.get_param(params, 'fldk')
        self.which_fold = maker.get_param(params, 'wfldk')
        self.params = params.get_copy()

    def get_training_names(self):
        source_training_names = self.source.get_training_indicies()
        return [source_training_names[i] for i in range(len(source_training_names)) if i % self.k == self.which_fold]

    def get_testing_names(self):
        source_training_names = self.source.get_training_names()
        return [source_training_names[i] for i in range(len(source_training_names)) if i % self.k != self.which_fold]


class theta_for_fixed_H_fixed_fold(object):

    def __init__(self, maker, params, recalculate):
        self.fold = maker.get_param(params, 'finst')
        self.H = maker.get_param(params, 'h')
        # call c++ code to get theta, or rather, submit a job.  without orchestra, can just wait for result.  


        
class results_for_fixed_H_fixed_fold(object):

    def __init__(self, maker, params, recalculate):
        self.fold = maker.get_param(params, 'finst')
        self.H = maker.get_param(params, 'h')
        theta = self.get_var_or_file(params, objects.theta_for_fixed_H_fixed_fold_wrapper, recalculate, True, False)


# has function that, given H and data/fold, performs a specified k fold CV.  stores these results
class results_for_fixed_H(object):

    def __init__(self, maker, params, recalculate):
        self.source = maker.get_param(params, 'dsc')
        self.k = maker.get_param(params, 'fldk')
        self.h = maker.get_param(params, 'h')
        self.params = params.get_copy()
        self.results = {}  # stores CV results as (fold, H) to result mapping
        
        for fold in self.source.get_fold_iter(self.k):
            self.results.add(
        

        
        # return data object
        source_indicies = self.source.get_indicies()
        training_indicies = [source_indicies[i] for i in range(len(source_indicies)) if i % self.k == self.which_fold]
        
        

class base_fold(keyed_object):

    # fold source will be a data object.  also specify fold number
    def __init__(self, maker, params, recalculate):
        self.fold_source = maker.get_param(params, 'fsc')
        self.k = maker.get_param(params, 'fldk')
        

class fold_iter(object):

    def __init__(self, the_data, k):
        self.the_data = the_data
        self.k = k
        self

    def __iter__(self):
        return self

    def next(self):
        


class data(keyed_object):

    def get_data_file(self):
        return self.data_file

    def get_indicies(self):
        

    def __init__(self, maker, params, recalculate):
        self.maker = maker
        self.data_file = maker.get_var_or_file(objects.mW, params, recalculate, True, True)

    def get_data(self):
        return self

    def get_fold_iter(self, k):
        return fold_iter(self, k)























all_indicies = range(num_samples)
F = get_folds(indicies, k)  # returns list of lists
performance_results = [] # dictionary of sample name to [ [positions], [marginals] ]

# testing/training should be specified as pdb_names
class fold(object):

    def __init__(self, all_indicies, k):
        pass

    # returns testing, training indicies
    def return_info(self):
        pass

# iterator that returns folds
class fold_iterator(object):

    def __init__(self, all_indicies, k):
        pass

    def __iter__(self):
        return self

    def next(self):
        return next_fold

class hyper_parameters(object):

    def __init__(self, names, values):
        pass

class data_set:

class 

class CV_results(object):

    def __init__(self, CV_params):

class CV_params:

    def __init__(self, data_set, k)

class performance_results_generator_specific_fold(object):

    def __init__(self, fold, hyperparameters):
        pass

    def 
    

# for a fixed hyperparameter, gets performance object for the result of CV on specified indicies.  results should include all indicies
class performance_results_generator_CV(object):

    def __init__(self, indicies, data_names, k, hyperparameters)


class hyper_parameter_iterator(object):

    def __init__(self, low_limits, high_limits, step_sizes, performance_logger):
        pass

    def __iter__(self):
        return self

    def next(self):
        pass


class performance_results_logger(


hp_iter = hyper_parameter_iterator(low_limits, high_limits, step_sizes)

for fold in F:
    get_performance_results(all_indicies, fold, data_set_name, performance_results)


def get_performance_results(all_indicies, fold, data_set_name, performance_results):
    
