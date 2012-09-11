from wc import *

import new_new_objects as objects

import wrapper
from wrapper_decorator import dec

from sklearn.cross_validation import StratifiedKFold
from sklearn import svm
import numpy
import global_stuff

# takes in data_list, num_folds, fold_number and returns the trained SVM model
class asW(wrapper.obj_wrapper, wrapper.experiment_results_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = False, always_recalculate = False, olb_obj = None):

        all_features, all_classes = self.get_var_or_file(objects.arW, params, recalculate, True, True, False)
        # get the indicies in training
        pdb.set_trace()
        y = numpy.array(all_classes)
        X = numpy.array(all_features)
        cv = list(enumerate(StratifiedKFold(y, k = self.get_param(params, 'nfld'))))
        which_fold = self.get_param(params, 'wfld')
        train_idx = cv[which_fold][1][0]
        svmC = self.get_param(params, 'svmC')
        classifier = svm.SVC(kernel = 'linear', probability = True, C = svmC)
        classifier.fit(X[train_idx], y[train_idx])
        return classifier
        

# takes in data_list, num_folds, fold_number, and return the roc curve input as object(but it will be filelized)
class atW(wrapper.mat_obj_wrapper, wrapper.experiment_results_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = True, to_filelize = True, always_recalculate = False, olb_obj = None):

        # get the model
        pdb.set_trace()
        model = self.get_var_or_file(asW, params, recalculate, True, False, False)
        all_features, all_classes = self.get_var_or_file(objects.arW, params, recalculate, True, True, False)
        X = numpy.array(all_features)
        y = numpy.array(all_classes)
        which_fold = self.get_param(params,'wfld')
        cv = list(enumerate(StratifiedKFold(y, k = self.get_param(params,'nfld'))))
        test_idx = cv[which_fold][1][1]
        model = self.get_var_or_file(asW, params, recalculate, True, False, False)
        probs = model.predict_proba(X[test_idx])[:,1]

        pdb.set_trace()
        # convert probs, true_classes back to regular lists so i can pass to
        y_test = y[test_idx]
        score = model.score(X[test_idx],y_test)
        self.set_param(params, 'obj_val', score)
        y_list = y_test.tolist()
        probs_list = probs.tolist()
#        pdb.set_trace()
        return global_stuff.get_transpose([y_list, probs_list])
