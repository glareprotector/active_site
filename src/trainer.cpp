/*
Provide:
1. function that takes in list of pdb_names, hyper_parameters, and trains, eventually returning a theta vector.  this will require (among other things) adding to model the ability to use those hyperparameters.  will also require allowing a constructor to specify explicitly the test/training names.  a minimizer object will be created, but report will not be called
2. function that, for a given list of pdb_names, hyper_parameters, theta, gives results on those pdb_names given the theta and hyper_parameters.  all pdb_names go into 'testing'
what would wrapper do?  it would call cpp code(which generates results file via object wrapper that pickles).  then would call the same wrapper since result is pickled, and return that.
alternatively, python function would call results wrapper, which in turn calls model code that calls report and returns results in c++ form, which would then be converted using typemap to python form
