
import _test

import wc

import new_new_objects as objects



import pdb

_test.init_crf()

from param import param


#import run_small_search
import helper

import sys
info_file = sys.argv[1]

relative_folder, the_params = helper.read_param(info_file)

the_params.set_param('tj',1)

results = wc.get_stuff(objects.ceW, the_params, False, True, True)

print 3
