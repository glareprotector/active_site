import pdb
pdb.set_trace()
import f
import new_new_objects as objects

import wrapper
from wrapper_decorator import dec
import wc
import global_stuff

import sys
info_file = sys.argv[1]

import helper
asdf, the_params = helper.read_param(info_file)


the_params.set_param('which_wrapperq', objects.fW)
wc.get_stuff(objects.abW, the_params, True, False, False)
