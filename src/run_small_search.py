import _test

import wc

import new_new_objects as objects
import helper
import sys
import pdb
import constants

from param import param

input: number of outer folds, which fold, total size to split hp_values into, which job it is
num_outer_fold = int(sys.argv[1])
which_outer_fold = int(sys.argv[2])
num_pieces = int(sys.argv[3])
which_piece = int(sys.argv[4])

#info_file = 'mf0'


info_file = sys.argv[5]
#info_file = constants.INFO_FOLDER + info_file

relative_folder, the_params = helper.read_param(info_file)



#num_outer_fold = 3
#which_outer_fold = 0
#num_pieces = 6
#which_piece = 0

the_params.set_param('tj', num_pieces)
the_params.set_param('wj', which_piece)


the_data = wc.get_stuff(objects.brW, the_params, True, False, False)
the_params.set_param('s', the_data)
the_params.set_param('m', num_outer_fold)
the_params.set_param('k', which_outer_fold)






the_fold = wc.get_stuff(objects.buW, the_params, True, False, False)

import cross_validation_pseudo as cv

the_params.set_param('f', the_fold)
asdf = wc.get_stuff(objects.cbW, the_params, False, False, False)
