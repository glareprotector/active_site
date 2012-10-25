import wc

import new_new_objects as objects
import global_stuff
import helper

import pdb
import helper
import sys
from param import param

file_location = sys.argv[1]


wrappers = [objects.bhW, objects.cfW]
to_pickle = [True, True]

folder_name, the_params = helper.read_param(file_location)


the_params.set_param('tj',1)
the_params.set_param('wj',0)

hp_stash = wc.get_stuff(objects.caW, the_params, False, False, False)




f = open('bin/get_features_serial', 'w', 0)

data_list = wc.get_stuff(objects.ciW, the_params, False, True, True)
works = []

all_keys = [ ['nvjd'],  ['wjd', 'wpw']]


for i in range(len(data_list)):

    pdb_name = data_list[i].pdb_name
    chain_letter = data_list[i].chain_letter
    start = data_list[i].start
    end = data_list[i].end
    k = 0




    for wrapper in wrappers:

        """
        pruned_hp_stash = set()
        pdb.set_trace()
        all_keys = wc.get_wrapper_instance(wrapper).all_keys_cache.get(False)
        """
        pruned_hp_stash = set()
        
        for hp in hp_stash:

            pruned_hp_stash.add(param.restriction(hp, all_keys[k]))


        j = 0
        for hp in pruned_hp_stash:
            print j, hp
            j += 1
            the_params.set_param('hp', hp)
            the_params.flatten_hp(wc.useless)
            
            try:
                
              
            


                the_params.set_param('p', pdb_name)
                the_params.set_param('c', chain_letter)
                the_params.set_param('st', start)
                the_params.set_param('en', end)
                print pdb_name, wrapper
                start_string = pdb_name + ' start ' + str(wrapper) + '\n'
                f.write(start_string)
                print start_string

                wc.get_stuff(wrapper, the_params, False, to_pickle[k], False, False)

                end_string = pdb_name + ' end ' + str(wrapper) + '\n'
                

            except Exception, err:

                print 'error in ', pdb_name, chain_letter, err
                f.write(pdb_name + ' error\n')
                raise
            else:
                works.append(pdb_name + ',' + chain_letter + ',' + str(start) + ',' + str(end))

        k += 1

