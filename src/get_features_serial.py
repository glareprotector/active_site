import wc

import new_new_objects as objects
import global_stuff
import helper

import pdb
import helper
import sys

file_location = sys.argv[1]


wrappers = [objects.bhW, objects.cfW, objects.oW]
to_pickle = [False, False, False]

folder_name, the_params = helper.read_param(file_location)


the_params.set_param('tj',1)
the_params.set_param('wj',0)
pdb.set_trace()
hp_stash = wc.get_stuff(objects.caW, the_params, False, False, False)

pdb.set_trace()


f = open('bin/get_features_serial', 'w', 0)

data_list = wc.get_stuff(objects.mW, the_params, False, True, True)
works = []
for i in range(len(data_list)):

    pdb_name = data_list[i].pdb_name
    chain_letter = data_list[i].chain_letter
    start = data_list[i].start
    end = data_list[i].end
    k = 0
    for hp in hp_stash:
        print k
        k += 1
        the_params.set_param('hp', hp)

        
        try:

            j = 0
            for wrapper in wrappers:


                the_params.set_param('p', pdb_name)
                the_params.set_param('c', chain_letter)
                the_params.set_param('st', start)
                the_params.set_param('en', end)
                print pdb_name, wrapper
                start_string = pdb_name + ' start ' + str(wrapper) + '\n'
                f.write(start_string)
                print start_string

                wc.get_stuff(wrapper, the_params, False, to_pickle[j], False, False)

                end_string = pdb_name + ' end ' + str(wrapper) + '\n'
                j += 1

        except Exception, err:

            print 'error in ', pdb_name, chain_letter, err
            f.write(pdb_name + ' error\n')
        else:
            works.append(pdb_name + ',' + chain_letter + ',' + str(start) + ',' + str(end))



