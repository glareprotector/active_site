import wc

import new_new_objects as objects
import global_stuff
import helper
from mpi4py import MPI
import pdb
import helper
import sys

file_location = sys.argv[1]


wrappers = [objects.bhW, objects.cfW]


folder_name, the_params = helper.read_param(file_location)

the_params.set_param('tj',1)
the_params.set_param('wj',0)
hp_stash = wc.get_stuff(objects.caW, the_params, False, False, False)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

f = open('bin/get_features_parallel_'+str(rank), 'w', 0)

data_list = wc.get_stuff(objects.ciW, the_params, global_stuff.recalculate, True, True)
works = []
for i in range(len(data_list)):
    if i % size == rank:
        pdb_name = data_list[i].pdb_name
        chain_letter = data_list[i].chain_letter
        start = data_list[i].start
        end = data_list[i].end
        try:

            for wrapper in wrappers:

                for hp in hp_stash:

                    the_params.set_param('p', pdb_name)
                    the_params.set_param('c', chain_letter)
                    the_params.set_param('st', start)
                    the_params.set_param('en', end)
                    the_params.set_param('hp', hp)
                    the_params.flatten_hp(wc.useless)
                    print rank, pdb_name, wrapper
                    start_string = str(rank) + ' ' + pdb_name + ' start ' + str(wrapper) + '\n'
                    f.write(start_string)
                    print start_string

                    wc.get_stuff(wrapper, the_params, False, True, True, False)

                    end_string = str(rank) + ' ' + pdb_name + ' end ' + str(wrapper) + '\n'

        except Exception, err:

            print 'error in ', pdb_name, chain_letter, err
            f.write(str(rank) + ' ' + pdb_name + ' error\n')
        else:
            works.append(pdb_name + ',' + chain_letter + ',' + str(start) + ',' + str(end))


f.write('before barrier')
comm.Barrier()
f.write('past barrier')
print rank, works
root_works = comm.gather(works, root=0)

if rank == 0:
    ans = []
    for elt in root_works:
        ans = ans + elt
    helper.write_mat([[x] for x in ans], 'catres_works.pdb_list')
        
