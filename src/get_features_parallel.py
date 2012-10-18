import wc

import new_new_objects as objects
import global_stuff
import helper
from mpi4py import MPI
import pdb
import helper
import sys

file_location = sys.argv[1]


folder_name, the_params = helper.read_param(file_location)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

f = open('bin/get_features_parallel_'+str(rank), 'w', 0)

data_list = wc.get_stuff(objects.mW, the_params, global_stuff.recalculate, True, True)
works = []
for i in range(len(data_list)):
    if i % size == rank:
        pdb_name = data_list[i].pdb_name
        chain_letter = data_list[i].chain_letter
        try:

            the_params.set_param('pdb_name', pdb_name)
            the_params.set_param('chain_letter', chain_letter)
            print rank, pdb_name, 'start jw'
            f.write(str(rank) + ' ' + pdb_name + ' start jw\n')
            wc.get_stuff(objects.jW, the_params, global_stuff.recalculate, True, True, False)
            print rank, pdb_name, 'start kw'
            f.write(str(rank) + ' ' + pdb_name + ' start kw\n')
            wc.get_stuff(objects.kW, the_params, global_stuff.recalculate, True, True, False)
            print rank, pdb_name, 'end kw'
            f.write(str(rank) + ' ' + pdb_name + ' end kw\n')
            works.append(pdb_name + '_' + chain_letter)
        except:
            print 'error in ', pdb_name, chain_letter
            f.write(str(rank) + ' ' + pdb_name + ' error\n')
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
        
