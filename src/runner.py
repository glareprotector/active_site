import subprocess

num_outer_folds = 3
num_pieces = 6

import pdb

for i in range(num_outer_folds):
    for j in range(num_pieces):


        log = 'bin/' + str(i) + '_' + str(j)
        error = log + 'err'
        f = open(log,'w')
        g = open(error,'w')
        
        subprocess.Popen(['python','run_small_search.py',str(num_outer_folds),str(i),str(num_pieces),str(j)], stdout=f,stderr=g)


