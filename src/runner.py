import subprocess

import sys

num_outer_folds = int(sys.argv[1])
num_pieces = int(sys.argv[2])

info_file = sys.argv[3]

#num_outer_folds = 3
#num_pieces = 200

import pdb

for i in range(num_outer_folds):
    for j in range(num_pieces):


        log = 'bin/' + str(i) + '_' + str(j)
        error = log + 'err'
        f = open(log,'w')
        g = open(error,'w')
        subprocess.call(['bsub', '-R', 'rusage[mem=700]', '-o', log, '-q', 'shared_12h', 'python', 'run_small_search.py', str(num_outer_folds), str(i), str(num_pieces), str(j), info_file])
        #subprocess.Popen(['python','run_small_search.py',str(num_outer_folds),str(i),str(num_pieces),str(j)], stdout=f,stderr=g)


